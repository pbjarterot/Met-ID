#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
import random
from pathlib import Path

import torch
import torch.nn as nn
import torch.nn.functional as F

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors
from torch_geometric.data import Batch, Data
from torch_geometric.loader import DataLoader
from torch_geometric.nn import global_add_pool
from torch_geometric.nn.models.visnet import (
    GatedEquivariantBlock,
    ViSNetBlock,
)
from tqdm.auto import tqdm
import sys
import json



def bundled_resource(filename: str) -> Path:
    """Resolve a file bundled by PyInstaller."""
    bundle_directory = getattr(sys, "_MEIPASS", None)

    if bundle_directory is not None:
        return Path(bundle_directory) / filename

    return Path(filename)


REQUIRED_COLUMNS = {"ccs", "m/z", "smiles", "charge"}
PREDICTION_REQUIRED_COLUMNS = {"smiles"}
DESCRIPTOR_NAMES = ["m/z", "charge", "estimated_polarizability_A3"]
DEFAULT_FALLBACK_ADDUCT = "[M+H]+"

ELECTRON_MASS_DA = 0.000548579909
PROTON_MASS_DA = 1.007276466621
SODIUM_CATION_MASS_DA = 22.989220702
POTASSIUM_CATION_MASS_DA = 38.963157906
ADDUCT_SPECS: dict[str, tuple[int, float, int]] = {
    "[M+H]+": (1, PROTON_MASS_DA, 1),
    "[M+Na]+": (1, SODIUM_CATION_MASS_DA, 1),
    "[M+K]+": (1, POTASSIUM_CATION_MASS_DA, 1),
    "[M]+": (1, -ELECTRON_MASS_DA, 1),
    "[M-H]-": (1, -PROTON_MASS_DA, -1),
    "[M]-": (1, ELECTRON_MASS_DA, -1),
    "[M+2H]2+": (1, 2.0 * PROTON_MASS_DA, 2),
    "[M-2H]2-": (1, -2.0 * PROTON_MASS_DA, -2),
    "[2M+H]+": (2, PROTON_MASS_DA, 1),
}

# Crippen.MolMR is a fragment-based estimate of molar refractivity in cm^3/mol.
# The Lorentz-Lorenz relation gives molecular polarizability volume as:
# alpha [A^3] = MR * 3e24 / (4*pi*N_A).
AVOGADRO_CONSTANT = 6.02214076e23
MOLAR_REFRACTIVITY_TO_POLARIZABILITY_A3 = (
    3.0e24 / (4.0 * math.pi * AVOGADRO_CONSTANT)
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=Path("ccs.csv"))
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--lr", type=float, default=2e-4)
    parser.add_argument("--test-fraction", type=float, default=0.20)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--checkpoint", type=Path, default=bundled_resource("trained_weights.pt"))
    parser.add_argument("--use-adduct", action="store_true", help="Train with a learned categorical adduct embedding.")
    parser.add_argument("--fallback-adduct", default=DEFAULT_FALLBACK_ADDUCT,
                        help=("Adduct used when prediction input is missing or unknown."
                              "For adduct-aware training it must occur in the training CSV."))
    parser.add_argument("--predict-csv", type=Path, default=None, help="Predict an unlabeled CSV instead of training a model.")
    parser.add_argument("--prediction-output", type=Path, default=Path("visnet_ccs_predictions.csv"), help="Output used with --predict-csv.")
    parser.add_argument("--predictions",type=Path, default=Path("visnet_ccs_test_predictions.csv"))
    parser.add_argument(
        "--sidecar-mode",
        action="store_true",
        help="Run as a persistent JSON-lines prediction process.",
    )
    return parser.parse_args()


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def estimate_polarizability_a3(mol: Chem.Mol) -> float:
    """Estimate molecular polarizability volume from Crippen MolMR."""
    molar_refractivity = float(Crippen.MolMR(mol))
    polarizability_a3 = (
        molar_refractivity * MOLAR_REFRACTIVITY_TO_POLARIZABILITY_A3
    )
    if not math.isfinite(polarizability_a3) or polarizability_a3 <= 0:
        raise ValueError(
            f"Invalid estimated polarizability: {polarizability_a3} A^3"
        )
    return polarizability_a3

def calculate_monoisotopic_mz(
    smiles: str,
    adduct: str,
) -> tuple[float, int]:
    """
    Calculate monoisotopic m/z from a parent SMILES and an adduct.

    Returns:
        (mz, charge)
    """
    smiles = str(smiles).strip()
    adduct = str(adduct).strip()

    if not smiles:
        raise ValueError("SMILES cannot be empty")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")

    if adduct not in ADDUCT_SPECS:
        raise ValueError(
            f"Unsupported adduct {adduct!r}. "
            f"Supported adducts: {sorted(ADDUCT_SPECS)}"
        )

    molecule_count, mass_shift, charge = ADDUCT_SPECS[adduct]

    neutral_monoisotopic_mass = float(Descriptors.ExactMolWt(mol))

    ion_mass = (
        molecule_count * neutral_monoisotopic_mass
        + mass_shift
    )

    mz = ion_mass / abs(charge)

    if not math.isfinite(mz) or mz <= 0:
        raise ValueError(
            f"Calculated invalid m/z={mz} for "
            f"SMILES={smiles!r}, adduct={adduct!r}"
        )

    return mz, charge

def smiles_to_3d_graph(
    smiles: str, seed: int
) -> tuple[torch.Tensor, torch.Tensor, float]:
    """Return atoms, one centered conformer, and polarizability in A^3."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit could not parse the SMILES")

    # Calculate from the sanitized molecule. Crippen.MolMR accounts for normal
    # implicit hydrogens and is independent of the generated conformer.
    polarizability_a3 = estimate_polarizability_a3(mol)

    # ViSNet can include hydrogens as ordinary atoms. This adds normal implicit
    # hydrogens; it does not add the extra proton represented by [M+H]+.
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.useRandomCoords = True
    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        raise ValueError("RDKit conformer embedding failed")

    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    else:
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)

    conformer = mol.GetConformer()
    pos = np.asarray(conformer.GetPositions(), dtype=np.float32)
    pos -= pos.mean(axis=0, keepdims=True)

    z = np.asarray([atom.GetAtomicNum() for atom in mol.GetAtoms()], dtype=np.int64)
    return torch.from_numpy(z), torch.from_numpy(pos), polarizability_a3


def clean_adduct(value: object, fallback_adduct: str) -> str:
    """Return a stripped adduct label, replacing empty/NaN values."""
    if pd.isna(value):
        return fallback_adduct
    adduct = str(value).strip()
    return adduct if adduct else fallback_adduct


def load_graphs(
    csv_path: Path,
    seed: int,
    use_adduct: bool,
    fallback_adduct: str,
) -> list[Data]:
    frame = pd.read_csv(csv_path)
    frame.columns = [str(column).strip() for column in frame.columns]
    required_columns = set(REQUIRED_COLUMNS)
    if use_adduct:
        required_columns.add("adduct")

    missing = required_columns.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing CSV columns: {sorted(missing)}")

    if "adduct" not in frame.columns:
        frame["adduct"] = "[M+H]+"#fallback_adduct

    frame = frame.dropna(subset=list(required_columns)).reset_index(drop=True)

    graphs: list[Data] = []
    rows = tqdm(
        frame.iterrows(),
        total=len(frame),
        desc="Generating training conformers",
        unit="mol",
        dynamic_ncols=True,
    )
    for row_index, row in rows:
        smiles = str(row["smiles"])
        try:
            z, pos, polarizability_a3 = smiles_to_3d_graph(
                smiles, seed + row_index
            )
        except Exception as exc:
            print(f"Skipping CSV row {row_index}: {exc}")
            continue

        ccs = float(row["ccs"])
        if not math.isfinite(ccs) or ccs <= 0:
            print(f"Skipping CSV row {row_index}: invalid CCS {ccs}")
            continue

        name = str(row.get("Molecule Name", f"row_{row_index}"))
        adduct = clean_adduct(row.get("adduct"), fallback_adduct)

        mz = float(row["m/z"])
        charge = float(row["charge"])

        # Shape [1, 3] makes PyG concatenate descriptors into [batch, 3].
        descriptors = torch.tensor(
            [[mz, charge, polarizability_a3]],
            dtype=torch.float32,
        )

        graphs.append(
            Data(
                z=z,
                pos=pos,
                descriptors=descriptors,
                y=torch.tensor([math.log(ccs)], dtype=torch.float32),
                row_index=torch.tensor([row_index], dtype=torch.long),
                input_mz=torch.tensor([mz], dtype=torch.float32),
                input_charge=torch.tensor([charge], dtype=torch.float32),
                input_polarizability=torch.tensor(
                    [polarizability_a3], dtype=torch.float32
                ),
                molecule_name=name,
                smiles=smiles,
                adduct=adduct,
            )
        )

    if len(graphs) < 4:
        raise ValueError(f"Only {len(graphs)} valid molecules were generated")

    return graphs


def load_prediction_graphs(
    csv_path: Path,
    seed: int,
    adduct_to_id: dict[str, int] | None = None,
    descriptor_dim: int = len(DESCRIPTOR_NAMES),
    fallback_adduct: str = DEFAULT_FALLBACK_ADDUCT,
) -> list[Data]:
    """Generate graphs for a CSV that does not need a CCS column."""
    frame = pd.read_csv(csv_path)
    frame.columns = [str(column).strip() for column in frame.columns]

    # Preserve supplied values. Defaults are used only when a column/value is
    # missing, rather than overwriting the prediction input.
    if "adduct" not in frame.columns:
        print(
            "Prediction CSV has no adduct column; "
            f"using fallback {fallback_adduct}."
        )
        frame["adduct"] = fallback_adduct

    frame["adduct"] = [
        clean_adduct(value, fallback_adduct) for value in frame["adduct"]
    ]

    required_columns = set(PREDICTION_REQUIRED_COLUMNS)

    missing = required_columns.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing prediction CSV columns: {sorted(missing)}")

    if descriptor_dim not in (2, 3):
        raise ValueError(
            f"Unsupported checkpoint descriptor_dim={descriptor_dim}; "
            "expected 2 (legacy) or 3 (with polarizability)."
        )

    frame = frame.dropna(subset=list(required_columns)).copy()
    if adduct_to_id is not None:
        unknown_adducts = sorted(
            set(frame["adduct"].astype(str)).difference(adduct_to_id)
        )
        if unknown_adducts:
            if fallback_adduct not in adduct_to_id:
                raise ValueError(
                    f"Unknown prediction adducts: {unknown_adducts}. "
                    f"Fallback {fallback_adduct!r} is not available in this "
                    f"checkpoint. Available adducts: {sorted(adduct_to_id)}"
                )
            print(
                f"Unknown prediction adducts {unknown_adducts}; "
                f"using fallback {fallback_adduct}."
            )
            frame.loc[
                ~frame["adduct"].astype(str).isin(adduct_to_id), "adduct"
            ] = fallback_adduct

    graphs: list[Data] = []
    rows = tqdm(
        frame.iterrows(),
        total=len(frame),
        desc="Generating prediction conformers",
        unit="mol",
        dynamic_ncols=True,
    )
    for row_index, row in rows:
        smiles = str(row["smiles"])
        try:
            z, pos, polarizability_a3 = smiles_to_3d_graph(
                smiles, seed + int(row_index)
            )
        except Exception as exc:
            print(f"Skipping prediction CSV row {row_index}: {exc}")
            continue

        adduct = clean_adduct(
            row.get("adduct"),
            fallback_adduct,
        )

        try:
            mz, charge = calculate_monoisotopic_mz(
                smiles=smiles,
                adduct=adduct,
            )
        except Exception as exc:
            print(
                f"Skipping prediction CSV row {row_index}: "
                f"could not calculate m/z: {exc}"
            )
            continue

        name_value = row.get("Molecule Name", row.get("name", f"row_{row_index}"))
        name = f"row_{row_index}" if pd.isna(name_value) else str(name_value)
        descriptor_values = [mz, charge]
        if descriptor_dim == 3:
            descriptor_values.append(polarizability_a3)
        graph = Data(
            z=z,
            pos=pos,
            descriptors=torch.tensor([descriptor_values], dtype=torch.float32),
            row_index=torch.tensor([row_index], dtype=torch.long),
            input_mz=torch.tensor([mz], dtype=torch.float32),
            input_charge=torch.tensor([charge], dtype=torch.float32),
            input_polarizability=torch.tensor(
                [polarizability_a3], dtype=torch.float32
            ),
            molecule_name=name,
            smiles=smiles,
            adduct=adduct,
        )
        if adduct_to_id is not None:
            graph.adduct_id = torch.tensor(
                [adduct_to_id[adduct]], dtype=torch.long
            )
        graphs.append(graph)

    if not graphs:
        raise ValueError("No valid prediction molecules were generated")
    return graphs


class ViSNetCCS(nn.Module):
    """One-conformer ViSNet with optional adduct conditioning."""

    def __init__(
        self,
        descriptor_dim: int,
        hidden: int = 64,
        num_adducts: int | None = None,
    ):
        super().__init__()

        self.representation = ViSNetBlock(
            hidden_channels=hidden,
            num_layers=3,
            num_heads=4,
            num_rbf=32,
            cutoff=5.0,
            max_num_neighbors=32,
            vertex=False,
        )
        self.atom_readout = GatedEquivariantBlock(
            hidden_channels=hidden,
            out_channels=hidden,
            intermediate_channels=hidden,
            scalar_activation=True,
        )

        self.geometry_mlp = nn.Sequential(
            nn.Linear(hidden + 1, hidden),
            nn.SiLU(),
            nn.LayerNorm(hidden),
        )
        self.descriptor_mlp = nn.Sequential(
            nn.Linear(descriptor_dim, 32),
            nn.SiLU(),
            nn.LayerNorm(32),
        )
        self.adduct_embedding = (
            nn.Embedding(num_adducts, 8) if num_adducts is not None else None
        )
        fusion_dim = hidden + 32 + (8 if self.adduct_embedding is not None else 0)
        self.regressor = nn.Sequential(
            nn.Linear(fusion_dim, 64),
            nn.SiLU(),
            nn.Dropout(0.10),
            nn.Linear(64, 1),
        )

    def forward(self, batch: Batch) -> torch.Tensor:
        scalar, vector = self.representation(batch.z, batch.pos, batch.batch)
        scalar, _ = self.atom_readout(scalar, vector)

        # Add pooling retains molecular-size information. Dividing by sqrt(N)
        # keeps its scale better behaved; log(N) explicitly supplies atom count.
        num_graphs = batch.descriptors.size(0)
        num_atoms = torch.bincount(
            batch.batch, minlength=num_graphs
        ).to(scalar.dtype).unsqueeze(1)
        pooled = global_add_pool(scalar, batch.batch) / num_atoms.sqrt().clamp_min(1)
        geometry = self.geometry_mlp(
            torch.cat([pooled, torch.log1p(num_atoms)], dim=1)
        )

        descriptors = self.descriptor_mlp(batch.descriptors)
        features = [geometry, descriptors]
        if self.adduct_embedding is not None:
            features.append(self.adduct_embedding(batch.adduct_id))
        fused = torch.cat(features, dim=1)
        return self.regressor(fused).squeeze(1)


class CcsPredictor:
    """Loaded ViSNet model for repeated single-molecule prediction."""

    def __init__(
        self,
        model: ViSNetCCS,
        device: torch.device,
        descriptor_mean: torch.Tensor,
        descriptor_std: torch.Tensor,
        log_ccs_mean: torch.Tensor,
        log_ccs_std: torch.Tensor,
        descriptor_dim: int,
        seed: int,
        adduct_to_id: dict[str, int] | None,
        fallback_adduct: str,
    ) -> None:
        self.model = model
        self.device = device

        # Descriptor normalization stays on CPU while constructing the graph.
        self.descriptor_mean = descriptor_mean.cpu()
        self.descriptor_std = descriptor_std.cpu().clamp_min(1e-8)

        self.log_ccs_mean = log_ccs_mean.to(device)
        self.log_ccs_std = log_ccs_std.to(device)

        self.descriptor_dim = descriptor_dim
        self.seed = seed
        self.adduct_to_id = adduct_to_id
        self.fallback_adduct = fallback_adduct

    def _resolve_adduct(self, adduct: str | None) -> str:
        """Clean the adduct and replace unknown values when necessary."""

        if adduct is None:
            cleaned_adduct = self.fallback_adduct
        else:
            cleaned_adduct = str(adduct).strip()

            if not cleaned_adduct:
                cleaned_adduct = self.fallback_adduct

        # For an adduct-agnostic checkpoint, the adduct is metadata only.
        if self.adduct_to_id is None:
            return cleaned_adduct

        if cleaned_adduct in self.adduct_to_id:
            return cleaned_adduct

        if self.fallback_adduct not in self.adduct_to_id:
            raise ValueError(
                f"Unknown adduct {cleaned_adduct!r}, and fallback adduct "
                f"{self.fallback_adduct!r} is not available in the checkpoint. "
                f"Available adducts: {sorted(self.adduct_to_id)}"
            )

        # Match the behavior of load_prediction_graphs().
        return self.fallback_adduct

    @torch.inference_mode()
    def predict(
        self,
        smiles: str,
        adduct: str | None = None,
    ) -> float:
        """Predict CCS after calculating monoisotopic m/z from SMILES."""

        smiles = str(smiles).strip()

        if not smiles:
            raise ValueError("SMILES cannot be empty")

        resolved_adduct = self._resolve_adduct(adduct)

        mz, charge = calculate_monoisotopic_mz(
            smiles=smiles,
            adduct=resolved_adduct,
        )

        z, pos, polarizability_a3 = smiles_to_3d_graph(
            smiles=smiles,
            seed=self.seed,
        )

        descriptor_values = [
            mz,
            float(charge),
        ]

        if self.descriptor_dim == 3:
            descriptor_values.append(polarizability_a3)

        descriptors = torch.tensor(
            [descriptor_values],
            dtype=torch.float32,
        )

        descriptors = (
            descriptors - self.descriptor_mean
        ) / self.descriptor_std

        graph = Data(
            z=z,
            pos=pos,
            descriptors=descriptors,
        )

        if self.adduct_to_id is not None:
            graph.adduct_id = torch.tensor(
                [self.adduct_to_id[resolved_adduct]],
                dtype=torch.long,
            )

        batch = Batch.from_data_list([graph]).to(self.device)

        prediction_z = self.model(batch)

        prediction_ccs = torch.exp(
            prediction_z * self.log_ccs_std
            + self.log_ccs_mean
        )

        result = float(prediction_ccs.squeeze().cpu())

        if not math.isfinite(result):
            raise RuntimeError("The model produced a non-finite CCS prediction")

        return result

def load_predictor(
    checkpoint_path: Path,
    seed: int,
) -> CcsPredictor:
    """Load the checkpoint and construct a reusable predictor."""

    checkpoint_path = Path(checkpoint_path)

    if not checkpoint_path.exists():
        raise FileNotFoundError(
            f"Checkpoint not found: {checkpoint_path}"
        )

    try:
        checkpoint = torch.load(
            checkpoint_path,
            map_location="cpu",
            weights_only=True,
        )
    except TypeError:
        # Compatibility with older PyTorch versions.
        checkpoint = torch.load(
            checkpoint_path,
            map_location="cpu",
        )

    required_keys = {
        "model_state_dict",
        "model_kwargs",
        "descriptor_mean",
        "descriptor_std",
        "log_ccs_mean",
        "log_ccs_std",
    }

    missing_keys = required_keys.difference(checkpoint)

    if missing_keys:
        raise ValueError(
            f"Checkpoint is missing keys: {sorted(missing_keys)}"
        )

    model_kwargs = dict(checkpoint["model_kwargs"])

    descriptor_dim = int(model_kwargs["descriptor_dim"])

    if descriptor_dim not in (2, 3):
        raise ValueError(
            f"Unsupported descriptor_dim={descriptor_dim}; "
            "expected 2 or 3"
        )

    uses_adduct = model_kwargs.get("num_adducts") is not None

    adduct_to_id: dict[str, int] | None = None

    if uses_adduct:
        stored_mapping = checkpoint.get("adduct_to_id")

        if stored_mapping is None:
            raise ValueError(
                "The checkpoint is adduct-aware but does not contain "
                "adduct_to_id"
            )

        adduct_to_id = {
            str(adduct): int(adduct_id)
            for adduct, adduct_id in stored_mapping.items()
        }

        expected_adducts = int(model_kwargs["num_adducts"])

        if len(adduct_to_id) != expected_adducts:
            raise ValueError(
                "Checkpoint adduct mapping size does not match "
                f"num_adducts: {len(adduct_to_id)} != "
                f"{expected_adducts}"
            )

    fallback_adduct = str(
        checkpoint.get(
            "fallback_adduct",
            DEFAULT_FALLBACK_ADDUCT,
        )
    ).strip()

    if not fallback_adduct:
        fallback_adduct = DEFAULT_FALLBACK_ADDUCT

    if (
        adduct_to_id is not None
        and fallback_adduct not in adduct_to_id
    ):
        raise ValueError(
            f"Fallback adduct {fallback_adduct!r} is absent from "
            f"the checkpoint mapping: {sorted(adduct_to_id)}"
        )

    descriptor_mean = torch.as_tensor(
        checkpoint["descriptor_mean"],
        dtype=torch.float32,
    ).reshape(1, -1)

    descriptor_std = torch.as_tensor(
        checkpoint["descriptor_std"],
        dtype=torch.float32,
    ).reshape(1, -1)

    if descriptor_mean.shape[1] != descriptor_dim:
        raise ValueError(
            "descriptor_mean dimension does not match model: "
            f"{descriptor_mean.shape[1]} != {descriptor_dim}"
        )

    if descriptor_std.shape[1] != descriptor_dim:
        raise ValueError(
            "descriptor_std dimension does not match model: "
            f"{descriptor_std.shape[1]} != {descriptor_dim}"
        )

    if not torch.isfinite(descriptor_mean).all():
        raise ValueError("descriptor_mean contains NaN or Inf")

    if not torch.isfinite(descriptor_std).all():
        raise ValueError("descriptor_std contains NaN or Inf")

    descriptor_std = descriptor_std.clamp_min(1e-8)

    log_ccs_mean = torch.as_tensor(
        checkpoint["log_ccs_mean"],
        dtype=torch.float32,
    ).reshape(())

    log_ccs_std = torch.as_tensor(
        checkpoint["log_ccs_std"],
        dtype=torch.float32,
    ).reshape(())

    if not torch.isfinite(log_ccs_mean):
        raise ValueError("log_ccs_mean is NaN or Inf")

    if not torch.isfinite(log_ccs_std) or log_ccs_std <= 0:
        raise ValueError("log_ccs_std must be positive and finite")

    device = torch.device(
        "cuda" if torch.cuda.is_available() else "cpu"
    )

    model = ViSNetCCS(**model_kwargs)
    model.load_state_dict(checkpoint["model_state_dict"])
    model = model.to(device).eval()

    return CcsPredictor(
        model=model,
        device=device,
        descriptor_mean=descriptor_mean,
        descriptor_std=descriptor_std,
        log_ccs_mean=log_ccs_mean,
        log_ccs_std=log_ccs_std,
        descriptor_dim=descriptor_dim,
        seed=seed,
        adduct_to_id=adduct_to_id,
        fallback_adduct=fallback_adduct,
    )

def split_indices(
    num_examples: int, test_fraction: float, seed: int
) -> tuple[list[int], list[int]]:
    if not 0.0 < test_fraction < 1.0:
        raise ValueError("--test-fraction must be between 0 and 1")

    generator = torch.Generator().manual_seed(seed)
    indices = torch.randperm(num_examples, generator=generator).tolist()
    num_test = max(1, round(num_examples * test_fraction))
    num_test = min(num_test, num_examples - 1)
    return indices[num_test:], indices[:num_test]


def split_indices_by_adduct(
    graphs: list[Data], test_fraction: float, seed: int
) -> tuple[list[int], list[int]]:
    """Split within each adduct so every learned adduct remains in training."""
    if not 0.0 < test_fraction < 1.0:
        raise ValueError("--test-fraction must be between 0 and 1")

    generator = torch.Generator().manual_seed(seed)
    groups: dict[str, list[int]] = {}
    for index, graph in enumerate(graphs):
        groups.setdefault(graph.adduct, []).append(index)

    train_indices: list[int] = []
    test_indices: list[int] = []
    for adduct, indices in sorted(groups.items()):
        order = torch.randperm(len(indices), generator=generator).tolist()
        shuffled = [indices[position] for position in order]

        if len(shuffled) == 1:
            print(
                f"Adduct {adduct!r} has one example; keeping it in training."
            )
            train_indices.extend(shuffled)
            continue

        num_test = max(1, round(len(shuffled) * test_fraction))
        num_test = min(num_test, len(shuffled) - 1)
        test_indices.extend(shuffled[:num_test])
        train_indices.extend(shuffled[num_test:])

    if not test_indices:
        raise ValueError(
            "Could not create a test split. At least one adduct needs two "
            "or more valid examples."
        )
    return train_indices, test_indices


def training_statistics(
    graphs: list[Data], train_indices: list[int]
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    train_descriptors = torch.cat(
        [graphs[index].descriptors for index in train_indices], dim=0
    )
    descriptor_mean = train_descriptors.mean(dim=0, keepdim=True)
    descriptor_std = train_descriptors.std(dim=0, keepdim=True).clamp_min(1e-8)

    train_y = torch.cat([graphs[index].y for index in train_indices])
    y_mean = train_y.mean()
    y_std = train_y.std().clamp_min(1e-8)
    return descriptor_mean, descriptor_std, y_mean, y_std


def train_model(
    model: nn.Module,
    loader: DataLoader,
    device: torch.device,
    y_mean: torch.Tensor,
    y_std: torch.Tensor,
    epochs: int,
    learning_rate: float,
) -> None:
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=learning_rate, weight_decay=1e-5
    )
    y_mean = y_mean.to(device)
    y_std = y_std.to(device)

    model.train()
    progress = tqdm(
        total=epochs * len(loader),
        desc="Training ViSNet",
        unit="batch",
        dynamic_ncols=True,
    )
    try:
        for epoch in range(1, epochs + 1):
            total_loss = 0.0
            total_examples = 0

            for batch in loader:
                batch = batch.to(device)
                target = (batch.y - y_mean) / y_std
                prediction = model(batch)
                loss = F.smooth_l1_loss(prediction, target, beta=0.5)

                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
                optimizer.step()

                total_loss += loss.item() * batch.num_graphs
                total_examples += batch.num_graphs
                running_loss = total_loss / total_examples

                progress.update(1)
                progress.set_postfix(
                    epoch=f"{epoch}/{epochs}",
                    loss=f"{running_loss:.5f}",
                    lr=f"{optimizer.param_groups[0]['lr']:.2e}",
                )
    finally:
        progress.close()


@torch.no_grad()
def evaluate(
    model: nn.Module,
    loader: DataLoader,
    device: torch.device,
    y_mean: torch.Tensor,
    y_std: torch.Tensor,
) -> pd.DataFrame:
    model.eval()
    rows: list[dict[str, object]] = []
    y_mean = y_mean.to(device)
    y_std = y_std.to(device)

    for batch in loader:
        batch = batch.to(device)
        prediction_z = model(batch)
        prediction_ccs = torch.exp(prediction_z * y_std + y_mean)
        true_ccs = torch.exp(batch.y)

        # String/list metadata stays on CPU when batch.to(device) is called.
        for index in range(batch.num_graphs):
            true_value = float(true_ccs[index].cpu())
            predicted_value = float(prediction_ccs[index].cpu())
            rows.append(
                {
                    "csv_row": int(batch.row_index[index].cpu()),
                    "Molecule Name": batch.molecule_name[index],
                    "smiles": batch.smiles[index],
                    "adduct": batch.adduct[index],
                    "m/z": float(batch.input_mz[index].cpu()),
                    "charge": float(batch.input_charge[index].cpu()),
                    "estimated_polarizability_A3": float(
                        batch.input_polarizability[index].cpu()
                    ),
                    "true_ccs": true_value,
                    "predicted_ccs": predicted_value,
                    "absolute_error": abs(predicted_value - true_value),
                    "relative_error_pct": (
                        100.0 * abs(predicted_value - true_value) / true_value
                    ),
                }
            )

    return pd.DataFrame(rows).sort_values("csv_row").reset_index(drop=True)


@torch.no_grad()
def predict_external_csv(
    csv_path: Path,
    checkpoint_path: Path,
    output_path: Path,
    batch_size: int,
    seed: int,
) -> None:
    """Load a trained checkpoint and predict CCS for a separate CSV."""
    if not checkpoint_path.exists():
        raise FileNotFoundError(f"Checkpoint not found: {checkpoint_path}")

    # weights_only is available in newer PyTorch versions. The fallback keeps
    # compatibility with older versions.
    try:
        checkpoint = torch.load(
            checkpoint_path,
            map_location="cpu",
            weights_only=True,
        )
    except TypeError:
        checkpoint = torch.load(checkpoint_path, map_location="cpu")

    model_kwargs = checkpoint["model_kwargs"]
    uses_adduct = model_kwargs.get("num_adducts") is not None
    adduct_to_id = checkpoint.get("adduct_to_id") if uses_adduct else None
    if uses_adduct and adduct_to_id is None:
        raise ValueError(
            "Checkpoint model_kwargs contains num_adducts, but the checkpoint "
            "does not contain the required adduct_to_id mapping."
        )

    fallback_adduct = str(
        checkpoint.get("fallback_adduct", DEFAULT_FALLBACK_ADDUCT)
    )
    model_type = "adduct-aware" if uses_adduct else "adduct-agnostic"
    print(f"checkpoint_model={model_type}")
    if uses_adduct:
        print(f"fallback_adduct={fallback_adduct}")
    descriptor_dim = int(model_kwargs["descriptor_dim"])
    graphs = load_prediction_graphs(
        csv_path,
        seed,
        adduct_to_id,
        descriptor_dim=descriptor_dim,
        fallback_adduct=fallback_adduct,
    )

    descriptor_mean = checkpoint["descriptor_mean"].cpu()
    descriptor_std = checkpoint["descriptor_std"].cpu().clamp_min(1e-8)
    for graph in graphs:
        graph.descriptors = (
            graph.descriptors - descriptor_mean
        ) / descriptor_std

    loader = DataLoader(graphs, batch_size=batch_size, shuffle=False)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = ViSNetCCS(**model_kwargs)
    model.load_state_dict(checkpoint["model_state_dict"])
    model = model.to(device).eval()

    y_mean = checkpoint["log_ccs_mean"].to(device)
    y_std = checkpoint["log_ccs_std"].to(device)
    output_rows: list[dict[str, object]] = []

    batches = tqdm(
        loader,
        desc="Predicting CCS",
        unit="batch",
        dynamic_ncols=True,
    )
    for batch in batches:
        batch = batch.to(device)
        prediction_z = model(batch)
        prediction_ccs = torch.exp(prediction_z * y_std + y_mean)

        for index in range(batch.num_graphs):
            output_rows.append(
                {
                    "csv_row": int(batch.row_index[index].cpu()),
                    "Molecule Name": batch.molecule_name[index],
                    "smiles": batch.smiles[index],
                    "adduct": batch.adduct[index],
                    "m/z": float(batch.input_mz[index].cpu()),
                    "charge": float(batch.input_charge[index].cpu()),
                    "estimated_polarizability_A3": float(
                        batch.input_polarizability[index].cpu()
                    ),
                    "predicted_ccs": float(prediction_ccs[index].cpu()),
                }
            )

    output = pd.DataFrame(output_rows).sort_values("csv_row").reset_index(drop=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, index=False)
    print(f"predicted={len(output)} device={device}")
    print(f"predictions={output_path.resolve()}")




def run_sidecar_mode(
    checkpoint_path: Path,
    seed: int,
) -> None:
    predictor = load_predictor(
        checkpoint_path=checkpoint_path,
        seed=seed,
    )

    print(
        f"CCS model ready on {predictor.device}",
        file=sys.stderr,
        flush=True,
    )

    for line in sys.stdin:
        line = line.strip()

        if not line:
            continue

        request_id = None

        try:
            request = json.loads(line)
            request_id = request.get("id")
            action = request.get("action")

            if action == "shutdown":
                print(
                    json.dumps({
                        "id": request_id,
                        "ok": True,
                    }),
                    flush=True,
                )
                break

            if action != "predict":
                raise ValueError(
                    f"Unsupported action: {action!r}"
                )

            charge_value = request.get("charge", 1.0)

            if charge_value is None:
                charge_value = 1.0

            smiles = str(request["smiles"])
            adduct = request.get("adduct")
            resolved_adduct = predictor._resolve_adduct(adduct)
            
            calculated_mz, calculated_charge = calculate_monoisotopic_mz(
                smiles=smiles,
                adduct=resolved_adduct,
            )

            predicted_ccs = predictor.predict(
                smiles=smiles,
                adduct=adduct,
            )

            

            print(
                json.dumps({
                    "id": request_id,
                    "ok": True,
                    "mz": calculated_mz,
                    "predicted_ccs": predicted_ccs,
                }),
                flush=True,
            )

        except Exception as error:
            print(
                json.dumps({
                    "id": request_id,
                    "ok": False,
                    "error": str(error),
                }),
                flush=True,
            )


def main() -> None:
    args = parse_args()
    set_seed(args.seed)

    if args.sidecar_mode:
        run_sidecar_mode(
            checkpoint_path=args.checkpoint,
            seed=args.seed,
        )
        return

    if args.predict_csv is not None:
        predict_external_csv(
            csv_path=args.predict_csv,
            checkpoint_path=args.checkpoint,
            output_path=args.prediction_output,
            batch_size=args.batch_size,
            seed=args.seed,
        )
        return

    fallback_adduct = str(args.fallback_adduct).strip()
    if not fallback_adduct:
        raise ValueError("--fallback-adduct cannot be empty")

    graphs = load_graphs(
        args.csv,
        args.seed,
        use_adduct=args.use_adduct,
        fallback_adduct=fallback_adduct,
    )

    adduct_to_id: dict[str, int] | None = None
    if args.use_adduct:
        observed_adducts = sorted({graph.adduct for graph in graphs})
        if fallback_adduct not in observed_adducts:
            raise ValueError(
                f"Fallback adduct {fallback_adduct!r} is absent from the "
                "training data. Add training examples for it or choose an "
                "observed value with --fallback-adduct."
            )
        adduct_to_id = {
            adduct: index for index, adduct in enumerate(observed_adducts)
        }
        for graph in graphs:
            graph.adduct_id = torch.tensor(
                [adduct_to_id[graph.adduct]], dtype=torch.long
            )
        train_indices, test_indices = split_indices_by_adduct(
            graphs, args.test_fraction, args.seed
        )
        print(f"adduct_to_id={adduct_to_id}")
    else:
        train_indices, test_indices = split_indices(
            len(graphs), args.test_fraction, args.seed
        )

    descriptor_mean, descriptor_std, y_mean, y_std = training_statistics(
        graphs, train_indices
    )
    for graph in graphs:
        graph.descriptors = (
            graph.descriptors - descriptor_mean
        ) / descriptor_std

    train_graphs = [graphs[index] for index in train_indices]
    test_graphs = [graphs[index] for index in test_indices]

    train_loader = DataLoader(
        train_graphs,
        batch_size=args.batch_size,
        shuffle=True,
    )
    test_loader = DataLoader(
        test_graphs,
        batch_size=args.batch_size,
        shuffle=False,
    )

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model_kwargs: dict[str, int] = {
        "descriptor_dim": len(DESCRIPTOR_NAMES),
        "hidden": 64,
    }
    if adduct_to_id is not None:
        model_kwargs["num_adducts"] = len(adduct_to_id)

    model = ViSNetCCS(**model_kwargs).to(device)

    print(
        f"device={device} valid={len(graphs)} "
        f"train={len(train_graphs)} test={len(test_graphs)}"
    )
    train_model(
        model=model,
        loader=train_loader,
        device=device,
        y_mean=y_mean,
        y_std=y_std,
        epochs=args.epochs,
        learning_rate=args.lr,
    )

    predictions = evaluate(model, test_loader, device, y_mean, y_std)
    args.predictions.parent.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(args.predictions, index=False)

    mae = predictions["absolute_error"].mean()
    median_relative_error = predictions["relative_error_pct"].median()
    print(f"test_MAE={mae:.3f} A^2")
    print(f"test_median_relative_error={median_relative_error:.2f}%")
    print(f"predictions={args.predictions.resolve()}")

    args.checkpoint.parent.mkdir(parents=True, exist_ok=True)
    checkpoint = {
        "model_version": (
            "visnet_ccs_adduct_aware_polarizability_v1"
            if args.use_adduct
            else "visnet_ccs_adduct_agnostic_polarizability_v1"
        ),
        "model_state_dict": model.state_dict(),
        "descriptor_names": DESCRIPTOR_NAMES,
        "descriptor_mean": descriptor_mean,
        "descriptor_std": descriptor_std,
        "log_ccs_mean": y_mean,
        "log_ccs_std": y_std,
        "model_kwargs": model_kwargs,
        "fallback_adduct": fallback_adduct,
    }
    if adduct_to_id is not None:
        checkpoint["adduct_to_id"] = adduct_to_id

    torch.save(checkpoint, args.checkpoint)
    print(f"checkpoint={args.checkpoint.resolve()}")


if __name__ == "__main__":
    main()
