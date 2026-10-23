"""Python client for an already installed Met-ID application."""
import json
import os
import platform
import shutil
import socket
import subprocess
import threading
import time
import uuid
from collections import deque
from pathlib import Path

import numpy as np
import pandas as pd
import requests


def _as_1d_list(value, name):
    values = np.asarray(value)
    if values.ndim > 1:
        raise ValueError(f"{name} must be a scalar or a one-dimensional sequence")
    return values.reshape(-1).tolist()


class MetIDServer:
    def __init__(self, project_dir=None, port=None, binary_name=None,
                 kill_conflicts=False, *, server_path=None, sidecar_path=None,
                 startup_timeout=180, auto_shutdown_after=3600):
        """Launch installed Met-ID; never invoke Cargo or rustc.

        Override discovery with server_path / sidecar_path or the environment
        variables METID_SERVER_PATH / METID_CCS_SIDECAR. project_dir is an
        optional installation directory or existing src-tauri directory.
        Occupied ports are rejected; unrelated processes are never killed.
        """
        self.instance_id = uuid.uuid4().hex[:6]
        self.process = None
        self._timer = None
        self._output = deque(maxlen=50)
        self.startup_timeout = startup_timeout
        self.auto_shutdown_after = auto_shutdown_after
        self.kill_conflicts = kill_conflicts  # retained for call compatibility
        self.project_dir = str(Path(project_dir).expanduser().resolve()) if project_dir else None
        self.port = port if port is not None else self._get_free_port()
        if not isinstance(self.port, int) or not 1 <= self.port <= 65535:
            raise ValueError("port must be an integer between 1 and 65535")
        self.base_url = f"http://127.0.0.1:{self.port}"
        self.binary_path = str(self._find_server(server_path, binary_name))
        self.binary_name = Path(self.binary_path).name
        self.ccs_sidecar_path = str(self._find_sidecar(sidecar_path))
        self.start()

    @staticmethod
    def _checked_file(value, label):
        path = Path(value).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"{label} not found: {path}")
        if os.name != "nt" and not os.access(path, os.X_OK):
            raise PermissionError(f"{label} is not executable: {path}")
        return path

    def _find_server(self, explicit, binary_name):
        override = explicit or os.environ.get("METID_SERVER_PATH")
        if override:
            return self._checked_file(override, "Met-ID executable")
        name = binary_name or "Met-ID"
        if os.name == "nt" and not name.lower().endswith(".exe"):
            name += ".exe"
        candidates = []
        if self.project_dir:
            root = Path(self.project_dir)
            candidates.extend([root / name, root / "target/release" / name,
                               root / "Contents/MacOS" / name])
        found = shutil.which(name)
        if found:
            candidates.append(Path(found))
        system = platform.system()
        if system == "Windows":
            for key in ("LOCALAPPDATA", "ProgramFiles", "ProgramFiles(x86)"):
                root = os.environ.get(key)
                if root:
                    for folder in ("Met-ID", "met-id", "Programs/Met-ID", "Programs/met-id"):
                        candidates.append(Path(root) / folder / name)
        elif system == "Darwin":
            for root in (Path("/Applications"), Path.home() / "Applications"):
                candidates.append(root / "Met-ID.app/Contents/MacOS" / name)
        else:
            candidates.extend([Path("/usr/bin") / name, Path("/usr/local/bin") / name,
                               Path.home() / ".local/bin" / name,
                               Path("/usr/bin/met-id")])
        for candidate in candidates:
            if candidate.is_file():
                return self._checked_file(candidate, "Met-ID executable")
        raise FileNotFoundError(
            "Could not locate installed Met-ID. Pass server_path=r'C:\\path\\Met-ID.exe' "
            "or set METID_SERVER_PATH. Custom installations and AppImages require an explicit path."
        )

    def _find_sidecar(self, explicit):
        override = (explicit or os.environ.get("METID_CCS_SIDECAR")
                    or os.environ.get("METID_CCS_SIDECAR_PATH"))
        if override:
            return self._checked_file(override, "CCS sidecar")
        binary_dir = Path(self.binary_path).parent
        roots = [binary_dir, binary_dir / "pyinstaller/dist",
                 binary_dir / "resources/pyinstaller/dist"]
        if self.project_dir:
            roots.append(Path(self.project_dir) / "pyinstaller/dist")
        if binary_dir.name == "MacOS" and binary_dir.parent.name == "Contents":
            roots.extend([binary_dir.parent / "Resources",
                          binary_dir.parent / "Resources/pyinstaller/dist"])
        if platform.system() == "Linux":
            for folder in ("met-id", "Met-ID", "com.farmbio.metid"):
                root = Path("/usr/lib") / folder
                roots.extend([root, root / "pyinstaller/dist"])
        suffix = ".exe" if os.name == "nt" else ""
        machine = platform.machine().lower()
        arch = "aarch64" if machine in ("arm64", "aarch64") else "x86_64"
        triple = {"Windows": f"{arch}-pc-windows-msvc",
                  "Darwin": f"{arch}-apple-darwin",
                  "Linux": f"{arch}-unknown-linux-gnu"}.get(platform.system())
        names = [f"prediction_model{suffix}"]
        if triple:
            names.append(f"prediction_model-{triple}{suffix}")
        for root in roots:
            for name in names:
                candidate = root / name
                if candidate.is_file():
                    return self._checked_file(candidate, "CCS sidecar")
        raise FileNotFoundError(
            "Could not locate prediction_model beside installed Met-ID. "
            "Pass sidecar_path or set METID_CCS_SIDECAR to its full executable path."
        )

    @staticmethod
    def _get_free_port():
        with socket.socket() as sock:
            sock.bind(("127.0.0.1", 0))
            return sock.getsockname()[1]

    def _is_port_in_use(self):
        with socket.socket() as sock:
            sock.settimeout(0.5)
            return sock.connect_ex(("127.0.0.1", self.port)) == 0

    def _print_stdout(self, stream):
        try:
            for line in stream:
                self._output.append(line.rstrip())
                print(f"[{self.instance_id} RUST] {line}", end="")
        finally:
            stream.close()

    def start(self, auto_shutdown_after=None):
        if self.process is not None and self.process.poll() is None:
            return self
        if self._is_port_in_use():
            raise RuntimeError(f"Port {self.port} is already in use. Choose another port.")
        if self._timer:
            self._timer.cancel()
        self._output.clear()
        env = os.environ.copy()
        env["METID_CCS_SIDECAR"] = self.ccs_sidecar_path
        self.process = subprocess.Popen(
            [self.binary_path, "--server", "--port", str(self.port)],
            cwd=str(Path(self.binary_path).parent), env=env,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace",
        )
        threading.Thread(target=self._print_stdout,
                         args=(self.process.stdout,), daemon=True).start()
        try:
            self._wait_for_server(self.startup_timeout)
        except Exception:
            self.stop()
            raise
        delay = self.auto_shutdown_after if auto_shutdown_after is None else auto_shutdown_after
        if delay is not None and delay > 0:
            self._timer = threading.Timer(delay, self.stop)
            self._timer.daemon = True
            self._timer.start()
        return self

    def _wait_for_server(self, timeout):
        deadline = time.monotonic() + timeout
        with requests.Session() as session:
            session.trust_env = False
            while time.monotonic() < deadline:
                code = self.process.poll()
                if code is not None:
                    raise RuntimeError(
                        f"Met-ID exited during startup (code {code}).\n"
                        + "\n".join(self._output)
                        + "\nCheck the Met-ID application log too; Windows GUI binaries may not emit console output."
                    )
                try:
                    response = session.get(f"{self.base_url}/warmup", timeout=(0.5, 1))
                    if response.status_code == 200 and response.text.strip() == "OK":
                        return
                except requests.RequestException:
                    pass
                time.sleep(0.25)
        raise TimeoutError(f"Met-ID was not ready after {timeout} seconds.\n" + "\n".join(self._output))

    def stop(self):
        if self._timer:
            self._timer.cancel()
            self._timer = None
        process = self.process
        if process is not None:
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait(timeout=5)
            self.process = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.stop()

    def _raise_with_rust_error(self, res):
        try:
            payload = res.json()
            message = (payload.get("error") or payload.get("message") or json.dumps(payload)
                       if isinstance(payload, dict) else json.dumps(payload))
        except ValueError:
            message = res.text
        raise MetIDApiError(res.status_code, message)

    def return_name(self, name:str, retries=3):
        for attempt in range(retries):
            try:
                res = requests.get(f"{self.base_url}/name", json=name)
                if not res.ok:
                    self._raise_with_rust_error(res)
                return res.json()
            except Exception as e:
                print(f"[{self.instance_id}] ❌ return_name failed: {e}")
                if attempt < retries - 1:
                    time.sleep(0.5)
        raise RuntimeError(f"[{self.instance_id}] Failed to return name after retries")

    def process_items(self, items, retries=3):
        for attempt in range(retries):
            try:
                res = requests.post(f"{self.base_url}/process", json=items)
                if not res.ok:
                    self._raise_with_rust_error(res)
                return res.json()
            except Exception as e:
                print(f"[{self.instance_id}] ❌ process_items failed: {e}")
                if attempt < retries - 1:
                    time.sleep(0.5)
        raise RuntimeError(f"[{self.instance_id}] Failed to process items after retries")
    
    def identify(self, 
                 metabolome: str = "HMDB(All)", 
                 matrix: str = "FMP-10", 
                 typ: list = "All", 
                 adducts: list = "All",
                 mass_error: str = 0,
                 mz = 421.1911,
                 mzwindow = "5",
                 ccs = -1,
                 retries=1
                 ):

        masses = _as_1d_list(mz, "mz")

        if ccs is None:
            ccss = [-1] * len(masses)
        else:
            ccss = _as_1d_list(ccs, "ccs")

            # Broadcast the no-CCS sentinel across all masses.
            if len(ccss) == 1 and float(ccss[0]) == -1:
                ccss = [-1] * len(masses)

        if len(ccss) != len(masses):
            raise ValueError("Provide one CCS value per m/z value")
        
        items = {
            "met": metabolome,
            "mat": matrix,
            "typ": typ, 
            "adducts": adducts,
            "mass_error": mass_error,
            "masses": masses,
            "ccss": ccss,
            "mzwindow": mzwindow,
        }
        for attempt in range(retries):
            try:
                res = requests.post(f"{self.base_url}/identify", json=items)
                if not res.ok:
                    self._raise_with_rust_error(res)    
                df = (
                    pd.DataFrame.from_dict(res.json(), orient="index")
                    .transpose()
                    .set_index("names")
                )

                return IdentificationDataFrame(df)
            
            except Exception as e:
                print(f"[{self.instance_id}] ❌ identify failed: {e}")
                if attempt < retries - 1:
                    time.sleep(0.5)
        self.stop()
        raise RuntimeError(f"[{self.instance_id}] Failed to process items after retries")
    
    def list_tables(self, retries=3):
        for attempt in range(retries):
            try: 
                res = requests.get(f"{self.base_url}/get_database", json="tables")
                res.raise_for_status()
                return res.json()["Tables"]
            except Exception as e:
                print(f"[{self.instance_id}] ❌ get_database failed: {e}")
                if attempt < retries - 1:
                    time.sleep(0.1)
        self.stop()
        raise RuntimeError(f"[{self.instance_id}] Failed to process items after retries")
    
    def list_metabolites(self, retries=3):
        for attempt in range(retries):
            try:
                res = requests.get(f"{self.base_url}/get_database", json="metabolites")
                res.raise_for_status()
                #return res.json()["Table"]["rows"]
                return pd.DataFrame(res.json()["Table"]["rows"], columns= res.json()["Table"]["col_names"])
            except Exception as e:
                print(f"[{self.instance_id}] get_database failed: {e}")
                if attempt < retries - 1:
                    time.sleep(0.1)
        self.stop()
        raise RuntimeError(f"[{self.instance_id}] Failed to process items after retries")
    
    def list_matrices(self, retries=3):
        for attempt in range(retries):
            try:
                res = requests.get(f"{self.base_url}/get_database", json="matrices")
                res.raise_for_status()
                #return res.json()["Table"]["rows"]
                return pd.DataFrame(res.json()["Table"]["rows"], columns= res.json()["Table"]["col_names"])
            except Exception as e:
                print(f"[{self.instance_id}] get_database failed: {e}")
                if attempt < retries - 1:
                    time.sleep(0.1)
        self.stop()
        raise RuntimeError(f"[{self.instance_id}] Failed to process items after retries")
    
    def list_functional_groups(self, retries=3):
        for attempt in range(retries):
            try:
                res = requests.get(f"{self.base_url}/get_database", json="functional_groups")
                res.raise_for_status()
                #return res.json()["Table"]["rows"]
                return pd.DataFrame(res.json()["Table"]["rows"], columns= res.json()["Table"]["col_names"])
            except Exception as e:
                print(f"[{self.instance_id}] get_database failed: {e}")
                if attempt < retries - 1:
                    time.sleep(0.1)
        self.stop()
        raise RuntimeError(f"[{self.instance_id}] Failed to process items after retries")

    def list_adducts(self, retries=3):
        for attempt in range(retries):
            try:
                res = requests.get(f"{self.base_url}/get_database", json="adducts")
                res.raise_for_status()
                #return res.json()["Table"]["rows"]
                return pd.DataFrame(res.json()["Table"]["rows"], columns= res.json()["Table"]["col_names"])
            except Exception as e:
                print(f"[{self.instance_id}] get_database failed: {e}")
                if attempt < retries - 1:
                    time.sleep(0.1)
        self.stop()
        raise RuntimeError(f"[{self.instance_id}] Failed to process items after retries")
    
    def search_by_name(self, name, ccs = False, retries=3):
        items = {
            "name": name,
            "ccs": ccs
        }
        for attempt in range(retries):
            try: 
                res = requests.get(f"{self.base_url}/search_by_name", json=items)
                res.raise_for_status()
                print(res)
                return pd.DataFrame(res.json()["Table"]["rows"], columns=res.json()["Table"]["col_names"]).set_index("name")
            except Exception as e:
                print(f"[{self.instance_id}] search_by_name failed: {e}")
                if attempt < retries - 1:
                    time.sleep(0.1)
        self.stop()
        raise RuntimeError(f"[{self.instance_id}] Failed to process items after retries")

    def read_imzml(self, path: str, n: int, tolerance: float, retries=3):
        payload = {
            "path": os.fspath(path),
            "n": n,
            "tolerance": tolerance,
        }

        for attempt in range(retries):
            try:
                res = requests.post(
                    f"{self.base_url}/read_imzml",
                    json=payload,
                )
                if not res.ok:
                    self._raise_with_rust_error(res)

                data = res.json()
                return pd.DataFrame({
                    "mz": data["mz"],
                    "intensity": data["intensity"],
                })

            except Exception as e:
                print(f"[{self.instance_id}] ❌ read_imzml failed: {e}")
                if attempt < retries - 1:
                    time.sleep(0.5)

        raise RuntimeError(
            f"[{self.instance_id}] Failed to read imzML after retries"
        )

    def predict_ccs(self, smiles, adduct=None, retries=3):
        def json_value(value):
            return value.tolist() if hasattr(value, "tolist") else value

        smiles_value = json_value(smiles)
        smiles_items = (
            smiles_value
            if isinstance(smiles_value, (list, tuple))
            else [smiles_value]
        )

        if not smiles_items:
            raise ValueError("Provide at least one SMILES")

        for value in smiles_items:
            if not isinstance(value, str) or not value.strip():
                raise ValueError("SMILES values cannot be empty")

        payload = {"smiles": smiles_value}

        if adduct is not None:
            adduct_value = json_value(adduct)
            adduct_items = (
                adduct_value
                if isinstance(adduct_value, (list, tuple))
                else [adduct_value]
            )

            allowed_adducts = {
                "[M+H]+",
                "[M]+",
                "[M+Na]+",
                "[M+K]+",
            }
            fallback_adduct = "[M+H]+"

            adduct_value = json_value(adduct)
            adduct_items = (
                list(adduct_value)
                if isinstance(adduct_value, (list, tuple))
                else [adduct_value]
            )

            normalized_adducts = []

            for value in adduct_items:
                requested = "" if value is None else str(value).strip()

                if not requested:
                    normalized_adducts.append(fallback_adduct)
                elif requested not in allowed_adducts:
                    print(
                        f"Adduct '{requested}' is not allowed; "
                        f"processing with fallback adduct '{fallback_adduct}'."
                        f"Allowed adducts are: {allowed_adducts}"
                    )
                    normalized_adducts.append(fallback_adduct)
                else:
                    normalized_adducts.append(requested)

            payload["adducts"] = (
                normalized_adducts
                if isinstance(adduct_value, (list, tuple))
                else normalized_adducts[0]
            )

        for attempt in range(retries):
            try:
                res = requests.post(
                    f"{self.base_url}/predict_ccs",
                    json=payload,
                )
                if not res.ok:
                    self._raise_with_rust_error(res)

                return pd.DataFrame(res.json())

            except MetIDApiError as e:
                if 400 <= e.status_code < 500:
                    # Invalid input: report it clearly and don't retry.
                    raise ValueError(e.message) from e

                print(f"[{self.instance_id}] ❌ CCS prediction failed: {e}")
                if attempt < retries - 1:
                    time.sleep(0.5)

            except Exception as e:
                print(f"[{self.instance_id}] ❌ CCS prediction failed: {e}")
                if attempt < retries - 1:
                    time.sleep(0.5)


        raise RuntimeError("Failed to predict CCS after retries; see preceding errors")


class IdentificationDataFrame(pd.DataFrame):
    @property
    def _constructor(self):
        return IdentificationDataFrame

    def as_dict(self, mass_column="oMass"):
        frame = (
            self.reset_index()
            if self.index.name == "names"
            else self.copy()
        )

        if mass_column not in frame.columns:
            raise ValueError(f"DataFrame has no {mass_column!r} column")

        result = {}

        for mass, group in frame.groupby(
            mass_column, sort=False, dropna=False
        ):
            try:
                mass_key = float(mass)
            except (TypeError, ValueError):
                mass_key = mass

            identifications = {}

            for number, (_, row) in enumerate(
                group.iterrows(), start=1
            ):
                name = str(row.get("names", "")).strip()
                if not name or name.lower() == "nan":
                    name = f"identification_{number}"

                key = name
                suffix = 2
                while key in identifications:
                    key = f"{name} ({suffix})"
                    suffix += 1

                identifications[key] = row.to_dict()

            result[mass_key] = identifications

        return result

class MetIDApiError(RuntimeError):
    def __init__(self, status_code, message):
        self.status_code = status_code
        self.message = message
        super().__init__(message)