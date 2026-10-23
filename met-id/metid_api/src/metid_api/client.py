import subprocess
import requests
import os
import platform
import time
import threading
import toml
import uuid
import socket
import psutil
import pandas as pd
import json

class MetIDServer:
    def __init__(self, project_dir=None, port=None, binary_name=None, kill_conflicts=True):
        if project_dir is None:
            # Default to the directory where this Python file is installed
            project_dir = os.path.dirname(os.path.abspath(__file__))
            project_dir = os.path.join(project_dir, "../src-tauri")
        self.project_dir = os.path.abspath(project_dir)
        self.port = port or self._get_free_port()
        self.binary_name = binary_name or self._detect_binary_name()
        self.binary_name += ".exe" if platform.system() == "Windows" else ""
        self.instance_id = uuid.uuid4().hex[:6]
        self.kill_conflicts = kill_conflicts
        self.base_url = f"http://127.0.0.1:{self.port}"
        self.process = None

        triple = self.get_rust_host_triple()
        suffix = ".exe" if platform.system() == "Windows" else ""

        #self.ccs_sidecar_path = os.path.abspath(
        #    os.path.join(
        #        self.project_dir,
        #        "bin",
        #        f"prediction_model-{triple}{suffix}",
        #    ),
        #)
        self.ccs_sidecar_path = os.path.abspath(
            os.path.join(
                self.project_dir,
                "pyinstaller",
                "dist",
                f"prediction_model-{triple}{suffix}",
            ),
        )

        self.ccs_sidecar_path = os.path.abspath(self.ccs_sidecar_path)

        self.binary_path = os.path.join(self.project_dir, "target", "release", self.binary_name)
        self.start()

    @staticmethod
    def get_rust_host_triple():
        result = subprocess.run(
            ["rustc", "-vV"],
            check=True,
            capture_output=True,
            text=True,
        )

        for line in result.stdout.splitlines():
            if line.startswith("host:"):
                return line.split(":", 1)[1].strip()

        raise RuntimeError("Could not find Rust host triple in `rustc -vV` output")
    
    def _detect_binary_name(self):
        cargo_toml = os.path.join(self.project_dir, "Cargo.toml")
        if not os.path.exists(cargo_toml):
            raise FileNotFoundError("Cargo.toml not found in project directory.")
        with open(cargo_toml, "r") as f:
            parsed = toml.load(f)
        return parsed["package"]["name"]

    def _get_free_port(self):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.bind(("127.0.0.1", 0))
            return s.getsockname()[1]

    def _is_port_in_use(self):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as sock:
            return sock.connect_ex(("127.0.0.1", self.port)) == 0

    def _kill_port(self):
        for conn in psutil.net_connections(kind='inet'):
            if conn.laddr.port == self.port and conn.status == psutil.CONN_LISTEN:
                print(f"💀 Killing process on port {self.port} (PID {conn.pid})")
                if conn.pid:
                    try:
                        psutil.Process(conn.pid).terminate()
                    except Exception as e:
                        print(f"⚠️ Could not kill process {conn.pid}: {e}")
                return

    def _kill_process_using_file(self, filepath):
        for proc in psutil.process_iter(['pid', 'open_files']):
            try:
                if any(f.path == filepath for f in (proc.info['open_files'] or [])):
                    print(f"💀 Killing process {proc.pid} locking {filepath}")
                    psutil.Process(proc.pid).terminate()
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass

    def _attempt_build(self, max_retries=2):
        for attempt in range(max_retries):
            print(f"🔨 Building release binary... (attempt {attempt + 1})")
            result = subprocess.run(["cargo", "build", "--release"], cwd=self.project_dir)
            if result.returncode == 0:
                return True
            print("❌ Build failed.")
            if attempt == 0:
                print("🔍 Trying to kill process that may be locking the binary...")
                self._kill_process_using_file(self.binary_path)
                time.sleep(1)
        return False

    def _print_stdout(self, stream):
        for line in stream:
            print(f"[{self.instance_id} RUST]", line, end="")

    def start(self, auto_shutdown_after=3600):
        if not os.path.isfile(self.binary_path):
            if not self._attempt_build():
                raise RuntimeError("❌ Cargo build failed after retry.")

        if not os.path.exists(self.binary_path):
            raise FileNotFoundError(f"❌ Binary not found after build: {self.binary_path}")

        if self._is_port_in_use():
            if self.kill_conflicts:
                self._kill_port()
                time.sleep(0.1)
            else:
                raise RuntimeError(f"❌ Port {self.port} is already in use!")

        print(f"🚀 Starting server [{self.instance_id}] on port {self.port}")
        env = os.environ.copy()
        env["METID_CCS_SIDECAR"] = self.ccs_sidecar_path

        print("CCS sidecar:", env["METID_CCS_SIDECAR"])
        print("Sidecar exists:", os.path.isfile(env["METID_CCS_SIDECAR"]))

        self.process = subprocess.Popen(
            [self.binary_path, "--server", "--port", str(self.port)],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            env=env,
        )

        threading.Thread(target=self._print_stdout, args=(self.process.stdout,), daemon=True).start()
        try:
            self._wait_for_server(timeout=60)
        except Exception:
            if self.process is not None:
                if self.process.poll() is None:
                    self.process.terminate()
                    try:
                        self.process.wait(timeout=5)
                    except subprocess.TimeoutExpired:
                        self.process.kill()
                        self.process.wait()
                self.process = None
            raise
        self._start_auto_shutdown_timer(timeout_seconds=auto_shutdown_after)

    def _wait_for_server(self, timeout=60):
        print(f"[{self.instance_id}] ⏳ Waiting for server...")
        deadline = time.monotonic() + timeout
        url = f"{self.base_url}/warmup"

        while time.monotonic() < deadline:
            if self.process is None:
                raise RuntimeError("Server process was not started.")

            exit_code = self.process.poll()
            if exit_code is not None:
                raise RuntimeError(
                    f"Met-ID server exited before becoming ready "
                    f"(exit code {exit_code})."
                )

            try:
                res = requests.get(url, timeout=(0.5, 1))
                if res.status_code == 200 and res.text.strip() == "OK":
                    print(f"[{self.instance_id}] ✅ Server is ready!")
                    return
            except requests.exceptions.RequestException:
                pass

            time.sleep(0.25)

        raise TimeoutError(
            f"[{self.instance_id}] Server did not become ready within {timeout} seconds."
        )
    
    def _start_auto_shutdown_timer(self, timeout_seconds=1800):  # default = 30 min
        def shutdown_later():
            print(f"[{self.instance_id}] ⏳ Auto-shutdown timer started ({timeout_seconds} seconds)")
            time.sleep(timeout_seconds)
            print(f"[{self.instance_id}] 🔔 Auto-shutdown timeout reached, stopping server.")
            self.stop()
        threading.Thread(target=shutdown_later, daemon=True).start()

    def _raise_with_rust_error(self, res: requests.Response):
        # Try JSON first
        ct = res.headers.get("content-type", "")
        if "application/json" in ct.lower():
            try:
                payload = res.json()
                # common patterns: {"error":"..."} or {"message":"..."}
                msg = payload.get("error") or payload.get("message") or json.dumps(payload)
            except Exception:
                msg = res.text
        else:
            msg = res.text

        raise MetIDApiError(res.status_code, msg)

    def stop(self):
        print(f"[{self.instance_id}] 🛑 Sending shutdown signal...")
        try:
            requests.post(f"{self.base_url}/shutdown")
        except requests.exceptions.RequestException as e:
            print(f"[{self.instance_id}] ⚠️ Shutdown failed: {e}")
        if self.process:
            self.process.terminate()
            self.process.wait()
            print(f"[{self.instance_id}] ✅ Server process terminated.")
            self.process = None

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

    def _as_1d_list(value, name):
        values = np.asarray(value)

        if values.ndim > 1:
            raise ValueError(f"{name} must be a scalar or a one-dimensional sequence")

        return values.reshape(-1).tolist()
    
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


        raise RuntimeError(f"Failed to predict CCS after retries: {e}")


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