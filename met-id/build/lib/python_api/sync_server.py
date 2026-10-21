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

        self.binary_path = os.path.join(self.project_dir, "target", "release", self.binary_name)

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

    def start(self, auto_shutdown_after=1800):
        if not os.path.exists(self.binary_path) or self.kill_conflicts:
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
        self.process = subprocess.Popen(
            [self.binary_path, "--server", "--port", str(self.port)],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )

        threading.Thread(target=self._print_stdout, args=(self.process.stdout,), daemon=True).start()
        self._wait_for_server()
        self._start_auto_shutdown_timer(timeout_seconds=auto_shutdown_after)

    def _wait_for_server(self, timeout=35):
        print(f"[{self.instance_id}] ⏳ Waiting for server...")
        start = time.time()
        url = f"{self.base_url}/warmup"

        while time.time() - start < timeout:
            try:
                res = requests.get(url)
                if res.status_code == 200:
                    print(f"[{self.instance_id}] ✅ Server is ready!")
                    return
            except requests.exceptions.RequestException:
                pass
            time.sleep(0.1)

        raise TimeoutError(f"[{self.instance_id}] ❌ Server did not start in time.")
    
    def _start_auto_shutdown_timer(self, timeout_seconds=1800):  # default = 30 min
        def shutdown_later():
            print(f"[{self.instance_id}] ⏳ Auto-shutdown timer started ({timeout_seconds} seconds)")
            time.sleep(timeout_seconds)
            print(f"[{self.instance_id}] 🔔 Auto-shutdown timeout reached, stopping server.")
            self.stop()
        threading.Thread(target=shutdown_later, daemon=True).start()

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
                res.raise_for_status()
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
                res.raise_for_status()
                return res.json()
            except Exception as e:
                print(f"[{self.instance_id}] ❌ process_items failed: {e}")
                if attempt < retries - 1:
                    time.sleep(0.5)
        raise RuntimeError(f"[{self.instance_id}] Failed to process items after retries")
    
    def identify(self, 
                 metabolome: str = "HMDB(All)", 
                 matrix: str = "FMP-10", 
                 typ: list = ["Endogenous", "Exogenous", "Unspecified"], 
                 adducts: list = ["Phenolic Hydroxyls"],
                 mass_error: str = "0",
                 masses: list = ["421.1911", "688.285"],
                 mzwindow: str = "50",
                 retries=3):
        items = {
            "met": metabolome,
            "mat": matrix,
            "typ": typ, 
            "adducts": adducts,
            "mass_error": mass_error,
            "masses": masses,
            "mzwindow": mzwindow,
        }
        for attempt in range(retries):
            try:
                res = requests.post(f"{self.base_url}/identify", json=items)
                res.raise_for_status()
                return pd.DataFrame.from_dict(res.json(), orient="index").transpose().set_index("names")
            except Exception as e:
                print(f"[{self.instance_id}] ❌ process_items failed: {e}")
                if attempt < retries - 1:
                    time.sleep(0.5)
        raise RuntimeError(f"[{self.instance_id}] Failed to process items after retries")
