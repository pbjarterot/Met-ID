from setuptools import setup
import os
from setuptools.command.install import install
import shutil

class CustomInstallCommand(install):
    def run(self):
        # Ensure the binary is executable and copied to the install dir
        print("\n\n\n\n", self.install_lib,  "\n \n ")
        print(os.path.dirname(os.path.realpath(__file__)))
        target_dir = os.path.join(self.install_lib, 'python_api', 'bin')
        os.makedirs(target_dir, exist_ok=True)
        shutil.copy('./src-tauri/target/release/Met-ID.exe', target_dir)
        os.chmod(os.path.join(target_dir, 'Met-ID.exe'), 0o755)
        install.run(self)

setup(
    name='met_id',
    version='0.1',
    packages=['python_api'],
    include_package_data=True,
    install_requires=[],
    cmdclass={'install': CustomInstallCommand},
    package_data={
        'python_api': ['bin/Met-ID.exe']},
    entry_points={
        'console_scripts': ['launch-tauri=python_api.sync_server:MetIDServer.start']
    }
)