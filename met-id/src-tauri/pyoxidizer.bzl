def make_python_distribution():
    return default_python_distribution()

def make_exe():
    print("Building exe...")
    dist = make_python_distribution()
    python_config = dist.make_python_interpreter_config()
    exe = dist.to_python_executable(
        config=python_config,
        name="met-id_python",
    )
    exe.add_in_memory_python_resources(dist.pip_install(["rdkit"]))
    print("Exe built successfully!")
    return exe

def make_embedded_resources():
    exe = make_exe()
    resources = FileManifest()
    resources.add_python_resource(".", exe)
    return resources

register_target("dist", make_python_distribution)
register_target("met-id_python", make_exe, default=True)
register_target("resources", make_embedded_resources)
resolve_target("met-id_python")