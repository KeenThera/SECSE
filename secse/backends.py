import json
import os
import shutil
from pathlib import Path
from tempfile import NamedTemporaryFile


NAMES = ("vina", "unidock", "autodock-gpu", "glide")
ENGINE_VARIABLES = {
    "vina": "VINA",
    "unidock": "UNIDOCK",
    "autodock-gpu": "AUTODOCK_GPU",
    "glide": "SCHRODINGER",
}
TARGET_SUFFIXES = {"vina": ".pdbqt", "unidock": ".pdbqt", "autodock-gpu": ".fld", "glide": ".zip"}


def registry_path():
    config = Path(os.environ.get("XDG_CONFIG_HOME") or Path.home() / ".config")
    return config / "secse" / "backends.json"


def load_registry():
    path = registry_path()
    try:
        with path.open(encoding="utf-8") as stream:
            registry = json.load(stream)
    except FileNotFoundError:
        return {"default": None, "backends": {}}
    if (not isinstance(registry, dict) or registry.get("default") not in (*NAMES, None)
            or not isinstance(registry.get("backends"), dict)
            or any(name not in NAMES or not isinstance(settings, dict)
                   or not all(key in {"path", "conda_env", "library_dir"} and isinstance(value, str)
                              for key, value in settings.items())
                   for name, settings in registry["backends"].items())
            or (registry["default"] and registry["default"] not in registry["backends"])):
        raise ValueError(f"Invalid backend registry: {path}")
    return registry


def save_registry(registry):
    path = registry_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    with NamedTemporaryFile(mode="w", encoding="utf-8", prefix=".backends-", dir=path.parent,
                            delete=False) as stream:
        temporary = Path(stream.name)
        try:
            json.dump(registry, stream, indent=2)
            stream.write("\n")
            stream.close()
            temporary.replace(path)
        finally:
            temporary.unlink(missing_ok=True)


def conda_command():
    return os.environ.get("CONDA_EXE") or shutil.which("conda") or "conda"


def engine_environment(name, settings):
    env = os.environ.copy()
    if "path" in settings:
        env[ENGINE_VARIABLES[name]] = settings["path"]
    if name == "unidock":
        for key in ("SECSE_UNIDOCK_CONDA_ENV", "SECSE_UNIDOCK_LIBRARY_DIR"):
            env.pop(key, None)
        if "conda_env" in settings:
            env["SECSE_UNIDOCK_CONDA_ENV"] = settings["conda_env"]
            env["SECSE_CONDA_EXE"] = conda_command()
        if "library_dir" in settings:
            env["SECSE_UNIDOCK_LIBRARY_DIR"] = settings["library_dir"]
    return env


def engine_commands(name, settings):
    env = engine_environment(name, settings)
    location = env.get(ENGINE_VARIABLES[name])
    if name == "unidock" and "conda_env" in settings:
        return [[conda_command(), "run", "--no-capture-output", "-n", settings["conda_env"], "unidock", "--version"]]
    if name in ("vina", "unidock"):
        return [[location or name, "--version"]]
    if not location:
        raise ValueError(f"Set up {name} with 'python -m secse backend add {name} --install-dir PATH'.")
    if name == "autodock-gpu":
        return [[str(Path(location) / "bin/autodock_gpu_128wi"), "--help"]]
    return [[str(Path(location) / "glide"), "-version"], [str(Path(location) / "ligprep"), "-h"]]


def probe_environment(name, settings):
    env = engine_environment(name, settings)
    if name == "unidock" and "library_dir" in settings:
        env["LD_LIBRARY_PATH"] = os.pathsep.join(filter(None, [settings["library_dir"], env.get("LD_LIBRARY_PATH")]))
    return env
