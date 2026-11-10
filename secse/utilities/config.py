import configparser
from pathlib import Path


SOURCE_DIR = Path(__file__).resolve().parents[1]


def read_config(path):
    config = configparser.ConfigParser()
    with open(path, encoding="utf-8") as stream:
        config.read_file(stream)
    return config


def docking_backend(name):
    # Keep the aliases accepted by the original CLI, in the same order.
    for backend in ("vina", "glide", "autodock-gpu", "unidock"):
        if backend in name.lower():
            return backend
    raise ValueError(f"Unsupported docking program: {name!r}")
