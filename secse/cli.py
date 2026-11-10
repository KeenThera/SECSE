import argparse
import configparser
import importlib
from itertools import islice
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import mkdtemp

from .backends import (
    NAMES, TARGET_SUFFIXES, engine_commands, engine_environment, load_registry,
    probe_environment, registry_path, save_registry,
)
from .utilities.config import docking_backend, read_config


DEMO_DIR = Path(__file__).resolve().parents[1] / "demo"
PDBQT_BACKENDS = ("vina", "unidock")


def validate_target(backend, target):
    path = Path(target).resolve(strict=True)
    if not path.is_file() or not path.stat().st_size or path.suffix.lower() != TARGET_SUFFIXES[backend]:
        raise ValueError(f"{backend} requires a {TARGET_SUFFIXES[backend]} target file: {path}")
    return path


def probe(command, env=None):
    try:
        result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=20)
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ValueError(f"Cannot start {command[0]}: {error}") from error
    output = (result.stdout + "\n" + result.stderr).strip()
    if result.returncode:
        reason = next(iter(output.splitlines()), "no diagnostic output")
        raise ValueError(f"{command[0]} exited with status {result.returncode}: {reason}")
    return output


def check_backend(name, settings):
    return [probe(command, probe_environment(name, settings)) for command in engine_commands(name, settings)]


def backend_version(name, settings):
    output = check_backend(name, settings)[0]
    return next(iter(output.splitlines()), "")


def choose_backend(explicit, registry):
    return explicit or registry["default"] or "vina"


def manage_backend(args, registry):
    registered = registry["backends"]
    if args.action == "list":
        if not registered:
            print("No backends registered. Vina is the default if it is available on PATH.")
        for name, settings in registered.items():
            default = " (default)" if registry["default"] == name else ""
            location = settings.get("path") or settings.get("conda_env") or "PATH / environment variables"
            print(f"{name}{default}: {location}")
        return 0
    name = args.backend
    if args.action == "default":
        if name not in registered:
            raise ValueError(f"Register {name} first with 'python -m secse backend add {name}'.")
        registry["default"] = name
    else:
        if args.executable and name not in PDBQT_BACKENDS:
            raise ValueError(f"Use --install-dir for {name}.")
        if args.install_dir and name in PDBQT_BACKENDS:
            raise ValueError(f"Use --executable for {name}.")
        if (args.conda_env or args.library_dir) and name != "unidock":
            raise ValueError("--conda-env and --library-dir are available for Uni-Dock.")
        if args.library_dir and not args.executable:
            raise ValueError("--library-dir requires --executable.")
        settings = {}
        if args.executable or args.install_dir:
            settings["path"] = str(Path(args.executable or args.install_dir).resolve(strict=True))
        if args.conda_env:
            settings["conda_env"] = args.conda_env
        if args.library_dir:
            library_dir = Path(args.library_dir).resolve(strict=True)
            if not library_dir.is_dir():
                raise ValueError(f"Not a library directory: {library_dir}")
            settings["library_dir"] = str(library_dir)
        check_backend(name, settings)
        registered[name] = settings
        if registry["default"] is None:
            registry["default"] = name
    save_registry(registry)
    if args.action == "add":
        print(f"Registered {name}.")
    print(f"Default backend: {registry['default']} ({registry_path()})")
    return 0


def doctor(args, registry):
    if args.all and args.target:
        raise ValueError("--target checks one backend. Select it with --backend.")
    names = list(registry["backends"]) if args.all else []
    names = names or [choose_backend(args.backend, registry)]
    failed = False

    def check(label, operation, hint):
        nonlocal failed
        try:
            result = operation()
        except Exception as error:
            print(f"FAIL {label}: {error}\n     {hint}")
            failed = True
        else:
            print(f"OK   {label}" + (f": {result}" if result else ""))

    for module in ("numpy", "pandas", "rdkit.Chem", "loguru", "pandarallel", "tqdm", "xlrd", "openbabel.pybel"):
        check(module, lambda module=module: getattr(importlib.import_module(module), "__version__", None),
              "Activate the SECSE environment or install requirements.txt into this Python.")

    def check_parallel():
        output = probe(["parallel", "--version"])
        if not output.startswith("GNU parallel"):
            raise ValueError("Expected GNU Parallel.")
        return output.splitlines()[0]

    check("GNU Parallel", check_parallel, "Install GNU Parallel in the SECSE environment.")
    check("Perl", lambda: probe(["perl", "-e", "print $^V"]), "Install Perl in the SECSE environment.")
    if any(name != "glide" for name in names):
        check("Open Babel CLI", lambda: probe(["obabel", "-V"]).splitlines()[0],
              "Install Open Babel and put obabel on PATH.")
    if args.prediction:
        check("Chemprop", lambda: probe(["chemprop", "--help"]) and "CLI starts",
              "Use Python 3.11+ and install 'chemprop~=2.1.0' with pip.")
    for name in names:
        check(name, lambda name=name: backend_version(name, registry["backends"].get(name, {})),
              f"See the {name} installation instructions in README.md.")
        if args.target:
            check("target", lambda: validate_target(name, args.target),
                  f"Supply a {TARGET_SUFFIXES[name]} target for {name}.")
    if not failed:
        print("Run a demo to verify docking, GPU access or license availability.")
    return int(failed)


def run_config(path, backend, registry):
    return subprocess.run(
        [sys.executable, "-m", "secse.run_secse", "--config", str(path)],
        env=engine_environment(backend, registry["backends"].get(backend, {})),
    ).returncode


def quick_run(args, registry):
    if getattr(args, "config", None):
        if any(value is not None for key, value in vars(args).items() if key not in {"command", "config"}):
            raise ValueError("--config cannot be combined with quick-run options.")
        path = Path(args.config).resolve(strict=True)
        config = read_config(path)
        backend = docking_backend(config.get("docking", "docking_program", fallback="vina"))
        return run_config(path, backend, registry)

    demo = args.command == "demo"
    backend = choose_backend(args.backend, registry)
    target = args.target
    if demo and backend in PDBQT_BACKENDS:
        target = target or DEMO_DIR / "PHGDH_6RJ3_for_vina.pdbqt"
    if not target:
        raise ValueError(f"Pass --target with a {TARGET_SUFFIXES[backend]} file for {backend}.")
    target = validate_target(backend, target)
    fragments = DEMO_DIR / "demo_1020.smi" if demo else Path(args.fragments).resolve(strict=True)
    if not fragments.is_file() or not fragments.stat().st_size:
        raise ValueError(f"Missing or empty fragments file: {fragments}")
    center = (20.9, -10.4, 3.0) if demo else args.center
    box = (15.0,) * 3 if demo else args.box_size or (20.0,) * 3
    if backend in PDBQT_BACKENDS and (not center or not all(map(math.isfinite, center))):
        raise ValueError(f"Pass --center X Y Z for {backend}.")
    if not all(math.isfinite(size) and size > 0 for size in box):
        raise ValueError("--box-size requires three positive values.")
    if not demo and backend not in PDBQT_BACKENDS and (args.center or args.box_size):
        raise ValueError(f"{backend} uses the prepared grid instead of --center or --box-size.")
    if args.cpu is not None and args.cpu < 1:
        raise ValueError("--cpu must be positive.")
    if args.generations is not None and args.generations < 0:
        raise ValueError("--generations must be zero or more.")
    check_backend(backend, registry["backends"].get(backend, {}))

    if args.workdir:
        workdir = Path(args.workdir).resolve()
        workdir.mkdir(parents=True, exist_ok=True)
        if any(workdir.iterdir()):
            raise ValueError(f"Choose an empty output directory: {workdir}")
    else:
        workdir = Path(mkdtemp(prefix="secse-phgdh-" if demo else "secse-run-", dir=Path.cwd()))
    if demo:
        with fragments.open(encoding="utf-8") as source:
            fragments = workdir / "phgdh-seeds.smi"
            fragments.write_text("".join(islice(source, 2)), encoding="utf-8")

    config = read_config(DEMO_DIR / "phgdh_demo_vina.ini")
    generations = args.generations if args.generations is not None else (0 if demo else 1)
    config["general"].update(
        project_code="PHG" if demo else "SECSE", workdir=escape_path(workdir),
        fragments=escape_path(fragments), num_gen=str(generations), start_gen="0", gpu="1",
        cpu=str(args.cpu or min(os.cpu_count() or 1, 4)),
    )
    config["docking"].update(docking_program=backend, target=escape_path(target))
    if backend in PDBQT_BACKENDS:
        for axis, coordinate, size in zip("xyz", center, box):
            config["docking"][axis] = str(coordinate)
            config["docking"][f"box_size_{axis}"] = str(size)
    config["prediction"]["mode"] = "0"
    path = workdir / "secse.ini"
    with path.open("x", encoding="utf-8") as stream:
        config.write(stream)
    print(f"Backend: {backend}\nResults: {workdir}\nConfig: {path}", flush=True)
    return run_config(path, backend, registry)


def escape_path(path):
    return str(path).replace("%", "%%")


def main(argv=None):
    parser = argparse.ArgumentParser(prog="python -m secse", description="Set up and run SECSE")
    subparsers = parser.add_subparsers(dest="command", required=True)
    doctor_parser = subparsers.add_parser("doctor", help="check the environment and docking tools")
    selection = doctor_parser.add_mutually_exclusive_group()
    selection.add_argument("--backend", choices=NAMES)
    selection.add_argument("--all", action="store_true")
    doctor_parser.add_argument("--target", type=Path)
    doctor_parser.add_argument("--prediction", action="store_true", help="also check Chemprop")

    backends = subparsers.add_parser("backend", help="register or select a docking backend")
    actions = backends.add_subparsers(dest="action", required=True)
    add = actions.add_parser("add", help="register a native docking installation")
    add.add_argument("backend", choices=NAMES)
    location = add.add_mutually_exclusive_group()
    location.add_argument("--executable", type=Path)
    location.add_argument("--install-dir", type=Path)
    location.add_argument("--conda-env", help="name of a separate native Uni-Dock environment")
    add.add_argument("--library-dir", type=Path, help="private runtime libraries for a Uni-Dock executable")
    actions.add_parser("list", help="list registered backends")
    default = actions.add_parser("default", help="change the default backend")
    default.add_argument("backend", choices=NAMES)

    demo = subparsers.add_parser("demo", help="run the two-fragment PHGDH example")
    demo.add_argument("example", choices=("phgdh",))
    run = subparsers.add_parser("run", help="dock your fragments or run an existing INI")
    run.add_argument("--config", type=Path)
    run.add_argument("--fragments", type=Path)
    run.add_argument("--center", type=float, nargs=3, metavar=("X", "Y", "Z"))
    run.add_argument("--box-size", type=float, nargs=3, metavar=("X", "Y", "Z"))
    for command in (demo, run):
        command.add_argument("--backend", choices=NAMES)
        command.add_argument("--target", type=Path)
        command.add_argument("--workdir", type=Path, help="empty directory for results and the generated INI")
        command.add_argument("--cpu", type=int)
        command.add_argument("--generations", type=int, help="number of growth generations (demo: 0, run: 1)")

    args = parser.parse_args(argv)
    if args.command == "run" and args.config is None and args.fragments is None:
        parser.error("run requires --fragments or --config")
    try:
        registry = load_registry()
        if args.command == "backend":
            return manage_backend(args, registry)
        if args.command == "doctor":
            return doctor(args, registry)
        return quick_run(args, registry)
    except (OSError, ValueError, configparser.Error) as error:
        print(f"SECSE: {error}", file=sys.stderr)
        return 1
