import configparser
import contextlib
import io
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from secse.backends import engine_environment, load_registry
from secse.cli import main


class QuickStartTests(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory(prefix="secse cli ")
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)
        env = patch.dict(os.environ, {"XDG_CONFIG_HOME": str(self.root / "config")})
        env.start()
        self.addCleanup(env.stop)

    def register(self, name):
        engine = self.root / f"{name} engine"
        engine.write_text(f"#!/bin/sh\nprintf '%s\\n' '{name} test'\n")
        engine.chmod(0o755)
        self.assertEqual(self.command("backend", "add", name, "--executable", str(engine))[0], 0)
        return engine

    def command(self, *args):
        output = io.StringIO()
        with contextlib.redirect_stdout(output), contextlib.redirect_stderr(output):
            status = main(args)
        return status, output.getvalue()

    def test_demo_uses_default_and_explicit_backend_without_changing_it(self):
        self.register("vina")
        self.register("unidock")
        with patch("secse.cli.run_config", return_value=0) as run:
            for backend, flags in (("vina", ()), ("unidock", ("--backend", "unidock"))):
                with self.subTest(backend=backend):
                    workdir = self.root / f"{backend} 100% results"
                    status, output = self.command("demo", "phgdh", "--workdir", str(workdir), *flags)
                    self.assertEqual(status, 0, output)
                    path, selected, _ = run.call_args.args
                    self.assertEqual(selected, backend)
                    config = configparser.ConfigParser()
                    config.read(path)
                    self.assertEqual(config["docking"]["docking_program"], backend)
                    self.assertEqual(config["general"]["workdir"], str(workdir))
                    self.assertEqual(config["prediction"]["mode"], "0")
                    seeds = Path(config["general"]["fragments"])
                    self.assertEqual([line.split()[1] for line in seeds.read_text().splitlines()], ["f1", "f2"])
        self.assertEqual(load_registry()["default"], "vina")
        self.assertEqual(self.command("backend", "default", "unidock")[0], 0)
        self.assertEqual(load_registry()["default"], "unidock")

    def test_doctor_checks_default_or_all_registered_engines(self):
        vina = self.register("vina")
        unidock = self.register("unidock")
        for flags, expected in (((), [str(vina)]), (("--all",), [str(vina), str(unidock)])):
            with patch("secse.cli.probe", return_value="GNU parallel test") as probe, \
                    patch("secse.cli.importlib.import_module"):
                self.assertEqual(self.command("doctor", *flags)[0], 0)
                engines = [call.args[0][0] for call in probe.call_args_list if call.args[0][0] in {str(vina), str(unidock)}]
                self.assertEqual(engines, expected)

    def test_quick_run_protects_existing_results(self):
        self.register("vina")
        workdir = self.root / "existing"
        workdir.mkdir()
        result = workdir / "selected.sdf"
        result.write_text("previous run")
        status, output = self.command("demo", "phgdh", "--workdir", str(workdir))
        self.assertEqual(status, 1)
        self.assertIn("Choose an empty output directory", output)
        self.assertEqual(result.read_text(), "previous run")

    def test_ini_backend_takes_precedence_over_default(self):
        self.register("vina")
        unidock = self.register("unidock")
        path = self.root / "existing.ini"
        path.write_text("[docking]\ndocking_program = UniDock\n")
        with patch("secse.cli.subprocess.run") as run:
            run.return_value.returncode = 0
            status, output = self.command("run", "--config", str(path))
        self.assertEqual(status, 0, output)
        self.assertEqual(run.call_args.kwargs["env"]["UNIDOCK"], str(unidock))

    def test_unidock_environment_and_library_are_scoped_to_the_engine(self):
        environment = engine_environment("unidock", {"conda_env": "gpu-tools"})
        self.assertEqual(environment["SECSE_UNIDOCK_CONDA_ENV"], "gpu-tools")
        self.assertNotIn("SECSE_UNIDOCK_LIBRARY_DIR", environment)
        environment = engine_environment("unidock", {"library_dir": "/opt/unidock/lib"})
        self.assertEqual(environment["SECSE_UNIDOCK_LIBRARY_DIR"], "/opt/unidock/lib")
        self.assertEqual(environment.get("LD_LIBRARY_PATH"), os.environ.get("LD_LIBRARY_PATH"))
