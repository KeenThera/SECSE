import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory
from textwrap import dedent
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd
from rdkit import Chem

SOURCE = Path(__file__).resolve().parents[1] / "secse"
sys.path.insert(0, str(SOURCE))

from evaluate.docking import dock_by_py_autodock_gpu, dock_by_py_vina, dock_by_unidock
from evaluate.glide_docking import dock_by_glide
from evaluate.sdf_output import validate_scored_sdf
from grow_processes import Grow
from report.grow_path import write_growth
from scoring.diversity_score import clustering
from scoring.docking_score_prediction import get_pre


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        directory = TemporaryDirectory(prefix="secse regression ")
        self.addCleanup(directory.cleanup)
        self.root = Path(directory.name)

    def executable(self, name, source):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("#!/usr/bin/env python3\n" + dedent(source))
        path.chmod(0o755)
        return path

    def sdf(self, path, field="docking score"):
        mol = Chem.MolFromSmiles("CCO")
        mol.SetProp("_Name", "seed")
        mol.SetProp(field, "-7.25")
        with Chem.SDWriter(str(path)) as writer:
            writer.write(mol)

    def test_exit_codes(self):
        for arguments, expected in [(["--unknown"], 2), (["--config", str(self.root / "absent.ini")], 1)]:
            result = subprocess.run([sys.executable, SOURCE / "run_secse.py", *arguments],
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, expected, result.stderr)

    def test_resume_reads_legacy_seed_order_and_preserves_scores(self):
        workdir = self.root / "generation_1"
        workdir.mkdir()
        scores = workdir / "docked_gen_1.csv"
        scores.write_text("smiles,id,docking score,le_ln\nCCO,001,-7,-4\nCC,002,-6,-3\n")
        (workdir / "best_fragment_gen_1.smi").write_text("smiles_gen_1\tid_gen_1\nCC\t002\nCCO\t001\n")
        (workdir / "filter.csv").write_text("id_gen_1\n001\n002\n")
        self.sdf(workdir / "docking_outputs_with_score.sdf")
        workflow = Grow(0, "unused", self.root, 10, "vina", "unused", 1, 0, "unused")
        with patch.object(workflow, "docking_sh") as dock, patch("grow_processes.write_growth"):
            workflow.grow()
        dock.assert_not_called()
        self.assertEqual(workflow.winner_df["id_gen_1"].tolist(), ["002", "001"])
        self.assertEqual(workflow.winner_df["docking score_gen_1"].tolist(), [-6, -7])

    @unittest.skipUnless(shutil.which("perl"), "Perl is needed for SDF selection")
    def test_generation_zero_produces_report(self):
        workdir = self.root / "generation_0"
        workdir.mkdir()
        (workdir / "docked_gen_0.csv").write_text("smiles,id,docking score,le_ln\nCCO,seed,-7.25,-4\n")
        self.sdf(workdir / "docking_outputs_with_score.sdf")
        config = self.root / "config file.ini"
        config.write_text(f"[general]\nworkdir = {self.root}\n")
        write_growth(config, 0, 0)
        report, = self.root.glob("*_with_grow_path.csv")
        self.assertEqual(pd.read_csv(report)["id"].tolist(), ["seed"])
        validate_scored_sdf(self.root / "selected.sdf")

    def test_glide_status_and_output_normalization(self):
        ligand = self.root / "seed input.smi"
        ligand.write_text("CCO\tseed\n")
        receptor = self.root / "target grid.zip"
        receptor.touch()
        sdf = self.root / "reference.sdf"
        self.sdf(sdf, "r_i_glide_gscore")
        self.executable("Glide tools/ligprep", "raise SystemExit(42)")
        self.executable("Glide tools/glide", "raise SystemExit(43)")
        with patch.dict(os.environ, {"SCHRODINGER": str(self.root / "Glide tools")}):
            with self.assertRaisesRegex(RuntimeError, "status 42"):
                dock_by_glide(self.root, ligand, receptor, 0, "SP", 1)
            self.executable("Glide tools/ligprep", "pass")
            with self.assertRaisesRegex(RuntimeError, "status 43"):
                dock_by_glide(self.root, ligand, receptor, 0, "SP", 1)
            self.executable("Glide tools/glide", f"""
                from shutil import copyfile
                copyfile({str(sdf)!r}, "glide_gen_0_lib.sdf")
            """)
            dock_by_glide(self.root, ligand, receptor, 0, "SP", 1)
        validate_scored_sdf(self.root / "docking_outputs_with_score.sdf")

    @unittest.skipUnless(shutil.which("parallel"), "GNU Parallel is needed for the backend wrappers")
    def test_open_source_backends_with_tool_stubs(self):
        ligand = self.root / "input with spaces.smi"
        ligand.write_text("CCO\tseed\n")
        receptor = self.root / "target with spaces.pdbqt"
        receptor.touch()
        template = self.root / "reference.sdf"
        self.sdf(template)
        prep = self.executable("ligprep stub", f"""
            from pathlib import Path
            from shutil import copyfile
            import sys
            assert Path(sys.argv[3]).read_text() == "CCO\\tseed\\n"
            folder = Path(sys.argv[2]) / "ligands_for_docking"
            copyfile({str(template)!r}, folder / "seed.sdf")
            (folder / "seed.pdbqt").touch()
        """)
        self.executable("bin/obabel", f"""
            from pathlib import Path
            import sys
            from rdkit import Chem
            from rdkit.Chem import AllChem
            args = sys.argv
            content = Path(args[args.index("-ipdbqt") + 1]).read_text()
            out = Path(args[args.index("-O") + 1])
            mol = Chem.SDMolSupplier({str(template)!r})[0]
            AllChem.EmbedMolecule(mol, randomSeed=7)
            out.with_name(out.stem + "1.pdb").write_text(content + Chem.MolToPDBBlock(mol))
        """)
        vina = self.executable("vina stub", """
            from pathlib import Path
            import sys
            args = sys.argv
            assert Path(args[args.index("--ligand") + 1]).is_file()
            assert "seed = 12345" in Path(args[args.index("--config") + 1]).read_text()
            Path(args[args.index("--out") + 1]).write_text("REMARK VINA RESULT: -7.25 0 0\\n")
        """)
        unidock = self.executable("unidock stub", """
            from pathlib import Path
            import sys
            args = sys.argv
            assert Path(args[args.index("--receptor") + 1]).is_file()
            index = Path(args[args.index("--ligand_index") + 1])
            assert index.read_text().split() == ["seed.pdbqt"]
            assert Path("seed.pdbqt").is_file()
            (Path(args[args.index("--dir") + 1]) / "seed_out.pdbqt").write_text("REMARK VINA RESULT: -7.25 0 0\\n")
        """)
        self.executable("gpu tools/bin/autodock_gpu_128wi", """
            from pathlib import Path
            import sys
            args = sys.argv
            assert Path(args[args.index("--ffile") + 1]).is_file()
            assert args[args.index("-D") + 1] == "1"
            output = Path(args[args.index("--resnam") + 1]).with_suffix(".dlg")
            output.write_text("DOCKED: USER    Estimated Free Energy of Binding    = -7.25 kcal/mol\\n")
        """)
        env = {"PATH": f"{self.root / 'bin'}:{os.environ['PATH']}", "SECSE_PYTHON": str(prep),
               "VINA": str(vina), "UNIDOCK": str(unidock), "AUTODOCK_GPU": str(self.root / "gpu tools")}
        with patch.dict(os.environ, env):
            for backend, args in [(dock_by_py_vina, (1, 0, 0, 0)), (dock_by_unidock, (1, 0, 0, 0)),
                                  (dock_by_py_autodock_gpu, (1, 1))]:
                with self.subTest(backend=backend.__name__):
                    workdir = self.root / backend.__name__
                    workdir.mkdir()
                    backend(workdir, ligand, receptor, *args)
                    sdf = workdir / "docking_outputs_with_score.sdf"
                    validate_scored_sdf(sdf)
                    self.assertEqual(Chem.SDMolSupplier(str(sdf))[0].GetProp("docking score"), "-7.25")

    def test_prediction_excludes_exact_ids_only(self):
        workdir = self.root / "generation_1"
        workdir.mkdir()
        (self.root / "prediction").mkdir()
        (workdir / "mols_for_docking.smi").write_text("CC\t001\n")
        (workdir / "filter.csv").write_text("smiles_gen_1,id_gen_1\nCC,001\nCCC,001-other\n")
        result = get_pre(self.root, 1)
        self.assertEqual(Path(result).read_text(), "001-other,CCC\n")

    @unittest.skipUnless(shutil.which("parallel"), "GNU Parallel is needed for Chemprop prediction")
    def test_chemprop_stages_checkpoint_and_selects_top_predictions(self):
        train = self.root / "train input.csv"
        train.write_text("Smiles,docking score\n" + "".join(f"{'C' * n},-5\n" for n in range(1, 9)))
        candidates = self.root / "predict input.csv"
        candidates.write_text("m1,CC\nm2,CCO\n")
        self.executable("bin/chemprop", """
            from pathlib import Path
            import os
            import sys
            args = sys.argv
            if args[1] == "train":
                assert os.environ["TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD"] == "1"
                assert args[args.index("--split-sizes") + 1:args.index("--split-sizes") + 4] == ["0.5", "0.25", "0.25"]
                directory = Path(args[args.index("--save-dir") + 1])
                assert not list(directory.iterdir())
                model = directory / "model_0"
                model.mkdir()
                (model / "best.pt").write_text("new model")
            else:
                assert Path(args[args.index("--model-paths") + 1]).read_text() == "new model"
                input = Path(args[args.index("--test-path") + 1])
                assert input.read_text().startswith("id,smiles\\n")
                count = len(input.read_text().splitlines()) - 1
                assert args[args.index("--batch-size") + 1] == ("1" if count == 1 else "64")
                output = Path(args[args.index("--preds-path") + 1])
                first = "m1,CC,-6\\n" if count == 2 else ""
                output.write_text("id,smiles,docking score\\n" + first + "m2,CCO,-8\\n")
        """)
        saved_model = self.root / "prediction/models/G1_seed22/model_0/best.pt"
        saved_model.parent.mkdir(parents=True)
        saved_model.write_text("old model")
        command = ["bash", SOURCE / "scoring/chemprop_pre.sh", self.root, train, candidates, "1", "1", "22"]
        with patch.dict(os.environ, {"PATH": f"{self.root / 'bin'}:{os.environ['PATH']}"}):
            subprocess.run(command, check=True)
            candidates.write_text("m2,CCO")
            subprocess.run(command, check=True)
        self.assertEqual(saved_model.read_text(), "new model")
        self.assertEqual((self.root / "generation_1_pre/mols_for_docking_pred.smi").read_text(), "CCO\tm2\n")

    def test_rejects_missing_score_and_truncated_sdf(self):
        path = self.root / "output.sdf"
        self.sdf(path, "unrelated")
        with self.assertRaisesRegex(ValueError, "unscored SDF record"):
            validate_scored_sdf(path)
        self.sdf(path)
        path.write_text(path.read_text().replace("$$$$\n", ""))
        with self.assertRaisesRegex(ValueError, "missing final SDF record delimiter"):
            validate_scored_sdf(path)

    def test_clustering_preserves_first_center_for_ties(self):
        mols = pd.DataFrame({"smiles": ["C", "N", "O", "O", "Cl"]})
        np.random.seed(13)
        with patch.object(pd.Series, "parallel_apply", pd.Series.apply, create=True), \
                patch("scoring.diversity_score.pandarallel.initialize"):
            result = clustering(mols, "smiles", 2, 1, k=2)
        np.random.seed(13)
        first = mols.sample(1).index[0]
        other = next(i for i in mols.index if mols.at[i, "smiles"] != mols.at[first, "smiles"])
        expected = [other if smi == mols.at[other, "smiles"] else first for smi in mols["smiles"]]
        self.assertEqual(result["cluster_center_gen_2"].tolist(), expected)


if __name__ == "__main__":
    unittest.main()
