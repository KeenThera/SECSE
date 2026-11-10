#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: docking.py
@time: 2021/9/6/11:22
"""
import argparse
import shutil
import re
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem
if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from evaluate.sdf_output import validate_scored_sdf
from utilities.config import SOURCE_DIR
from utilities.function_helper import run_command


VINA_SHELL = SOURCE_DIR / "evaluate" / "ligprep_vina_parallel.sh"
AUTODOCK_GPU_SHELL = SOURCE_DIR / "evaluate" / "ligprep_autodock_gpu.sh"
UNIDOCK_SHELL = SOURCE_DIR / "evaluate" / "ligprep_unidock.sh"


def dock_by_py_vina(workdir, smi, receptor, cpu_num, x, y, z, box_size_x=20, box_size_y=20, box_size_z=20):
    cmd = [VINA_SHELL, workdir, smi, receptor, x, y, z, box_size_x, box_size_y, box_size_z, cpu_num]
    run_command(cmd)
    merged_sdf(workdir, 0)


def dock_by_py_autodock_gpu(workdir, smi, receptor, cpu_num, gpu_num):
    cmd = [AUTODOCK_GPU_SHELL, workdir, smi, receptor, cpu_num, gpu_num]
    run_command(cmd)
    merged_sdf(workdir, 1)


def dock_by_unidock(workdir, smi, receptor, cpu_num, x, y, z, box_size_x=20, box_size_y=20, box_size_z=20):
    cmd = [UNIDOCK_SHELL, workdir, smi, receptor, x, y, z, box_size_x, box_size_y, box_size_z, cpu_num]
    run_command(cmd)
    for pose in (Path(workdir) / "pdb_files").glob("*_out-dp*.pdb"):
        pose.rename(pose.with_name(pose.name.replace("_out-dp", "-dp")))
    merged_sdf(workdir, 2)


def merged_sdf(workdir, program):
    workdir = Path(workdir)
    check_mols(workdir, program)
    name = "docking_outputs_with_score.sdf"
    with TemporaryDirectory(prefix=".secse-docking-", dir=workdir) as temp_dir:
        staged = Path(temp_dir) / name
        with staged.open("wb") as destination:
            for path in (workdir / "sdf_files").glob("*.sdf"):
                with path.open("rb") as source:
                    shutil.copyfileobj(source, destination)
        validate_scored_sdf(staged)
        staged.replace(workdir / name)
    for name in ("pdb_files", "ligands_for_docking", "docking_poses", "docking_split"):
        shutil.rmtree(workdir / name)


def check_mols(workdir, program):
    workdir = Path(workdir)
    for pdb_path in (workdir / "pdb_files").glob("*.pdb"):
        raw_id = pdb_path.stem.rsplit("-dp", 1)[0]
        template = workdir / "ligands_for_docking" / f"{raw_id}.sdf"
        raw_mol = Chem.SDMolSupplier(str(template))[0]
        mol = AllChem.MolFromPDBFile(str(pdb_path), removeHs=True)
        if mol is None or raw_mol is None:
            logger.warning("Cannot read pose or ligand template: {}", pdb_path)
            continue
        try:
            mol = AllChem.AssignBondOrdersFromTemplate(raw_mol, mol)
        except ValueError as error:
            logger.warning("Cannot assign bond orders for {}: {}", pdb_path, error)
            continue
        score = None
        with pdb_path.open() as source:
            for line in source:
                if program in (0, 2) and line.startswith("REMARK VINA RESULT:"):
                    score = line.split(":", 1)[1].strip().split()[0]
                    break
                if program == 1:
                    match = re.search(r"REMARK\s+([-+0-9.eE]+)\s*kcal", line)
                    if match:
                        score = match.group(1)
                        break
        if score is None:
            raise ValueError(f"No docking score in {pdb_path}")
        mol = Chem.AddHs(mol, addCoords=True)
        sdf_path = workdir / "sdf_files" / pdb_path.with_suffix(".sdf").name
        Chem.MolToMolFile(mol, str(sdf_path))
        with sdf_path.open("a") as destination:
            destination.write(f"> <docking score>\n{score}\n\n$$$$\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run Open-source Docking Program for SMILES Format.")
    parser.add_argument("program", type=str.lower, choices=("vina", "autodock-gpu", "unidock"))
    parser.add_argument("workdir", help="Workdir")
    parser.add_argument("mols_smi", help="Seed fragments")
    parser.add_argument("receptor", help="Target PDBQT")

    parser.add_argument("cpu_num", help="Number of CPU cores")

    parser.add_argument("--gpu_num", help="Number of GPUs", type=int, default=1)
    parser.add_argument("--x", help="Docking box x", type=float)
    parser.add_argument("--y", help="Docking box y", type=float)
    parser.add_argument("--z", help="Docking box z", type=float)

    parser.add_argument("--box_size_x", help="Docking box size x, default 20", type=float, default=20)
    parser.add_argument("--box_size_y", help="Docking box size y, default 20", type=float, default=20)
    parser.add_argument("--box_size_z", help="Docking box size z, default 20", type=float, default=20)

    args = parser.parse_args()
    if args.program in ("vina", "unidock"):
        if None in (args.x, args.y, args.z):
            parser.error("--x, --y and --z are required for vina and unidock")
        dock = dock_by_py_vina if args.program == "vina" else dock_by_unidock
        dock(args.workdir, args.mols_smi, args.receptor, args.cpu_num, args.x, args.y, args.z,
                        args.box_size_x, args.box_size_y, args.box_size_z)
    elif args.program == "autodock-gpu":
        logger.info("Docking by Autodock-GPU with {} CPUs and {} GPUs...".format(args.cpu_num, args.gpu_num))
        dock_by_py_autodock_gpu(args.workdir, args.mols_smi, args.receptor, args.cpu_num, args.gpu_num)
