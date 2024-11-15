#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: docking.py
@time: 2021/9/6/11:22
"""
import argparse
import os
import shutil
import glob
import sys
from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem
from utilities.function_helper import shell_cmd_execute

sys.path.append(os.getenv("SECSE"))

VINA_SHELL = os.path.join(os.getenv("SECSE"), "evaluate", "ligprep_vina_parallel.sh")
AUTODOCK_GPU_SHELL = os.path.join(os.getenv("SECSE"), "evaluate", "ligprep_autodock_gpu.sh")
UNIDOCK_SHELL = os.path.join(os.getenv("SECSE"), "evaluate", "ligprep_unidock.sh")


def dock_by_py_vina(workdir, smi, receptor, cpu_num, x, y, z, box_size_x=20, box_size_y=20, box_size_z=20):
    cmd = list(map(str, [VINA_SHELL, workdir, smi, receptor, x, y, z, box_size_x, box_size_y, box_size_z, cpu_num]))
    shell_cmd_execute(cmd)
    merged_sdf(workdir, 0)


def dock_by_py_autodock_gpu(workdir, smi, receptor, cpu_num, gpu_num):
    cmd = list(map(str, [AUTODOCK_GPU_SHELL, workdir, smi, receptor, cpu_num, gpu_num]))
    shell_cmd_execute(cmd)
    merged_sdf(workdir, 1)


def dock_by_unidock(workdir, smi, receptor, cpu_num, x, y, z, box_size_x=20, box_size_y=20, box_size_z=20):
    if not os.environ.get("UNIDOCK"):
        os.environ["UNIDOCK"] = "unidock"
    cmd = list(map(str, [UNIDOCK_SHELL, workdir, smi, receptor, x, y, z, box_size_x, box_size_y, box_size_z, cpu_num]))
    shell_cmd_execute(cmd)
    for res_file in glob.glob(os.path.join(workdir, "pdb_files", "*.pdb")):
        new_name = os.path.basename(res_file).replace("_out", "")
        os.rename(res_file, os.path.join(workdir, "pdb_files", new_name))
    merged_sdf(workdir, 2)


def merged_sdf(workdir, program):
    # modify output sdf
    check_mols(workdir, program)
    out_sdf = os.path.join(workdir, "docking_outputs_with_score.sdf")
    cmd_cat = ["find", os.path.join(workdir, "sdf_files"), "-name \"*sdf\" | xargs -n 100 cat >", out_sdf]
    shell_cmd_execute(cmd_cat)
    # remove temporary files
    shutil.rmtree(os.path.join(workdir, "pdb_files"))
    shutil.rmtree(os.path.join(workdir, "ligands_for_docking"))
    shutil.rmtree(os.path.join(workdir, "docking_poses"))
    shutil.rmtree(os.path.join(workdir, "docking_split"))


def check_mols(workdir, program):
    files = os.listdir(os.path.join(workdir, "pdb_files"))
    for i in files:
        raw_id = i.rsplit("-dp", 1)[0]
        pdb_path = os.path.join(workdir, "pdb_files", i)
        sdf_path = os.path.join(workdir, "sdf_files", i.replace("pdb", "sdf"))
        raw_mol = Chem.SDMolSupplier(os.path.join(workdir, "ligands_for_docking", raw_id + ".sdf"))[0]
        mol = AllChem.MolFromPDBFile(pdb_path, removeHs=True)
        if mol:
            try:
                new = AllChem.AssignBondOrdersFromTemplate(raw_mol, mol)
            except ValueError:
                logger.error("Failed check: ", i)
                continue
            new = Chem.AddHs(new, addCoords=True)
            Chem.MolToMolFile(new, sdf_path)
            if program == 0 or program == 2:
                with open(pdb_path, "r") as pdb:
                    for line in pdb.readlines():
                        if line.startswith("REMARK VINA RESULT"):
                            score = line.split(":")[1][:10].replace(" ", "")
                            with open(sdf_path, "a") as sdf:
                                newline = "\n".join(["> <docking score>", score, "\n$$$$\n"])
                                sdf.write(newline)
            elif program == 1:
                with open(pdb_path, "r") as pdb:
                    for line in pdb.readlines():
                        if "kcal" in line:
                            score = line.split("kcal")[0].replace(" ", "")[6:]
                            with open(sdf_path, "a") as sdf:
                                newline = "\n".join(["> <docking score>", score, "\n$$$$\n"])
                                sdf.write(newline)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run Open-source Docking Program for SMILES Format.")
    parser.add_argument("program", help="Name of docking program, input vina or autodock-gpu", type=str)
    parser.add_argument("workdir", help="Workdir")
    parser.add_argument("mols_smi", help="Seed fragments")
    parser.add_argument("receptor", help="Target PDBQT")

    parser.add_argument("cpu_num", help="Number of CPU cores")

    parser.add_argument("--gpu_num", help="Number of GPUs")
    parser.add_argument("--x", help="Docking box x", type=float)
    parser.add_argument("--y", help="Docking box y", type=float)
    parser.add_argument("--z", help="Docking box z", type=float)

    parser.add_argument("--box_size_x", help="Docking box size x, default 20", type=float, default=20)
    parser.add_argument("--box_size_y", help="Docking box size y, default 20", type=float, default=20)
    parser.add_argument("--box_size_z", help="Docking box size z, default 20", type=float, default=20)

    args = parser.parse_args()
    if args.program == "vina":
        logger.info("Docking by Autodock Vina with {} CPUs...".format(args.cpu_num))
        dock_by_py_vina(args.workdir, args.mols_smi, args.receptor, args.cpu_num, args.x, args.y, args.z,
                        args.box_size_x, args.box_size_y, args.box_size_z)
    elif args.program == "autodock-gpu":
        logger.info("Docking by Autodock-GPU with {} CPUs and {} GPUs...".format(args.cpu_num, args.gpu_num))
        dock_by_py_autodock_gpu(args.workdir, args.mols_smi, args.receptor, args.cpu_num, args.gpu_num)
    else:
        logger.error("Please choose a docking program.")
