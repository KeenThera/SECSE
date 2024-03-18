#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Liu Shien
@file: ligprep.py
@time: 2021/4/1/16:28
@modify: 2022/3/1/12:04
@modify: 2023/5/5/14:22
"""
import argparse
import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolAlign
from openbabel import pybel
from openbabel import openbabel as ob

sys.path.append(os.getenv("SECSE"))
from utilities.wash_mol import charge_mol


def setero(mol, onlyUnassigned=True):
    if onlyUnassigned:
        opts = StereoEnumerationOptions(tryEmbedding=True)
    else:
        opts = StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
    isomers = tuple(EnumerateStereoisomers(mol, options=opts))
    res = []
    if len(isomers) > 1:
        for idx, tmp in enumerate(isomers):
            name = tmp.GetProp("_Name") + "-CC" + str(idx)
            tmp.SetProp("_Name", name)
            res.append(tmp)
        return res
    else:
        return list(isomers)


def tau(mol, can=True):
    params = rdMolStandardize.CleanupParameters()
    params.tautomerRemoveSp3Stereo = False
    params.tautomerRemoveBondStereo = False
    params.maxTautomers = 1000
    params.maxTransforms = 10000
    enumerator = rdMolStandardize.TautomerEnumerator(params)
    try:
        canon = enumerator.Canonicalize(mol)
    except Exception as e:
        print(e)
        return [mol]

    if can:
        return [canon]
    csmi = Chem.MolToSmiles(canon)
    res = [canon]
    tauts = enumerator.Enumerate(mol)
    smis = [Chem.MolToSmiles(x) for x in tauts]
    stpl = sorted((x, y) for x, y in zip(smis, tauts) if x != csmi)
    res += [y for x, y in stpl]

    new = []
    for idx, tmp in enumerate(res):
        name = tmp.GetProp("_Name") + "-CT" + str(idx)
        tmp.SetProp("_Name", name)
        new.append(tmp)

    return new


def to_3D(mol):
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True, maxAttempts=10000,
                          useRandomCoords=True)
    if mol.GetNumConformers() > 0:
        AllChem.UFFOptimizeMolecule(mol, 200, 10.0, -1)
        return mol
    else:
        return None


def gen_minimized_3D(path, rdmol, numConformer=1, rms_cutoff=1, addH=True):
    name = rdmol.GetProp("_Name")
    sdf_path = os.path.join(path, name + ".sdf")
    writer = Chem.SDWriter(sdf_path)
    if addH:
        rdmol = Chem.AddHs(rdmol, addCoords=True)

    param = rdDistGeom.ETKDGv2()
    param.pruneRmsThresh = rms_cutoff
    cids = rdDistGeom.EmbedMultipleConfs(rdmol, 50, param)
    mp = AllChem.MMFFGetMoleculeProperties(rdmol, mmffVariant='MMFF94s')
    AllChem.MMFFOptimizeMoleculeConfs(rdmol, numThreads=0, mmffVariant='MMFF94s')
    res = []
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(rdmol, mp, confId=cid)
        # ff.Initialize()
        ff.Minimize()
        e = ff.CalcEnergy()
        res.append((cid, e))
    sorted_res = sorted(res, key=lambda x: x[1])
    rdMolAlign.AlignMolConformers(rdmol)
    if len(sorted_res) > numConformer:
        selected = numConformer
    else:
        selected = len(sorted_res)
    # new = Chem.Mol(rdmol)
    # new.RemoveAllConformers()
    # min_conf = rdmol.GetConformer(sorted_res[0][0])
    # new.AddConformer(min_conf)
    for i in range(selected):
        cid = sorted_res[i][0]
        writer.write(rdmol, cid)
    writer.close()

    return sdf_path


def ionization(smi_string):
    return charge_mol(smi_string)


def sdf2pdbqt(sdf_path):
    path = os.path.dirname(sdf_path)
    name = os.path.basename(sdf_path).split(".")[0]
    num = 0
    for mol in pybel.readfile("sdf", sdf_path):
        mol.write("pdbqt", "{}.pdbqt".format(os.path.join(path, name)), overwrite=True)
        num += 1

    return num == 1


class LigPrep:
    def __init__(self, infile, workdir):
        self.infile = infile
        self.workdir = workdir
        self.mol_dict = {}

    def parse_infile(self):
        with open(self.infile, "r") as inf:
            for line in inf:
                tmp = line.strip().split()
                if len(tmp) < 2:
                    continue
                smi = tmp[0]
                id1 = tmp[1]
                smi = ionization(smi)

                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    continue
                mol.SetProp("_Name", id1)
                self.mol_dict[id1] = mol

    def process(self, des):
        dirc_name = "ligands_for_" + des
        path = os.path.join(self.workdir, dirc_name)
        os.makedirs(path, exist_ok=True)

        self.parse_infile()
        for gid in self.mol_dict:
            mol = self.mol_dict[gid]
            mystereo = setero(mol)

            mytau = []
            for stereo in mystereo:
                tmp = tau(stereo)
                mytau += tmp

            for newmol in mytau:
                if newmol is not None:
                    try:
                        if des == 'docking':
                            sdf_path = gen_minimized_3D(path, newmol)
                            sdf2pdbqt(sdf_path)
                        if des == 'shape':
                            gen_minimized_3D(path, newmol, 10)
                    except Exception as e:
                        print(e)
                        continue


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="LigPrep @dalong")
    parser.add_argument("workdir", help="Workdir")
    parser.add_argument("mols_smi", help="Seed fragments")
    parser.add_argument("--mode",
                        help="1: prepare pdbqt file for docking input; 2: prepare sdf file for shape based screening.",
                        type=int, default=1)

    args = parser.parse_args()
    lig = LigPrep(args.mols_smi, args.workdir)
    lig.parse_infile()
    if args.mode == 1:
        lig.process(des="docking")
    elif args.mode == 2:
        lig.process(des="shape")
