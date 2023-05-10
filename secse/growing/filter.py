#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: filter.py 
@time: 2020/11/16/13:14
"""
import argparse
import os
import sys
import time
import configparser

sys.path.append(os.getenv("SECSE"))

import rdkit.Chem as Chem
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt, CalcNumHBD, CalcNumHBA
from rdkit.Chem import Descriptors, AllChem
import json

from utilities.ring_tool import RingSystems
from utilities.substructure_filter import StructureFilter
from utilities.wash_mol import wash_mol, neutralize, charge_mol, get_rotatable_bound_num, get_rigid_body_num


class Filter:
    def __init__(self, gen, config_path):
        self.gen = int(gen)
        self.input_smiles = None
        self.mol = None
        self.pains_smarts = None
        self.strutFilter = StructureFilter()

        config = configparser.ConfigParser()
        config.read(config_path)
        print(config_path)
        self.MW = config.getfloat("properties", "MW")
        self.logP_lower = config.getfloat("properties", "logP_lower")
        self.logP_upper = config.getfloat("properties", "logP_upper")
        self.chiral_center = config.getint("properties", "chiral_center")
        self.heteroatom_ratio = config.getfloat("properties", "heteroatom_ratio")
        self.rotatable_bound_num = config.getint("properties", "rotatable_bound_num")
        self.rigid_body_num = config.getint("properties", "rigid_body_num")

    def load_mol(self, input_smiles):
        self.clean()
        self.input_smiles = input_smiles
        self.mol = Chem.MolFromSmiles(self.input_smiles)

        # uncharged each atom
        if self.input_smiles.count("-") + self.input_smiles.count("+") > 0:
            self.mol, self.input_smiles = neutralize(self.input_smiles)

        if self.mol is None:
            self.input_smiles = wash_mol(self.input_smiles)
            self.mol = Chem.MolFromSmiles(self.input_smiles)
            if self.mol is None:
                self.input_smiles = "C"
                self.mol = Chem.MolFromSmiles(self.input_smiles)

    def clean(self):
        self.input_smiles = None
        self.mol = None

    def lipinski_filter(self):
        mol = Chem.MolFromSmiles(self.input_smiles)
        violation_counter = 0

        if CalcExactMolWt(mol) >= 500:
            violation_counter += 1

        if CalcNumHBD(mol) >= 5:
            violation_counter += 1

        if CalcNumHBA(mol) >= 10:
            violation_counter += 1

        if Descriptors.MolLogP(mol) > 5:
            violation_counter += 1

        return violation_counter < 2

    def pp_filter(self):
        """
        property filter
        """
        mw = CalcExactMolWt(self.mol)

        if mw > self.MW:
            yield "MW"
        if self.gen > 3:
            if 81 > mw:
                yield "MW"
        if CalcNumHBD(self.mol) > 5:
            yield "HBD"
        if CalcNumHBA(self.mol) > 10:
            yield "HBA"

        logp = Descriptors.MolLogP(self.mol)
        if logp < self.logP_lower or logp > self.logP_upper:
            yield "cLogP"

        if Descriptors.TPSA(self.mol) > 200:
            yield "TPSA"
        if get_rotatable_bound_num(self.mol) > self.rotatable_bound_num:
            # rotatable bound customized @dalong
            yield "Rotatable Bound"
        if get_rigid_body_num(self.mol) > self.rigid_body_num:
            # rotatable bound customized @dalong
            yield "Rigid Body"
        yield "PASS"

    def load_pains_filter(self):
        # read smarts for pains
        with open(os.path.join(os.getenv("SECSE"), 'growing/pains_smarts.json')) as f:
            data = json.load(f)
        pains_smarts = dict((k, Chem.MolFromSmarts(v)) for k, v in data.items())
        self.pains_smarts = pains_smarts

    def alert_filter(self):
        self.load_pains_filter()
        for name in self.pains_smarts:
            sma = self.pains_smarts[name]
            if self.mol.HasSubstructMatch(sma):
                yield "PAINS"
        yield "PASS"

    def element_filter(self):
        f_count = self.input_smiles.count("F")
        br_count = self.input_smiles.count("Br")
        cl_count = self.input_smiles.count("Cl")
        i_count = self.input_smiles.count("I")
        s_count = self.input_smiles.count("S") + self.input_smiles.count("s")
        p_count = self.input_smiles.count("P")
        if not all([f_count <= 5, br_count < 3, cl_count <= 3, i_count <= 1, s_count <= 2, p_count <= 1]):
            yield "element"
        yield "PASS"

    def substructure_filter(self):
        # self.element_filter()
        yield self.strutFilter.sfilter(self.mol)

    def ring_system_filter(self):
        ring_sys = RingSystems(self.mol)
        if ring_sys.ring_check():
            yield "PASS"
        yield "RS"

    def custom_filter(self):
        # add Chiral center filter, cycle size less than 7, remove 3 continues hetero-atom
        chiral_tags = Chem.FindMolChiralCenters(self.mol, includeUnassigned=True, useLegacyImplementation=True)
        # the maximum number of chiral center <= 3
        if len(chiral_tags) > self.chiral_center:
            yield "CC"

        chiral_atom_list = set([x[0] for x in chiral_tags])
        rings = self.mol.GetRingInfo().AtomRings()

        if rings:
            # the maximum of ring size <= 7
            max_ring_size = max([len(x) for x in rings])
            if max_ring_size > 7:
                yield "max ring size >7"

            if len(chiral_tags) == 3:
                # 3 CCs should not in the same ring
                for ring in rings:
                    if len(set(ring).intersection(chiral_atom_list)) >= 3:
                        yield "chiral center in one ring >2"
        yield "PASS"

    def heteroatom_filter(self):
        hetero_ratio = Chem.rdMolDescriptors.CalcNumHeteroatoms(self.mol) / self.mol.GetNumHeavyAtoms()
        if hetero_ratio > self.heteroatom_ratio:
            yield "heteroatom_ratio"
        else:
            yield "PASS"

    def charge_filter(self):
        negative_charge = Chem.MolFromSmarts("[*-1]")
        positive_charge = Chem.MolFromSmarts("[*+1]")
        charged_smi = charge_mol(self.input_smiles)
        mol = Chem.MolFromSmiles(charged_smi)
        nc = len(mol.GetSubstructMatches(negative_charge))
        pc = len(mol.GetSubstructMatches(positive_charge))
        npc = nc + pc
        if npc <= 1:
            yield "PASS"
        elif npc == 2:
            if nc <= 1:
                yield "PASS"
            else:
                yield "Charge"
        else:
            yield "Charge"

    def similarity_filter(self):
        fp = AllChem.GetMorganFingerprintAsBitVect(self.mol, 2, 512)


def mol_filter(molfilter: Filter, smi):
    molfilter.load_mol(smi)
    pass_filter = [molfilter.pp_filter(),
                   molfilter.custom_filter(),
                   molfilter.charge_filter(),
                   molfilter.heteroatom_filter(),
                   molfilter.substructure_filter(),
                   molfilter.ring_system_filter(),
                   molfilter.alert_filter(),
                   ]
    for i in pass_filter:
        res = next(i)
        if res != "PASS":
            return res
    return "PASS"


def file_filter(file_path, workdir, gen, config):
    molsfilter = Filter(gen, config)
    with open(file_path, "r") as inf:
        with open(os.path.join(workdir, "filter_flag", os.path.basename(file_path)), "w") as outf:
            for line in inf.readlines():
                line = line.strip()
                smi = line.split(",")[-5]
                flag = mol_filter(molsfilter, smi)
                new_line = line + "," + flag + "\n"
                outf.write(new_line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Filter per file")
    parser.add_argument("file_path", help="File path")
    parser.add_argument("workdir", help="Workdir")
    parser.add_argument("gen", help="generation number")
    parser.add_argument("config", help="Configuration file")
    args = parser.parse_args()
    time1 = time.time()
    file_filter(args.file_path, args.workdir, args.gen, args.config)

    time2 = time.time()
    # mfilter = Filter()
    # mfilter.load_mol("C12CCCC3(CCCCC3)C1C4C5C(CC(C(C6CCCC7C6C8C9C(C%10CCC9C%10)C7C8)CCC%11)C%11C5)C2C4")
    # print(next(mfilter.ring_system_filter()))
