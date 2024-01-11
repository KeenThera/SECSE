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
import rdkit.Chem as Chem
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt, CalcNumHBD, CalcNumHBA, CalcNumRotatableBonds
from rdkit.Chem import Descriptors, AllChem
from rdkit.Chem import QED
from rdkit.Chem import RDConfig
import json

sys.path.append(os.getenv("SECSE"))
from utilities.ring_tool import RingSystems
from utilities.substructure_filter import StructureFilter
from utilities.wash_mol import wash_mol, neutralize, charge_mol, get_keen_rotatable_bound_num, get_rigid_body_num
from utilities.open_filter import user_filter

sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer


class Filter:
    def __init__(self, gen, config_path):

        self.gen = int(gen)
        self.input_smiles = None
        self.mol = None
        self.pains_smarts = None

        config = configparser.ConfigParser()
        config.read(config_path)

        substructure_filter_file = config.get("properties", "substructure_filter")
        if substructure_filter_file == "0":
            self.strutFilter = StructureFilter()
        else:
            # print("Use additional substructure filter patters.")
            self.strutFilter = StructureFilter(substructure_filter_file)

        self.MW = config.getfloat("properties", "MW")
        self.logP_lower = config.getfloat("properties", "logP_lower")
        self.logP_upper = config.getfloat("properties", "logP_upper")
        self.chiral_center = config.getint("properties", "chiral_center")
        self.heteroatom_ratio = config.getfloat("properties", "heteroatom_ratio")
        self.rdkit_rotatable_bound_num = config.getint("properties", "rdkit_rotatable_bound_num")
        self.keen_rotatable_bound_num = config.getint("properties", "keen_rotatable_bound_num")
        self.rigid_body_num = config.getint("properties", "rigid_body_num")
        self.hbd = config.getint("properties", "HBD")
        self.hba = config.getint("properties", "HBA")
        self.tpsa = config.getfloat("properties", "TPSA")
        self.lipinski_violation = config.getint("properties", "lipinski_violation")
        self.qed = config.getfloat("properties", "qed")
        self.max_ring_size = config.getint("properties", "max_ring_size")
        self.max_ring_system_size = config.getint("properties", "max_ring_system_size")
        self.ring_system_count = config.getint("properties", "ring_system_count")
        self.bridged_site_count = config.getint("properties", "bridged_site_count")
        self.spiro_site_count = config.getint("properties", "spiro_site_count")
        self.fused_site_count = config.getint("properties", "fused_site_count")
        self.rdkit_sa_score = config.getint("properties", "rdkit_sa_score")

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

    def pp_filter(self):
        """
        property filter
        """
        violation_counter = 0

        mw = CalcExactMolWt(self.mol)
        if mw > self.MW:
            yield "MW"
        if mw > 500:
            violation_counter += 1
        if self.gen > 3:
            if 81 > mw:
                yield "MW"

        mol_hbd = CalcNumHBD(self.mol)
        if mol_hbd > self.hbd:
            yield "HBD"
        if mol_hbd > 5:
            violation_counter += 1

        mol_hba = CalcNumHBA(self.mol)
        if mol_hba > self.hba:
            yield "HBA"
        if mol_hba > 10:
            violation_counter += 1

        logp = Descriptors.MolLogP(self.mol)
        if logp < self.logP_lower or logp > self.logP_upper:
            yield "cLogP"
        if logp > 5:
            violation_counter += 1

        if violation_counter > self.lipinski_violation:
            yield "Lipinski Violation"

        if Descriptors.TPSA(self.mol) > self.tpsa:
            yield "TPSA"

        if CalcNumRotatableBonds(self.mol) > self.rdkit_rotatable_bound_num:
            yield "RDKit Rotatable Bonds"

        if get_keen_rotatable_bound_num(self.mol) > self.keen_rotatable_bound_num:
            # rotatable bound customized @dalong
            yield "Keen Rotatable Bounds"
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

    def substructure_filter(self):
        yield self.strutFilter.sfilter(self.mol)

    def ring_system_filter(self):
        ring_sys = RingSystems(self.mol)
        if ring_sys.ring_check(self.max_ring_system_size, self.bridged_site_count, self.spiro_site_count,
                               self.fused_site_count, self.ring_system_count):
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
            mol_max_ring_size = max([len(x) for x in rings])
            if mol_max_ring_size > self.max_ring_size:
                yield "max ring size"

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
        if mol is None:
            mol = self.mol
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

    def QED_filter(self):
        if QED.qed(self.mol) >= self.qed:
            yield "PASS"
        else:
            yield "QED"

    def SA_filter(self):
        sa_score = sascorer.calculateScore(self.mol)
        if sa_score <= self.rdkit_sa_score:
            yield "PASS"
        else:
            yield "SA score"

    def my_filter(self):
        if self.gen > 3:
            tag = user_filter(self.mol)
            if tag:
                yield "PASS"
            else:
                yield "CUSTOM"
        else:
            yield "PASS"


def mol_filter(molfilter: Filter, smi):
    molfilter.load_mol(smi)
    pass_filter = [molfilter.pp_filter(),
                   molfilter.custom_filter(),
                   molfilter.charge_filter(),
                   molfilter.heteroatom_filter(),
                   molfilter.substructure_filter(),
                   molfilter.ring_system_filter(),
                   molfilter.alert_filter(),
                   molfilter.QED_filter(),
                   molfilter.SA_filter(),
                   molfilter.my_filter()
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
