#!/usr/bin/env python
# -*- coding:utf-8 _*-
"""
@author: Lu Chong
@file: wash_mol.py
@time: 2021/02/08/14:13
"""

import random
from openbabel import openbabel
from openbabel import pybel
from rdkit import Chem


def wash_mol(smi):
    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetInAndOutFormats("smi", "can")
    ob_mol = openbabel.OBMol()
    ob_conversion.ReadString(ob_mol, smi)
    ob_conversion.Convert()
    res = ob_conversion.WriteString(ob_mol).strip()
    return res


def retreat_aromatic_nitrogen(smi):
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    mol.UpdatePropertyCache()
    Chem.GetSymmSSSR(mol)
    ri = mol.GetRingInfo()
    aromatic_n_atoms = mol.GetSubstructMatches(Chem.MolFromSmarts('[nr5]'))
    res = set()
    for ring in ri.AtomRings():
        n_at_ring = set()
        for n_atom in aromatic_n_atoms:
            tmp = set(n_atom).intersection(set(ring))
            if tmp:
                n_at_ring = n_at_ring.union(n_atom)
        if n_at_ring:
            res.add(random.choice(list(n_at_ring)))
    for index in res:
        atom = mol.GetAtomWithIdx(index)
        atom.SetNumExplicitHs(1)

    return Chem.MolToSmiles(mol)


def neutralize(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        smi = wash_mol(smi)
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return "C"
    new_mol = neutralize_atoms(mol)
    return new_mol, Chem.MolToSmiles(new_mol)


def neutralize_atoms(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol


def charge_mol(smi):
    mol = pybel.readstring("smi", smi)
    mol.removeh()
    mol.OBMol.AddHydrogens(False, True, 7.4)
    # mol.OBMol.CorrectForPH(7.4)
    charged_smi = mol.write("can", None, overwrite=False).strip()
    return charged_smi


def radical_filter(smi):
    mol = Chem.MolFromSmiles(smi)
    for a in mol.GetAtoms():
        if a.GetNumRadicalElectrons() == 1:
            return False
    return True


def get_bridged_atoms(mol):
    ri = mol.GetRingInfo()
    bond_rings = ri.BondRings()
    bridged_atoms = set()

    for i in range(len(bond_rings)):
        bond_ring_i = set(bond_rings[i])
        for j in range(i):
            bond_ring_j = set(bond_rings[j])
            common_bonds = bond_ring_i.intersection(bond_ring_j)

            if len(common_bonds) > 1:
                atoms = [0] * len(mol.GetAtoms())
                bridged_unit = ()
                for b in common_bonds:
                    atoms[mol.GetBondWithIdx(b).GetBeginAtomIdx()] += 1
                    atoms[mol.GetBondWithIdx(b).GetEndAtomIdx()] += 1
                for idx in range(len(atoms)):
                    if atoms[idx] == 1:
                        bridged_unit += (idx,)
                bridged_atoms.add(bridged_unit)
    return bridged_atoms


def get_rotatable_bound_num(mol):
    rb_smarts = Chem.MolFromSmarts(
        '[C^3!D1;!$(C(F)(F)F)]-!@[!Br!F!Cl!I!H3&!$(*#*)!D1;!$([!Br!F!Cl!I](F)(F)F)]')
    # sma = '[C^3!D1;!$(C(F)(F)F);!R;!$(C=O(N));!$(NC(=O));!$(C(=O)O);!$(C(=O)O)]-!@[!Br!F!Cl!I!H3&!$(*#*)!D1;!$([!Br!F!Cl!I](F)(F)F);!R;!$(C=O([N,O]));!$(NC(=O));!$(C(=O)O)]'
    return len((mol.GetSubstructMatches(rb_smarts)))


def get_rigid_body_num(mol):
    pattern = "[C^3!D1;!$(C(F)(F)F);!R;!$(C=O(N));!$(NC(=O));!$(C(=O)O);!$(C(=O)O)]-!@[!Br!F!Cl!I!H3&!$(*#*)!D1;!$([!Br!F!Cl!I](F)(F)F);!R;!$(C=O([N,O]));!$(NC(=O));!$(C(=O)O)]"
    rb = Chem.MolFromSmarts(pattern)
    return len((mol.GetSubstructMatches(rb)))
