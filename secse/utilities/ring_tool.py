#!/usr/bin/env python  
# -*- coding:utf-8 _*-
"""
@author: Lu Chong
@file: ring tool.py
@time: 2021/02/07/14:17
"""
from rdkit import Chem


def ring_site_count(ring_atoms, systems):
    site_count = [-1]  # add -1 in case no ring site
    for ring_s in systems:
        ring_s = set(ring_s)
        count = 0
        for site in ring_atoms:
            site = set(site)
            if ring_s.intersection(site):
                count += 1
        site_count.append(count)
    return site_count


class RingSystems(object):
    def __init__(self, mol):
        self.mol = mol
        self.ri = self.mol.GetRingInfo()
        self.atom_rings = self.ri.AtomRings()
        self.bond_rings = self.ri.BondRings()
        self.systems = self.ring_systems()

    def ring_systems(self):
        systems = []
        for ring in self.atom_rings:
            ringAts = set(ring)
            nSystems = []
            for system in systems:
                nInCommon = len(ringAts.intersection(system))
                if nInCommon:
                    ringAts = ringAts.union(system)
                else:
                    nSystems.append(system)
            nSystems.append(ringAts)
            systems = nSystems
        return systems

    # ring size of each ring system
    def ring_systems_size(self):
        ring_sys_size = []
        for ring_s in self.systems:
            ring_s = set(ring_s)
            size = 0
            for ring in self.atom_rings:
                ring = set(ring)
                if ring_s.intersection(ring):
                    size += 1
            ring_sys_size.append(size)
        return ring_sys_size

    def get_spiro_atoms(self):
        spiro = []
        spiro_atoms = set()
        for i in range(len(self.atom_rings)):
            atom_ring_i = set(self.atom_rings[i])
            for j in range(i):
                atom_ring_j = set(self.atom_rings[j])
                common_atoms = atom_ring_i.intersection(atom_ring_j)
                if len(common_atoms) == 1:
                    atoms = [0] * len(self.mol.GetAtoms())
                    for a in common_atoms:
                        atoms[a] += 1

                    for idx in range(len(atoms)):
                        if atoms[idx] == 1:
                            spiro = (idx,)
                    spiro_atoms.add(spiro)
        return spiro_atoms

    def get_fused_atoms(self):
        fused_atoms = set()

        for i in range(len(self.bond_rings)):
            bond_ring_i = set(self.bond_rings[i])
            for j in range(i):
                bond_ring_j = set(self.bond_rings[j])
                common_bonds = bond_ring_i.intersection(bond_ring_j)
                if len(common_bonds) == 1:
                    atoms = [0] * len(self.mol.GetAtoms())
                    fused_unit = ()

                    for b in common_bonds:
                        atoms[self.mol.GetBondWithIdx(b).GetBeginAtomIdx()] += 1
                        atoms[self.mol.GetBondWithIdx(b).GetEndAtomIdx()] += 1
                    for idx in range(len(atoms)):
                        if atoms[idx] == 1:
                            fused_unit += (idx,)
                    fused_atoms.add(fused_unit)

        return fused_atoms

    def get_bridged_atoms(self):
        bridged_atoms = set()

        for i in range(len(self.bond_rings)):
            bond_ring_i = set(self.bond_rings[i])
            for j in range(i):
                bond_ring_j = set(self.bond_rings[j])
                common_bonds = bond_ring_i.intersection(bond_ring_j)

                if len(common_bonds) > 1:
                    atoms = [0] * len(self.mol.GetAtoms())
                    bridged_unit = ()
                    for b in common_bonds:
                        atoms[self.mol.GetBondWithIdx(b).GetBeginAtomIdx()] += 1
                        atoms[self.mol.GetBondWithIdx(b).GetEndAtomIdx()] += 1
                    for idx in range(len(atoms)):
                        if atoms[idx] == 1:
                            bridged_unit += (idx,)
                    bridged_atoms.add(bridged_unit)
        return bridged_atoms

    def spiro_site_count(self):
        return ring_site_count(self.get_spiro_atoms(), self.systems)

    def bridged_site_count(self):
        return ring_site_count(self.get_bridged_atoms(), self.systems)

    def fused_site_count(self):
        return ring_site_count(self.get_fused_atoms(), self.systems)

    def ring_system_count_filter(self, num=4):
        return len(self.systems) <= num

    def largest_ring_system_size_filter(self, num=3):
        return max(self.ring_systems_size() + [-1]) <= num

    def largest_spiro_site_filter(self, num=1):
        return max(self.spiro_site_count()) <= num

    def largest_fused_site_filter(self, num=3):
        return max(self.fused_site_count()) <= num

    def largest_bridged_site_filter(self, num=2):
        return max(self.bridged_site_count()) <= num

    def bridged_atom_is_aromatic_filter(self):
        bridged_atoms = self.get_bridged_atoms()
        for atom_cubic in bridged_atoms:
            for atom_idx in atom_cubic:
                atom = self.mol.GetAtomWithIdx(atom_idx)
                if atom.GetIsAromatic():
                    return False
        return True

    def ring_check(self, rssc, bsc, ssc, fsc, rsc):
        return all([self.largest_ring_system_size_filter(rssc),
                    self.largest_bridged_site_filter(bsc),
                    self.largest_spiro_site_filter(ssc),
                    self.largest_fused_site_filter(fsc),
                    self.ring_system_count_filter(rsc),
                    self.bridged_atom_is_aromatic_filter()])


if __name__ == '__main__':
    mol = Chem.MolFromSmiles("C1(C2)CC(CC3)NC3CC2C1")
    ringcheck = RingSystems(mol)
    print("Not Pass Filter" if not ringcheck.ring_check() else "Pass Filter")
    print(ringcheck.ring_systems_size())
