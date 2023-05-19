import os
import sys

sys.path.append(os.getenv("SECSE"))
import copy
import sqlite3
import pandas as pd
import rdkit
from pandarallel import pandarallel
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from utilities.wash_mol import get_bridged_atoms, neutralize_atoms
from utilities.load_rules import json_to_DB
from utilities.function_helper import shell_cmd_execute

rdkit.RDLogger.DisableLog("rdApp.*")

RULE_DB = os.path.join(os.getenv("SECSE"), "growing/mutation/rules_demo.db")


class Mutation:

    def __init__(self, num, workdir, rule_db=RULE_DB):
        # self.load_reaction()
        self.workdir = workdir
        self.rule_db = rule_db
        # self.load_buildingblock(num=num)
        self.rules_dict = {}
        self.load_common_rules()
        self.load_spacer_rings_rules()

        # drop unwanted rules where Priority < 0
        self.rules_dict = {k: v for k, v in self.rules_dict.items() if int(v[1]) > 0}
        self.out_product_smiles = []
        self.input_smiles = None
        self.mol = None

    def load_common_rules(self, tables=None):
        if tables is None:
            tables = ['B-001',
                      'G-001', 'G-003', 'G-004', 'G-005', 'G-006', 'G-007',
                      'M-001', 'M-002', 'M-003', 'M-004', 'M-005', 'M-006', 'M-007', 'M-008', 'M-009', 'M-010'
                      ]
        rules_dict = {}
        for table in tables:
            try:
                sql = 'select * from "{0}"'.format(table)
                conn = sqlite3.connect(self.rule_db)
                conn.row_factory = sqlite3.Row
                c = conn.cursor()
                c.execute(sql)
                rs = c.fetchall()
                for row in rs:
                    row = dict(row)
                    rules_dict[row["Rule ID"]] = (rdChemReactions.ReactionFromSmarts(row["SMARTS"]), row['Priority'])
            except sqlite3.OperationalError:
                print("No rule class: ", table)
                pass
        self.rules_dict.update(rules_dict)

    def load_spacer_rings_rules(self):
        rules_dict = {}
        try:
            sql = 'select * from "{}"'.format("G-002")
            conn = sqlite3.connect(self.rule_db)
            conn.row_factory = sqlite3.Row
            c = conn.cursor()
            c.execute(sql)
            rs = c.fetchall()
            for row in rs:
                row = dict(row)
                pri = int(row['Spacer Priority']) * int(row['Ring Priority'])
                rules_dict[row["Rule ID"]] = (rdChemReactions.ReactionFromSmarts(row["SMARTS"]), str(pri))
            self.rules_dict.update(rules_dict)
        except sqlite3.OperationalError:
            print("No rule class: G-002")

    # set smiles
    def load_mol(self, input_smiles):
        self.clean()
        self.input_smiles = input_smiles
        # uncharged each atom
        self.mol = Chem.MolFromSmiles(self.input_smiles)
        assert self.mol, "Can not read smiles"
        if self.input_smiles.count("-") + self.input_smiles.count("+") > 0:
            self.mol = neutralize_atoms(self.mol)
            # self.input_smiles = Chem.MolToSmiles(self.mol)

    def reaction(self, rxn, react, item, partner, priority):
        try:
            products = rxn.RunReactants(react)
            uniq = set()
            for mol_tuple in products:
                Chem.SanitizeMol(mol_tuple[0])
                # enumerator = rdMolStandardize.TautomerEnumerator()
                # canon = enumerator.Canonicalize(mol_tuple[0])
                # smi = Chem.MolToSmiles(Chem.RemoveHs(canon), isomericSmiles=True, kekuleSmiles=False)
                smi = Chem.MolToSmiles(Chem.RemoveHs(mol_tuple[0]), isomericSmiles=True, kekuleSmiles=False)
                uniq.add(smi)
            for smi in uniq:
                self.out_product_smiles.append((smi, item, partner, priority))
        except Exception as e:
            # print(e)
            pass

    # add 2021.1.7
    # modify 2021.01.14
    def single_point_mutate(self):
        mol = self.spiro_atom_label()
        for item in self.rules_dict:
            rxn = self.rules_dict[item][0]
            priority = self.rules_dict[item][1]
            if mol.HasSubstructMatch(rxn.GetReactantTemplate(0)):
                self.reaction(rxn, (mol,), item, "", priority)
        self.protected_atom_label_remove()
        return self.out_product_smiles

    def spiro_atom_label(self):
        mol = copy.deepcopy(self.mol)
        ri = mol.GetRingInfo()

        # spiro_ sma = '[*r3,*r4,*r5,*r6;R2X4$([*,*,*,*](@[r3,r4,r5,r6,r7])(@[r3,r4,r5,r6,r7])(@[r3,r4,r5,r6,
        # r7])@[r3,r4,r5,r6,r7])]'
        spiro_sma = '[x4]'
        spiro_atoms = mol.GetSubstructMatches(Chem.MolFromSmarts(spiro_sma))

        res = set()
        for ring in ri.AtomRings():
            for spi in spiro_atoms:
                tmp = set(spi).intersection(set(ring))
                if tmp:
                    res = res.union(ring)

        for index in res:
            mol.GetAtomWithIdx(index).SetProp('_protected', '1')

        self.mol = mol
        return mol

    def bridged_atom_label(self):

        mol = self.mol
        brigded_atoms = get_bridged_atoms(mol)
        ri = mol.GetRingInfo()
        res = set()
        for ring in ri.AtomRings():
            for bri in brigded_atoms:
                tmp = set(bri).intersection(set(ring))
                if tmp:
                    res = res.union(ring)

        for index in res:
            mol.GetAtomWithIdx(index).SetProp('_protected', '1')

        self.mol = mol
        return mol

    def protected_atom_label_remove(self):
        mol = self.mol
        for idx in range(len(mol.GetAtoms())):
            if mol.GetAtomWithIdx(idx).HasProp('_protected'):
                mol.GetAtomWithIdx(idx).ClearProp('_protected')
        self.mol = mol
        return mol

    def clean(self):
        self.input_smiles = None
        self.out_product_smiles = []


def mutation_df(df: pd.DataFrame, workdir, cpu_num, gen=1, rule_db=None, project_code="GEN"):
    workdir = os.path.join(workdir, "generation_" + str(gen))

    if rule_db is None:
        mutation = Mutation(5000, workdir)
    else:
        mutation = Mutation(5000, workdir, rule_db=rule_db)

    def mutation_per_row(mut: Mutation, smi):
        # mutation for each seed molecule
        try:
            mut.load_mol(smi)
        except AssertionError:
            return None
        mut.single_point_mutate()
        return mut.out_product_smiles

    mut_df = df.copy()
    if mut_df.shape[0] == 1:
        mut_df["smiles_gen_" + str(gen)] = mut_df["smiles_gen_" + str(gen - 1)].apply(
            lambda x: mutation_per_row(mutation, x))
    else:
        pandarallel.initialize(verbose=0, nb_workers=cpu_num)
        mut_df["smiles_gen_" + str(gen)] = mut_df["smiles_gen_" + str(gen - 1)].parallel_apply(
            lambda x: mutation_per_row(mutation, x))
    mut_df = mut_df.dropna(subset=["smiles_gen_" + str(gen)]).reset_index(drop=True)
    n = 1
    mut_path = os.path.join(workdir, "mutation")
    with open(mut_path + ".raw", "w") as f:
        header = list(mut_df.columns[:-1]) + ["smiles_gen_" + str(gen), "id_gen_" + str(gen),
                                              "reaction_id_gen_" + str(gen), "partner_gen_" + str(gen),
                                              "priority_gen_" + str(gen)]
        for i in mut_df.values.tolist():
            last_gen_info = list(map(str, i[:-1]))
            # keep parent mol
            f.write(",".join(last_gen_info + [last_gen_info[0], last_gen_info[1].split("-dp")[0].split("-C")[0],
                                              "Na-Na-Na", "", "3"]) + "\n")
            # write mutation mols
            for info in i[-1]:
                info = list(map(str, info))
                new_line = last_gen_info + [info[0]] + [
                    project_code.upper() + "_" + str(gen) + "_M_" + str(n).zfill(9)] + info[1:]
                f.write(",".join(new_line) + "\n")
                n += 1
    # drop duplicates product smiles by awk
    cmd_dedup = ["awk -F',' '!seen[$(NF-4)]++'", mut_path + ".raw ", ">", mut_path + ".csv"]
    shell_cmd_execute(cmd_dedup)

    return header
