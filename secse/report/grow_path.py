#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: grow_path.py 
@time: 2021/01/19/13:42
"""
import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt
from pandarallel import pandarallel
import configparser

sys.path.append(os.getenv("SECSE"))
from scoring.ranking import read_dock_file
from utilities.function_helper import shell_cmd_execute

pandarallel.initialize(verbose=0)


def cal_mutation_dic(workdir, max_gen):
    mut_dic_all = dict()

    while max_gen > 0:
        mut_file = os.path.join(workdir, "generation_" + str(max_gen), "filter.csv")
        print(mut_file)
        with open(mut_file, "r") as f:
            lines = f.readlines()
        lines = [i.strip().split(",") for i in lines]
        mut_dic = {i[-6].split("-dp")[0].split("-C")[0]: [i[0], i[1].split("-dp")[0].split("-C")[0], i[-5], i[-4]] for i
                   in lines}
        mut_dic_all["gen" + str(max_gen)] = mut_dic
        max_gen -= 1
    return mut_dic_all


def merge_multi_generation(workdir, max_gen, file_path, dl_mode, config_path):
    df_lst = [pd.read_csv(os.path.join(workdir, "generation_" + str(i),
                                       "docked_gen_" + str(i) + ".csv")) for i in range(1, max_gen + 1)]
    if dl_mode == 2:
        dl_df = read_dock_file(os.path.join(workdir, "generation_{}_pre".format(max_gen),
                                            "docking_outputs_with_score.sdf"))
        dl_df["le_ln"] = dl_df.apply(
            lambda x: x["docking score"] / Chem.MolFromSmiles(x["smiles"]).GetNumHeavyAtoms(),
            axis=1)
        dl_df.columns = [i.lower() for i in list(dl_df.columns)]
        dl_df = dl_df.drop(columns=["molecule"])
        dl_df = dl_df.reindex(columns=df_lst[0].columns)
        config = configparser.ConfigParser()
        config.read(config_path)
        score_cutoff = config.getfloat("deep learning", "dl_score_cutoff")
        dl_df = dl_df[dl_df["docking score"] < score_cutoff]
        df_lst.append(dl_df)

    final_df = pd.concat(df_lst, axis=0).drop_duplicates(subset=["smiles"])
    final_df.to_csv(file_path, index=False)
    return final_df


def grow_path(mut_dic_all, mut_id):
    mut_id = mut_id.split("-dp")[0].split("-C")[0]
    try:
        gen_mol = int(mut_id.split("_")[-3])
    except IndexError:
        print(mut_id)
        return None
    mut_info_lst = []

    while gen_mol > 0:
        mut_info = mut_dic_all["gen" + str(gen_mol)][mut_id]
        if "." in mut_info[2]:
            gen_mol -= 1
            continue
        mut_info_lst.append(mut_info)
        mut_id = mut_info[1]
        gen_mol -= 1
    return mut_info_lst


def add_prop(merged_df_path):
    merged_df = pd.read_csv(merged_df_path)
    raw_cols = list(merged_df.columns)
    merged_df["mol"] = merged_df["smiles"].apply(Chem.MolFromSmiles)
    # check charge
    merged_df["charge flag"] = merged_df["mol"].apply(charge_filter)
    merged_df = merged_df[merged_df["charge flag"]]
    # add MW, logP
    merged_df["MW"] = merged_df["mol"].apply(CalcExactMolWt)
    merged_df["LogP"] = merged_df["mol"].apply(Descriptors.MolLogP)
    new_cols = ["smiles", "MW", "LogP"] + raw_cols[1:]
    return merged_df[new_cols]


def charge_filter(mol):
    negative_charge = Chem.MolFromSmarts("[*-1]")
    positive_charge = Chem.MolFromSmarts("[*+1]")
    nc = len(mol.GetSubstructMatches(negative_charge))
    pc = len(mol.GetSubstructMatches(positive_charge))
    npc = nc + pc
    if npc <= 1:
        return True
    elif npc == 2:
        if nc <= 1:
            return True
    return False


def grep_sdf(workdir, merge_file):
    merged_sdf = os.path.join(workdir, "merged_all.sdf")
    selected_sdf = os.path.join(workdir, "selected.sdf")
    # merge all sdf
    cmd_merge = ["find", workdir, "-name \"docking_outputs_with_score.sdf\" | xargs cat >", merged_sdf]
    shell_cmd_execute(cmd_merge)
    # create ids
    df = pd.read_csv(merge_file)
    ids = list(set(df["id"].apply(lambda x: x.split("-dp")[0])))
    # subset sdf
    trak = False
    with open(merged_sdf, "r") as sdf:
        with open(selected_sdf, "w") as sel_sdf:
            for line in sdf.readlines():
                if trak is False:
                    if line.strip() in ids:
                        sel_sdf.write(line)
                        trak = True
                elif trak:
                    sel_sdf.write(line)
                    if "$$$$" in line:
                        trak = False


def write_growth(max_gen: int, workdir: str, dl_mode: int, config_path: str):
    now = str(int(time.time()))
    file_path = os.path.join(workdir, "merged_docked_best_" + now + ".csv")
    merge_multi_generation(workdir, max_gen, file_path, dl_mode, config_path)

    new_file = file_path.replace(".csv", "_tmp.csv")
    final_file = file_path.replace(".csv", "_with_grow_path.csv")

    mut_dic_all = cal_mutation_dic(workdir, max_gen)
    with open(file_path, 'r') as raw:
        header = raw.readline().strip().split(",")
        path_header = list(zip(["smi_gen_", "id_gen_", "rxn_gen_", "partner_gen_"] * max_gen,
                               np.repeat(list(range(max_gen)), 4).astype(str)))
        header += ["".join(i) for i in path_header]
        new_header = ",".join(header) + "\n"
        with open(new_file, "w") as new:
            new.write(new_header)
            for line in raw.readlines():
                line = line.strip().split(",")
                mol_id = line[1]
                # find grow path per line
                mut_info_lst = grow_path(mut_dic_all, mol_id)
                if mut_info_lst is None:
                    continue
                mut_info_lst.reverse()
                mut_info_lst = list(np.concatenate(mut_info_lst))
                new_line = ",".join(line + mut_info_lst)

                # fill empty columns
                cols = new_header.count(",")
                new_line += "," * (cols - new_line.count(",")) + "\n"
                new.write(new_line)

    grow_df = add_prop(new_file)
    grow_df.to_csv(final_file, index=False)
    grep_sdf(workdir, final_file)
    print("\n", "*" * 100)
    print("Output file: ", final_file)
    print("*" * 100)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="SCESE -- find path")
    parser.add_argument("max_gen", help="Max number of generation.", type=int)
    parser.add_argument("workdir", help="Workdir")
    parser.add_argument("dl_mode",
                        help="Mode of deep learning modeling, 0: not use, 1: modeling per generation, 2: modeling overall after all the generation")
    parser.add_argument("config_path", help="config file path", type=str)
    args = parser.parse_args()
    write_growth(args.max_gen, args.workdir, args.dl_mode, args.config_path)
