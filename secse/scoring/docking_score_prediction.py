#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: docking_score_prediction.py
@time: 2021/10/27/14:26
"""
import argparse

from openbabel import openbabel
import pandas as pd
import os
import rdkit
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem.MolStandardize import rdMolStandardize
from tqdm import tqdm
from utilities.function_helper import shell_cmd_execute

rdkit.RDLogger.DisableLog("rdApp.*")


def get_train(sdf, dock):
    g = PandasTools.LoadSDF(sdf, molColName='Molecule')

    g_smi = pd.read_csv(dock, sep="\t", header=None)
    g_smi.columns = ["Smiles", "ID"]
    g_smi = g_smi.drop_duplicates(subset="ID")
    g_smi = g_smi.set_index("ID")

    g = g[["ID", "Molecule", "docking score"]]
    g["docking score"] = g["docking score"].astype(float)
    g = g.sort_values("docking score", ascending=True)

    g["Smiles"] = g["ID"].apply(lambda x: g_smi.loc[x.rsplit("-C", 1)[0]][0])
    g_new = g.sort_values(by="docking score", ascending=True).drop_duplicates(subset="Smiles", keep="first")

    smi = g_new["Smiles"].apply(lambda x: neutralize(x))
    g_new["Smiles"] = smi
    g_new = g_new.drop_duplicates(subset="Smiles", keep="first")
    return g_new


def get_pre(workdir, max_gen, get_all=False):
    pre_dir = os.path.join(workdir, "prediction")
    if get_all:
        pre_raw = os.path.join(pre_dir, "all_G" + str(max_gen) + "_for_pre.raw")
        pre_file = os.path.join(pre_dir, "all_G" + str(max_gen) + "_for_pre.csv")

        cmd_cat = ["find", workdir, "-name \"filter.csv\" |xargs awk -F, 'NR>1{{print $(NF-5)\",\"$(NF-6)}}' >",
                   pre_raw]
        shell_cmd_execute(cmd_cat)
        cmd_dedup = ["awk -F',' '!seen[$2]++'", pre_raw, ">", pre_file]
        shell_cmd_execute(cmd_dedup)

        drop_mols = os.path.join(pre_dir, "drop_ids.txt")
        mols_id_cat = ["find", workdir, "-name \"mols_for_docking.smi\" |xargs cut -f2  >", drop_mols]
        shell_cmd_execute(mols_id_cat)
        final_file = os.path.join(pre_dir, "all_G" + str(max_gen) + "_for_pre_uniq.csv")
    else:
        pre_file = os.path.join(pre_dir, "gen_" + str(max_gen) + "_for_pre.csv")
        cmd_cp = ["awk -F, 'NR>1{{print $(NF-5)\",\"$(NF-6)}}'",
                  os.path.join(workdir, "generation_" + str(max_gen), "filter.csv"), ">", pre_file]
        shell_cmd_execute(cmd_cp)

        drop_mols = os.path.join(pre_dir, "drop_ids_{}.txt".format(max_gen))
        mols_id_cat = ["cut -f2", os.path.join(workdir, "generation_" + str(max_gen), "mols_for_docking.smi"), ">",
                       drop_mols]
        shell_cmd_execute(mols_id_cat)
        final_file = os.path.join(pre_dir, "gen_" + str(max_gen) + "_for_pre_uniq.csv")

    try:
        cmd_drop = ["grep -wvf", drop_mols, pre_file, ">", final_file]
        shell_cmd_execute(cmd_drop)
    except:
        final_file = None
    return final_file


def neutralize(smi):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        smi = wash_mol(smi)
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return "C"
    uc = rdMolStandardize.Uncharger()
    return Chem.MolToSmiles(uc.uncharge(mol))


def wash_mol(smi):
    ob_conversion = openbabel.OBConversion()
    ob_conversion.SetInAndOutFormats("smi", "can")
    ob_mol = openbabel.OBMol()
    ob_conversion.ReadString(ob_mol, smi)
    ob_conversion.Convert()
    res = ob_conversion.WriteString(ob_mol).strip()
    return res


def prepare_files(max_gen, workdir, dl_mode):
    pre_dir = os.path.join(workdir, "prediction")
    os.makedirs(pre_dir, exist_ok=True)

    def pre_train_per_gen(gen):
        sdf = os.path.join(workdir, "generation_{}/docking_outputs_with_score.sdf".format(gen))
        dock = os.path.join(workdir, "generation_{}/mols_for_docking.smi".format(gen))
        df_train = get_train(sdf, dock)[['Smiles', 'docking score']]
        # write per generation
        df_train.to_csv(os.path.join(pre_dir, "train_G{}.csv".format(gen)), index=False)
        return df_train

    if dl_mode == 1:
        # prepare current generation data
        pre_train_per_gen(max_gen)
        train = os.path.join(pre_dir, "train_G{}.csv".format(max_gen))
        pre = get_pre(workdir, max_gen, False)
        return train, pre

    elif dl_mode == 2:
        # prepare files for all the generation and merge together
        cum_path = os.path.join(pre_dir, "train_G" + str(max_gen) + "_all.csv")
        df_lst = []
        for i in tqdm(range(1, max_gen + 1)):
            df = pre_train_per_gen(i)
            # write cumulative dataframe
            df_lst.append(df)

        df_all = pd.concat(df_lst, axis=0).sort_values(
            by="docking score", ascending=True).drop_duplicates(subset="Smiles", keep="first")
        df_all.to_csv(cum_path, index=False)
        pre = get_pre(workdir, max_gen, True)
        return cum_path, pre


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="SCESE -- Prepare Data for Deep Learning")
    parser.add_argument("max_gen", help="Max number of generation.", type=int)
    parser.add_argument("workdir", help="Workdir")
    parser.add_argument("dl_mode",
                        help="Mode of deep learning modeling, 1: modeling per generation, 2: modeling overall after all the generation",
                        type=int, default=0)
    args = parser.parse_args()
    prepare_files(args.max_gen, args.workdir, args.dl_mode)
