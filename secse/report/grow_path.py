#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: grow_path.py 
@time: 2021/01/19/13:42
"""
import os
import sys
from pathlib import Path
import time
import argparse
import csv
import shutil
from tempfile import TemporaryDirectory
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt
from loguru import logger

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scoring.ranking import read_dock_file
from utilities.config import read_config
from utilities.function_helper import run_command

SELECT_SDF_SHELL = Path(__file__).with_name("filter_sdf_by_titles.pl")


def cal_mutation_dic(workdir, max_gen):
    mut_dic_all = dict()

    while max_gen > 0:
        mut_file = os.path.join(workdir, "generation_" + str(max_gen), "filter.csv")
        logger.info(mut_file)
        with open(mut_file, newline="") as source:
            rows = csv.DictReader(source)
            mut_dic = {
                row[f"id_gen_{max_gen}"].split("-dp")[0].split("-C")[0]: [
                    row[f"smiles_gen_{max_gen - 1}"],
                    row[f"id_gen_{max_gen - 1}"].split("-dp")[0].split("-C")[0],
                    row[f"reaction_id_gen_{max_gen}"],
                    row[f"partner_gen_{max_gen}"],
                ]
                for row in rows
            }
        mut_dic_all["gen" + str(max_gen)] = mut_dic
        max_gen -= 1
    return mut_dic_all


def merge_multi_generation(workdir, max_gen, file_path, dl_mode, config_path):
    df_lst = [pd.read_csv(os.path.join(workdir, "generation_" + str(i),
                                       "docked_gen_" + str(i) + ".csv"), dtype={"id": str})
              for i in range(min(1, max_gen), max_gen + 1)]
    if dl_mode == 2 and max_gen > 0:
        dl_df = read_dock_file(os.path.join(workdir, "generation_{}_pre".format(max_gen),
                                            "docking_outputs_with_score.sdf"))
        dl_df["le_ln"] = dl_df.apply(
            lambda x: x["docking score"] / Chem.MolFromSmiles(x["smiles"]).GetNumHeavyAtoms(),
            axis=1)
        dl_df.columns = [i.lower() for i in list(dl_df.columns)]
        dl_df = dl_df.drop(columns=["molecule"])
        dl_df = dl_df.reindex(columns=df_lst[0].columns)
        config = read_config(config_path)
        score_cutoff = config.getfloat("prediction", "dl_score_cutoff")
        dl_df = dl_df[dl_df["docking score"] < score_cutoff]
        df_lst.append(dl_df)

    final_df = pd.concat(df_lst, axis=0).drop_duplicates(subset=["smiles"])
    final_df.to_csv(file_path, index=False)
    return final_df


def grow_path(mut_dic_all, mut_id):
    mut_id = mut_id.split("-dp")[0].split("-C")[0]
    try:
        gen_mol = int(mut_id.split("_")[-3])
    except (IndexError, ValueError):
        logger.debug("Skipping molecule with no generation in its ID: {}", mut_id)
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
    merged_df = pd.read_csv(merged_df_path, dtype={"id": str})
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


def grep_sdf(workdir, merge_file, generations=None):
    workdir = Path(workdir)
    name = "docking_outputs_with_score.sdf"
    sources = (workdir.rglob(name) if generations is None
               else (workdir / generation / name for generation in generations))
    df = pd.read_csv(merge_file, dtype={"id": str})
    ids = dict.fromkeys(value.split("-dp")[0] for value in df["id"])
    with TemporaryDirectory(prefix=".secse-report-", dir=workdir) as temp_dir:
        temp = Path(temp_dir)
        merged_sdf = temp / "merged_all.sdf"
        with merged_sdf.open("wb") as destination:
            for path in sources:
                with path.open("rb") as source:
                    shutil.copyfileobj(source, destination)
        ids_txt = temp / "seleted_ids.txt"
        ids_txt.write_text("".join(f"{title}\n" for title in ids))
        selected_sdf = temp / "selected.sdf"
        run_command(["perl", SELECT_SDF_SHELL, merged_sdf, ids_txt, selected_sdf])
        selected_sdf.replace(workdir / selected_sdf.name)
        ids_txt.replace(workdir / ids_txt.name)


def write_growth(config_path: str, max_gen: int, dl_mode: int):
    config = read_config(config_path)
    workdir = config.get("general", "workdir")
    now = str(int(time.time()))
    file_path = os.path.join(workdir, "merged_docked_best_" + now + ".csv")
    merge_multi_generation(workdir, max_gen, file_path, dl_mode, config_path)

    new_file = file_path.replace(".csv", "_tmp.csv")
    final_file = file_path.replace(".csv", "_with_grow_path.csv")

    mut_dic_all = cal_mutation_dic(workdir, max_gen)
    with open(file_path, newline="") as raw, open(new_file, "w", newline="") as new:
        reader = csv.DictReader(raw)
        path_header = [f"{field}_gen_{gen}" for gen in range(max_gen)
                       for field in ("smi", "id", "rxn", "partner")]
        writer = csv.DictWriter(new, fieldnames=reader.fieldnames + path_header)
        writer.writeheader()
        for row in reader:
            path = grow_path(mut_dic_all, row["id"]) if max_gen else []
            if path is not None:
                fields = (value for step in reversed(path) for value in step)
                writer.writerow({**row, **dict(zip(path_header, fields))})

    grow_df = add_prop(new_file)
    grow_df.to_csv(final_file, index=False)
    generations = [f"generation_{gen}" for gen in range(max_gen + 1)]
    if dl_mode == 2 and max_gen > 0:
        generations.append(f"generation_{max_gen}_pre")
    grep_sdf(workdir, final_file, generations)
    # logger.info("\n", "*" * 100)
    logger.info(f"Output file: {final_file}")
    # logger.info("*" * 100)

    # remove temporary files
    os.remove(file_path)
    os.remove(new_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="SCESE -- find path")
    parser.add_argument("config_path", help="config file path", type=str)
    parser.add_argument("max_gen", help="Max number of generation.", type=int)
    parser.add_argument("dl_mode",
                        type=int, choices=(0, 1, 2), help="Deep learning prediction mode")
    args = parser.parse_args()
    write_growth(args.config_path, args.max_gen, args.dl_mode)
