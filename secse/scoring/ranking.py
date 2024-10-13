#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: ranking.py 
@time: 2020/11/04/13:35
"""
import pandas as pd
from rdkit.Chem import PandasTools
from scoring.diversity_score import *
import numpy as np
import os
import configparser
from pandarallel import pandarallel

pandarallel.initialize(verbose=0)
rdkit.RDLogger.DisableLog("rdApp.*")


def read_dock_file(sdf):
    # assign new id for duplicates, with suffix -1, -2, ...
    sdf_df = PandasTools.LoadSDF(sdf, smilesName='smiles', molColName='Molecule')[
        ["ID", "Molecule", "smiles", "docking score"]]
    sdf_df["docking score"] = sdf_df["docking score"].astype(float)
    sdf_df = sdf_df.sort_values(by="docking score", ascending=True)
    name_groups = sdf_df.groupby("ID")["ID"]
    suffix = name_groups.cumcount() + 1
    repeats = name_groups.transform("size")
    sdf_df["ID"] = np.where(repeats > 1, sdf_df['ID'] + "-dp" + suffix.map(str), sdf_df["ID"])
    return sdf_df


def clean_id(raw_id, gen):
    new_id = raw_id
    # if "GEN_" + str(gen) in raw_id:
    if "-C" in new_id:
        new_id = new_id.rsplit("-C", 1)[0]
    # elif new_id.count("dp") > 1:
    #     new_id = new_id.rsplit("-dp", 1)[0]
    # elif new_id.count("-C") > 1:
    #     new_id = new_id.rsplit("-C", 1)[0]
    return new_id


class Ranking(object):
    def __init__(self, sdf, gen, config_file):
        self.sdf = sdf
        self.gen = gen

        config = configparser.ConfigParser()
        config.read(config_file)
        self.docking_score_cutoff = config.getfloat("docking", "score_cutoff")
        self.RMSD = config.getfloat("docking", "rmsd")
        self.delta_docking_score = config.getfloat("docking", "delta_score")

        self.docked_df = pd.DataFrame(None)
        self.diff = None
        self.score_min = None
        self.winner = None
        self.final_df = None
        self.keep_mols = None

        self.load_sdf()
        self.ranking_flag = True
        if self.gen > 0:
            if self.filter_rmsd_docking_score():
                self.cal_le_rank()
            else:
                self.ranking_flag = False
                print("No molecule left, stopping generation.")
        elif self.gen == 0:
            self.cal_le_rank()

        self.size = min(config.getint("general", "seed_per_gen"), self.docked_df.shape[0])

    def load_sdf(self):
        raw_df = PandasTools.LoadSDF(self.sdf, smilesName='smiles', molColName='Molecule')[
            ["ID", "Molecule", "smiles", "docking score"]]
        raw_df["docking score"] = raw_df["docking score"].astype(float)
        raw_df = raw_df.sort_values(by="docking score", ascending=True)

        raw_df.columns = [i.lower() for i in list(raw_df.columns)]

        self.docked_df = raw_df[["smiles", "id", "docking score", "molecule"]].copy()
        # assign new id for duplicates, with suffix -1, -2, ...
        name_groups = self.docked_df.groupby("id")["id"]
        suffix = name_groups.cumcount() + 1
        repeats = name_groups.transform("size")
        self.docked_df["id_raw"] = self.docked_df["id"].copy()
        self.docked_df["id"] = np.where(repeats > 1, self.docked_df['id'] + "-dp" + suffix.map(str),
                                        self.docked_df["id"])

        print("{} cmpds after evaluate".format(self.docked_df.shape[0]))

    def load_parents_sdf(self):
        gen = str(self.gen - 1)
        read_dock_file(os.path.join(os.path.dirname(os.path.dirname(self.sdf)), "generation_" + gen,
                                    "docking_outputs_with_score.sdf"))

    def mols_score_below_cutoff(self):
        self.docking_score_cutoff = min(self.docking_score_cutoff,
                                        self.docked_df["docking score"].astype(float).quantile(0.01))
        print("The evaluate score cutoff is: {}".format(self.docking_score_cutoff))
        self.keep_mols = self.docked_df[self.docked_df["docking score"].astype(float) <= self.docking_score_cutoff]
        self.final_df = pd.concat([self.keep_mols, self.winner]).drop_duplicates(subset="id")
        cols = list(self.final_df.columns)
        cols = [i + "_gen_" + str(self.gen) for i in cols]
        self.final_df.columns = cols
        print("{} final seeds.".format(self.final_df.shape[0]))

    def filter_rmsd_docking_score(self):
        last_sdf = self.sdf.replace("generation_" + str(self.gen), "generation_" + str(self.gen - 1))
        last_df = read_dock_file(last_sdf).set_index("ID")
        mut_df = pd.read_csv(os.path.join(os.path.dirname(self.sdf), "filter.csv"), low_memory=False)
        parent_dic = dict(zip(mut_df["id_gen_" + str(self.gen)], zip(mut_df["id_gen_" + str(self.gen - 1)],
                                                                     mut_df["type"])))
        self.docked_df["id_find_parent"] = self.docked_df["id_raw"].apply(lambda x: clean_id(x, self.gen))

        # calculate RMSD: parent from last generation
        def cal_rmsd_docked(row):
            # do not care rmsd for the first generation
            if self.gen == 1:
                return -1
            # do not care rmsd except for Grow type
            if "G" not in parent_dic[row["id_find_parent"]][1]:
                return -2
            return cal_rmsd(last_df.loc[parent_dic[row["id_find_parent"]][0]]["Molecule"], row["molecule"])

        # calculate RMSD only for Type Grow mutation, assign -1 for other mutation
        self.docked_df["rmsd"] = self.docked_df.apply(cal_rmsd_docked, axis=1)
        # calculate change of evaluate score after growing
        self.docked_df["delta_docking_score"] = self.docked_df.apply(lambda x: float(x["docking score"]) - float(
            last_df.loc[parent_dic[x["id_find_parent"]][0]]["docking score"]), axis=1)

        # keep same binding mode (RMSD < 2A and delta evaluate score < -0.3) or
        # find a better binding mode (delta evaluate score < -1.2kcal )
        print("{} cmpds before RMSD/Docking Score filter".format(self.docked_df.shape[0]))
        self.docked_df = self.docked_df[(self.docked_df["delta_docking_score"] <= self.delta_docking_score) | (
                (self.docked_df["rmsd"] <= self.RMSD) & (self.docked_df["delta_docking_score"] <= -0.2))]
        rest_cmpds = self.docked_df.shape[0]
        print("{} cmpds after RMSD/Docking Score filter".format(rest_cmpds))
        if rest_cmpds == 0:
            return False
        return True

    def cal_le_rank(self):
        # calculate ln LE and fitness rank
        self.docked_df["le_ln"] = self.docked_df.apply(
            lambda x: x["docking score"] / (1 + np.log(x["molecule"].GetNumHeavyAtoms())),
            axis=1)
        self.diff = self.docked_df["le_ln"].max() - self.docked_df["le_ln"].min()
        self.score_min = self.docked_df["le_ln"].min()
        self.docked_df["fitness"] = 1 - ((self.docked_df["le_ln"] - self.score_min) / self.diff)
        self.docked_df["fitness"] = self.docked_df["fitness"].fillna(-1)
        self.docked_df["fitness_rank"] = self.docked_df["fitness"].rank(ascending=False)
        self.docked_df["fitness_rank"] = self.docked_df["fitness_rank"].fillna(-1)
        # drop molecule columns
        self.docked_df = self.docked_df.drop(columns=["molecule", "id_raw"])

    def roulette_selection(self):
        self.winner = self.docked_df.sample(n=self.size, weights="fitness")

    def tournament_selection(self):
        # random sample 3 molecules the one with smallest evaluate score win, repeat until get 20% of original data
        win_lst = []
        if self.size == 1:
            self.winner = self.docked_df.copy()
        pool = self.docked_df.copy()
        for i in range(int(self.size)):
            winner = pool.sample(min(10, pool.shape[0])).nsmallest(1, "le_ln", keep="first")
            win_lst.append(winner)
            pool = pool.drop(winner.index)

        self.winner = pd.concat(win_lst)
