#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: sampling.py
@time: 2022/2/8/10:25
"""
import os
import pandas as pd
from loguru import logger

from scoring.diversity_score import cal_morgan_fp, tanimoto_smi


def sample_by_rule_weight(gen, filter_df, workdir_now):
    if "G-002" in list(filter_df["type"]):
        # control ratio of G-002 mutation
        spacer_df = filter_df[filter_df["type"] == "G-002"]

        common_df = filter_df.drop(spacer_df.index, axis=0)
        # control ratio of ring with spacer based on different stage
        if gen <= 3:
            spacer_ratio = 0.3
        elif gen <= 7:
            spacer_ratio = 0.1
        else:
            spacer_ratio = 0.01
        sample_size = min(filter_df.shape[0], 500000)

        spacer_df = spacer_df.sample(min(int(sample_size * spacer_ratio), spacer_df.shape[0]),
                                     replace=False,
                                     weights="priority_gen_" + str(gen))

        common_df = common_df.sample(min(int(sample_size * (1 - spacer_ratio)), common_df.shape[0]),
                                     replace=False,
                                     weights="priority_gen_" + str(gen))
        sampled_df = pd.concat([spacer_df, common_df], axis=0)
        sampled_df.to_csv(os.path.join(workdir_now, "sampled.csv"), index=False)
    else:
        logger.error("No cmpds generated from ring with spacer in the generation!")
        sampled_df = filter_df.sample(min(filter_df.shape[0], 500000), replace=False,
                                      weights="priority_gen_" + str(gen))
        sampled_df.to_csv(os.path.join(workdir_now, "sampled.csv"), index=False)

    return sampled_df


def sample_by_similarity(gen, filter_df, workdir_now, num_per_gen,
                         ref_smi="O=C(C1=CC=C(C(C)NC(C2=CC(C3=CC=CC=C3)=NN2C)=O)C=C1)O"):
    ref_fp = cal_morgan_fp(ref_smi)
    filter_df["similarity"] = filter_df["smiles_gen_" + str(gen)].apply(
        lambda x: tanimoto_smi(cal_morgan_fp(x), ref_fp))
    sampled_df = filter_df.nlargest(num_per_gen, columns="similarity")
    sampled_df.to_csv(os.path.join(workdir_now, "sampled.csv"), index=False)
    return sampled_df
