#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: diversity_score.py 
@time: 2020/11/18/9:47
"""
import math
import numpy as np
import pandas as pd
import rdkit
from rdkit.Chem import AllChem, rdFMCS, rdShapeHelpers
from rdkit import Chem, DataStructs
from pandarallel import pandarallel
from loguru import logger

def cal_morgan_fp(smi):
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        mol = Chem.MolFromSmiles("C")
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 512)
    return fp


def tanimoto_smi(fp1, fp2):
    return rdkit.DataStructs.cDataStructs.TanimotoSimilarity(fp1, fp2)


def tanimoto_shape(ref, mol):
    return 1 - rdShapeHelpers.ShapeTanimotoDist(ref, mol)


def protrude_shape(ref, mol):
    return 1 - rdShapeHelpers.ShapeProtrudeDist(ref, mol)


def clustering(df: pd.DataFrame, smi, gen, cpu_num, k=500):
    df = df.reset_index(drop=True)
    center_column = f"cluster_center_gen_{gen}"
    similarity_column = f"cluster_center_dis_gen_{gen}"
    if df.empty:
        return df.assign(**{center_column: pd.Series(dtype=int), similarity_column: pd.Series(dtype=float)})
    if k < 1:
        raise ValueError("At least one cluster center is required")
    pandarallel.initialize(verbose=0, nb_workers=cpu_num)
    fingerprints = df[smi].parallel_apply(cal_morgan_fp).tolist()
    center = df.sample(1).index[0]
    cumulative_similarity = np.zeros(len(df))
    best_similarity = np.full(len(df), -np.inf)
    best_center = np.empty(len(df), dtype=int)
    excluded = 999999999
    for _ in range(k):
        similarity = np.asarray(DataStructs.BulkTanimotoSimilarity(fingerprints[center], fingerprints))
        # A tie stays with the first center, matching the previous nlargest(1) assignment.
        improved = similarity > best_similarity
        best_center[improved] = center
        best_similarity[improved] = similarity[improved]
        similarity[similarity >= 0.6] = excluded
        cumulative_similarity += similarity
        center = np.argmin(cumulative_similarity)
        if cumulative_similarity[center] >= excluded:
            break
    df[center_column] = best_center
    df[similarity_column] = best_similarity
    return df


def cal_rmsd(parent, c):
    mcs = rdFMCS.FindMCS([parent, c], threshold=1, completeRingsOnly=True, ringMatchesRingOnly=True,
                         bondCompare=rdFMCS.BondCompare.CompareOrderExact,

                         timeout=1).queryMol
    if mcs is None:  # no common substructure
        return -2
    p_match = parent.GetSubstructMatch(mcs)
    c_match = c.GetSubstructMatch(mcs)

    delta2 = 0.0
    for pi, ci in zip(p_match, c_match):
        d = (parent.GetConformer().GetAtomPosition(pi) - c.GetConformer().GetAtomPosition(ci)).LengthSq()
        delta2 += d
    return math.sqrt(delta2 / len(p_match))
