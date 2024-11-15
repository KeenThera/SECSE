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
from rdkit import Chem
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
    pandarallel.initialize(verbose=0, nb_workers=cpu_num)
    df["fp2"] = df[smi].parallel_apply(cal_morgan_fp)
    df = df.dropna(subset=["fp2"])
    c = df["fp2"].sample(1)
    c_next = c.index[0]
    c = c.iloc[0]
    c_lst = []
    dis = np.zeros(df.shape[0])
    dis_dic = dict()
    for i in range(k):
        new_dis = np.array(df["fp2"].apply(lambda x: tanimoto_smi(c, x)))
        dis_dic[c_next] = new_dis.copy()
        # mask mols with similarity larger than 0.6, those mols with not be consider as cluster center in next loops
        new_dis[new_dis >= 0.6] = 999999999
        dis += new_dis
        if np.min(dis) >= 999999999:
            break
        else:
            c_next = np.argmin(dis)
            c = df["fp2"].iloc[c_next]
            c_lst.append(c_next)

    df_cluster = pd.DataFrame(dis_dic)
    df["cluster_center_gen_" + str(gen)] = df_cluster.parallel_apply(lambda x: x.nlargest(1).index[0], axis=1)
    df["cluster_center_dis_gen_" + str(gen)] = df_cluster.parallel_apply(lambda x: x.nlargest(1).iloc[0], axis=1)
    df = df.drop(columns="fp2")
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
