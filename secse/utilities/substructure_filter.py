#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: substructure_filter.py 
@time: 2021/02/08/14:13
"""
import os
import pandas as pd
from rdkit import Chem

FILTER_FILE = "Structure Filter_20211015_v1.12.xls"


class StructureFilter:
    def __init__(self):
        df = pd.read_excel(os.path.join(os.getenv("SECSE"), "utilities", FILTER_FILE),
                           usecols=["Pattern", "ID", "Max"]).dropna()
        df = df.set_index("ID")
        df["Pattern_sma"] = df["Pattern"].apply(lambda x: Chem.MolFromSmarts(x))
        self.fdic = df[["Pattern_sma", "Max"]].T.to_dict()

    def sfilter(self, mol):
        for k, v in self.fdic.items():
            pattern = v["Pattern_sma"]
            if int(v["Max"]) == 0:
                if mol.HasSubstructMatch(pattern):
                    return k
            else:
                mts = mol.GetSubstructMatches(pattern)
                if len(mts) > int(v['Max']):
                    return k
        return "PASS"

    def sfilter_all(self, mol):
        res = []
        for k, v in self.fdic.items():
            pattern = v["Pattern_sma"]
            if int(v["Max"]) == 0:
                if mol.HasSubstructMatch(pattern):
                    res.append(k)
            else:
                mts = mol.GetSubstructMatches(pattern)
                if len(mts) > int(v['Max']):
                    res.append(k)
        if len(res) == 0:
            return "PASS"
        else:
            return res


if __name__ == '__main__':
    sf = StructureFilter()
    tmol = Chem.MolFromSmiles("CC(Cc1ncccn1)(c2ncccc2)C")
    print(sf.sfilter(tmol))
