#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: proprep.py
@time: 2021/9/10/14:14

prepare the protein file (pdbqt format)
"""
import os
import subprocess
from loguru import logger
from biopandas.pdb import PandasPdb


def clean(code, chain=None):
    ppdb = PandasPdb().fetch_pdb(code)
    if chain is not None:
        ppdb.df['ATOM'] = ppdb.df['ATOM'][ppdb.df['ATOM'].chain_id == chain]
    name = code + "_clean.pdb"
    ppdb.to_pdb(path=name,
                records=['ATOM', 'OTHERS'],
                gz=False,
                append_newline=True)

    ADFRsuit = r"C:\Program Files (x86)\ADFRsuite-1.0\bin"
    prepare_ligand = "prepare_receptor.bat"
    exe = os.path.join(ADFRsuit, prepare_ligand)
    p = subprocess.Popen([exe, "-r", name,
                          "-A", "hydrogens", '-w'], stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdout_data, stderr_data) = p.communicate()


def boxinfo(code, resn, extend=6):
    ppdb = PandasPdb().fetch_pdb(code)
    df_het = ppdb.df['HETATM'][ppdb.df['HETATM'].residue_name == resn]
    x_center = df_het.x_coord.mean()
    y_center = df_het.y_coord.mean()
    z_center = df_het.z_coord.mean()

    x_max = df_het.x_coord.max() + extend
    x_min = df_het.x_coord.min() - extend
    y_max = df_het.y_coord.max() + extend
    y_min = df_het.y_coord.min() - extend
    z_max = df_het.z_coord.max() + extend
    z_min = df_het.z_coord.min() - extend

    x_size = x_max - x_min
    y_size = y_max - y_min
    z_size = z_max - z_min

    return x_center, y_center, z_center, x_size, y_size, z_size
