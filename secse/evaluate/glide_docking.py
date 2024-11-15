#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: glide_docking.py
@time: 2021/11/19/10:05
"""
import os
from loguru import logger
from utilities.function_helper import shell_cmd_execute

GLIDE_SHELL = os.path.join(os.getenv("SECSE"), "evaluate", "ligprep_glide.sh")


def dock_by_glide(workdir, mols_smi, target, gen, dock_mode, cpu_num):
    ligprep_glide = [GLIDE_SHELL, mols_smi, workdir, target, str(gen), dock_mode, str(cpu_num)]
    shell_cmd_execute(ligprep_glide)
    glide_out = os.path.join(workdir, "glide_gen_{}_lib.sdf".format(gen))
    sdf_path = os.path.join(workdir, "docking_outputs_with_score.sdf")
    write_score = False
    pass_line = 0
    with open(glide_out, "r") as glide:
        with open(sdf_path, "w") as sdf:
            for line in glide.readlines():
                if line.startswith("> <r_i_glide_gscore>"):
                    # write docking score
                    write_score = True
                    continue
                elif write_score:
                    score = line.strip()
                    newline = "> <docking score>\n{}\n".format(score)
                    write_score = False
                elif line.startswith("> <"):
                    # drop other fields
                    pass_line = 2
                    continue
                elif pass_line > 0:
                    pass_line -= 1
                    continue
                else:
                    newline = line
                sdf.write(newline)
