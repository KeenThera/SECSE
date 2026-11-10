#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: glide_docking.py
@time: 2021/11/19/10:05
"""
import re
from pathlib import Path
from tempfile import TemporaryDirectory

from evaluate.sdf_output import validate_scored_sdf
from utilities.config import SOURCE_DIR
from utilities.function_helper import run_command

GLIDE_SHELL = SOURCE_DIR / "evaluate" / "ligprep_glide.sh"
PROPERTY = re.compile(r"^>\s*<([^>]+)>")


def dock_by_glide(workdir, mols_smi, target, gen, dock_mode, cpu_num):
    run_command([GLIDE_SHELL, mols_smi, workdir, target, gen, dock_mode, cpu_num])
    workdir = Path(workdir)
    glide_out = workdir / f"glide_gen_{gen}_lib.sdf"
    output = workdir / "docking_outputs_with_score.sdf"
    with TemporaryDirectory(prefix=".secse-docking-", dir=workdir) as temp_dir:
        staged = Path(temp_dir) / output.name
        with glide_out.open() as source, staged.open("w") as destination:
            skip_property = False
            for line in source:
                field = PROPERTY.match(line)
                if field:
                    skip_property = field.group(1) != "r_i_glide_gscore"
                    if not skip_property:
                        destination.write("> <docking score>\n")
                elif not skip_property or line.strip() == "$$$$":
                    destination.write(line)
                    skip_property = False
                elif not line.strip():
                    skip_property = False
        validate_scored_sdf(staged)
        staged.replace(output)
