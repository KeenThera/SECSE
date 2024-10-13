#!/usr/bin/env python
# -*- coding:utf-8 _*-
"""
@author: Lu Chong
@file: run_secse.py
@time: 2020/11/02/13:35
"""
import argparse
import time
import configparser

from grow_processes import Grow
from report.grow_path import write_growth


def main():
    parser = argparse.ArgumentParser(description="SECSE")

    parser.add_argument("--config", help="path of config file", default=False)
    args = parser.parse_args()

    try:
        # config file given
        config = configparser.ConfigParser()
        config.read(args.config)
        num_gen = config.getint("general", "num_gen")
        mols_smi = config.get("general", "fragments")
        workdir = config.get("general", "workdir")
        num_per_gen = config.getint("general", "num_per_gen")
        start_gen = config.getint("general", "start_gen")
        docking_program = config.get("docking", "docking_program")
        cpu_num = config.getint("general", "cpu")
        gpu_num = config.getint("general", "gpu")
        rule_db = config.get("general", "rule_db")
        project_code = config.get("general", "project_code")

        receptor = config.get("docking", "target")
        dl_mode = config.getint("prediction", "mode")
        if "vina" in docking_program.lower() or "unidock" in docking_program.lower():
            x = config.getfloat("docking", "x")
            y = config.getfloat("docking", "y")
            z = config.getfloat("docking", "z")
            box_size_x = config.getfloat("docking", "box_size_x")
            box_size_y = config.getfloat("docking", "box_size_y")
            box_size_z = config.getfloat("docking", "box_size_z")

    except Exception as e:
        print(e)
        print("Please check your input arguments.")
        return None

    if "vina" in docking_program.lower():
        workflow = Grow(num_gen, mols_smi, workdir, num_per_gen, docking_program, receptor, start_gen, dl_mode,
                        args.config, cpu_num=cpu_num, rule_db=rule_db, project_code=project_code, x=x, y=y, z=z,
                        box_size_x=box_size_x, box_size_y=box_size_y, box_size_z=box_size_z)
    elif "glide" in docking_program.lower():
        workflow = Grow(num_gen, mols_smi, workdir, num_per_gen, docking_program, receptor, start_gen, dl_mode,
                        args.config, cpu_num=cpu_num, rule_db=rule_db, project_code=project_code)
    elif "autodock-gpu" in docking_program.lower():
        workflow = Grow(num_gen, mols_smi, workdir, num_per_gen, docking_program, receptor, start_gen, dl_mode,
                        args.config, cpu_num=cpu_num, gpu_num=gpu_num, rule_db=rule_db, project_code=project_code)
    elif "unidock" in docking_program.lower():
        workflow = Grow(num_gen, mols_smi, workdir, num_per_gen, docking_program, receptor, start_gen, dl_mode,
                        args.config, cpu_num=cpu_num, rule_db=rule_db, project_code=project_code, x=x, y=y, z=z,
                        box_size_x=box_size_x, box_size_y=box_size_y, box_size_z=box_size_z)
    else:
        print("Please check your input docking program argument.")
        return None
    workflow.grow()


if __name__ == '__main__':
    time1 = time.time()
    print(
        "\n",
        "*" * 88, "\n",
        "      ____    _____    ____   ____    _____ \n",
        "     / ___|  | ____|  / ___| / ___|  | ____|\n",
        "     \\___ \\  |  _|   | |     \\___ \\  |  _|  \n",
        "      ___) | | |___  | |___   ___) | | |___ \n",
        "     |____/  |_____|  \\____| |____/  |_____| v1.3")

    try:
        main()
    except SystemExit as err:
        print(err)

    time2 = time.time()
    print("Time consumption (total): {} hours".format(round((time2 - time1) / 3600, 2)))
    print("*" * 88)
