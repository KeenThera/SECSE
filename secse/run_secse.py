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
        num_gen = config.getint("DEFAULT", "num_gen")
        mols_smi = config.get("DEFAULT", "fragments")
        workdir = config.get("DEFAULT", "workdir")
        num_per_gen = config.getint("DEFAULT", "num_per_gen")
        start_gen = config.getint("DEFAULT", "start_gen")
        docking_program = config.get("DEFAULT", "docking_program")

        receptor = config.get("docking", "target")
        dl_mode = config.getint("deep learning", "mode")
        if "vina" in docking_program.lower():
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
        workflow = Grow(num_gen, mols_smi, workdir, num_per_gen, docking_program, receptor,
                        start_gen, dl_mode, args.config, x, y, z, box_size_x, box_size_y, box_size_z)
    else:
        workflow = Grow(num_gen, mols_smi, workdir, num_per_gen, docking_program, receptor, start_gen, dl_mode,
                        args.config)
    workflow.grow()
    write_growth(num_gen, workdir, dl_mode, args.config)


if __name__ == '__main__':
    time1 = time.time()
    print(
        "\n",
        "*" * 88, "\n",
        "      ____    _____    ____   ____    _____ \n",
        "     / ___|  | ____|  / ___| / ___|  | ____|\n",
        "     \\___ \\  |  _|   | |     \\___ \\  |  _|  \n",
        "      ___) | | |___  | |___   ___) | | |___ \n",
        "     |____/  |_____|  \\____| |____/  |_____|")

    main()
    time2 = time.time()
    print("Time consumption (total): {} hours".format(round((time2 - time1) / 3600, 2)))
    print("*" * 88)
