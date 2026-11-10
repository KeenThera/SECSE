#!/usr/bin/env python
# -*- coding:utf-8 _*-
"""
@author: Lu Chong
@file: run_secse.py
@time: 2020/11/02/13:35
"""
import argparse
import time
from loguru import logger
from pathlib import Path

if __package__:
    # The legacy standalone scripts import modules from the source directory.
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))

from utilities.config import docking_backend, read_config


def setup_logger(project_code, work_directory):
    work_directory = Path(work_directory)
    log_file_path = work_directory / f'{project_code}.log'
    error_file_path = work_directory / f'{project_code}_error.log'

    logger.add(log_file_path, rotation="10 MB", backtrace=True, diagnose=False, level="INFO", mode='a')
    logger.add(error_file_path, rotation="5 MB", backtrace=True, diagnose=False, level="ERROR", mode='a')
    return logger


def main(argv=None):
    parser = argparse.ArgumentParser(description="SECSE")
    parser.add_argument("--config", help="path of config file", required=True)
    args = parser.parse_args(argv)
    started = time.monotonic()

    try:
        config_path = Path(args.config).resolve()
        config = read_config(config_path)
        general = config["general"]
        docking = config["docking"]
        prediction = config["prediction"]
        backend = docking_backend(docking.get("docking_program", "vina"))
        num_gen = int(general["num_gen"])
        start_gen = general.getint("start_gen", 0)
        dl_mode = prediction.getint("mode", 0)
        num_per_gen = int(general["num_per_gen"])
        cpu_num = int(general["cpu"])
        gpu_num = general.getint("gpu", 1)
        if num_gen < 0 or start_gen < 0 or num_per_gen < 1 or cpu_num < 0:
            raise ValueError("Check num_gen, start_gen, num_per_gen and cpu in [general]")
        if general.getint("seed_per_gen", 1000) < 1 or dl_mode not in (0, 1, 2):
            raise ValueError("seed_per_gen must be positive and prediction mode must be 0, 1 or 2")

        workdir = str(Path(general["workdir"]).resolve())
        project_code = general["project_code"]
        fragments = general.get("fragments", "")
        receptor = str(Path(docking["target"]).resolve())
        if not Path(receptor).is_file():
            raise FileNotFoundError(f"Docking target: {receptor}")
        if start_gen == 0 and (not fragments or not Path(fragments).is_file()):
            raise FileNotFoundError(f"Input fragments: {fragments!r}")
        box = {}
        if backend in ("vina", "unidock"):
            box = {axis: float(docking[axis]) for axis in ("x", "y", "z")}
            box.update({f"box_size_{axis}": docking.getfloat(f"box_size_{axis}", 20) for axis in ("x", "y", "z")})

        from grow_processes import Grow

        workflow = Grow(
            num_gen, str(Path(fragments).resolve()) if fragments else "", workdir,
            num_per_gen, backend, receptor, start_gen, dl_mode, str(config_path),
            cpu_num=cpu_num, gpu_num=gpu_num, rule_db=general.get("rule_db", "0"),
            project_code=project_code, **box,
        )
        setup_logger(project_code, workdir)
        workflow.grow()
    except KeyboardInterrupt:
        logger.info("Program interrupted by user")
        return 130
    except Exception as error:
        logger.opt(exception=True).error("SECSE failed: {}", error)
        return 1
    finally:
        logger.info("Time consumption (total): {:.2f} hours", (time.monotonic() - started) / 3600)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
