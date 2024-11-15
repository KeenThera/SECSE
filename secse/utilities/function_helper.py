#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: function_helper.py
@time: 2022/10/13/16:36
"""
import subprocess
from loguru import logger


def shell_cmd_execute(cmd_lst, capture_mode="all"):
    cmd = " ".join(cmd_lst)
    logger.info(f"Executing command:\n{cmd}")

    try:
        # Set subprocess options based on the capture_mode
        if capture_mode == "all":
            result = subprocess.run(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, shell=True, check=True
            )
            if len(result.stdout) > 0:
                logger.info("Command output:\n" + result.stdout)
            return result.stdout

        elif capture_mode == "error":
            result = subprocess.run(
                cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True, shell=True, check=True
            )
            logger.error("Captured stderr:\n" + result.stderr)
            return result.stderr

        elif capture_mode == 0:
            subprocess.run(cmd, shell=True, check=True)
            return None

        else:
            raise ValueError("Invalid capture_mode. Use 'all', 'error', or 0.")

    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with return code {e.returncode}.")
        if capture_mode in {"all", "error"}:
            logger.error("Captured error:\n" + e.output if e.output else "No error captured.")
        raise Exception(f"Error executing command: {cmd}") from e
