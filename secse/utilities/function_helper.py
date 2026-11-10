#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: function_helper.py
@time: 2022/10/13/16:36
"""
import subprocess
import shlex
import os
import sys
from loguru import logger


def run_command(args, capture_mode="all", *, cwd=None, env=None):
    args = [str(arg) for arg in args]
    logger.info("Executing command: {}", shlex.join(args))
    capture = {
        "all": (subprocess.PIPE, subprocess.STDOUT),
        "error": (subprocess.DEVNULL, subprocess.PIPE),
        0: (None, None),
    }
    stdout, stderr = capture[capture_mode]
    child_env = dict(os.environ if env is None else env)
    child_env.setdefault("SECSE_PYTHON", sys.executable)
    result = subprocess.run(args, cwd=cwd, env=child_env, stdout=stdout, stderr=stderr, text=True)
    output = result.stdout if capture_mode == "all" else result.stderr
    if result.returncode:
        detail = f"\n{output.rstrip()}" if output else ""
        raise RuntimeError(f"Command exited with status {result.returncode}: {shlex.join(args)}{detail}")
    if output:
        logger.info("Command output:\n{}", output.rstrip())
    return output


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
