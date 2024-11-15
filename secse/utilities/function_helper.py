#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: function_helper.py
@time: 2022/10/13/16:36
"""
import os
import subprocess
import sys
from loguru import logger


def shell_cmd_execute(cmd_lst):
    cmd = " ".join(cmd_lst)
    logger.info(cmd)
    try:
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        logger.error(e.output.decode())
        raise Exception("Error executing command: {}".format(cmd))
