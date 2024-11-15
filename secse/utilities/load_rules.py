#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: load_rules.py
@time: 2022/2/28/09:52
"""

import sqlite3
import pandas as pd
from loguru import logger


def json_to_DB(in_json, out_db_path):
    df = pd.read_json(in_json)
    conn = sqlite3.connect(out_db_path)
    try:
        df.to_sql("G-001", conn)
    except Exception as e:
        logger.error(e)
    conn.close()
