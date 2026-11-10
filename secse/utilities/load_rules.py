#!/usr/bin/env python  
# -*- coding:utf-8 _*-
""" 
@author: Lu Chong
@file: load_rules.py
@time: 2022/2/28/09:52
"""

import sqlite3
from pathlib import Path
from contextlib import closing
from tempfile import TemporaryDirectory
import pandas as pd


def json_to_DB(in_json, out_db_path):
    df = pd.read_json(in_json)
    required = {"Rule ID", "SMARTS", "Priority"}
    if df.empty or not required.issubset(df.columns):
        raise ValueError(f"Invalid rules JSON {in_json}: expected Rule ID, SMARTS and Priority")
    output = Path(out_db_path)
    with TemporaryDirectory(prefix=".secse-rules-", dir=output.parent) as directory:
        temporary = Path(directory) / output.name
        with closing(sqlite3.connect(temporary)) as connection:
            df.to_sql("G-001", connection)
        temporary.replace(output)
