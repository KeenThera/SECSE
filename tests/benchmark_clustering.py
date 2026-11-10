"""Measure the assignment step with identical fingerprints and a fixed first center."""

import argparse
import hashlib
from pathlib import Path
import random
import subprocess
import sys
from time import perf_counter
from unittest.mock import patch

import numpy as np
import pandas as pd
from rdkit import DataStructs

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "secse"))

from scoring import diversity_score


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-ref", help="Git revision containing the original clustering implementation")
    args = parser.parse_args()
    implementation = vars(diversity_score)
    if args.baseline_ref:
        source = subprocess.check_output(
            ["git", "show", f"{args.baseline_ref}:secse/scoring/diversity_score.py"], cwd=ROOT, text=True,
        )
        implementation = {}
        exec(source, implementation)

    randomizer = random.Random(67)
    fingerprints = []
    for _ in range(30_000):
        vector = DataStructs.ExplicitBitVect(512)
        vector.SetBitsFromList(randomizer.sample(range(512), 30))
        fingerprints.append(vector)
    frame = pd.DataFrame({"fp": range(len(fingerprints))})
    np.random.seed(42)
    with patch.dict(implementation, {"cal_morgan_fp": fingerprints.__getitem__}), \
            patch.object(implementation["pandarallel"], "initialize"), \
            patch.object(pd.Series, "parallel_apply", pd.Series.apply, create=True), \
            patch.object(pd.DataFrame, "parallel_apply", pd.DataFrame.apply, create=True):
        started = perf_counter()
        result = implementation["clustering"](frame, "fp", 1, 1, k=60)
        elapsed = perf_counter() - started
    digest = hashlib.sha256(result.to_csv(index=False).encode()).hexdigest()
    print(f"seconds={elapsed:.3f} centers={result['cluster_center_gen_1'].nunique()} sha256={digest}")


if __name__ == "__main__":
    main()
