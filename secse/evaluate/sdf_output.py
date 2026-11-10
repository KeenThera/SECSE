import math
import os

from rdkit import Chem


def validate_scored_sdf(path):
    count = 0
    with open(path, "rb") as source:
        size = source.seek(0, os.SEEK_END)
        source.seek(-min(size, 4096), os.SEEK_END)
        if not source.read().rstrip().endswith(b"$$$$"):
            raise ValueError(f"{path}: missing final SDF record delimiter")
        source.seek(0)
        for mol in Chem.ForwardSDMolSupplier(source):
            count += 1
            if mol is None or not mol.HasProp("docking score"):
                raise ValueError(f"{path}: invalid or unscored SDF record {count}")
            try:
                score = float(mol.GetProp("docking score"))
            except ValueError as error:
                raise ValueError(f"{path}: invalid docking score on record {count}") from error
            if not math.isfinite(score):
                raise ValueError(f"{path}: non-finite docking score on record {count}")
    if not count:
        raise ValueError(f"{path}: no scored docking poses were produced")
