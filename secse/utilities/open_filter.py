# we use logs for generated molecule filter

"""
The user can define their own filter function as needed.
The input parameter of the function is an rdkit mol object,
and the return value is a boolean. If the molecule is needed, return true;
if it is not needed, return false.
The user can modify this Python script file according to their own requirements.

The following code is just an example.
LogS = 0.26 -  0.74 LogP - 0.0066 MW + 0.0034 RB - 0.42 AP
ref :https://practicalcheminformatics.blogspot.com/2023/06/
getting-real-with-molecular-property.html

"""
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski
from loguru import logger


def user_filter(mol):
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    rotors = Lipinski.NumRotatableBonds(mol)
    ap = len(mol.GetSubstructMatches(Chem.MolFromSmarts("a"))) / mol.GetNumAtoms()
    intercept = 0.16
    coef = {"logp": -0.63, "mw": -0.0062, "rotors": 0.066, "ap": -0.74}
    esol = intercept + coef["logp"] * logp + coef["mw"] * mw + coef["rotors"] * rotors + coef["ap"] * ap

    if esol <= -4.5:
        return True
    else:
        return False
