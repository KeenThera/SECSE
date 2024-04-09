#!/tools/miniconda3/envs/cadd/bin/python
'''
This script will calculate ligand efficiency of a user specified property in a SDF file, bin the list by the property, select the rows with highest LE and save into a new SDF file
## Author: zhentgpicasa@gmail.com
## Revision History:
- 2024/3/17
  - Fist version
'''

import argparse
import os.path
import pandas as pd
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
import os
import time
import pandas as pd
from rdkit.Chem import PandasTools
import numpy as np


startTime = time.time()
'''
if RDkit is installed as an virtual environment other than 'base' environment, within Jupyter some paths are not included, and thus rdkit will not be imported, so the missing path need to be added manually
Zhenting has tracked this bug on 2/21/2020.
'''
pythonPath = os.__file__.split("lib")[0]
os.environ['PATH'] = os.environ[
    'PATH'] + os.pathsep + pythonPath + r'Library\bin' + os.pathsep

# import click

parser = argparse.ArgumentParser(
    description='Calculate ligand efficiency of a user specified property in a SDF file, bin the list by the property, select the rows with highest LE and save into a new SDF file')
parser.add_argument('-i', required=True, help='SDF input file')
parser.add_argument('-o', required=True, help='SDF output file')
parser.add_argument('-p', required=False,
                    help='Property for ligand efficiency calculation', default='docking score')
parser.add_argument('-d', default='ID', required=False,
                    help='Molecule ID column name')
parser.add_argument('-b', type=int, default=100,
                    required=False, help='Bin count')
args = parser.parse_args()

prop4LE = args.p
idCol = args.d
sdfFile = args.i
outputSdfFile = args.o
binCount = args.b


def workflow():

    # Read the SDF file into a DataFrame
    df = PandasTools.LoadSDF(sdfFile, removeHs=False)
    if idCol not in df.columns:
        print('Molecule ID column name is not detected', df.columns)
        quit()
    if prop4LE not in df.columns:
        print('Column name of the property for ligand efficiency calculation is not detected', df.columns)
        quit()
    try:  # Set data type to float
        df[prop4LE] = df[prop4LE].astype(float)
    except:
        '''Do nothing'''

    # Calculate the heavy atom count for each molecule
    df['HeavyAtomCount'] = df['ROMol'].apply(lambda x: x.GetNumHeavyAtoms())
    df['LE'] = df[prop4LE]/df['HeavyAtomCount']

    # Sort by prop4LE and remove duplicated rows
    df.sort_values([prop4LE], inplace=True, ascending=[True])
    df.drop_duplicates([idCol], inplace=True)

    # Calculate the range and step size
    min_value = df[prop4LE].min()
    max_value = df[prop4LE].max()
    range_of_values = max_value - min_value
    step_size = range_of_values / binCount

    # Create an array of bins
    bins = list(np.arange(min_value, max_value + step_size, step_size))

    # Bin the data into intervals
    df['bin'] = pd.cut(df[prop4LE], bins=bins, right=False)

    # Sort by LE
    df.sort_values('LE', inplace=True)
    # Keep the row with minimum LE in each bin
    resultDf = df.drop_duplicates(['bin']).copy()
    # Sort the result dataframe by prop4LE
    resultDf.sort_values([prop4LE], inplace=True)
    # Write the output SDF file
    PandasTools.WriteSDF(resultDf, outputSdfFile,
                         molColName='ROMol', properties=list(resultDf))

    print('The script took {:.2f} second!'.format(time.time() - startTime))


if __name__ == '__main__':
    workflow()
