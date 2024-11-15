import sqlite3
import pandas as pd
import json
import re
import os
from rdkit import Chem
from rdkit.Chem import Descriptors, rdChemReactions
from rdkit.Chem.rdMolDescriptors import CalcExactMolWt, CalcFractionCSP3, CalcNumRings
from rdkit.Chem import FindMolChiralCenters
from loguru import logger


def read_excel(filename, sheet_name):
    """Read an Excel sheet into a DataFrame."""
    return pd.read_excel(filename, sheet_name=sheet_name)


def write_to_json(df, filename):
    """Write a DataFrame to a JSON file."""
    df.to_json(filename, orient='records', force_ascii=False, indent=4)


def write_to_sqlite(df, table, db_path):
    """Write a DataFrame to an SQLite database."""
    with sqlite3.connect(db_path) as conn:
        df.to_sql(table, conn, if_exists='replace', index=False)


def test_rxn(sma):
    """Test if a SMARTS string can be converted to an RDKit reaction."""
    try:
        rdChemReactions.ReactionFromSmarts(sma)
    except Exception as e:
        logger.error(f"Error processing SMARTS: {sma}\n{e}")


def add_prop(df, ref_smi):
    mol = Chem.MolFromSmiles(ref_smi)
    mol_weight_ref = CalcExactMolWt(mol)
    fsp3_ref = CalcFractionCSP3(mol)
    ring_num_ref = CalcNumRings(mol)
    logp_ref = Descriptors.MolLogP(mol)
    chiral_num_ref = len(FindMolChiralCenters(mol, includeUnassigned=True))

    for index, row in df.iterrows():
        sma = row['SMARTS']
        try:
            rxn = rdChemReactions.ReactionFromSmarts(sma)
            products = rxn.RunReactants((mol,))
            new = Chem.MolFromSmiles(Chem.MolToSmiles(products[0][0]))
            mol_weight = CalcExactMolWt(new)
            df.at[index, 'ΔMW'] = mol_weight - mol_weight_ref
            fsp3 = CalcFractionCSP3(new)
            df.at[index, 'ΔFsp3'] = fsp3 - fsp3_ref
            ring_num = CalcNumRings(new)
            df.at[index, 'ΔNR'] = ring_num - ring_num_ref
            logp = Descriptors.MolLogP(mol)
            df.at[index, 'ΔlogP'] = logp - logp_ref
            chiral_num = len(FindMolChiralCenters(new, includeUnassigned=True))
            df.at[index, 'ΔNCC'] = chiral_num - chiral_num_ref
        except Exception as e:
            logger.error(e)
            logger.error(sma)
    return df


def process_sheet(sheet_df, ref_smi):
    """Process a sheet DataFrame to calculate properties and test reactions."""
    for sma in sheet_df['SMARTS']:
        test_rxn(sma)
    add_prop(sheet_df, ref_smi)
    return sheet_df


def main(excel_filename, output_type):
    """Main function to convert Excel to DB or JSON based on user input."""
    pattern = r'^[A-Za-z]-\d{3}$'
    collect_df = []
    db_name = f"{os.path.splitext(excel_filename)[0]}.db" if output_type == 'db' else None

    xls = pd.ExcelFile(excel_filename)
    for sheet_name in xls.sheet_names:
        if re.match(pattern, sheet_name):
            logger.info(f"Processing sheet: {sheet_name}")
            sheet_df = read_excel(excel_filename, sheet_name)
            # sheet_df = process_sheet(sheet_df, ref_smi='YourReferenceSMILES')

            if output_type == 'db':
                write_to_sqlite(sheet_df, sheet_name, db_name)
            else:
                if sheet_name == "G-002":
                    new_df = sheet_df[['Rule ID', 'SMARTS', 'Spacer Priority', 'Ring Priority']]
                else:
                    new_df = sheet_df[['Rule ID', 'SMARTS', 'Priority']]  # Adjust columns as needed
                collect_df.append(new_df)
            logger.info(f"Finished processing {output_type} for {sheet_name}")

    if output_type == 'json':
        combined_df = pd.concat(collect_df, ignore_index=True)
        output_filename = f"{os.path.splitext(excel_filename)[0]}.json"
        write_to_json(combined_df, output_filename)
        logger.info("All sheets processed and JSON file created.")


if __name__ == "__main__":
    excel_filename = input("Enter the Excel file name: ")
    output_type = input("Choose the output type (db for database, json for JSON file): ")
    if output_type not in ['db', 'json']:
        logger.error("Invalid output type. Please choose 'db' for database or 'json' for JSON file.")
    else:
        main(excel_filename, output_type)
