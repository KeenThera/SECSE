import os
import subprocess
import pandas as pd
# os.environ["BABEL_DATADIR"] = r"C:\Users\shien\miniconda3\share\openbabel"
from plip.basic import config
from plip.structure.preparation import PDBComplex
from openbabel import pybel


def checkFilesExists(filename):
    if os.path.isfile(filename):
        return True
    else:
        return False


def complex_analysis(pdb_file: str):
    config.NOHYDRO = True
    pdb_complex = PDBComplex()
    pdb_complex.load_pdb(pdb_file)
    pdb_complex.analyze()
    my_lig_set = {}
    for ligand in pdb_complex.ligands:
        if ligand.hetid in my_lig_set.keys():
            my_lig_set[ligand.hetid] += ligand.members
        else:
            my_lig_set[ligand.hetid] = ligand.members
    return pdb_complex, my_lig_set


def rename(title):
    new = title.split("/")
    new = new[:5]
    out = "/".join(new)

    return out


def get_atom_contact_area(tablefile, skip_none_contact=False):
    with open(tablefile) as f:
        l = f.readline()
        columns = l.strip().split("\t")
    df = pd.read_csv(tablefile, sep="\t", index_col=0, skipinitialspace=True,
                     skip_blank_lines=True, usecols=columns)

    if skip_none_contact:
        columns = df.columns
        df = df.drop([col for col in columns if df[col].sum() == 0.0], axis=1)
        indexes = df.index
        df = df.drop([idx for idx in indexes if df.loc[idx].sum() == 0.0])

    # rename = lambda a: "{}/{}/{}/{}/{}".format(*a.split("/"))
    df.columns = [rename(c) for c in df.columns]
    iname = df.index.name
    df.index = [rename(i) for i in df.index]
    df.index.name = iname

    return df


def get_atom_bsa_vs_asa(tablefile, atmasafile, skip_none_contact=True):
    contacts_df = get_atom_contact_area(tablefile, skip_none_contact)
    atom_bsa_asa_a = {}
    atom_bsa_asa_b = {}

    atmasa = pd.read_csv(atmasafile, sep='\t')
    atmasa_by_id = atmasa.set_index('ID')
    atom_total_asa = {str(atnum): atmasa_by_id['total_ASA'][atnum] for atnum in atmasa_by_id.index}
    atmasa_by_resnum = atmasa.groupby(['resnum']).sum()
    res_total_asa = {str(resnum): atmasa_by_resnum['total_ASA'][resnum] for resnum in atmasa_by_resnum.index}

    atom_bsa_a = {atom.split('/')[-1]: contacts_df[atom].sum() \
                  for atom in contacts_df.columns}
    atom_bsa_b = {atom.split('/')[-1]: contacts_df.loc[atom].sum() \
                  for atom in contacts_df.index}

    atom_total_asa_a = {atom: atom_bsa_a[atom] + atom_total_asa[atom] for atom in atom_bsa_a}
    atom_total_asa_b = {atom: atom_bsa_b[atom] + atom_total_asa[atom] for atom in atom_bsa_b}

    for atom in atom_total_asa_a:
        atom_bsa_asa_a[atom] = atom_bsa_a[atom] / atom_total_asa_a[atom]
    for atom in atom_total_asa_b:
        atom_bsa_asa_b[atom] = atom_bsa_b[atom] / atom_total_asa_b[atom]

    return atom_bsa_asa_a, atom_bsa_asa_b


def get_contact_atom_bsaasa(basename, skip_none_contact=True):
    tablefile = basename + ".PROTEIN_vs_LIGAND.by_atom.tsv"
    atmasafile = basename + ".atmasa"
    atom_bsa_asa_a, atom_bsa_asa_b = get_atom_bsa_vs_asa(tablefile, atmasafile, skip_none_contact)
    df = get_atom_contact_area(tablefile, skip_none_contact)
    x = list(df.columns)
    y = list(atom_bsa_asa_a[e.split('/')[-1]] * 100 for e in df.columns)
    my_dict = {k: v for k, v in zip(x, y)}
    return my_dict


def hbond_info(pdbid):
    complex_name = pdbid + ".pdb"
    pdb_complex, ligs = complex_analysis(complex_name)
    protein_atom_index = []
    hbond_selected_list = {}
    for lig in ligs:
        tmp = list(ligs[lig][0])
        mem = ":".join([str(x) for x in tmp])
        interaction = pdb_complex.interaction_sets[mem]
        all_hbonds = interaction.hbonds_ldon + interaction.hbonds_pdon
        for hbond in all_hbonds:
            if hbond.protisdon:
                protein_atom_index.append(hbond.d_orig_idx)
                hbond_selected_list[hbond.d_orig_idx] = hbond
            else:
                protein_atom_index.append(hbond.a_orig_idx)
                hbond_selected_list[hbond.a_orig_idx] = hbond
    return protein_atom_index, hbond_selected_list


def BSA_filter(atoms, hbond_bsa):
    for info in hbond_bsa:
        atom = info.split("/")[-1]
        if int(atom) in atoms:
            atom_bsa = hbond_bsa[info]
            if atom_bsa <= 90:
                return False

    return True


def gen_complex(pdbid, protein_pdbqt_filename, ligand_pdbqt_filename):
    protein = next(pybel.readfile("pdbqt", protein_pdbqt_filename))
    ligand = next(pybel.readfile("pdb", ligand_pdbqt_filename))
    pro = protein.OBMol
    lig = ligand.OBMol
    pro += lig
    pro = pybel.Molecule(pro)
    pro.addh()
    pro.write("pdb", pdbid + ".pdb", overwrite=True)
    if checkFilesExists(pdbid + ".pdb"):
        return True
    return False


class BsaFilter:
    def __init__(self, pdbid, protein_filename, ligand_filename):
        self.pdbid = pdbid
        self.protein_filename = protein_filename
        self.ligand_filename = ligand_filename
        self.atmasa_file = None
        self.interaction_table_file = None
        if gen_complex(self.pdbid, self.protein_filename, self.ligand_filename):
            self.complex_filename = pdbid + ".pdb"

    def calculate_atom_sa(self):
        cmd_string = "dr_sasa -m 1 -i {0} protein ligand".format(self.complex_filename)
        print(cmd_string)
        try:
            p = subprocess.Popen(cmd_string, shell=True)
            p.communicate()
            return True
        except subprocess.CalledProcessError as e:
            print(f"Error: {e}")
            return False

    def check_key_file(self):
        basename = self.complex_filename.split(".")[0]
        # check if key files are available
        interaction_table_file = basename + ".PROTEIN_vs_LIGAND.by_atom.tsv"
        atmasa_file = basename + ".atmasa"
        if checkFilesExists(interaction_table_file) and checkFilesExists(atmasa_file):
            self.interaction_table_file = interaction_table_file
            self.atmasa_file = atmasa_file
            return True
        else:
            return False

    def bsa_filter(self):
        self.calculate_atom_sa()
        if self.check_key_file():
            hbond_bsa = get_contact_atom_bsaasa(self.pdbid)
            p_atom_list, b = hbond_info(self.pdbid)
            if BSA_filter(p_atom_list, hbond_bsa):
                return True
        return False


def main():
    myBSA = BsaFilter("3og7", "protein_2.pdb", "ligand.pdb")
    print(myBSA.bsa_filter())


if __name__ == "__main__":
    main()
