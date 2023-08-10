from plip.basic import config
from plip.structure.preparation import PDBComplex
from openbabel import pybel
import os

#os.environ['BABEL_DATADIR'] = r"C:\Users\shien\miniconda3\share\openbabel"


def complex_analysis(pdb_file):
    config.NOHYDRO = True
    pdb_complex = PDBComplex()
    pdb_complex.load_pdb(pdb_file)
    pdb_complex.analyze()
    my_lig_set = {}
    for ligand in pdb_complex.ligands:
        my_lig_set[ligand.hetid] = ligand.members
    return pdb_complex, my_lig_set


class HbondFilter:
    def __init__(self, rawid, ligand_filename, protein_filename):
        self.rawid = rawid
        self.ligand_filename = ligand_filename
        self.protein_filename = protein_filename
        self.ligand = next(pybel.readfile("pdb", self.ligand_filename))
        self.protein = next(pybel.readfile("pdbqt", self.protein_filename))
        self.protein.OBMol += self.ligand.OBMol
        self.complex_filename = self.rawid + "_complex.pdb"

    def generateComplexPDBFile(self):
        myComplex = pybel.Molecule(self.protein.OBMol)
        myComplex.write("pdb", self.complex_filename, overwrite=True)

    def complexAnalysis(self) -> list:
        self.generateComplexPDBFile()
        pdbcomplex, ligset = complex_analysis(self.complex_filename)
        lig_qualified_list = []
        hbond_selected_list = []
        for lig in ligset:
            hbonds_qualified_list = []
            tmp = list(ligset[lig][0])
            mem = ":".join([str(x) for x in tmp])
            interaction = pdbcomplex.interaction_sets[mem]
            all_hbonds = interaction.hbonds_ldon + interaction.hbonds_pdon
            for hb in all_hbonds:
                if not hb.sidechain and hb.type == "strong":
                    hbonds_qualified_list.append(hb)
            if len(hbonds_qualified_list) >= 1:
                lig_qualified_list.append(lig)
                hbond_selected_list = hbond_selected_list + hbonds_qualified_list

        if len(lig_qualified_list) >= 1:
            return hbond_selected_list
        else:
            return []

    def hbondFilter(self):
        selected_hbond_list = self.complexAnalysis()
        if len(selected_hbond_list) >= 1:
            return True
        else:
            return False


def main():
    protein = "protein.pdbqt"
    ligand = "ligand.pdb"
    myhbondfiler = HbondFilter("3og7", ligand, protein)
    print(myhbondfiler.hbondFilter())


if __name__ == '__main__':
    main()
