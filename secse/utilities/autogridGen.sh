#!/bin/bash

# Zhenting Gao
# Command
# - autogridGen.sh pro.pdbqt grid.gpf
# Update
# - 2023/5/16
#  - This script is created for AutoDock Grid generation

# Parameters
pdbqtFile=$1
gridInputFile=$2 #gpf file
# Please download autogrid4 from https://autodock.scripps.edu/download-autodock4/
autogrid='/tools/docking/autodock/cpu/autogrid4'

if [ ! -f ${autogrid} ]; then
    echo ${autogrid}" is needed but does not exist!"
    echo "
    - Please download autogrid4 from https://autodock.scripps.edu/download-autodock4/
    - Modify this script at line 14 to set the correct path of autogrid4
    "
    exit
fi

if [ ! $# -eq 2 ]; then #Test input parameter
    echo 'autogridGen.sh protein.pdbqt gpfFile'
    echo
    echo "grid.gpf.example is created for your reference."
    echo "
npts 70 70 70
spacing 0.375
gridcenter    17.510   29.510   32.520
" > grid.gpf.example
    cat grid.gpf.example
    exit
fi
npts=$(grep npts ${gridInputFile})
spacing=$(grep spacing ${gridInputFile})
gridcenter=$(grep gridcenter ${gridInputFile})
echo $pdbqtFile
prefix=$(basename ${pdbqtFile} | sed -e 's/.pdbqt$//')
gpfPrefix=$(basename ${gridInputFile} | sed -e 's/.gpf$//')
cat >${gpfPrefix}_production.gpf <<EOF
${npts}
gridfld ${prefix}.maps.fld # grid_data_file
${spacing}
receptor_types A C HD N NA OA SA        # receptor atom types
ligand_types A Br C Cl F HD N NA OA S SA # ligand atom types
receptor ${pdbqtFile}   # macromolecule
${gridcenter}
smooth 0.5                           # store minimum energy w/in rad(A)
map ${prefix}.A.map
map ${prefix}.Br.map
map ${prefix}.C.map
map ${prefix}.Cl.map
map ${prefix}.F.map
map ${prefix}.HD.map
map ${prefix}.N.map
map ${prefix}.NA.map
map ${prefix}.OA.map
map ${prefix}.S.map
map ${prefix}.SA.map
elecmap ${prefix}.e.map    # electrostatic potential map
dsolvmap ${prefix}.d.map              # desolvation potential map
dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant
EOF

${autogrid} -p ${gpfPrefix}_production.gpf -l ${gpfPrefix}_production.glg
