#! /bin/bash
# -*- coding:utf-8 _*-
# @author: Yannan Yuan
# @file: ligprep_unidock.sh
# @time: 2024/3/11/17:00

SECONDS=0
workdir=${1}
smi=${2}
receptor=${3}
x=${4}
y=${5}
z=${6}
box_size_x=${7}
box_size_y=${8}
box_size_z=${9}
cpu_num=${10}
script=$SECSE/evaluate/ligprep.py
split_dir=$workdir/docking_split
docking_dir=$workdir/docking_poses
lig_dir=$workdir/ligands_for_docking
pdb_dir=$workdir/pdb_files
sdf_dir=$workdir/sdf_files
conf=$workdir/vina_config.txt
cd "$workdir" || exit
create_clean_directory() {
  dir_name=$1
  if [ -d "$dir_name" ]; then
    echo "Directory $dir_name already exists, removing $dir_name ..."
    rm -rf "$dir_name"
  fi
  if mkdir "$dir_name"; then
    return 0
  else
    echo "Creating directory failed: $dir_name"
    return 1
  fi
}
for dir in "$split_dir" "$docking_dir" "$lig_dir" "$pdb_dir" "$sdf_dir"; do
  create_clean_directory "$dir"
done
# split by line
split -l 100 -d "$smi" "$split_dir"/part --additional-suffix ".smi"

# run ligprep
cd "$split_dir" || exit
find . -name "*smi" | parallel --jobs "$cpu_num" python "$script" "$workdir"

# run unidock
files=ligand_index.txt
cd "$lig_dir" || exit
for i in *pdbqt; do
  echo "$lig_dir/$i"
done >$files

$UNIDOCK --receptor $receptor --ligand_index $files --dir $docking_dir \
    --center_x $x --center_y $y --center_z $z \
    --size_x $box_size_x --size_y $box_size_y --size_z $box_size_z \
    --exhaustiveness 128 --max_step 20 --refine_step 3 \
    --num_modes 3 --energy_range 3 --verbosity 2 >/dev/null
rm $files

find "$docking_dir" -name "*pdbqt" | parallel --jobs "$cpu_num" obabel -ipdbqt {} -O "$pdb_dir"/{/.}-dp.pdb -m &>/dev/null

duration=$SECONDS
echo "Docking runtime: $((duration / 60)) minutes $((duration % 60)) seconds."
