#! /bin/bash
# -*- coding:utf-8 _*-
# @author: Lu Chong
# @file: ligprep_vina_parallel.sh
# @time: 2021/9/8/09:52

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
files=$RANDOM
script=$SECSE/evaluate/ligprep.py
split_dir=$workdir/docking_split
vina_dir=$workdir/vina_poses
lig_dir=$workdir/ligands_for_vina
pdb_dir=$workdir/pdb_files
sdf_dir=$workdir/sdf_files
conf=$workdir/vina_config.txt
cd "$workdir" || exit
mkdir -p "$split_dir" "$vina_dir" "$lig_dir" "$pdb_dir" "$sdf_dir"
# split by line
split -l 100 -d "$smi" "$split_dir"/part --additional-suffix ".smi"

# run ligprep
cd "$split_dir" || exit
find . -name "*smi" | parallel --jobs "$cpu_num" --bar python "$script" "$workdir"

# write vina config file
cat >"$conf" <<EOF
receptor = $receptor
center_x =  $x
center_y =  $y
center_z =  $z

size_x = $box_size_x
size_y = $box_size_y
size_z = $box_size_z

seed = 12345
cpu = 1
num_modes = 3
energy_range = 3
exhaustiveness = 16
verbosity = 0
EOF

# run vina
cd "$lig_dir" || exit
for i in *pdbqt; do
  echo "$lig_dir/$i;$vina_dir/$i"
done >$files

# ignore Vina stdout
parallel --jobs "$cpu_num" --bar -I {} -a ${files} -C ";" "$VINA" --config "$conf" --ligand {1} --out {2} >/dev/null
rm $files

find "$vina_dir" -name "*pdbqt" | parallel --jobs "$cpu_num" obabel -ipdbqt {} -O "$pdb_dir"/{/.}-dp.pdb -m &>/dev/null

duration=$SECONDS
echo "Docking runtime: $((duration / 60)) minutes $((duration % 60)) seconds."
