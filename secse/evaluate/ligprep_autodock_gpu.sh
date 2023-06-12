#! /bin/bash
# -*- coding:utf-8 _*-
# @author: Lu Chong
# @file: ligprep_autodock_gpu.sh
# @time: 2022/2/16/15:25

SECONDS=0
workdir=${1}
smi=${2}
receptor=${3}
cpu_num=${4}
gpu_num=${5}

files=tmp.txt
script=$SECSE/evaluate/ligprep.py
split_dir=$workdir/docking_split
docking_dir=$workdir/docking_poses
lig_dir=$workdir/ligands_for_docking
pdb_dir=$workdir/pdb_files
sdf_dir=$workdir/sdf_files

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
find . -name "*smi" | parallel --jobs "$cpu_num" --bar python "$script" "$workdir"

# run autdock gpu
cd "$lig_dir" || exit
for i in *pdbqt; do
  echo "$lig_dir/$i;$docking_dir/${i%.*}"
done >$files

parallel --jobs "$gpu_num" --bar -I {} -a ${files} -C ";" "$AUTODOCK_GPU/bin/autodock_gpu_128wi" --ffile "$receptor" --lfile {1} --resnam {2} --seed 12345 -D '$(({%}))' -x 0 -n 3 # >/dev/null
#rm $files

# covert dlg file to pdb
cd "$docking_dir" || exit
find . -name "*.dlg" | parallel "grep '^DOCKED' {} >{.}.tmp"
find . -name "*.tmp" | parallel "cut -c9- {} >{.}.pdbqt"
rm ./*.tmp

sed -e "s/USER    Estimated Free Energy of Binding    =/REMARK/g" -i *pdbqt
find "$docking_dir" -name "*pdbqt" | parallel --jobs "$cpu_num" obabel -ipdbqt {} -O "$pdb_dir"/{/.}-dp.pdb -m &>/dev/null

duration=$SECONDS
echo "Docking runtime: $((duration / 60)) minutes $((duration % 60)) seconds."
