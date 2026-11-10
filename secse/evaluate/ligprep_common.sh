#!/bin/bash

resolve_program() {
  local program
  program=$(command -v -- "$1") || {
    echo "Executable not found: $1" >&2
    return 127
  }
  realpath -- "$program"
}

prepare_docking() {
  workdir=$(realpath -- "$1")
  smi=$(realpath -e -- "$2")
  receptor=$(realpath -e -- "$3")
  cpu_num=$4
  split_dir=$workdir/docking_split
  docking_dir=$workdir/docking_poses
  lig_dir=$workdir/ligands_for_docking
  pdb_dir=$workdir/pdb_files
  sdf_dir=$workdir/sdf_files

  for directory in "$split_dir" "$docking_dir" "$lig_dir" "$pdb_dir" "$sdf_dir"; do
    rm -rf -- "$directory"
    mkdir -- "$directory"
  done

  split -l 100 -d --additional-suffix .smi -- "$smi" "$split_dir/part"
  script=$script_dir/ligprep.py
  export script workdir pdb_dir docking_dir
  export SECSE_PYTHON=${SECSE_PYTHON:-python}
  find "$split_dir" -name '*.smi' -print0 |
    parallel -0 --jobs "$cpu_num" '"$SECSE_PYTHON" "$script" "$workdir" {}'

  ligands=("$lig_dir"/*.pdbqt)
  if ((${#ligands[@]} == 0)); then
    echo "Ligand preparation produced no PDBQT files in $lig_dir" >&2
    return 1
  fi
}

convert_poses() {
  find "$docking_dir" -name '*.pdbqt' -print0 |
    parallel -0 --jobs "$cpu_num" 'obabel -ipdbqt {} -O "$pdb_dir"/{/.}-dp.pdb -m' >/dev/null
}
