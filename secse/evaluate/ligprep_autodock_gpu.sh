#!/bin/bash
set -euo pipefail
shopt -s nullglob
SECONDS=0
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/ligprep_common.sh"

gpu_executable=$(resolve_program "${AUTODOCK_GPU:?Set AUTODOCK_GPU to the installation directory}/bin/autodock_gpu_128wi")
prepare_docking "$1" "$2" "$3" "$4"
gpu_num=$5
export gpu_executable receptor

printf '%s\0' "${ligands[@]}" |
  parallel -0 --jobs "$gpu_num" '"$gpu_executable" --ffile "$receptor" --lfile {} --resnam "$docking_dir"/{/.} --seed 12345 -D $(({%})) -x 0 -n 3'

# Extract only docked records and retain the existing score remark.
cd "$docking_dir"
for dlg in *.dlg; do
  grep '^DOCKED' "$dlg" | cut -c9- |
    sed 's/USER    Estimated Free Energy of Binding    =/REMARK/g' >"${dlg%.dlg}.pdbqt"
done
convert_poses
echo "Docking runtime: $((SECONDS / 60)) minutes $((SECONDS % 60)) seconds."
