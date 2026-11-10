#!/bin/bash
set -euo pipefail
shopt -s nullglob
SECONDS=0
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/ligprep_common.sh"

if [[ -n "${SECSE_UNIDOCK_CONDA_ENV:-}" ]]; then
  conda=$(resolve_program "${SECSE_CONDA_EXE:-conda}")
  engine=("$conda" run --no-capture-output -n "$SECSE_UNIDOCK_CONDA_ENV" unidock)
else
  engine=("$(resolve_program "${UNIDOCK:-unidock}")")
fi
if [[ -n "${SECSE_UNIDOCK_LIBRARY_DIR:-}" ]]; then
  engine=(env "LD_LIBRARY_PATH=$SECSE_UNIDOCK_LIBRARY_DIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" "${engine[@]}")
fi
prepare_docking "$1" "$2" "$3" "${10}"
files=$lig_dir/ligand_index.txt
cd "$lig_dir"
# Uni-Dock reads the index as whitespace-separated paths.
for ligand in "${ligands[@]}"; do
  printf '%s\n' "${ligand##*/}"
done >"$files"

"${engine[@]}" --receptor "$receptor" --ligand_index "$files" --dir "$docking_dir" \
    --center_x "$4" --center_y "$5" --center_z "$6" \
    --size_x "$7" --size_y "$8" --size_z "$9" \
    --exhaustiveness 128 --max_step 20 --refine_step 3 \
    --num_modes 3 --energy_range 3 --verbosity 2 >/dev/null
rm -- "$files"
convert_poses
echo "Docking runtime: $((SECONDS / 60)) minutes $((SECONDS % 60)) seconds."
