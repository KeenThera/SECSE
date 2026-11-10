#!/bin/bash
set -euo pipefail
shopt -s nullglob
SECONDS=0
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$script_dir/ligprep_common.sh"

vina=$(resolve_program "${VINA:-vina}")
prepare_docking "$1" "$2" "$3" "${10}"
x=$4
y=$5
z=$6
box_size_x=$7
box_size_y=$8
box_size_z=$9
conf=$workdir/vina_config.txt
export vina conf

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

printf '%s\0' "${ligands[@]}" |
  parallel -0 --jobs "$cpu_num" '"$vina" --config "$conf" --ligand {} --out "$docking_dir"/{/}' >/dev/null
convert_poses
echo "Docking runtime: $((SECONDS / 60)) minutes $((SECONDS % 60)) seconds."
