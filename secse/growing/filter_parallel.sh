#!/bin/bash
set -euo pipefail
shopt -s nullglob
SECONDS=0
workdir=$(realpath -e -- "$1")
gen=$2
config=$(realpath -e -- "$3")
cpu_num=$4
script=$(cd -- "$(dirname -- "$0")" && pwd)/filter.py
export SECSE_PYTHON=${SECSE_PYTHON:-python}
export script workdir gen config

mkdir -p -- "$workdir/filter_flag"
files=("$workdir/generation_split_by_seed"/*.csv)
flags=()
for file in "${files[@]}"; do
  flag="$workdir/filter_flag/${file##*/}"
  rm -f -- "$flag"
  flags+=("$flag")
done
if ((${#files[@]} > 0)); then
  printf '%s\0' "${files[@]}" |
    parallel -0 --jobs "$cpu_num" '"$SECSE_PYTHON" "$script" {} "$workdir" "$gen" "$config"'
fi

: >"$workdir/filter.csv"
if ((${#flags[@]} > 0)); then
  printf '%s\0' "${flags[@]}" | xargs -0 -r awk -F, '$NF == "PASS"' >>"$workdir/filter.csv"
fi
rm -rf -- "$workdir/generation_split_by_seed"
rm -f -- "$workdir/mutation.csv" "$workdir/mutation.raw" "$workdir/generation.raw"
echo "Filter runtime: $((SECONDS / 60)) minutes $((SECONDS % 60)) seconds."
