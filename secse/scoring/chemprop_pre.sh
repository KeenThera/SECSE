#!/bin/bash
set -euo pipefail
shopt -s nullglob
# -*- coding:utf-8 _*-
# @author: Lu Chong
# @file: chemprop_pre.sh
# @time: 2021/10/27/16:32
workdir=$(realpath -e -- "$1")
train=$(realpath -e -- "$2")
pre=$(realpath -e -- "$3")
max_gen=${4}
num_output=${5}
seed=${6}
model_dir=$workdir/prediction/models
if ! [[ $max_gen =~ ^[0-9]+$ && $seed =~ ^[0-9]+$ && $num_output =~ ^[1-9][0-9]*$ ]]; then
  echo "Expected numeric generation, seed and positive prediction count" >&2
  exit 2
fi

mkdir -p "$model_dir"

# Chemprop 2.1 drops a singleton batch even in validation and prediction.
training_rows=$(awk 'END { print NR - 1 }' "$train")
split_args=()
if ((training_rows < 8)); then
  echo "Chemprop requires at least eight scored molecules for this workflow" >&2
  exit 2
elif ((training_rows < 20)); then
  split_args=(--split-sizes 0.5 0.25 0.25)
fi

# Chemprop 2.1 reloads its own Lightning checkpoint during training.
# Use a fresh directory before permitting PyTorch to load that checkpoint.
model_name=G"${max_gen}"_seed"${seed}"
model=$(mktemp -d "$model_dir/.$model_name.XXXXXX")
trap 'rm -rf -- "$model"' EXIT
TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1 chemprop train --data-path "$train" --task-type regression --save-dir \
  "$model" --data-seed "$seed" --show-individual-scores --split-type random "${split_args[@]}" -qq

# split files and prediction with CPU Parallelization
split_dir=$workdir/prediction/pre_split_$max_gen
pre_dir="$workdir/prediction/pre_dir_$max_gen"
mkdir -p -- "$split_dir" "$pre_dir"
rm -f -- "$split_dir"/part*.csv "$pre_dir"/part*.csv
split -l 1000 -d --additional-suffix .csv -- "$pre" "$split_dir/part"
parts=("$split_dir"/part*.csv)
if ((${#parts[@]} == 0)); then
  echo "No input for Chemprop prediction: $pre" >&2
  exit 1
fi
sed -i "1i\\id,smiles" "${parts[@]}"
checkpoint="$model/model_0/best.pt"
test -s "$checkpoint"
export pre_dir checkpoint
printf '%s\0' "${parts[@]}" |
  parallel -0 '
    rows=$(awk "END { print NR }" {})
    count=$((rows - 1))
    batch=64
    if [ "$count" -eq 1 ]; then
      batch=1
    elif [ "$((count % 64))" -eq 1 ]; then
      batch=63
    fi
    chemprop predict --test-path {} --preds-path "$pre_dir"/{/} --smiles-columns smiles --model-paths "$checkpoint" --accelerator cpu --batch-size "$batch" -qq
  '
for part in "${parts[@]}"; do
  output="$pre_dir/${part##*/}"
  if ! test -s "$output" || [ "$(awk 'END { print NR }' "$part")" -ne "$(awk 'END { print NR }' "$output")" ]; then
    echo "Chemprop prediction output does not match input rows for $part" >&2
    exit 1
  fi
done

cd "$workdir/prediction"
tail -n +2 -q "$pre_dir"/part*.csv >pre_G"$max_gen".csv

# fetch top predicted compounds
sort -nk3 -t, pre_G"$max_gen".csv >pre_G"$max_gen"_sorted.csv
echo "id,smiles,pred score" >pre_G"$max_gen".csv
head -n "$num_output" pre_G"$max_gen"_sorted.csv >>pre_G"$max_gen".csv
#rm ../pre_G"$max_gen"_sorted.csv

# write mols for next round of docking
pre_docking_dir=$workdir/generation_"$max_gen"_pre
mkdir -p "$pre_docking_dir"
tail -n+2 pre_G"$max_gen".csv | awk -F, '{print $2"\t"$1}' >"$pre_docking_dir"/mols_for_docking_pred.smi
test -s "$pre_docking_dir/mols_for_docking_pred.smi"

rm -rf -- "$model_dir/$model_name"
mv -- "$model" "$model_dir/$model_name"
