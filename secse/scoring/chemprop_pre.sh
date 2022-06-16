#! /bin/bash
# -*- coding:utf-8 _*-
# @author: Lu Chong
# @file: chemprop_pre.sh
# @time: 2021/10/27/16:32
workdir=${1}
train=${2}
pre=${3}
max_gen=${4}
num_output=${5}
seed=${6}
model_dir=$workdir/prediction/models/
files=tmp.txt

mkdir -p "$model_dir"

# all data
model="$model_dir"/G"$max_gen"_seed"$seed"
chemprop_train --data_path "$train" --dataset_type regression --save_dir \
  "$model" --seed "$seed" --save_smiles --save_preds --show_individual_scores \
  --extra_metrics {r2,mae,mse} --split_type random

# split files and prediction with CPU Parallelization
split_dir=$workdir/prediction/pre_split_$max_gen
mkdir -p "$split_dir"
split -l 1000 -d "$pre" "$split_dir"/part --additional-suffix ".csv"

pre_dir="$workdir"/prediction/pre_dir_$max_gen
cd "$split_dir" || exit
# add header
sed -i "1i\\id,smiles" part*.csv
for i in *.csv; do
  echo "$split_dir/$i;$pre_dir/$i"
done >$files

# run chemprop_predict
parallel --bar -I {} -a ${files} -C ";" chemprop_predict --test_path {1} --preds_path {2} --smiles_columns smiles --checkpoint_path "$model"/fold_0/model_0/model.pt --no_cuda

# merge prediction
cd "$workdir"/prediction || exit
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
