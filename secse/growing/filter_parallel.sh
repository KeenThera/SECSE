#! /bin/bash
# @author: Lu Chong
# @file: filter_parallel.sh
# @time: 2021/ 03/03/9:26

SECONDS=0
workdir=${1}
gen=${2}
config=${3}
cpu_num=${4}
script=$SECSE/growing/filter.py
files=tmp.txt
cd "${workdir}"/generation_split_by_seed || exit
for i in *.csv; do
  echo "$i;$workdir;$gen;$config"
done >$files

mkdir -p ../filter_flag
# filter default
parallel --jobs "$cpu_num" --bar -I {} -a ${files} -C ";" python "$script"
rm $files
cd "${workdir}"/filter_flag || exit
for i in *.csv; do
  echo "$i" | parallel grep PASS
done >"${workdir}"/filter.csv
cd "${workdir}" || exit
rm -r filter_flag/
rm -r generation_split_by_seed/ mutation.csv mutation.raw generation.raw
duration=$SECONDS
echo "Filter runtime: $((duration / 60)) minutes $((duration % 60)) seconds."
