#!/bin/bash
mols=${1}
workdir=${2}
target=${3}
generation=${4}
docking_precision=${5}
#docking_precision=SP
#docking_precision=XP
#docking_precision=HTVS
ligprep_in=ligprep_gen_$generation.inp
glide_in=glide_gen_$generation.in
glide_mae=ligprep_gen_$generation.maegz
cpu=$(grep -c ^processor /proc/cpuinfo)
#

cd "$workdir" || exit

# LigPreparation
echo "Run ligprep ..."

cat >"$ligprep_in" <<EOF
INPUT_FILE_NAME $mols
OUT_MAE $glide_mae
FORCE_FIELD 16
IONIZATION  2
PH  7.2
PH_THRESHOLD  1.0
EPIK  yes
DETERMINE_CHIRALITIES no
IGNORE_CHIRALITIES  no
NUM_STEREOISOMERS 8
EOF

"${SCHRODINGER}/ligprep" -inp $ligprep_in -HOST "localhost:$cpu" -TMPDIR "$workdir" -WAIT

# generate glide input file
cat >"$glide_in" <<EOF
FORCEFIELD  OPLS3e
GRIDFILE  $target
LIGANDFILE  $glide_mae
POSTDOCK_NPOSE  3
PRECISION $docking_precision
POSE_OUTTYPE  ligandlib_sd
COMPRESS_POSES  FALSE
EOF

if [ "$generation" -le 1 ]; then
  # add parameter for fragments docking, add constrains
  cat >>"$glide_in" <<EOF
EXPANDED_SAMPLING   True
MAXKEEP   50000
MAXREF   1200
SCORING_CUTOFF   500.0
EOF
fi

echo "Run glide ..."

# Docking
"${SCHRODINGER}/glide" "$glide_in" -OVERWRITE -adjust -HOST "localhost:$cpu" -TMPDIR "$workdir" -WAIT

echo "Finshing docking of genration $generation !"
