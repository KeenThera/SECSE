#!/bin/bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: ./setup.sh BACKEND [--env NAME] [--install-dir PATH]

Create a new native Conda environment for SECSE.
  vina           Install Vina and SECSE dependencies (CPU).
  unidock        Install native Uni-Dock 1.2.0 in a separate engine environment.
  autodock-gpu   Connect an existing installation using --install-dir.
  glide          Connect an existing licensed installation using --install-dir.

The default environment name is secse. Add more backends to an existing
environment with 'python -m secse backend add'.
EOF
}

fail() { echo "SECSE setup: $*" >&2; exit 2; }

if [[ ${1:-} == "--help" || ${1:-} == "-h" ]]; then
  usage
  exit 0
fi
backend=${1:-}
[[ -n "$backend" ]] || { usage; exit 2; }
shift
environment=secse
install_dir=
while (($#)); do
  case "$1" in
    --env) [[ $# -ge 2 && -n "$2" ]] || fail "--env requires a name"; environment=$2 ;;
    --install-dir) [[ $# -ge 2 && -n "$2" ]] || fail "--install-dir requires a path"; install_dir=$2 ;;
    *) fail "unknown option $1" ;;
  esac
  shift 2
done

packages=()
backend_args=()
engine_environment=
case "$backend" in
  vina) packages+=('vina>=1.2.5,<2') ;;
  unidock) engine_environment="$environment-unidock" ;;
  autodock-gpu)
    [[ -n "$install_dir" && -x "$install_dir/bin/autodock_gpu_128wi" ]] ||
      fail "pass --install-dir pointing to the AutoDock GPU installation"
    ;;
  glide)
    [[ -n "$install_dir" && -x "$install_dir/glide" && -x "$install_dir/ligprep" ]] ||
      fail "pass --install-dir pointing to your Schrödinger installation"
    ;;
  *) fail "unknown backend $backend" ;;
esac
if [[ -n "$install_dir" ]]; then
  [[ "$backend" == "autodock-gpu" || "$backend" == "glide" ]] ||
    fail "--install-dir is for AutoDock GPU or Glide"
  backend_args+=(--install-dir "$(cd -- "$install_dir" && pwd)")
fi

[[ $(uname -s) == Linux ]] || fail "this setup supports Linux"
[[ $environment =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]] || fail "invalid environment name"
conda=${CONDA_EXE:-conda}
command -v -- "$conda" >/dev/null || fail "install Conda first and put conda on PATH"
if "$conda" run -n "$environment" true >/dev/null 2>&1; then
  fail "environment $environment already exists. Activate it to add a backend, or choose another --env name"
fi

cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
if [[ -n "$engine_environment" ]]; then
  if "$conda" run -n "$engine_environment" true >/dev/null 2>&1; then
    fail "engine environment $engine_environment already exists. Choose another --env name or register it with backend add"
  fi
  "$conda" create --yes --name "$engine_environment" --override-channels -c conda-forge 'unidock=1.2.0'
  backend_args+=(--conda-env "$engine_environment")
fi
# BioPandas' Python metadata requires looseversion 1.1.2.
"$conda" create --yes --name "$environment" --override-channels -c conda-forge \
  python=3.12 pip 'numpy>=1.26,<3' 'pandas>=2.2,<3' 'rdkit>=2024.9' 'pandarallel>=1.6.5,<1.7' \
  'openbabel>=3.1.1,<3.2' 'xlrd>=2.0.1,<2.1' 'tqdm>=4.67,<5' 'biopandas>=0.5.1,<0.6' \
  'looseversion=1.1.2' 'loguru>=0.7,<1' parallel perl "${packages[@]}"
"$conda" run --no-capture-output -n "$environment" python -s -m pip check
if [[ "$backend" == "vina" ]]; then
  "$conda" run --no-capture-output -n "$environment" bash -c \
    'exec python -m secse backend add "$1" --executable "$CONDA_PREFIX/bin/$1"' _ "$backend"
else
  "$conda" run --no-capture-output -n "$environment" python -m secse backend add "$backend" "${backend_args[@]}"
fi
"$conda" run --no-capture-output -n "$environment" python -m secse doctor --backend "$backend"
echo "Activate with: conda activate $environment"
if [[ "$backend" == "vina" || "$backend" == "unidock" ]]; then
  echo "Try: python -m secse demo phgdh --backend $backend"
elif [[ "$backend" == "autodock-gpu" ]]; then
  echo "Try: python -m secse demo phgdh --backend $backend --target /path/to/phgdh.maps.fld"
else
  echo "Try: python -m secse demo phgdh --backend $backend --target /path/to/phgdh-grid.zip"
fi
