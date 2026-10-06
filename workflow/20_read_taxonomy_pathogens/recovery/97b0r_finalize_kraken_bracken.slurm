#!/bin/bash

#SBATCH --job-name=97b0r_k2_finalize
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=97b0r_k2_finalize_%j.out
#SBATCH --error=97b0r_k2_finalize_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

SCRIPT="${SCRIPT_DIR}/97b0r_finalize_kraken_bracken.sh"

"$SCRIPT"

echo
echo "============================================================"
echo "97B0R FINALIZO CORRECTAMENTE"
echo "KRAKEN2 + BRACKEN + PLUSPF VALIDADOS"
echo "ES SEGURO SALIR"
echo "============================================================"

