#!/bin/bash

#SBATCH --job-name=97b0_k2_prepare
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=97b0_k2_prepare_%j.out
#SBATCH --error=97b0_k2_prepare_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

SCRIPT="${SCRIPT_DIR}/97b0_prepare_kraken_bracken.sh"

"$SCRIPT"

echo
echo "============================================================"
echo "97B0 FINALIZO CORRECTAMENTE"
echo "KRAKEN2 + BRACKEN + PLUSPF PREPARADOS"
echo "ES SEGURO SALIR"
echo "============================================================"

