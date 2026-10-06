#!/bin/bash

#SBATCH --job-name=97c2a_val_audit
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:20:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=97c2a_val_audit_%j.out
#SBATCH --error=97c2a_val_audit_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

SCRIPT="${SCRIPT_DIR}/97c2a_pathogen_validation_audit.sh"

"$SCRIPT"

echo
echo "============================================================"
echo "97C2A FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
