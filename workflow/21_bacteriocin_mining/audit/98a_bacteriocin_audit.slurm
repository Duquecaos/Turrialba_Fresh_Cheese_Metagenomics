#!/bin/bash

#SBATCH --job-name=98a_bact_audit
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=98a_bact_audit_%j.out
#SBATCH --error=98a_bact_audit_%j.err

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

SCRIPT="${SCRIPT_DIR}/98a_bacteriocin_audit.sh"

"$SCRIPT"

echo
echo "============================================================"
echo "98A FINALIZO CORRECTAMENTE"
echo "ES SEGURO SALIR"
echo "============================================================"
