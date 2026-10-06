#!/bin/bash

#SBATCH --job-name=93e2b_AH
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=93e2b_AH_%A_%a.out
#SBATCH --error=93e2b_AH_%A_%a.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"

SIF="${ROOT}/software/genomic_islands/alien-hunter-1.7.7.sif"

BASE="${ROOT}/95_mobilome/93E2_alienhunter_final18"

MANIFEST="${BASE}/93E2A_manifest.tsv"

OUTROOT="${BASE}/predictions"

IDX="${SLURM_ARRAY_TASK_ID}"

mkdir -p "$OUTROOT"

export LC_ALL=C
export LANG=C
export LANGUAGE=C


echo "======================================================================"
echo "93E2B - ALIENHUNTER"
echo "Array index: ${IDX}"
echo "Job: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Inicio: $(date)"
echo "======================================================================"


# ======================================================================
# 1. Validar entradas
# ======================================================================

[[ -x "$APPTAINER" ]] || {
    echo "ERROR: Apptainer no disponible" >&2
    exit 2
}

[[ -s "$SIF" ]] || {
    echo "ERROR: falta contenedor $SIF" >&2
    exit 3
}

[[ -s "$MANIFEST" ]] || {
    echo "ERROR: falta manifest $MANIFEST" >&2
    exit 4
}


# ======================================================================
# 2. Recuperar exactamente la fila del array
# ======================================================================

ROWFILE="${OUTROOT}/.manifest_row_${IDX}.tmp"

python3 - \
    "$MANIFEST" \
    "$IDX" \
    "$ROWFILE" <<'PY'

import csv
import sys
from pathlib import Path

manifest = Path(sys.argv[1])
idx = int(sys.argv[2])
outfile = Path(sys.argv[3])

with manifest.open() as fh:
    rows = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )

if len(rows) != 271:
    raise RuntimeError(
        f"Manifest esperado=271 observado={len(rows)}"
    )

matches = [
    r
    for r in rows
    if int(r["array_index"]) == idx
]

if len(matches) != 1:
    raise RuntimeError(
        f"array_index={idx}: "
        f"esperaba 1 fila y observé {len(matches)}"
    )

r = matches[0]

fields = [
    "array_index",
    "input_id",
    "coassembly",
    "MAG",
    "contig",
    "contig_length_bp",
    "input_fasta",
    "sequence_sha256",
]

with outfile.open("w") as fh:
    fh.write(
        "\t".join(
            r[x]
            for x in fields
        )
        + "\n"
    )

PY


IFS=$'\t' read -r \
    ARRAY_INDEX \
    INPUT_ID \
    COASSEMBLY \
    MAG \
    CONTIG \
    CONTIG_LENGTH \
    INPUT_FASTA \
    SEQ_SHA \
    < "$ROWFILE"

rm -f "$ROWFILE"


[[ "$ARRAY_INDEX" == "$IDX" ]] || {
    echo "ERROR: array index inconsistente" >&2
    exit 5
}

[[ -s "$INPUT_FASTA" ]] || {
    echo "ERROR: falta input FASTA $INPUT_FASTA" >&2
    exit 6
}


TASKOUT="${OUTROOT}/${INPUT_ID}"

mkdir -p "$TASKOUT"

OUTPUT="${TASKOUT}/${INPUT_ID}.alienhunter.embl"


echo
echo "Input ID:    $INPUT_ID"
echo "MAG:         $MAG"
echo "Coassembly:  $COASSEMBLY"
echo "Contig:      $CONTIG"
echo "Length:      $CONTIG_LENGTH"
echo "Input FASTA: $INPUT_FASTA"
echo "Output:      $OUTPUT"


# ======================================================================
# 3. QC del FASTA individual
# ======================================================================

python3 - \
    "$INPUT_FASTA" \
    "$CONTIG" \
    "$CONTIG_LENGTH" <<'PY'

import sys
from pathlib import Path

fasta = Path(sys.argv[1])
expected_contig = sys.argv[2]
expected_length = int(sys.argv[3])

headers = []
seq = []

with fasta.open() as fh:

    for line in fh:

        if line.startswith(">"):
            headers.append(
                line[1:].strip().split()[0]
            )

        else:
            seq.append(
                line.strip()
            )

if len(headers) != 1:
    raise RuntimeError(
        f"Esperaba 1 secuencia; observé {len(headers)}"
    )

if headers[0] != expected_contig:
    raise RuntimeError(
        f"Header={headers[0]} "
        f"esperado={expected_contig}"
    )

observed_length = len(
    "".join(seq)
)

if observed_length != expected_length:
    raise RuntimeError(
        f"Longitud observada={observed_length}, "
        f"esperada={expected_length}"
    )

print(
    f"INPUT_QC=PASS "
    f"contig={headers[0]} "
    f"length={observed_length}"
)

PY


# ======================================================================
# 4. Confirmar visibilidad en Apptainer
# ======================================================================

"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$SIF" \
    test -r "$INPUT_FASTA" || {
        echo "ERROR: contenedor no puede leer $INPUT_FASTA" >&2
        exit 7
    }


# ======================================================================
# 5. Ejecutar AlienHunter
# ======================================================================

echo
echo "Ejecutando AlienHunter..."

"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$SIF" \
    alien_hunter \
    "$INPUT_FASTA" \
    "$OUTPUT" \
    > "${TASKOUT}/${INPUT_ID}.alienhunter.stdout" \
    2> "${TASKOUT}/${INPUT_ID}.alienhunter.stderr"


# AlienHunter puede legítimamente producir cero regiones.
# El archivo debe existir, pero no exigimos que tenga predicciones.

[[ -e "$OUTPUT" ]] || {
    echo "ERROR: AlienHunter no creó $OUTPUT" >&2
    exit 8
}


# ======================================================================
# 6. Auditoría básica del output
# ======================================================================

python3 - \
    "$OUTPUT" \
    "$TASKOUT" \
    "$INPUT_ID" \
    "$IDX" \
    "$MAG" \
    "$COASSEMBLY" \
    "$CONTIG" \
    "$CONTIG_LENGTH" <<'PY'

import re
import sys
from pathlib import Path

output = Path(sys.argv[1])
taskout = Path(sys.argv[2])

input_id = sys.argv[3]
array_index = sys.argv[4]
mag = sys.argv[5]
coassembly = sys.argv[6]
contig = sys.argv[7]
contig_length = int(sys.argv[8])

text = output.read_text(
    errors="replace"
)

# AlienHunter entrega feature-table estilo EMBL.
# Este conteo es solamente QC preliminar;
# el parseo formal de coordenadas será 93E2C.

feature_lines = []

for line in text.splitlines():

    if not line.startswith("FT"):
        continue

    # línea de feature con key + localización,
    # no qualifier.
    if re.match(
        r"^FT\s+\S+\s+\S+",
        line
    ):
        feature_lines.append(line)


summary = (
    taskout
    / f"{input_id}.task_summary.tsv"
)

with summary.open("w") as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        f"array_index\t{array_index}\n"
    )

    fh.write(
        f"input_id\t{input_id}\n"
    )

    fh.write(
        f"coassembly\t{coassembly}\n"
    )

    fh.write(
        f"MAG\t{mag}\n"
    )

    fh.write(
        f"contig\t{contig}\n"
    )

    fh.write(
        f"contig_length_bp\t{contig_length}\n"
    )

    fh.write(
        f"alienhunter_output_bytes\t"
        f"{output.stat().st_size}\n"
    )

    fh.write(
        f"feature_lines_preliminary\t"
        f"{len(feature_lines)}\n"
    )

    fh.write(
        "prediction_interpretation\t"
        "composition_based_candidate_GI_region\n"
    )

    fh.write(
        "zero_predictions_equals_GI_absence\tNO\n"
    )


print(
    f"ALIENHUNTER_OUTPUT_BYTES="
    f"{output.stat().st_size}"
)

print(
    f"PRELIMINARY_FEATURE_LINES="
    f"{len(feature_lines)}"
)

PY


# ======================================================================
# 7. Provenance
# ======================================================================

sha256sum "$OUTPUT" \
    > "${TASKOUT}/${INPUT_ID}.alienhunter.sha256.tsv"


{
    echo -e "field\tvalue"
    echo -e "tool\tAlienHunter"
    echo -e "release\t1.7"
    echo -e "container\t${SIF}"
    echo -e "container_sha256\t542b00d424230e06551f585d968eddeeaa82f98f1b1fb16cfbac8c1536f492e4"
    echo -e "mode\tdefault_IVOM_prediction"
    echo -e "change_point_HMM_c_option_used\tNO"
    echo -e "minimum_input_contig_bp\t50000"
    echo -e "prediction_equals_demonstrated_HGT\tNO"
    echo -e "negative_result_equals_GI_absence\tNO"
} > "${TASKOUT}/${INPUT_ID}.metadata.tsv"


printf \
"COMPLETED\t%s\tjob=%s_%s\tinput_id=%s\tMAG=%s\tcontig=%s\n" \
"$(date -Iseconds)" \
"${SLURM_ARRAY_JOB_ID}" \
"${SLURM_ARRAY_TASK_ID}" \
"${INPUT_ID}" \
"${MAG}" \
"${CONTIG}" \
> "${TASKOUT}/93E2B_COMPLETE.ok"


echo
echo "======================================================================"
echo "93E2B ${INPUT_ID} FINALIZO CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "======================================================================"

