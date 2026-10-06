#!/bin/bash

#SBATCH --job-name=93d1_integron
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=23:00:00
#SBATCH --array=0-5%3

#SBATCH --mail-type=END,FAIL

#SBATCH --output=93d1_integron_%A_%a.out
#SBATCH --error=93d1_integron_%A_%a.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"

IF_SIF="${ROOT}/software/mobilome/integron_finder-2.0.6.sif"

OUTROOT="${ROOT}/95_mobilome/93D1_integronfinder"

GROUPS_LIST=(
    L1
    L2
    L3
    M1
    M2
    M3
)

GROUP="${GROUPS_LIST[$SLURM_ARRAY_TASK_ID]}"

INPUT="${ROOT}/05_assembly/${GROUP}/final.contigs.fa"

OUT="${OUTROOT}/${GROUP}"

mkdir -p "$OUT"

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"


echo "======================================================================"
echo "93D1 - IntegronFinder"
echo "Grupo: $GROUP"
echo "Input: $INPUT"
echo "Inicio: $(date)"
echo "Job: ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "======================================================================"


# ======================================================================
# 1. Validación de inputs
# ======================================================================

[[ -x "$APPTAINER" ]] || {
    echo "ERROR: Apptainer no disponible" >&2
    exit 2
}

[[ -s "$IF_SIF" ]] || {
    echo "ERROR: falta $IF_SIF" >&2
    exit 3
}

[[ -s "$INPUT" ]] || {
    echo "ERROR: falta $INPUT" >&2
    exit 4
}


# ======================================================================
# 2. Metadata del input
# ======================================================================

python3 - "$INPUT" "$OUT" "$GROUP" <<'PY'
import sys
from pathlib import Path

fasta = Path(sys.argv[1])
out = Path(sys.argv[2])
group = sys.argv[3]

n = 0
bp = 0
min_len = None
max_len = 0

current = 0

with fasta.open() as fh:

    for line in fh:

        if line.startswith(">"):

            if current:

                min_len = (
                    current
                    if min_len is None
                    else min(
                        min_len,
                        current
                    )
                )

                max_len = max(
                    max_len,
                    current
                )

            n += 1
            current = 0

        else:

            x = len(
                line.strip()
            )

            bp += x
            current += x


if current:

    min_len = (
        current
        if min_len is None
        else min(
            min_len,
            current
        )
    )

    max_len = max(
        max_len,
        current
    )


with (
    out
    / "93D1_input_metadata.tsv"
).open("w") as fh:

    fh.write(
        "coassembly\tcontigs\ttotal_bp\t"
        "min_contig_bp\tmax_contig_bp\n"
    )

    fh.write(
        f"{group}\t{n}\t{bp}\t"
        f"{min_len or 0}\t{max_len}\n"
    )

PY


# ======================================================================
# 3. IntegronFinder
#
# --local-max:
#   búsqueda más sensible de attC.
#
# --promoter-attI:
#   buscar promotores conocidos y sitios attI.
#
# --lin:
#   los contigs de un ensamblaje se consideran fragmentos lineales.
#
# NO --func-annot:
#   la carga AMR se integrará después con AMRFinderPlus.
# ======================================================================

echo
echo "Ejecutando IntegronFinder..."

# Verificar explícitamente que el contenedor vea el scratch
"$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$IF_SIF" \
    test -r "$INPUT" || {
        echo "ERROR: el contenedor no puede leer $INPUT" >&2
        exit 5
    }

echo "Container input visibility: PASS"

srun --exact \
    -n 1 \
    -c "$SLURM_CPUS_PER_TASK" \
    "$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$IF_SIF" \
    integron_finder \
    "$INPUT" \
    --local-max \
    --promoter-attI \
    --linear \
    --cpu "$SLURM_CPUS_PER_TASK" \
    --outdir "$OUT" \
    > "${OUT}/93D1_integronfinder.stdout" \
    2> "${OUT}/93D1_integronfinder.stderr"


# ======================================================================
# 4. Localizar outputs oficiales
# ======================================================================

mapfile -t INTEGRON_FILES < <(
    find "$OUT" \
        -type f \
        -name '*.integrons' \
        -print \
    | sort
)

mapfile -t SUMMARY_FILES < <(
    find "$OUT" \
        -type f \
        -name '*.summary' \
        -print \
    | sort
)


if [[ "${#INTEGRON_FILES[@]}" -ne 1 ]]
then
    echo \
    "ERROR: esperaba 1 .integrons; encontré ${#INTEGRON_FILES[@]}" \
    >&2

    printf '%s\n' \
        "${INTEGRON_FILES[@]:-}" \
        >&2

    exit 5
fi


if [[ "${#SUMMARY_FILES[@]}" -ne 1 ]]
then
    echo \
    "ERROR: esperaba 1 .summary; encontré ${#SUMMARY_FILES[@]}" \
    >&2

    printf '%s\n' \
        "${SUMMARY_FILES[@]:-}" \
        >&2

    exit 6
fi


INTEGRONS="${INTEGRON_FILES[0]}"
SUMMARY="${SUMMARY_FILES[0]}"


printf \
"file_type\tpath\nintegrons\t%s\nsummary\t%s\n" \
"$INTEGRONS" \
"$SUMMARY" \
> "${OUT}/93D1_result_files.tsv"


# ======================================================================
# 5. Parseo robusto de resultados
# ======================================================================

python3 - \
    "$INTEGRONS" \
    "$SUMMARY" \
    "$OUT" \
    "$GROUP" <<'PY'

import csv
import sys
from collections import Counter
from pathlib import Path


integrons_file = Path(sys.argv[1])
summary_file = Path(sys.argv[2])
out = Path(sys.argv[3])
group = sys.argv[4]


def read_if_tsv(path):

    lines = []

    with path.open(
        errors="replace"
    ) as fh:

        for line in fh:

            if not line.strip():
                continue

            if line.startswith("#"):
                continue

            lines.append(
                line.rstrip("\n")
            )


    if not lines:
        return [], []


    reader = csv.DictReader(
        lines,
        delimiter="\t"
    )

    rows = list(reader)

    return (
        reader.fieldnames or [],
        rows,
    )


int_fields, int_rows = read_if_tsv(
    integrons_file
)

sum_fields, sum_rows = read_if_tsv(
    summary_file
)


# ==================================================================
# Summary counts
# ==================================================================

complete = 0
in0 = 0
calin = 0


for row in sum_rows:

    def as_int(field):

        value = str(
            row.get(
                field,
                "0"
            )
        ).strip()

        try:
            return int(
                float(
                    value or 0
                )
            )

        except ValueError:
            return 0


    complete += as_int(
        "complete"
    )

    in0 += as_int(
        "In0"
    )

    calin += as_int(
        "CALIN"
    )


# Independent cross-check from .integrons
integron_keys = set()

type_counts_from_integrons = Counter()

replicons_with_elements = set()

type_elt_counts = Counter()

annotation_counts = Counter()


for row in int_rows:

    rid = str(
        row.get(
            "ID_replicon",
            ""
        )
    ).strip()

    iid = str(
        row.get(
            "ID_integron",
            ""
        )
    ).strip()

    typ = str(
        row.get(
            "type",
            ""
        )
    ).strip()

    type_elt = str(
        row.get(
            "type_elt",
            ""
        )
    ).strip()

    annotation = str(
        row.get(
            "annotation",
            ""
        )
    ).strip()


    if rid:
        replicons_with_elements.add(
            rid
        )


    if rid or iid:

        integron_keys.add(
            (
                rid,
                iid,
            )
        )


    if rid or iid:

        if typ:

            type_counts_from_integrons[
                (
                    rid,
                    iid,
                    typ,
                )
            ] += 1


    if type_elt:

        type_elt_counts[
            type_elt
        ] += 1


    if annotation:

        annotation_counts[
            annotation
        ] += 1


unique_by_type = Counter()

for rid, iid, typ in (
    type_counts_from_integrons
):

    unique_by_type[
        typ
    ] += 1


# ==================================================================
# Validation of type counts
# ==================================================================

for label, expected in [
    (
        "complete",
        complete,
    ),
    (
        "In0",
        in0,
    ),
    (
        "CALIN",
        calin,
    ),
]:

    observed = unique_by_type.get(
        label,
        0
    )

    if observed != expected:

        raise RuntimeError(
            f"{group}: summary/integrons mismatch "
            f"{label}: summary={expected}, "
            f"integrons={observed}"
        )


# ==================================================================
# Output group summary
# ==================================================================

metrics = [
    (
        "coassembly",
        group,
    ),

    (
        "integron_element_rows",
        len(
            int_rows
        ),
    ),

    (
        "integron_structures_total",
        (
            complete
            + in0
            + calin
        ),
    ),

    (
        "complete_integrons",
        complete,
    ),

    (
        "In0_elements",
        in0,
    ),

    (
        "CALIN_elements",
        calin,
    ),

    (
        "contigs_with_integron_elements",
        len(
            replicons_with_elements
        ),
    ),
]


with (
    out
    / "93D1_group_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for key, value in metrics:

        fh.write(
            f"{key}\t{value}\n"
        )


# ==================================================================
# Type-element counts
# ==================================================================

with (
    out
    / "93D1_type_elt_counts.tsv"
).open("w") as fh:

    fh.write(
        "type_elt\tcount\n"
    )

    for key, value in sorted(
        type_elt_counts.items(),
        key=lambda x: (
            -x[1],
            x[0],
        )
    ):

        fh.write(
            f"{key}\t{value}\n"
        )


# ==================================================================
# Annotation counts
# ==================================================================

with (
    out
    / "93D1_annotation_counts.tsv"
).open("w") as fh:

    fh.write(
        "annotation\tcount\n"
    )

    for key, value in sorted(
        annotation_counts.items(),
        key=lambda x: (
            -x[1],
            x[0],
        )
    ):

        fh.write(
            f"{key}\t{value}\n"
        )


# ==================================================================
# Schema
# ==================================================================

with (
    out
    / "93D1_output_schema.tsv"
).open("w") as fh:

    fh.write(
        "file_type\tcolumns\n"
    )

    fh.write(
        "integrons\t"
        + ";".join(
            int_fields
        )
        + "\n"
    )

    fh.write(
        "summary\t"
        + ";".join(
            sum_fields
        )
        + "\n"
    )


# ==================================================================
# Preview first 10 integron-element rows
# ==================================================================

preview = (
    out
    / "93D1_integrons_preview.tsv"
)


with preview.open(
    "w",
    newline=""
) as fh:

    if int_fields:

        writer = csv.DictWriter(
            fh,
            fieldnames=int_fields,
            delimiter="\t"
        )

        writer.writeheader()

        writer.writerows(
            int_rows[:10]
        )

    else:

        fh.write(
            "NO_INTEGRON_ROWS\n"
        )


print(
    f"{group}: "
    f"complete={complete}, "
    f"In0={in0}, "
    f"CALIN={calin}, "
    f"elements={len(int_rows)}"
)

PY


# ======================================================================
# 6. Metadata/provenance
# ======================================================================

{
    echo -e "field\tvalue"

    echo -e \
"tool\tIntegronFinder"

    echo -e \
"tool_version\t2.0.6"

    echo -e \
"input\t${INPUT}"

    echo -e \
"coassembly\t${GROUP}"

    echo -e \
"local_max\tYES"

    echo -e \
"promoter_attI_search\tYES"

    echo -e \
"default_topology_for_all_contigs\tlinear"

    echo -e \
"functional_annotation_inside_IntegronFinder\tNO"

    echo -e \
"AMR_integration_planned\tAMRFinderPlus_community_loci"

    echo -e \
"CALIN_threshold\tIntegronFinder_default_2_attC"

    echo -e \
"presence_or_mobility_claim_from_IntegronFinder_alone\tNO"

} > "${OUT}/93D1_run_metadata.tsv"


# ======================================================================
# 7. Checksums
# ======================================================================

sha256sum \
    "$INTEGRONS" \
    "$SUMMARY" \
    > "${OUT}/93D1_primary_results_sha256.tsv"


# ======================================================================
# 8. Final marker
# ======================================================================

for F in \
    "${OUT}/93D1_group_summary.tsv" \
    "${OUT}/93D1_type_elt_counts.tsv" \
    "${OUT}/93D1_annotation_counts.tsv" \
    "${OUT}/93D1_output_schema.tsv" \
    "${OUT}/93D1_run_metadata.tsv" \
    "${OUT}/93D1_primary_results_sha256.tsv"
do

    [[ -s "$F" ]] || {
        echo "ERROR: falta/vacío $F" >&2
        exit 7
    }

done


printf \
"COMPLETED\t%s\tjob=%s_%s\tgroup=%s\n" \
"$(date -Iseconds)" \
"$SLURM_ARRAY_JOB_ID" \
"$SLURM_ARRAY_TASK_ID" \
"$GROUP" \
> "${OUT}/93D1_COMPLETE.ok"


echo
echo "======================================================================"
echo "RESULTADO $GROUP"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93D1_group_summary.tsv"

echo
echo "======================================================================"
echo "TYPE_ELT"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93D1_type_elt_counts.tsv"

echo
echo "======================================================================"
echo "93D1 ${GROUP} FINALIZO CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "======================================================================"

