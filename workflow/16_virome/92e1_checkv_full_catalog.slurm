#!/bin/bash

#SBATCH --job-name=92e1_checkv
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=23:30:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=92e1_checkv_%j.out
#SBATCH --error=92e1_checkv_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"

SIF="${ROOT}/software/checkv/checkv-1.1.1.sif"
DB="${ROOT}/databases/checkv/checkv-db-v1.5"

INPUT="${ROOT}/86_viral_catalog_preQC/92D_exact_representatives.fna"
CLUSTERS="${ROOT}/86_viral_catalog_preQC/92D_exact_sequence_clusters.tsv"

OUT="${ROOT}/88_checkv_viral_catalog"
RUN="${OUT}/checkv_exact_representatives"

mkdir -p "$OUT"

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

echo "======================================================================"
echo "PASO 92E1 - CheckV catálogo viral completo"
echo "Inicio: $(date)"
echo "Host: $(hostname)"
echo "Job: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "======================================================================"

# ======================================================================
# 1. VALIDACIONES
# ======================================================================

[[ -x "$APPTAINER" ]] || {
    echo "ERROR: falta Apptainer." >&2
    exit 2
}

[[ -s "$SIF" ]] || {
    echo "ERROR: falta CheckV SIF: $SIF" >&2
    exit 3
}

[[ -d "$DB" ]] || {
    echo "ERROR: falta CheckV DB: $DB" >&2
    exit 4
}

[[ -s "$INPUT" ]] || {
    echo "ERROR: falta FASTA viral: $INPUT" >&2
    exit 5
}

[[ -s "$CLUSTERS" ]] || {
    echo "ERROR: falta tabla de clusters: $CLUSTERS" >&2
    exit 6
}


# ======================================================================
# 2. VALIDAR INPUT
# ======================================================================

echo
echo "======================================================================"
echo "1. VALIDACION INPUT"
echo "======================================================================"

python3 - "$INPUT" "$CLUSTERS" <<'PY'

import csv
import sys
from pathlib import Path

fasta = Path(sys.argv[1])
clusters = Path(sys.argv[2])

ids = []
bp = 0

with fasta.open() as fh:
    for line in fh:

        if line.startswith(">"):
            ids.append(
                line[1:].strip().split()[0]
            )
        else:
            bp += len(line.strip())

if len(ids) != len(set(ids)):
    raise RuntimeError(
        "Hay IDs duplicados en FASTA."
    )

with clusters.open() as fh:
    rows = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )

cluster_ids = {
    x["exact_cluster_id"]
    for x in rows
}

print("FASTA_sequences:", len(ids))
print("FASTA_bp:", bp)
print("cluster_rows:", len(rows))
print("unique_cluster_ids:", len(cluster_ids))

if len(ids) != 6751:
    raise RuntimeError(
        f"Esperaba 6751 secuencias; hay {len(ids)}"
    )

if len(rows) != 6751:
    raise RuntimeError(
        f"Esperaba 6751 clusters; hay {len(rows)}"
    )

if set(ids) != cluster_ids:
    raise RuntimeError(
        "Los IDs del FASTA no coinciden con los clusters 92D."
    )

print("INPUT_QC: PASS")

PY


# ======================================================================
# 3. EJECUTAR CHECKV
# ======================================================================

echo
echo "======================================================================"
echo "2. CHECKV END-TO-END"
echo "======================================================================"

QUALITY="${RUN}/quality_summary.tsv"

if [[ -s "$QUALITY" ]]
then

    EXISTING_ROWS="$(
        awk 'END{print NR-1}' "$QUALITY"
    )"

    if [[ "$EXISTING_ROWS" -eq 6751 ]]
    then
        echo "Salida CheckV completa preexistente."
        echo "No se recalculará."
    else
        echo "Salida incompleta detectada (${EXISTING_ROWS} filas)."
        echo "Eliminando ejecución incompleta."
        rm -rf "$RUN"
    fi

fi


if [[ ! -s "$QUALITY" ]]
then

    rm -rf "$RUN"

    srun --exact \
        -n 1 \
        -c "$SLURM_CPUS_PER_TASK" \
        "$APPTAINER" exec \
        --bind "${ROOT}:/work" \
        --env CHECKVDB="/work/databases/checkv/checkv-db-v1.5" \
        "$SIF" \
        checkv end_to_end \
            /work/86_viral_catalog_preQC/92D_exact_representatives.fna \
            /work/88_checkv_viral_catalog/checkv_exact_representatives \
            -t "$SLURM_CPUS_PER_TASK"

fi


# ======================================================================
# 4. VALIDAR CHECKV
# ======================================================================

echo
echo "======================================================================"
echo "3. VALIDACION CHECKV"
echo "======================================================================"

for F in \
    "${RUN}/quality_summary.tsv" \
    "${RUN}/contamination.tsv" \
    "${RUN}/completeness.tsv"
do

    [[ -s "$F" ]] || {
        echo "ERROR: falta salida CheckV: $F" >&2
        exit 7
    }

    echo "OK: $(basename "$F")"

done


NQUALITY="$(
    awk 'END{print NR-1}' \
        "${RUN}/quality_summary.tsv"
)"

echo "quality_summary_records=${NQUALITY}"

[[ "$NQUALITY" -eq 6751 ]] || {
    echo "ERROR: se esperaban 6751 registros CheckV." >&2
    exit 8
}


# ======================================================================
# 5. INTEGRAR CHECKV + geNomad 92D
# ======================================================================

echo
echo "======================================================================"
echo "4. INTEGRACION"
echo "======================================================================"

srun --exact \
    -n 1 \
    -c 1 \
    python3 - <<'PY'

import csv
import math
import os
import re
import statistics
from collections import Counter
from pathlib import Path

root = Path(
    f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"
)

catalog = (
    root
    / "86_viral_catalog_preQC"
    / "92D_exact_sequence_clusters.tsv"
)

checkv_dir = (
    root
    / "88_checkv_viral_catalog"
    / "checkv_exact_representatives"
)

quality_file = (
    checkv_dir
    / "quality_summary.tsv"
)

out = (
    root
    / "88_checkv_viral_catalog"
)


def clean(x):
    if x is None:
        return ""

    x = str(x).strip()

    if x.upper() in {
        "NA",
        "NAN",
        "NONE",
        ""
    }:
        return ""

    return x


def num(x):
    x = clean(x)

    if not x:
        return None

    try:
        return float(x)
    except ValueError:
        return None


def slug(x):

    x = clean(x).lower()

    x = re.sub(
        r"[^a-z0-9]+",
        "_",
        x
    )

    return x.strip("_") or "blank"


# ----------------------------------------------------------
# 92D
# ----------------------------------------------------------

with catalog.open() as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    clusters = {
        x["exact_cluster_id"]: x
        for x in r
    }


# ----------------------------------------------------------
# CheckV
# ----------------------------------------------------------

with quality_file.open() as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    qrows = list(r)


if len(qrows) != 6751:
    raise RuntimeError(
        f"quality_summary tiene {len(qrows)} filas"
    )


checkv = {
    x["contig_id"]: x
    for x in qrows
}

if len(checkv) != 6751:
    raise RuntimeError(
        "IDs duplicados en quality_summary."
    )


if set(clusters) != set(checkv):

    missing_checkv = (
        set(clusters)
        - set(checkv)
    )

    missing_catalog = (
        set(checkv)
        - set(clusters)
    )

    raise RuntimeError(
        "IDs 92D / CheckV no coinciden. "
        f"missing_checkv={len(missing_checkv)}, "
        f"missing_catalog={len(missing_catalog)}"
    )


# ----------------------------------------------------------
# Tabla integrada
# ----------------------------------------------------------

integrated = []

for cid in sorted(clusters):

    c = clusters[cid]
    q = checkv[cid]

    comp = num(
        q.get("completeness")
    )

    contam = num(
        q.get("contamination")
    )

    ge_provirus = (
        "Provirus"
        in {
            x.strip()
            for x in c[
                "topologies"
            ].split(";")
            if x.strip()
        }
    )

    cv_provirus = (
        clean(
            q.get("provirus")
        ).lower()
        == "yes"
    )

    if ge_provirus and cv_provirus:
        concordance = (
            "both_provirus"
        )

    elif ge_provirus:
        concordance = (
            "genomad_only_provirus"
        )

    elif cv_provirus:
        concordance = (
            "checkv_only_provirus"
        )

    else:
        concordance = (
            "neither_provirus"
        )

    integrated.append({
        "exact_cluster_id":
            cid,

        "sequence_sha256":
            c["sequence_sha256"],

        "sequence_length":
            c["length"],

        "n_exact_members":
            c["n_members"],

        "n_coassemblies":
            c["n_coassemblies"],

        "coassemblies":
            c["coassemblies"],

        "n_producers":
            c["n_producers"],

        "producers":
            c["producers"],

        "genomad_topologies":
            c["topologies"],

        "genomad_terminal_taxa":
            c["terminal_taxa"],

        "genomad_max_virus_score":
            c["max_virus_score"],

        "genomad_min_fdr":
            c["min_fdr"],

        "genomad_max_n_hallmarks":
            c["max_n_hallmarks"],

        "checkv_contig_length":
            q["contig_length"],

        "checkv_provirus":
            q["provirus"],

        "checkv_proviral_length":
            clean(
                q["proviral_length"]
            ),

        "checkv_gene_count":
            q["gene_count"],

        "checkv_viral_genes":
            q["viral_genes"],

        "checkv_host_genes":
            q["host_genes"],

        "checkv_quality":
            q["checkv_quality"],

        "miuvig_quality":
            q["miuvig_quality"],

        "completeness":
            (
                comp
                if comp is not None
                else ""
            ),

        "completeness_method":
            clean(
                q[
                    "completeness_method"
                ]
            ),

        "contamination":
            (
                contam
                if contam is not None
                else ""
            ),

        "kmer_freq":
            clean(
                q["kmer_freq"]
            ),

        "warnings":
            clean(
                q["warnings"]
            ),

        "provirus_concordance":
            concordance,
    })


outfile = (
    out
    / "92E1_checkv_integrated.tsv"
)

with outfile.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            integrated[0].keys()
        ),
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(integrated)


# ----------------------------------------------------------
# Métricas globales
# ----------------------------------------------------------

quality_counts = Counter(
    x["checkv_quality"]
    for x in integrated
)

miuvig_counts = Counter(
    x["miuvig_quality"]
    for x in integrated
)

methods = Counter(
    x["completeness_method"]
    if x["completeness_method"]
    else "NA"
    for x in integrated
)

concordance_counts = Counter(
    x["provirus_concordance"]
    for x in integrated
)


comps = [
    float(x["completeness"])
    for x in integrated
    if x["completeness"] != ""
]

contams = [
    float(x["contamination"])
    for x in integrated
    if x["contamination"] != ""
]


summary = [
    (
        "exact_viral_sequences",
        len(integrated)
    ),

    (
        "completeness_nonNA",
        len(comps)
    ),

    (
        "completeness_NA",
        len(integrated) - len(comps)
    ),

    (
        "completeness_ge_90",
        sum(x >= 90 for x in comps)
    ),

    (
        "completeness_ge_75",
        sum(x >= 75 for x in comps)
    ),

    (
        "completeness_ge_50",
        sum(x >= 50 for x in comps)
    ),

    (
        "completeness_lt_50",
        sum(x < 50 for x in comps)
    ),

    (
        "median_completeness_nonNA",
        (
            statistics.median(comps)
            if comps
            else "NA"
        )
    ),

    (
        "contamination_nonNA",
        len(contams)
    ),

    (
        "contamination_eq_0",
        sum(x == 0 for x in contams)
    ),

    (
        "contamination_gt_0",
        sum(x > 0 for x in contams)
    ),

    (
        "contamination_gt_5",
        sum(x > 5 for x in contams)
    ),

    (
        "contamination_gt_10",
        sum(x > 10 for x in contams)
    ),

    (
        "contamination_gt_50",
        sum(x > 50 for x in contams)
    ),

    (
        "median_contamination",
        (
            statistics.median(contams)
            if contams
            else "NA"
        )
    ),

    (
        "warnings_nonempty",
        sum(
            bool(x["warnings"])
            for x in integrated
        )
    ),

    (
        "warning_no_viral_genes",
        sum(
            "no viral genes detected"
            in x["warnings"].lower()
            for x in integrated
        )
    ),
]


for key, value in sorted(
    quality_counts.items()
):
    summary.append(
        (
            f"checkv_quality_{slug(key)}",
            value
        )
    )


for key, value in sorted(
    miuvig_counts.items()
):
    summary.append(
        (
            f"miuvig_quality_{slug(key)}",
            value
        )
    )


for key, value in sorted(
    concordance_counts.items()
):
    summary.append(
        (
            f"provirus_{slug(key)}",
            value
        )
    )


summary_file = (
    out
    / "92E1_checkv_global_summary.tsv"
)

with summary_file.open("w") as fh:

    fh.write("metric\tvalue\n")

    for k, v in summary:
        fh.write(
            f"{k}\t{v}\n"
        )


# ----------------------------------------------------------
# Métodos de completitud
# ----------------------------------------------------------

method_file = (
    out
    / "92E1_completeness_methods.tsv"
)

with method_file.open("w") as fh:

    fh.write(
        "completeness_method\tcount\n"
    )

    for method, count in sorted(
        methods.items(),
        key=lambda x: (
            -x[1],
            x[0]
        )
    ):

        fh.write(
            f"{method}\t{count}\n"
        )


# ----------------------------------------------------------
# Calidad según contexto geNomad
# ----------------------------------------------------------

context_rows = []

for context in [
    "PROVIRUS",
    "NON_PROVIRUS",
]:

    if context == "PROVIRUS":

        rows = [
            x
            for x in integrated
            if "Provirus"
            in x["genomad_topologies"].split(";")
        ]

    else:

        rows = [
            x
            for x in integrated
            if "Provirus"
            not in x["genomad_topologies"].split(";")
        ]

    cvals = [
        float(x["completeness"])
        for x in rows
        if x["completeness"] != ""
    ]

    xvals = [
        float(x["contamination"])
        for x in rows
        if x["contamination"] != ""
    ]

    context_rows.append({
        "genomad_context":
            context,

        "sequences":
            len(rows),

        "checkv_provirus_yes":
            sum(
                x["checkv_provirus"].lower()
                == "yes"
                for x in rows
            ),

        "completeness_nonNA":
            len(cvals),

        "median_completeness":
            (
                statistics.median(cvals)
                if cvals
                else ""
            ),

        "completeness_ge_90":
            sum(
                x >= 90
                for x in cvals
            ),

        "completeness_ge_50":
            sum(
                x >= 50
                for x in cvals
            ),

        "contamination_gt_0":
            sum(
                x > 0
                for x in xvals
            ),

        "contamination_gt_5":
            sum(
                x > 5
                for x in xvals
            ),

        "contamination_gt_10":
            sum(
                x > 10
                for x in xvals
            ),
    })


context_file = (
    out
    / "92E1_quality_by_genomad_context.tsv"
)

with context_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            context_rows[0].keys()
        ),
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(context_rows)


# ----------------------------------------------------------
# Warning counts
# ----------------------------------------------------------

warning_counts = Counter()

for x in integrated:

    w = x["warnings"]

    if not w:
        continue

    for item in re.split(
        r"[;,]",
        w
    ):

        item = item.strip()

        if item:
            warning_counts[item] += 1


warning_file = (
    out
    / "92E1_warning_counts.tsv"
)

with warning_file.open("w") as fh:

    fh.write("warning\tcount\n")

    for warning, count in sorted(
        warning_counts.items(),
        key=lambda x: (
            -x[1],
            x[0]
        )
    ):

        fh.write(
            f"{warning}\t{count}\n"
        )


# ----------------------------------------------------------
# Print
# ----------------------------------------------------------

print()
print("=" * 78)
print("CHECKV GLOBAL")
print("=" * 78)

for k, v in summary:
    print(k, v, sep="\t")

print()
print("=" * 78)
print("COMPLETENESS METHODS")
print("=" * 78)

for k, v in sorted(
    methods.items(),
    key=lambda x: (
        -x[1],
        x[0]
    )
):
    print(k, v, sep="\t")

print()
print("=" * 78)
print("CONTEXT")
print("=" * 78)

for x in context_rows:
    print(
        x["genomad_context"],
        x["sequences"],
        x["checkv_provirus_yes"],
        x["completeness_nonNA"],
        x["median_completeness"],
        x["completeness_ge_90"],
        x["completeness_ge_50"],
        x["contamination_gt_5"],
        x["contamination_gt_10"],
        sep="\t"
    )

PY


# ======================================================================
# 6. METADATA
# ======================================================================

{
    echo -e "field\tvalue"
    echo -e "analysis\tCheckV_full_exact_viral_catalog"
    echo -e "input\t86_viral_catalog_preQC/92D_exact_representatives.fna"
    echo -e "input_sequences\t6751"
    echo -e "sequence_level\texact_sequence_representatives"
    echo -e "vOTU_clustering\tNOT_PERFORMED"
    echo -e "sequence_trimming\tNOT_PERFORMED"
    echo -e "CheckV_version\t1.1.1"
    echo -e "CheckV_DB\tcheckv-db-v1.5"
    echo -e "threads\t${SLURM_CPUS_PER_TASK}"
    echo -e "date\t$(date -Iseconds)"
} > "${OUT}/92E1_run_metadata.tsv"


# ======================================================================
# 7. VALIDACION FINAL
# ======================================================================

echo
echo "======================================================================"
echo "5. VALIDACION FINAL"
echo "======================================================================"

for F in \
    "${OUT}/92E1_checkv_integrated.tsv" \
    "${OUT}/92E1_checkv_global_summary.tsv" \
    "${OUT}/92E1_completeness_methods.tsv" \
    "${OUT}/92E1_quality_by_genomad_context.tsv" \
    "${OUT}/92E1_warning_counts.tsv" \
    "${OUT}/92E1_run_metadata.tsv"
do

    [[ -s "$F" ]] || {
        echo "ERROR: salida ausente/vacía: $F" >&2
        exit 9
    }

    echo "OK: $(basename "$F")"

done


echo
echo "======================================================================"
echo "RESUMEN GLOBAL CHECKV"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92E1_checkv_global_summary.tsv"


echo
echo "======================================================================"
echo "CALIDAD SEGUN CONTEXTO geNomad"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92E1_quality_by_genomad_context.tsv"


echo
echo "======================================================================"
echo "METODOS DE COMPLETITUD"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92E1_completeness_methods.tsv"


echo
echo "======================================================================"
echo "WARNINGS"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92E1_warning_counts.tsv" \
    | head -n 30


echo
echo "======================================================================"
echo "PASO 92E1 FINALIZÓ CORRECTAMENTE."
echo "NO se eliminaron ni recortaron secuencias."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "======================================================================"

