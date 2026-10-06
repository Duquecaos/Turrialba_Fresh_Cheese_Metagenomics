#!/bin/bash

#SBATCH --job-name=92e3_cleancheck
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=92e3_cleancheck_%j.out
#SBATCH --error=92e3_cleancheck_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPT="${APPT:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"

CHECKV_SIF="${ROOT}/software/checkv/checkv-1.1.1.sif"
CHECKV_DB="${ROOT}/databases/checkv/checkv-db-v1.5"

BLAST_SIF="${BLAST_SIF:-/opt/ohpc/pub/containers/BIO/blast-2.16.0.sif}"

INPUT="${ROOT}/89_viral_catalog_hostclean/92E2_hostclean_exact_representatives.fna"

OUT="${ROOT}/90_viral_hostclean_validation"
CVOUT="${OUT}/checkv_hostclean"

mkdir -p "$OUT"

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

echo "======================================================================"
echo "PASO 92E3 - VALIDACION POST-LIMPIEZA + AUDITORIA vOTU"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "======================================================================"

for F in \
    "$CHECKV_SIF" \
    "$BLAST_SIF" \
    "$INPUT"
do
    [[ -s "$F" ]] || {
        echo "ERROR: falta $F" >&2
        exit 2
    }
done

[[ -d "$CHECKV_DB" ]] || {
    echo "ERROR: falta CheckV DB." >&2
    exit 3
}


# ======================================================================
# 1. AUDITAR INPUT
# ======================================================================

echo
echo "======================================================================"
echo "1. INPUT VCLEAN"
echo "======================================================================"

NSEQ="$(
    grep -c '^>' "$INPUT"
)"

echo "sequences=${NSEQ}"

[[ "$NSEQ" -eq 6748 ]] || {
    echo "ERROR: esperaba 6748 secuencias." >&2
    exit 4
}


# ======================================================================
# 2. RE-EJECUTAR CHECKV SOBRE SECUENCIAS LIMPIAS
# ======================================================================

echo
echo "======================================================================"
echo "2. CHECKV POST-LIMPIEZA"
echo "======================================================================"

rm -rf "$CVOUT"

srun --exact \
    -n 1 \
    -c "$SLURM_CPUS_PER_TASK" \
    "$APPT" exec \
    --bind "${ROOT}:/work" \
    --env CHECKVDB="/work/databases/checkv/checkv-db-v1.5" \
    "$CHECKV_SIF" \
    checkv end_to_end \
        /work/89_viral_catalog_hostclean/92E2_hostclean_exact_representatives.fna \
        /work/90_viral_hostclean_validation/checkv_hostclean \
        -t "$SLURM_CPUS_PER_TASK"


QUALITY="${CVOUT}/quality_summary.tsv"

[[ -s "$QUALITY" ]] || {
    echo "ERROR: CheckV no generó quality_summary.tsv" >&2
    exit 5
}

NQ="$(
    awk 'END{print NR-1}' "$QUALITY"
)"

echo "CheckV records=${NQ}"

[[ "$NQ" -eq 6748 ]] || {
    echo "ERROR: CheckV produjo ${NQ}, esperaba 6748." >&2
    exit 6
}


# ======================================================================
# 3. RESUMEN POST-CLEAN
# ======================================================================

srun --exact -n 1 -c 1 python3 - <<'PY'

import csv
import os
import statistics
from collections import Counter
from pathlib import Path

root = Path(
    f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"
)

qfile = (
    root
    / "90_viral_hostclean_validation"
    / "checkv_hostclean"
    / "quality_summary.tsv"
)

out = (
    root
    / "90_viral_hostclean_validation"
)

with qfile.open() as fh:
    rows = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )

if len(rows) != 6748:
    raise RuntimeError(
        f"Esperaba 6748 filas; hay {len(rows)}"
    )


def num(x):
    x = str(x).strip()

    if x in {"", "NA", "nan", "None"}:
        return None

    return float(x)


comp = [
    num(x["completeness"])
    for x in rows
]

comp = [
    x for x in comp
    if x is not None
]

contam = [
    num(x["contamination"])
    for x in rows
]

contam = [
    x for x in contam
    if x is not None
]

quality = Counter(
    x["checkv_quality"]
    for x in rows
)

miuvig = Counter(
    x["miuvig_quality"]
    for x in rows
)

provirus = Counter(
    x["provirus"]
    for x in rows
)

warnings = Counter()

for row in rows:

    w = str(
        row.get("warnings", "")
    ).strip()

    if w in {"", "NA"}:
        continue

    for item in w.split(";"):
        item = item.strip()

        if item:
            warnings[item] += 1


summary = [
    ("hostclean_sequences", len(rows)),
    ("completeness_nonNA", len(comp)),
    ("completeness_NA", len(rows)-len(comp)),

    (
        "median_completeness_nonNA",
        statistics.median(comp)
        if comp else "NA"
    ),

    (
        "completeness_ge_90",
        sum(x >= 90 for x in comp)
    ),

    (
        "completeness_ge_75",
        sum(x >= 75 for x in comp)
    ),

    (
        "completeness_ge_50",
        sum(x >= 50 for x in comp)
    ),

    (
        "contamination_gt_0",
        sum(x > 0 for x in contam)
    ),

    (
        "contamination_gt_5",
        sum(x > 5 for x in contam)
    ),

    (
        "contamination_gt_10",
        sum(x > 10 for x in contam)
    ),

    (
        "checkv_provirus_yes",
        provirus.get("Yes", 0)
    ),
]


for key, n in sorted(
    quality.items()
):
    summary.append(
        (
            "checkv_quality_"
            + key.lower().replace(
                "-", "_"
            ).replace(" ", "_"),
            n
        )
    )


for key, n in sorted(
    miuvig.items()
):
    summary.append(
        (
            "miuvig_quality_"
            + key.lower().replace(
                "-", "_"
            ).replace(" ", "_"),
            n
        )
    )


summary_file = (
    out
    / "92E3_postclean_checkv_summary.tsv"
)

with summary_file.open("w") as fh:

    fh.write("metric\tvalue\n")

    for k, v in summary:
        fh.write(f"{k}\t{v}\n")


warning_file = (
    out
    / "92E3_postclean_warnings.tsv"
)

with warning_file.open("w") as fh:

    fh.write("warning\tcount\n")

    for k, v in sorted(
        warnings.items(),
        key=lambda x: (
            -x[1],
            x[0]
        )
    ):
        fh.write(f"{k}\t{v}\n")


print("===== POST-CLEAN CHECKV =====")

for k, v in summary:
    print(k, v, sep="\t")

PY


# ======================================================================
# 4. AUDITAR BLAST
# ======================================================================

echo
echo "======================================================================"
echo "3. BLAST PARA vOTU"
echo "======================================================================"

"$APPT" exec \
    "$BLAST_SIF" \
    blastn -version \
    | head -n 3

"$APPT" exec \
    "$BLAST_SIF" \
    makeblastdb -version \
    | head -n 3


# ======================================================================
# 5. LOCALIZAR anicalc.py / aniclust.py
# ======================================================================

echo
echo "======================================================================"
echo "4. SCRIPTS CHECKV DE CLUSTERING"
echo "======================================================================"

AUDIT="${OUT}/92E3_votu_tool_audit.tsv"

echo -e \
"tool\tstatus\tpath\tsha256" \
> "$AUDIT"


for TOOL in anicalc.py aniclust.py
do

    PATH_INSIDE="$(
        "$APPT" exec \
            "$CHECKV_SIF" \
            sh -lc "
                command -v ${TOOL} 2>/dev/null \
                || find /opt /usr /app \
                    -type f \
                    -name '${TOOL}' \
                    2>/dev/null \
                    | head -n 1
            "
    )"

    if [[ -n "$PATH_INSIDE" ]]
    then

        HASH="$(
            "$APPT" exec \
                "$CHECKV_SIF" \
                sh -lc \
                "sha256sum '$PATH_INSIDE' | awk '{print \$1}'"
        )"

        echo -e \
"${TOOL}\tFOUND\t${PATH_INSIDE}\t${HASH}" \
            | tee -a "$AUDIT"

        echo
        echo "===== ${TOOL} HELP ====="

        "$APPT" exec \
            "$CHECKV_SIF" \
            "$PATH_INSIDE" \
            -h \
            2>&1 \
            | head -n 80 \
            || true

    else

        echo -e \
"${TOOL}\tNOT_FOUND\tNA\tNA" \
            | tee -a "$AUDIT"

    fi

done


# ======================================================================
# 6. SMOKE TEST DE CLUSTERING CON 50 VCLEAN
# ======================================================================

echo
echo "======================================================================"
echo "5. SMOKE TEST vOTU"
echo "======================================================================"

ANICALC_PATH="$(
    awk -F '\t' \
        '$1=="anicalc.py" && $2=="FOUND"{print $3}' \
        "$AUDIT"
)"

ANICLUST_PATH="$(
    awk -F '\t' \
        '$1=="aniclust.py" && $2=="FOUND"{print $3}' \
        "$AUDIT"
)"


if [[ -n "$ANICALC_PATH" && -n "$ANICLUST_PATH" ]]
then

    SMALL="${OUT}/first50_vclean.fna"
    DBPREFIX="${OUT}/smoke_db"
    BLASTOUT="${OUT}/smoke_allvsall.tsv"
    ANIOUT="${OUT}/smoke_ani.tsv"
    CLUSTOUT="${OUT}/smoke_votu_clusters.tsv"

    python3 - \
        "$INPUT" \
        "$SMALL" <<'PY'

import sys

src, dst = sys.argv[1:3]

n = 0

with open(src) as inp, open(dst, "w") as out:

    for line in inp:

        if line.startswith(">"):

            n += 1

            if n > 50:
                break

        out.write(line)

print(
    "smoke_sequences",
    min(n, 50)
)

PY

    rm -f \
        "${DBPREFIX}".*

    "$APPT" exec \
        --bind "${ROOT}:/work" \
        "$BLAST_SIF" \
        makeblastdb \
            -in /work/90_viral_hostclean_validation/first50_vclean.fna \
            -dbtype nucl \
            -out /work/90_viral_hostclean_validation/smoke_db

    "$APPT" exec \
        --bind "${ROOT}:/work" \
        "$BLAST_SIF" \
        blastn \
            -task megablast \
            -query /work/90_viral_hostclean_validation/first50_vclean.fna \
            -db /work/90_viral_hostclean_validation/smoke_db \
            -outfmt "6 std qlen slen" \
            -max_target_seqs 10000 \
            -num_threads "$SLURM_CPUS_PER_TASK" \
            -out /work/90_viral_hostclean_validation/smoke_allvsall.tsv

    "$APPT" exec \
        --bind "${ROOT}:/work" \
        "$CHECKV_SIF" \
        "$ANICALC_PATH" \
            -i /work/90_viral_hostclean_validation/smoke_allvsall.tsv \
            -o /work/90_viral_hostclean_validation/smoke_ani.tsv

    "$APPT" exec \
        --bind "${ROOT}:/work" \
        "$CHECKV_SIF" \
        "$ANICLUST_PATH" \
            --fna /work/90_viral_hostclean_validation/first50_vclean.fna \
            --ani /work/90_viral_hostclean_validation/smoke_ani.tsv \
            --out /work/90_viral_hostclean_validation/smoke_votu_clusters.tsv \
            --min_ani 95 \
            --min_tcov 85 \
            --min_qcov 0

    [[ -s "$CLUSTOUT" ]] || {
        echo "ERROR: smoke clustering no produjo clusters." >&2
        exit 7
    }

    echo
    echo "===== SMOKE CLUSTERS ====="

    cat "$CLUSTOUT"

else

    echo "WARNING: scripts de clustering no encontrados."
    echo "El CheckV post-clean sí queda válido."
    echo "No se ejecutó smoke test de vOTU."

fi


# ======================================================================
# 7. METADATA
# ======================================================================

{
    echo -e "field\tvalue"
    echo -e "analysis\tpostclean_CheckV_and_vOTU_tool_audit"
    echo -e "input\t89_viral_catalog_hostclean/92E2_hostclean_exact_representatives.fna"
    echo -e "input_sequences\t6748"
    echo -e "CheckV_version\t1.1.1"
    echo -e "CheckV_DB\tcheckv-db-v1.5"
    echo -e "BLAST_container\tblast-2.16.0.sif"
    echo -e "vOTU_target_ANI\t95"
    echo -e "vOTU_target_AF_shorter\t85"
    echo -e "production_vOTU_clustering\tNOT_YET_PERFORMED"
    echo -e "date\t$(date -Iseconds)"
} > "${OUT}/92E3_run_metadata.tsv"


echo
echo "======================================================================"
echo "RESUMEN CHECKV POST-CLEAN"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92E3_postclean_checkv_summary.tsv"


echo
echo "======================================================================"
echo "AUDITORIA HERRAMIENTAS vOTU"
echo "======================================================================"

column -t -s $'\t' \
    "$AUDIT"


echo
echo "======================================================================"
echo "PASO 92E3 FINALIZÓ CORRECTAMENTE."
echo "Todavía NO se generaron los vOTUs definitivos."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "======================================================================"

