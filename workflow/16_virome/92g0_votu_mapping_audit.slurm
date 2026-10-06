#!/bin/bash

#SBATCH --job-name=92g0_votu_audit
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=92g0_votu_audit_%j.out
#SBATCH --error=92g0_votu_audit_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"
BOWTIE2_SIF="${BOWTIE2_SIF:-/opt/ohpc/pub/containers/BIO/bowtie2-2.5.4.sif}"

SOURCE_REF="${ROOT}/91_votu_catalog/92F_votu_representatives.fna"

OUT="${ROOT}/92_votu_abundance"
REFDIR="${OUT}/reference"
INDEXDIR="${OUT}/index"
AUDITDIR="${OUT}/audit"

REF="${REFDIR}/votu_4491.fna"
INDEX="${INDEXDIR}/votu_4491"

mkdir -p \
    "$REFDIR" \
    "$INDEXDIR" \
    "$AUDITDIR"

echo "======================================================================"
echo "PASO 92G0 - AUDITORIA PRE-MAPEO DE vOTUs"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "======================================================================"


# ======================================================================
# 1. VALIDACIONES
# ======================================================================

[[ -x "$APPTAINER" ]] || {
    echo "ERROR: falta Apptainer." >&2
    exit 2
}

[[ -s "$BOWTIE2_SIF" ]] || {
    echo "ERROR: falta Bowtie2 SIF." >&2
    exit 3
}

[[ -s "$SOURCE_REF" ]] || {
    echo "ERROR: falta catálogo vOTU." >&2
    exit 4
}


# ======================================================================
# 2. CONGELAR REFERENCIA PARA CUANTIFICACION
# ======================================================================

echo
echo "======================================================================"
echo "1. REFERENCIA vOTU"
echo "======================================================================"

cp -f \
    "$SOURCE_REF" \
    "$REF"

sha256sum "$REF" \
    > "${REFDIR}/votu_4491.sha256"


python3 - "$REF" "${REFDIR}/votu_4491_stats.tsv" <<'PY'

import statistics
import sys
from pathlib import Path

fasta = Path(sys.argv[1])
out = Path(sys.argv[2])

lengths = []
names = []
name = None
length = 0

with fasta.open() as fh:

    for line in fh:

        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):

            if name is not None:
                names.append(name)
                lengths.append(length)

            name = line[1:].split()[0]
            length = 0

        else:
            length += len(line)

    if name is not None:
        names.append(name)
        lengths.append(length)


if len(names) != 4491:
    raise RuntimeError(
        f"Esperaba 4491 vOTUs; encontré {len(names)}"
    )

if len(names) != len(set(names)):
    raise RuntimeError(
        "Hay IDs duplicados en referencia."
    )

if set(names) != {
    f"VOTU{i:05d}"
    for i in range(1, 4492)
}:
    raise RuntimeError(
        "IDs VOTU no corresponden a VOTU00001..VOTU04491"
    )


with out.open("w") as fh:

    fh.write("metric\tvalue\n")
    fh.write(f"sequences\t{len(names)}\n")
    fh.write(f"total_bp\t{sum(lengths)}\n")
    fh.write(f"min_length\t{min(lengths)}\n")
    fh.write(
        f"median_length\t{statistics.median(lengths)}\n"
    )
    fh.write(f"max_length\t{max(lengths)}\n")
    fh.write(
        f"length_ge_2500\t"
        f"{sum(x >= 2500 for x in lengths)}\n"
    )
    fh.write(
        f"length_ge_5000\t"
        f"{sum(x >= 5000 for x in lengths)}\n"
    )
    fh.write(
        f"length_ge_10000\t"
        f"{sum(x >= 10000 for x in lengths)}\n"
    )
    fh.write(
        f"length_ge_20000\t"
        f"{sum(x >= 20000 for x in lengths)}\n"
    )
    fh.write(
        f"length_ge_50000\t"
        f"{sum(x >= 50000 for x in lengths)}\n"
    )


print("reference_sequences:", len(names))
print("total_bp:", sum(lengths))
print("median_length:", statistics.median(lengths))
print("max_length:", max(lengths))

PY


column -t -s $'\t' \
    "${REFDIR}/votu_4491_stats.tsv"


# ======================================================================
# 3. CONSTRUIR INDICE BOWTIE2
# ======================================================================

echo
echo "======================================================================"
echo "2. INDICE BOWTIE2"
echo "======================================================================"

index_complete () {

    local prefix="$1"

    for suffix in \
        1 \
        2 \
        3 \
        4 \
        rev.1 \
        rev.2
    do

        if [[ ! -s "${prefix}.${suffix}.bt2" \
           && ! -s "${prefix}.${suffix}.bt2l" ]]
        then
            return 1
        fi

    done

    return 0
}


if index_complete "$INDEX"
then

    echo "Índice completo preexistente."
    echo "No se reconstruirá."

else

    echo "Construyendo índice..."

    rm -f \
        "${INDEX}".*.bt2 \
        "${INDEX}".*.bt2l

    "$APPTAINER" exec \
        --bind "${ROOT}:/work" \
        "$BOWTIE2_SIF" \
        bowtie2-build \
            /work/92_votu_abundance/reference/votu_4491.fna \
            /work/92_votu_abundance/index/votu_4491

fi


index_complete "$INDEX" || {
    echo "ERROR: índice Bowtie2 incompleto." >&2
    exit 5
}


echo
echo "Bowtie2:"

"$APPTAINER" exec \
    "$BOWTIE2_SIF" \
    bowtie2 --version \
    | sed -n '1,3p'


echo
echo "Secuencias reconocidas por bowtie2-inspect:"

NINDEX="$(
    "$APPTAINER" exec \
        --bind "${ROOT}:/work" \
        "$BOWTIE2_SIF" \
        bowtie2-inspect \
            -n \
            /work/92_votu_abundance/index/votu_4491 \
        | wc -l
)"

echo "index_sequences=${NINDEX}"

[[ "$NINDEX" -eq 4491 ]] || {
    echo "ERROR: índice contiene ${NINDEX}, esperaba 4491." >&2
    exit 6
}


# ======================================================================
# 4. AUDITAR LOS 18 PARES HOST-FREE
# ======================================================================

echo
echo "======================================================================"
echo "3. AUDITORIA DE READS HOST-FREE"
echo "======================================================================"

python3 - <<'PY'

import os
import re
from pathlib import Path

root = Path(
    f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"
)

base = root / "04_host_removed"

out = (
    root
    / "92_votu_abundance"
    / "audit"
    / "92G0_hostfree_pair_audit.tsv"
)

samples = [
    "L1_1", "L1_2", "L1_3",
    "L2_1", "L2_2", "L2_3",
    "L3_1", "L3_2", "L3_3",
    "M1_1", "M1_2", "M1_3",
    "M2_1", "M2_2", "M2_3",
    "M3_1", "M3_2", "M3_3",
]


def classify(sample, path):

    name = path.name.lower()

    # Quitamos el ID de muestra para evitar que,
    # por ejemplo, L1_1 sea interpretado como read 1.
    sample_lower = sample.lower()

    if sample_lower in name:
        name = name.replace(
            sample_lower,
            "",
            1
        )

    patterns_r1 = [
        r"(?:^|[._-])r1(?:[._-]|$)",
        r"(?:^|[._-])read1(?:[._-]|$)",
        r"(?:^|[._-])1(?:[._-]|$)",
    ]

    patterns_r2 = [
        r"(?:^|[._-])r2(?:[._-]|$)",
        r"(?:^|[._-])read2(?:[._-]|$)",
        r"(?:^|[._-])2(?:[._-]|$)",
    ]

    is_r1 = any(
        re.search(x, name)
        for x in patterns_r1
    )

    is_r2 = any(
        re.search(x, name)
        for x in patterns_r2
    )

    if is_r1 and not is_r2:
        return "R1"

    if is_r2 and not is_r1:
        return "R2"

    return "UNCLASSIFIED"


rows = []

for sample in samples:

    directory = base / sample

    candidates = []

    if directory.is_dir():

        candidates = sorted(
            p
            for p in directory.rglob("*.fastq.gz")
            if "hostfree" in p.name.lower()
        )

    roles = {
        "R1": [],
        "R2": [],
        "UNCLASSIFIED": [],
    }

    for p in candidates:

        roles[
            classify(
                sample,
                p
            )
        ].append(p)


    if (
        len(roles["R1"]) == 1
        and len(roles["R2"]) == 1
        and len(roles["UNCLASSIFIED"]) == 0
    ):

        status = "READY"

        r1 = roles["R1"][0]
        r2 = roles["R2"][0]

    else:

        status = "REVIEW"

        r1 = (
            roles["R1"][0]
            if len(roles["R1"]) == 1
            else None
        )

        r2 = (
            roles["R2"][0]
            if len(roles["R2"]) == 1
            else None
        )


    def pstr(p):
        return str(p) if p else ""


    def gib(p):

        if not p:
            return ""

        return f"{p.stat().st_size / (1024**3):.3f}"


    rows.append({
        "sample":
            sample,

        "sample_dir_exists":
            "YES"
            if directory.is_dir()
            else "NO",

        "hostfree_candidates":
            len(candidates),

        "status":
            status,

        "R1":
            pstr(r1),

        "R1_GiB":
            gib(r1),

        "R2":
            pstr(r2),

        "R2_GiB":
            gib(r2),

        "unclassified_count":
            len(
                roles[
                    "UNCLASSIFIED"
                ]
            ),

        "candidate_files":
            ";".join(
                str(x)
                for x in candidates
            ),
    })


fields = list(
    rows[0].keys()
)

with out.open("w") as fh:

    fh.write(
        "\t".join(fields)
        + "\n"
    )

    for row in rows:

        fh.write(
            "\t".join(
                str(
                    row[x]
                )
                for x in fields
            )
            + "\n"
        )


ready = sum(
    x["status"] == "READY"
    for x in rows
)

review = len(rows) - ready

print("samples_expected:", len(rows))
print("samples_READY:", ready)
print("samples_REVIEW:", review)

PY


echo
echo "===== PARES DETECTADOS ====="

column -t -s $'\t' \
    "${AUDITDIR}/92G0_hostfree_pair_audit.tsv"


# ======================================================================
# 5. RESUMEN
# ======================================================================

READY="$(
    awk -F '\t' \
        'NR>1 && $4=="READY"{n++} END{print n+0}' \
        "${AUDITDIR}/92G0_hostfree_pair_audit.tsv"
)"

REVIEW="$(
    awk -F '\t' \
        'NR>1 && $4=="REVIEW"{n++} END{print n+0}' \
        "${AUDITDIR}/92G0_hostfree_pair_audit.tsv"
)"


{
    echo -e "metric\tvalue"
    echo -e "reference_vOTUs\t4491"
    echo -e "bowtie2_index_sequences\t${NINDEX}"
    echo -e "samples_expected\t18"
    echo -e "samples_READY\t${READY}"
    echo -e "samples_REVIEW\t${REVIEW}"
} > "${AUDITDIR}/92G0_summary.tsv"


{
    echo -e "field\tvalue"
    echo -e "analysis\tvOTU_mapping_reference_and_read_audit"
    echo -e "reference\t91_votu_catalog/92F_votu_representatives.fna"
    echo -e "reference_sequences\t4491"
    echo -e "bowtie2\t2.5.4"
    echo -e "mapping_performed\tNO"
    echo -e "date\t$(date -Iseconds)"
} > "${AUDITDIR}/92G0_run_metadata.tsv"


echo
echo "======================================================================"
echo "RESUMEN 92G0"
echo "======================================================================"

column -t -s $'\t' \
    "${AUDITDIR}/92G0_summary.tsv"


echo
echo "======================================================================"
echo "PASO 92G0 FINALIZÓ CORRECTAMENTE."
echo "Todavía NO se mapearon reads."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "======================================================================"

