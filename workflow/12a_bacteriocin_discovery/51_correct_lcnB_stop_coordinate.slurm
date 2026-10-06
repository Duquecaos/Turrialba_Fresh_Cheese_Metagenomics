#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

#SBATCH --job-name=lcnB_stop_correct
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"
OUT="${ROOT}/33_lcnB_stop_corrected"

ANTISMASH="${ANTISMASH:-/opt/ohpc/pub/containers/BIO/antismash-7.1.0.sif}"
SAMTOOLS="${SAMTOOLS:-/opt/ohpc/pub/containers/BIO/samtools-1.21.sif}"

mkdir -p "${OUT}"

echo "============================================================"
echo "CORRECCIÓN COORDENADA lcnB-LIKE"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "============================================================"

# ============================================================
# LOCALIZAR REGIÓN antiSMASH EN EL CONTIG COMPLETO
# ============================================================

srun --exact -n 1 -c 1 \
apptainer exec \
    --cleanenv \
    --bind "${ROOT}:/project:ro" \
    --bind "${OUT}:/out" \
    "${ANTISMASH}" \
python3 - <<'PY'

from Bio import SeqIO
from pathlib import Path

ROOT = Path("/project")
OUT = Path("/out")

MAG = (
    ROOT /
    "22_final_representative_mags/fastas/"
    "M2__M2_maxbin2.004.fa"
)

GBK = (
    ROOT /
    "26_antismash_rep18/results/"
    "M2__M2_maxbin2.004/"
    "k141_64662.region001.gbk"
)

TARGET = "k141_64662"

# Coordenadas dentro del archivo region001.gbk
REGION_STOP_START = 4284
REGION_STOP_END = 4286

full_record = None

for rec in SeqIO.parse(str(MAG), "fasta"):
    if rec.id == TARGET:
        full_record = rec
        break

if full_record is None:
    raise SystemExit(
        f"ERROR: no se encontró {TARGET} en {MAG}"
    )

region_record = SeqIO.read(
    str(GBK),
    "genbank"
)

full = str(full_record.seq).upper()
region = str(region_record.seq).upper()

print("===== LONGITUDES =====")
print("Contig completo:", len(full))
print("Región antiSMASH:", len(region))
print("Diferencia:", len(full) - len(region))

# Buscar región exacta dentro del contig completo
positions = []

start = 0

while True:
    pos = full.find(region, start)

    if pos == -1:
        break

    positions.append(pos)
    start = pos + 1

print()
print("===== LOCALIZACIÓN DE LA REGIÓN =====")
print("Matches exactos:", len(positions))

if len(positions) != 1:
    raise SystemExit(
        "ERROR: la región antiSMASH no se encontró "
        "exactamente una vez en el contig completo."
    )

offset0 = positions[0]

full_region_start = offset0 + 1
full_region_end = offset0 + len(region)

print(
    "Región antiSMASH corresponde a contig completo:",
    f"{full_region_start}-{full_region_end}"
)

# Convertir coordenada 1-based de región a contig completo
full_stop_start = (
    offset0 + REGION_STOP_START
)

full_stop_end = (
    offset0 + REGION_STOP_END
)

codon = full[
    full_stop_start - 1:
    full_stop_end
]

print()
print("===== STOP CORREGIDO =====")
print(
    "Coordenada región:",
    f"{REGION_STOP_START}-{REGION_STOP_END}"
)
print(
    "Coordenada contig completo:",
    f"{full_stop_start}-{full_stop_end}"
)
print(
    "Codón contig completo:",
    codon
)

if codon != "TAA":
    raise SystemExit(
        f"ERROR: el codón corregido no es TAA: {codon}"
    )

context_start = max(
    1,
    full_stop_start - 15
)

context_end = min(
    len(full),
    full_stop_end + 15
)

print(
    "Contexto:",
    full[
        context_start - 1:
        context_end
    ]
)

# Guardar contig completo
with (
    OUT / "k141_64662.full.fasta"
).open("w") as fh:

    fh.write(
        f">{TARGET}\n"
        f"{full}\n"
    )

with (
    OUT / "corrected_stop_coordinates.tsv"
).open("w") as fh:

    fh.write(
        "contig\t"
        "full_contig_length\t"
        "antismash_region_length\t"
        "antismash_region_start_full\t"
        "antismash_region_end_full\t"
        "stop_start_region\t"
        "stop_end_region\t"
        "stop_start_full\t"
        "stop_end_full\t"
        "codon\n"
    )

    fh.write(
        f"{TARGET}\t"
        f"{len(full)}\t"
        f"{len(region)}\t"
        f"{full_region_start}\t"
        f"{full_region_end}\t"
        f"{REGION_STOP_START}\t"
        f"{REGION_STOP_END}\t"
        f"{full_stop_start}\t"
        f"{full_stop_end}\t"
        f"{codon}\n"
    )

with (
    OUT / "stop_position.txt"
).open("w") as fh:
    fh.write(str(full_stop_start) + "\n")

PY


# ============================================================
# INDEXAR CONTIG COMPLETO
# ============================================================

apptainer exec \
    --cleanenv \
    --bind "${OUT}:/out" \
    "${SAMTOOLS}" \
    samtools faidx \
    /out/k141_64662.full.fasta


POS=$(cat "${OUT}/stop_position.txt")

echo
echo "Posición real primera base del STOP: ${POS}"


# ============================================================
# VALIDAR LONGITUD EN BAMs Y HACER PILEUP CORRECTO
# ============================================================

printf \
"sample\tcontig\tposition\treference\tdepth\tA\tC\tG\tT\tT_frequency\tC_frequency\n" \
> "${OUT}/corrected_stop_pileup.tsv"

for SAMPLE in M2_1 M2_2 M2_3; do

    BAM="${ROOT}/06_mapping/M2/${SAMPLE}/${SAMPLE}.sorted.bam"

    INTERNAL="/project/${BAM#${ROOT}/}"

    echo
    echo "===== ${SAMPLE} ====="

    IDX=$(
        apptainer exec \
            --cleanenv \
            --bind "${ROOT}:/project:ro" \
            "${SAMTOOLS}" \
            samtools idxstats \
            "${INTERNAL}" \
        | awk '$1=="k141_64662"{print}'
    )

    echo "idxstats: ${IDX}"

    BAM_LEN=$(echo "${IDX}" | awk '{print $2}')

    if [[ "${BAM_LEN}" -ne 13106 ]]; then
        echo \
        "ERROR: longitud BAM inesperada: ${BAM_LEN}" >&2
        exit 1
    fi

    apptainer exec \
        --cleanenv \
        --bind "${ROOT}:/project:ro" \
        --bind "${OUT}:/out:ro" \
        "${SAMTOOLS}" \
        samtools mpileup \
            -B \
            -q 20 \
            -Q 20 \
            -d 1000000 \
            -f /out/k141_64662.full.fasta \
            -r "k141_64662:${POS}-${POS}" \
            "${INTERNAL}" \
        > "${OUT}/${SAMPLE}.corrected_stop.pileup"

done


# ============================================================
# CONTAR BASES
# ============================================================

srun --exact -n 1 -c 1 \
apptainer exec \
    --cleanenv \
    --bind "${OUT}:/out" \
    "${ANTISMASH}" \
python3 - <<'PY'

from pathlib import Path
import csv

OUT = Path("/out")


def parse_bases(text, ref):

    counts = {
        "A": 0,
        "C": 0,
        "G": 0,
        "T": 0,
    }

    i = 0

    while i < len(text):

        c = text[i]

        if c in ".,":
            counts[ref] += 1
            i += 1
            continue

        if c.upper() in counts:
            counts[c.upper()] += 1
            i += 1
            continue

        if c == "^":
            i += 2
            continue

        if c == "$":
            i += 1
            continue

        if c in "+-":

            i += 1
            digits = ""

            while (
                i < len(text)
                and text[i].isdigit()
            ):
                digits += text[i]
                i += 1

            if digits:
                i += int(digits)

            continue

        i += 1

    return counts


rows = []

for sample in [
    "M2_1",
    "M2_2",
    "M2_3"
]:

    path = (
        OUT /
        f"{sample}.corrected_stop.pileup"
    )

    if not path.exists():
        continue

    lines = [
        x.rstrip("\n")
        for x in path.open()
        if x.strip()
    ]

    if not lines:

        rows.append({
            "sample": sample,
            "contig": "k141_64662",
            "position": "",
            "reference": "",
            "depth": 0,
            "A": 0,
            "C": 0,
            "G": 0,
            "T": 0,
            "T_frequency": "",
            "C_frequency": "",
        })

        continue

    f = lines[0].split("\t")

    contig = f[0]
    pos = int(f[1])
    ref = f[2].upper()
    depth = int(f[3])
    bases = f[4]

    c = parse_bases(
        bases,
        ref
    )

    informative = (
        c["T"] + c["C"]
    )

    rows.append({
        "sample": sample,
        "contig": contig,
        "position": pos,
        "reference": ref,
        "depth": depth,
        "A": c["A"],
        "C": c["C"],
        "G": c["G"],
        "T": c["T"],
        "T_frequency":
            round(
                c["T"] / informative,
                6
            )
            if informative else "",
        "C_frequency":
            round(
                c["C"] / informative,
                6
            )
            if informative else "",
    })


fields = [
    "sample",
    "contig",
    "position",
    "reference",
    "depth",
    "A",
    "C",
    "G",
    "T",
    "T_frequency",
    "C_frequency",
]

with (
    OUT /
    "corrected_stop_pileup.tsv"
).open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)

PY


echo
echo "============================================================"
echo "COORDENADA CORREGIDA"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/corrected_stop_coordinates.tsv"

echo
echo "============================================================"
echo "PILEUP CORREGIDO"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/corrected_stop_pileup.tsv"

echo
echo "============================================================"
echo "VALIDACIÓN lcnB STOP CORREGIDA COMPLETADA"
echo "Resultados: ${OUT}"
echo "Fin: $(date)"
echo "============================================================"
