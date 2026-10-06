#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00

#SBATCH --job-name=lcn_read_validation
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"
OUT="${ROOT}/32_lactococcin_read_validation"

ANTISMASH="${ANTISMASH:-/opt/ohpc/pub/containers/BIO/antismash-7.1.0.sif}"

SAMTOOLS=$(find "${CONTAINER_ROOT:-/opt/ohpc/pub/containers/BIO}" \
    -maxdepth 1 \
    -type f \
    -name 'samtools*.sif' \
    -print \
    -quit)

mkdir -p "${OUT}"

echo "============================================================"
echo "VALIDACIÓN POR READS - LACTOCOCCIN"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "============================================================"

if [[ -z "${SAMTOOLS}" ]]; then
    echo "ERROR: no se encontró contenedor samtools." >&2
    exit 1
fi

echo "Samtools container:"
echo "${SAMTOOLS}"

# ============================================================
# BAMs
# ============================================================

M2_BAMS=(
    "${ROOT}/06_mapping/M2/M2_1/M2_1.sorted.bam"
    "${ROOT}/06_mapping/M2/M2_2/M2_2.sorted.bam"
    "${ROOT}/06_mapping/M2/M2_3/M2_3.sorted.bam"
)

L2_BAMS=(
    "${ROOT}/06_mapping/L2/L2_1/L2_1.sorted.bam"
    "${ROOT}/06_mapping/L2/L2_2/L2_2.sorted.bam"
    "${ROOT}/06_mapping/L2/L2_3/L2_3.sorted.bam"
)

for BAM in "${M2_BAMS[@]}" "${L2_BAMS[@]}"; do
    if [[ ! -s "${BAM}" ]]; then
        echo "ERROR: BAM inexistente: ${BAM}" >&2
        exit 1
    fi
done


# ============================================================
# EXTRAER k141_64662 DESDE antiSMASH
# ============================================================

echo
echo "============================================================"
echo "EXTRAYENDO REFERENCIA k141_64662"
echo "============================================================"

apptainer exec \
    --cleanenv \
    --bind "${ROOT}:/project:ro" \
    --bind "${OUT}:/out" \
    "${ANTISMASH}" \
python3 - <<'PY'

from Bio import SeqIO
from pathlib import Path

gbk = Path(
    "/project/26_antismash_rep18/results/"
    "M2__M2_maxbin2.004/"
    "k141_64662.region001.gbk"
)

record = SeqIO.read(
    str(gbk),
    "genbank"
)

if record.id != "k141_64662":
    raise SystemExit(
        f"ERROR: contig inesperado: {record.id}"
    )

if len(record.seq) != 10975:
    raise SystemExit(
        f"ERROR: longitud inesperada: {len(record.seq)}"
    )

out = Path(
    "/out/k141_64662.fasta"
)

with out.open("w") as fh:
    fh.write(
        f">{record.id}\n"
        f"{str(record.seq)}\n"
    )

seq = str(record.seq)

codon = seq[4283:4286]

print(
    "Contig:",
    record.id
)

print(
    "Longitud:",
    len(seq)
)

print(
    "Codón ensamblado 4284-4286:",
    codon
)

print(
    "Contexto 4275-4295:",
    seq[4274:4295]
)

if codon != "TAA":
    raise SystemExit(
        f"ERROR: se esperaba TAA y se encontró {codon}"
    )

PY


# ============================================================
# INDEXAR REFERENCIA
# ============================================================

apptainer exec \
    --cleanenv \
    --bind "${OUT}:/out" \
    "${SAMTOOLS}" \
    samtools faidx \
    /out/k141_64662.fasta


# ============================================================
# COMPROBAR QUE k141_64662 ESTÁ EN LOS BAM M2
# ============================================================

echo
echo "============================================================"
echo "VALIDANDO BAMs M2"
echo "============================================================"

for BAM in "${M2_BAMS[@]}"; do

    SAMPLE=$(basename "${BAM}" .sorted.bam)

    echo
    echo "${SAMPLE}"

    INTERNAL="/project/${BAM#${ROOT}/}"

    LINE=$(
        apptainer exec \
            --cleanenv \
            --bind "${ROOT}:/project:ro" \
            "${SAMTOOLS}" \
            samtools idxstats \
            "${INTERNAL}" \
        | awk '$1=="k141_64662"{print}'
    )

    if [[ -z "${LINE}" ]]; then
        echo "ERROR: k141_64662 no está en ${SAMPLE}" >&2
        exit 1
    fi

    echo "${LINE}"

done


# ============================================================
# PILEUP DEL CODÓN STOP
# ============================================================

echo
echo "============================================================"
echo "PILEUP TAA 4284-4286"
echo "MAPQ >=20; BaseQ >=20"
echo "============================================================"

for BAM in "${M2_BAMS[@]}"; do

    SAMPLE=$(basename "${BAM}" .sorted.bam)
    INTERNAL="/project/${BAM#${ROOT}/}"

    echo "Procesando ${SAMPLE}"

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
            -f /out/k141_64662.fasta \
            -r k141_64662:4284-4286 \
            "${INTERNAL}" \
        > "${OUT}/${SAMPLE}.stop_codon.pileup"

done


# ============================================================
# PARSEAR BASES DEL PILEUP
# ============================================================

apptainer exec \
    --cleanenv \
    --bind "${OUT}:/out" \
    "${ANTISMASH}" \
python3 - <<'PY'

from pathlib import Path
import csv

OUT = Path("/out")


def parse_pileup_bases(bases, ref):

    counts = {
        "A": 0,
        "C": 0,
        "G": 0,
        "T": 0,
        "N": 0,
        "DEL": 0,
    }

    i = 0

    while i < len(bases):

        c = bases[i]

        if c in ".,":
            counts[ref.upper()] += 1
            i += 1
            continue

        if c in "ACGTNacgtn":
            counts[c.upper()] += 1
            i += 1
            continue

        if c == "^":
            # siguiente carácter = calidad de mapeo codificada
            i += 2
            continue

        if c == "$":
            i += 1
            continue

        if c in "+-":

            i += 1

            digits = ""

            while (
                i < len(bases)
                and bases[i].isdigit()
            ):
                digits += bases[i]
                i += 1

            if digits:
                length = int(digits)
                i += length

            continue

        if c == "*":
            counts["DEL"] += 1
            i += 1
            continue

        # símbolos como < >
        i += 1

    return counts


rows = []

for path in sorted(
    OUT.glob("M2_*.stop_codon.pileup")
):

    sample = path.name.split(".")[0]

    observed_positions = set()

    with path.open() as fh:

        for line in fh:

            fields = line.rstrip().split("\t")

            if len(fields) < 5:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[2].upper()
            raw_depth = int(fields[3])
            bases = fields[4]

            counts = parse_pileup_bases(
                bases,
                ref
            )

            usable = sum(
                counts.values()
            )

            observed_positions.add(pos)

            row = {
                "sample": sample,
                "contig": chrom,
                "position": pos,
                "reference": ref,
                "mpileup_depth":
                    raw_depth,
                "count_A":
                    counts["A"],
                "count_C":
                    counts["C"],
                "count_G":
                    counts["G"],
                "count_T":
                    counts["T"],
                "count_N":
                    counts["N"],
                "count_DEL":
                    counts["DEL"],
                "usable_bases":
                    usable,
            }

            for base in "ACGT":

                row[
                    f"freq_{base}"
                ] = (
                    round(
                        counts[base] /
                        usable,
                        6
                    )
                    if usable
                    else 0
                )

            rows.append(row)

    missing = {
        4284, 4285, 4286
    } - observed_positions

    if missing:
        print(
            f"ADVERTENCIA {sample}: "
            f"sin cobertura en {sorted(missing)}"
        )


fields = [
    "sample",
    "contig",
    "position",
    "reference",
    "mpileup_depth",
    "usable_bases",
    "count_A",
    "count_C",
    "count_G",
    "count_T",
    "count_N",
    "count_DEL",
    "freq_A",
    "freq_C",
    "freq_G",
    "freq_T",
]

with (
    OUT / "stop_codon_base_counts.tsv"
).open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)


# Resumen específico de la primera base:
# TAA -> T en 4284
# CAA/CAG -> C en 4284

summary = []

for sample in [
    "M2_1",
    "M2_2",
    "M2_3"
]:

    sr = {
        x["position"]: x
        for x in rows
        if x["sample"] == sample
    }

    p = sr.get(4284)

    if p is None:
        continue

    informative = (
        p["count_T"]
        + p["count_C"]
    )

    summary.append({
        "sample": sample,
        "position": 4284,
        "assembled_base": "T",
        "TAA_T_reads":
            p["count_T"],
        "Q_rescue_C_reads":
            p["count_C"],
        "informative_T_plus_C":
            informative,
        "T_frequency":
            round(
                p["count_T"] /
                informative,
                6
            )
            if informative else 0,
        "C_frequency":
            round(
                p["count_C"] /
                informative,
                6
            )
            if informative else 0,
    })


fields2 = [
    "sample",
    "position",
    "assembled_base",
    "TAA_T_reads",
    "Q_rescue_C_reads",
    "informative_T_plus_C",
    "T_frequency",
    "C_frequency",
]

with (
    OUT / "stop_codon_T_vs_C.tsv"
).open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields2,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(summary)

PY


# ============================================================
# CREAR BED DE LOS DOS MAGS
# ============================================================

echo
echo "============================================================"
echo "CREANDO BEDs DE MAGS"
echo "============================================================"

apptainer exec \
    --cleanenv \
    --bind "${ROOT}:/project:ro" \
    --bind "${OUT}:/out" \
    "${ANTISMASH}" \
python3 - <<'PY'

from Bio import SeqIO
from pathlib import Path

pairs = [
    (
        "/project/22_final_representative_mags/fastas/"
        "M2__M2_maxbin2.004.fa",
        "/out/M2__M2_maxbin2.004.bed",
        "k141_64662",
    ),
    (
        "/project/22_final_representative_mags/fastas/"
        "L2__L2_maxbin2.011_sub.fa",
        "/out/L2__L2_maxbin2.011_sub.bed",
        "k141_100305",
    ),
]

for fasta, bed, target in pairs:

    records = list(
        SeqIO.parse(
            fasta,
            "fasta"
        )
    )

    if not records:
        raise RuntimeError(
            f"FASTA vacío: {fasta}"
        )

    names = {
        r.id
        for r in records
    }

    if target not in names:
        raise RuntimeError(
            f"{target} no está en {fasta}"
        )

    with open(bed, "w") as fh:

        for record in records:

            fh.write(
                f"{record.id}\t"
                f"0\t"
                f"{len(record.seq)}\n"
            )

    print(
        Path(fasta).name,
        "contigs=",
        len(records),
        "target=",
        target
    )

PY


# ============================================================
# COBERTURA POR CONTIG EN CADA MAG
# ============================================================

COV="${OUT}/bin_contig_coverage.tsv"

printf \
"MAG\tsample\tcontig\tlength_bp\tmean_depth\n" \
> "${COV}"


calculate_coverage () {

    MAG="$1"
    BED="$2"
    BAM="$3"

    SAMPLE=$(basename "${BAM}" .sorted.bam)
    INTERNAL="/project/${BAM#${ROOT}/}"

    echo "Cobertura: ${MAG} / ${SAMPLE}"

    apptainer exec \
        --cleanenv \
        --bind "${ROOT}:/project:ro" \
        --bind "${OUT}:/out:ro" \
        "${SAMTOOLS}" \
        samtools depth \
            -aa \
            -b "/out/${BED}" \
            "${INTERNAL}" \
    | awk \
        -v mag="${MAG}" \
        -v sample="${SAMPLE}" \
        'BEGIN{OFS="\t"}
         {
            sum[$1]+=$3;
            n[$1]++
         }
         END{
            for (c in n) {
                print mag,sample,c,n[c],sum[c]/n[c]
            }
         }' \
    >> "${COV}"
}


echo
echo "============================================================"
echo "COBERTURA M2 - MAG Lactococcus lactis"
echo "============================================================"

for BAM in "${M2_BAMS[@]}"; do

    calculate_coverage \
        "M2__M2_maxbin2.004" \
        "M2__M2_maxbin2.004.bed" \
        "${BAM}"

done


echo
echo "============================================================"
echo "COBERTURA L2 - MAG Atlantibacter"
echo "============================================================"

for BAM in "${L2_BAMS[@]}"; do

    calculate_coverage \
        "L2__L2_maxbin2.011_sub" \
        "L2__L2_maxbin2.011_sub.bed" \
        "${BAM}"

done


# ============================================================
# RESUMIR COBERTURA DEL CONTIG OBJETIVO VS SU MAG
# ============================================================

apptainer exec \
    --cleanenv \
    --bind "${OUT}:/out" \
    "${ANTISMASH}" \
python3 - <<'PY'

import csv
import statistics
from collections import defaultdict
from pathlib import Path

OUT = Path("/out")

targets = {
    "M2__M2_maxbin2.004":
        "k141_64662",

    "L2__L2_maxbin2.011_sub":
        "k141_100305",
}

data = defaultdict(list)

with (
    OUT / "bin_contig_coverage.tsv"
).open(newline="") as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in reader:

        row["length_bp"] = int(
            row["length_bp"]
        )

        row["mean_depth"] = float(
            row["mean_depth"]
        )

        data[
            (
                row["MAG"],
                row["sample"],
            )
        ].append(row)


summary = []

for (mag, sample), rows in sorted(
    data.items()
):

    target = targets[mag]

    target_rows = [
        r for r in rows
        if r["contig"] == target
    ]

    if len(target_rows) != 1:
        raise RuntimeError(
            f"{mag}/{sample}: "
            f"target {target} encontrado "
            f"{len(target_rows)} veces"
        )

    tr = target_rows[0]

    means = [
        r["mean_depth"]
        for r in rows
    ]

    total_bp = sum(
        r["length_bp"]
        for r in rows
    )

    weighted_mean = sum(
        r["mean_depth"]
        * r["length_bp"]
        for r in rows
    ) / total_bp

    median_mean = statistics.median(
        means
    )

    target_depth = tr["mean_depth"]

    sorted_means = sorted(
        means
    )

    rank = sum(
        x <= target_depth
        for x in sorted_means
    )

    percentile = (
        100 * rank /
        len(sorted_means)
    )

    summary.append({
        "MAG": mag,
        "sample": sample,
        "target_contig":
            target,
        "n_MAG_contigs":
            len(rows),
        "MAG_total_bp":
            total_bp,
        "target_length_bp":
            tr["length_bp"],
        "target_mean_depth":
            round(
                target_depth,
                4
            ),
        "MAG_weighted_mean_depth":
            round(
                weighted_mean,
                4
            ),
        "MAG_median_contig_depth":
            round(
                median_mean,
                4
            ),
        "target_over_weighted_mean":
            round(
                target_depth /
                weighted_mean,
                4
            )
            if weighted_mean > 0
            else "",
        "target_over_median":
            round(
                target_depth /
                median_mean,
                4
            )
            if median_mean > 0
            else "",
        "target_depth_percentile":
            round(
                percentile,
                2
            ),
    })


fields = [
    "MAG",
    "sample",
    "target_contig",
    "n_MAG_contigs",
    "MAG_total_bp",
    "target_length_bp",
    "target_mean_depth",
    "MAG_weighted_mean_depth",
    "MAG_median_contig_depth",
    "target_over_weighted_mean",
    "target_over_median",
    "target_depth_percentile",
]

with (
    OUT / "target_contig_coverage_summary.tsv"
).open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(summary)

PY


# ============================================================
# MOSTRAR RESULTADOS
# ============================================================

echo
echo "============================================================"
echo "BASES OBSERVADAS EN EL CODÓN STOP"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/stop_codon_base_counts.tsv"

echo
echo "============================================================"
echo "T vs C EN POSICIÓN 4284"
echo "T = mantiene TAA"
echo "C = rescataría el codón hacia CAA/CAG"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/stop_codon_T_vs_C.tsv"

echo
echo "============================================================"
echo "COBERTURA DEL CONTIG OBJETIVO VS MAG"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/target_contig_coverage_summary.tsv"

echo
echo "============================================================"
echo "VALIDACIÓN POR READS COMPLETADA"
echo "Resultados: ${OUT}"
echo "Fin: $(date)"
echo "============================================================"
