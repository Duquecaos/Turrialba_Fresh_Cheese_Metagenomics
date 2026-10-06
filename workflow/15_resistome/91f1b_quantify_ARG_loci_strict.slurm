#!/bin/bash

#SBATCH --job-name=91f1b_argcov
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=08:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=91f1b_argcov_%j.out
#SBATCH --error=91f1b_argcov_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

LOCI="${ROOT}/80_resistome_community_consolidated/91E_community_AMR_loci.tsv"

OUT="${ROOT}/81b_resistome_ARG_abundance_strict"

SAMTOOLS_SIF="${SAMTOOLS_SIF:-/opt/ohpc/pub/containers/BIO/samtools-1.21.sif}"
APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"

mkdir -p "$OUT"

echo "======================================================================"
echo "PASO 91F1B - CUANTIFICACION ESTRICTA DE LOCI ARG"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "======================================================================"

[[ -s "$LOCI" ]] || {
    echo "ERROR: falta tabla de loci: $LOCI" >&2
    exit 2
}

[[ -s "$SAMTOOLS_SIF" ]] || {
    echo "ERROR: falta contenedor samtools." >&2
    exit 3
}

# ----------------------------------------------------------------------
# Comprobar BAMs
# ----------------------------------------------------------------------

echo
echo "======================================================================"
echo "1. VALIDANDO LOS 18 BAM"
echo "======================================================================"

for G in L1 L2 L3 M1 M2 M3
do
    for T in 1 2 3
    do
        SAMPLE="${G}_${T}"
        BAM="${ROOT}/06_mapping/${G}/${SAMPLE}/${SAMPLE}.sorted.bam"

        [[ -s "$BAM" ]] || {
            echo "ERROR: falta BAM: $BAM" >&2
            exit 4
        }

        "$APPTAINER" exec \
            --bind "${ROOT}:/work" \
            "$SAMTOOLS_SIF" \
            samtools quickcheck \
            "/work/06_mapping/${G}/${SAMPLE}/${SAMPLE}.sorted.bam" \
        || {
            echo "ERROR: BAM inválido: $BAM" >&2
            exit 5
        }

        echo "OK	${SAMPLE}"
    done
done

# ----------------------------------------------------------------------
# Convertir loci a BED
#
# BED es 0-based, half-open.
# ----------------------------------------------------------------------

echo
echo "======================================================================"
echo "2. GENERANDO BED POR COENSAMBLAJE"
echo "======================================================================"

python3 - <<'PY'
import csv
import os
from pathlib import Path

root = Path(
    f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"
)

src = (
    root
    / "80_resistome_community_consolidated"
    / "91E_community_AMR_loci.tsv"
)

out = (
    root
    / "81_resistome_ARG_abundance"
)

handles = {}

try:
    for g in ["L1","L2","L3","M1","M2","M3"]:
        handles[g] = (
            out / f"{g}.ARG_loci.bed"
        ).open("w")

    with src.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")

        n = 0

        for row in r:
            group = row["coassembly"]

            if group not in handles:
                raise RuntimeError(
                    f"Coassembly inesperado: {group}"
                )

            start = int(row["start"])
            stop = int(row["stop"])

            lo = min(start, stop)
            hi = max(start, stop)

            # BED: start-1, stop exclusivo
            handles[group].write(
                "\t".join([
                    row["contig"],
                    str(lo - 1),
                    str(hi),
                    row["community_locus_id"],
                ])
                + "\n"
            )

            n += 1

finally:
    for fh in handles.values():
        fh.close()

print("ARG loci written:", n)

for g in ["L1","L2","L3","M1","M2","M3"]:
    p = out / f"{g}.ARG_loci.bed"

    with p.open() as fh:
        n = sum(1 for _ in fh)

    print(g, n, sep="\t")
PY

# ----------------------------------------------------------------------
# Información de normalización por BAM
# ----------------------------------------------------------------------

echo
echo "======================================================================"
echo "3. ESTADISTICAS DE LOS BAM"
echo "======================================================================"

BAMSTATS="${OUT}/91F_BAM_normalization.tsv"

echo -e \
"sample\tcoassembly\tproducer\tbiological_unit\tweek\teligible_primary_mapped_records" \
> "$BAMSTATS"

for G in L1 L2 L3 M1 M2 M3
do
    for T in 1 2 3
    do
        SAMPLE="${G}_${T}"

        BAM_HOST="${ROOT}/06_mapping/${G}/${SAMPLE}/${SAMPLE}.sorted.bam"
        BAM_CONT="/work/06_mapping/${G}/${SAMPLE}/${SAMPLE}.sorted.bam"

        # -F 3844 elimina:
        # 4    unmapped
        # 256  secondary
        # 2048 supplementary
        N=$(
            "$APPTAINER" exec \
                --bind "${ROOT}:/work" \
                "$SAMTOOLS_SIF" \
                samtools view \
                    -c \
                    -F 3844 \
                    "$BAM_CONT"
        )

        PRODUCER="${G:0:1}"
        BIO="${G:1}"

        # _1 = semana 0
        # _2 = semana 1
        # _3 = semana 2
        WEEK=$((T - 1))

        echo -e \
"${SAMPLE}\t${G}\t${PRODUCER}\t${BIO}\t${WEEK}\t${N}" \
>> "$BAMSTATS"

        echo "${SAMPLE}: primary mapped records = ${N}"
    done
done

# ----------------------------------------------------------------------
# Cuantificación
# ----------------------------------------------------------------------

echo
echo "======================================================================"
echo "4. COBERTURA POR LOCUS"
echo "======================================================================"

RAW="${OUT}/91F_ARG_coverage_long.tsv"

echo -e \
"community_locus_id\tcoassembly\tsample\tproducer\tbiological_unit\tweek\tcontig\tstart\tstop\tlength_bp\tdepth_mode\tbreadth_1x\tbreadth_5x\tbreadth_10x\tmean_depth\tmedian_depth\tmax_depth\teligible_primary_mapped_records\tmean_depth_per_million_mapped" \
> "$RAW"

for G in L1 L2 L3 M1 M2 M3
do

    BED="${OUT}/${G}.ARG_loci.bed"

    for T in 1 2 3
    do
        SAMPLE="${G}_${T}"

        BAM_CONT="/work/06_mapping/${G}/${SAMPLE}/${SAMPLE}.sorted.bam"

        PRODUCER="${G:0:1}"
        BIO="${G:1}"
        WEEK=$((T - 1))

        TOTAL=$(
            awk -F'\t' \
                -v S="$SAMPLE" \
                '$1==S {print $6}' \
                "$BAMSTATS"
        )

        echo
        echo ">>> ${SAMPLE}"

        SAMPLE="$SAMPLE" \
        GROUP="$G" \
        PRODUCER="$PRODUCER" \
        BIO="$BIO" \
        WEEK="$WEEK" \
        TOTAL="$TOTAL" \
        ROOT="$ROOT" \
        OUT="$OUT" \
        APPTAINER="$APPTAINER" \
        SAMTOOLS_SIF="$SAMTOOLS_SIF" \
        python3 - <<'PY'
import csv
import os
import statistics
import subprocess
from pathlib import Path

root = Path(os.environ["ROOT"])
out = Path(os.environ["OUT"])

group = os.environ["GROUP"]
sample = os.environ["SAMPLE"]
producer = os.environ["PRODUCER"]
bio = os.environ["BIO"]
week = os.environ["WEEK"]

total = int(os.environ["TOTAL"])

apptainer = os.environ["APPTAINER"]
sif = os.environ["SAMTOOLS_SIF"]

loci_file = (
    root
    / "80_resistome_community_consolidated"
    / "91E_community_AMR_loci.tsv"
)

bam = (
    f"/work/06_mapping/{group}/"
    f"{sample}/{sample}.sorted.bam"
)

rows = []

with loci_file.open() as fh:
    r = csv.DictReader(fh, delimiter="\t")

    rows = [
        x for x in r
        if x["coassembly"] == group
    ]

output_rows = []

# ----------------------------------------------------------
# Se calculan dos modos:
#
# ALL_PRIMARY:
#   alineamientos primarios; conserva multimapping.
#
# MAPQ20:
#   soporte más conservador; mapping quality >=20.
#
# No usamos uno para reemplazar al otro.
# ----------------------------------------------------------

for row in rows:

    contig = row["contig"]

    start = min(
        int(row["start"]),
        int(row["stop"])
    )

    stop = max(
        int(row["start"]),
        int(row["stop"])
    )

    length = stop - start + 1

    region = (
        f"{contig}:{start}-{stop}"
    )

    for mode, mapq in [
        ("ALL_PRIMARY", 0),
        ("MAPQ20", 20),
    ]:

        cmd = [
            apptainer,
            "exec",
            "--bind",
            f"{root}:/work",
            sif,
            "samtools",
            "depth",
            "-a",
            "-s",
            "-G",
            "UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY",
            "-Q",
            str(mapq),
            "-r",
            region,
            bam,
        ]

        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )

        depths = [
            int(line.split("\t")[2])
            for line in result.stdout.splitlines()
            if line.strip()
        ]

        # samtools depth -a con una región debe
        # devolver las posiciones de la región,
        # pero rellenamos defensivamente.
        if len(depths) < length:
            depths.extend(
                [0] * (length - len(depths))
            )

        if len(depths) > length:
            raise RuntimeError(
                f"{sample} {region}: "
                f"{len(depths)} posiciones > {length}"
            )

        covered1 = sum(
            d >= 1 for d in depths
        )

        covered5 = sum(
            d >= 5 for d in depths
        )

        covered10 = sum(
            d >= 10 for d in depths
        )

        mean_depth = (
            sum(depths) / length
            if length
            else 0
        )

        median_depth = (
            statistics.median(depths)
            if depths
            else 0
        )

        max_depth = (
            max(depths)
            if depths
            else 0
        )

        breadth1 = covered1 / length
        breadth5 = covered5 / length
        breadth10 = covered10 / length

        normalized = (
            mean_depth
            / (total / 1_000_000)
            if total > 0
            else 0
        )

        output_rows.append([
            row["community_locus_id"],
            group,
            sample,
            producer,
            bio,
            week,
            contig,
            start,
            stop,
            length,
            mode,
            breadth1,
            breadth5,
            breadth10,
            mean_depth,
            median_depth,
            max_depth,
            total,
            normalized,
        ])

outfile = (
    out
    / f"{sample}.ARG_coverage.tsv"
)

header = [
    "community_locus_id",
    "coassembly",
    "sample",
    "producer",
    "biological_unit",
    "week",
    "contig",
    "start",
    "stop",
    "length_bp",
    "depth_mode",
    "breadth_1x",
    "breadth_5x",
    "breadth_10x",
    "mean_depth",
    "median_depth",
    "max_depth",
    "eligible_primary_mapped_records",
    "mean_depth_per_million_mapped",
]

with outfile.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(header)
    w.writerows(output_rows)

print(
    sample,
    "loci=",
    len(rows),
    "coverage_rows=",
    len(output_rows),
)
PY

    done
done

# ----------------------------------------------------------------------
# Consolidar las 18 tablas
# ----------------------------------------------------------------------

python3 - <<'PY'
import csv
import os
from pathlib import Path

root = Path(
    f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"
)

out = root / "81_resistome_ARG_abundance"

samples = [
    f"{g}_{t}"
    for g in ["L1","L2","L3","M1","M2","M3"]
    for t in [1,2,3]
]

dest = out / "91F_ARG_coverage_long.tsv"

header = None
rows = []

for sample in samples:

    f = out / f"{sample}.ARG_coverage.tsv"

    if not f.exists():
        raise RuntimeError(
            f"Falta {f}"
        )

    with f.open() as fh:
        r = csv.reader(fh, delimiter="\t")

        this_header = next(r)

        if header is None:
            header = this_header

        elif this_header != header:
            raise RuntimeError(
                f"Headers incompatibles en {sample}"
            )

        rows.extend(r)

with dest.open("w", newline="") as fh:

    w = csv.writer(
        fh,
        delimiter="\t"
    )

    w.writerow(header)
    w.writerows(rows)

print("consolidated_rows:", len(rows))
print("output:", dest)
PY

echo
echo "======================================================================"
echo "5. INTEGRAR ANOTACION + COBERTURA"
echo "======================================================================"

python3 - <<'PY'
import csv
import os
from pathlib import Path
from collections import Counter

root = Path(
    f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"
)

base = root / "81_resistome_ARG_abundance"

loci_file = (
    root
    / "80_resistome_community_consolidated"
    / "91E_community_AMR_loci.tsv"
)

coverage_file = (
    base
    / "91F_ARG_coverage_long.tsv"
)

out_file = (
    base
    / "91F_ARG_annotation_coverage.tsv"
)

summary_file = (
    base
    / "91F_coverage_summary.tsv"
)

with loci_file.open() as fh:
    r = csv.DictReader(fh, delimiter="\t")
    loci = {
        x["community_locus_id"]: x
        for x in r
    }

with coverage_file.open() as fh:
    r = csv.DictReader(fh, delimiter="\t")
    cov = list(r)

coverage_fields = list(r.fieldnames)

annotation_fields = [
    "element_symbol",
    "element_name",
    "scope",
    "class",
    "subclass",
    "method",
    "coverage_reference",
    "identity_reference",
    "closest_reference_accession",
    "HMM_accession",
    "protein_sha256",
    "nucleotide_sha256",
    "representation_vs_final18",
]

fields = (
    coverage_fields
    + annotation_fields
    + [
        "structural_flag",
    ]
)

result = []

for row in cov:

    locus = loci[
        row["community_locus_id"]
    ]

    method = locus["method"]

    if (
        "PARTIAL" in method
        or "INTERNAL_STOP" in method
    ):
        structural = (
            "partial_or_interrupted_detection"
        )
    else:
        structural = (
            "no_partial_or_internal_stop_flag"
        )

    x = dict(row)

    for field in annotation_fields:
        x[field] = locus[field]

    x["structural_flag"] = structural

    result.append(x)

with out_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(result)

# ----------------------------------------------------------
# Resumen descriptivo, SIN imponer presencia/ausencia final
# ----------------------------------------------------------

metrics = []

for mode in [
    "ALL_PRIMARY",
    "MAPQ20",
]:

    subset = [
        x for x in result
        if x["depth_mode"] == mode
    ]

    metrics.extend([
        (
            f"{mode}_locus_sample_records",
            len(subset)
        ),
        (
            f"{mode}_breadth_ge_0.50",
            sum(
                float(x["breadth_1x"]) >= 0.50
                for x in subset
            )
        ),
        (
            f"{mode}_breadth_ge_0.80",
            sum(
                float(x["breadth_1x"]) >= 0.80
                for x in subset
            )
        ),
        (
            f"{mode}_breadth_ge_0.95",
            sum(
                float(x["breadth_1x"]) >= 0.95
                for x in subset
            )
        ),
        (
            f"{mode}_breadth5_ge_0.80",
            sum(
                float(x["breadth_5x"]) >= 0.80
                for x in subset
            )
        ),
    ])

with summary_file.open("w") as fh:

    fh.write("metric\tvalue\n")

    for k, v in metrics:
        fh.write(
            f"{k}\t{v}\n"
        )

print("integrated_records:", len(result))

for k, v in metrics:
    print(k, v, sep="\t")
PY

# ----------------------------------------------------------------------
# Provenance
# ----------------------------------------------------------------------

{
    echo -e "field\tvalue"
    echo -e "analysis\tARG_locus_coverage_existing_coassembly_BAMs_strict"
    echo -e "mapping_source\t06_mapping"
    echo -e "mapping_reference\t05_assembly/{coassembly}/final.contigs.fa"
    echo -e "samtools\t1.21"
    echo -e "depth_modes\tALL_PRIMARY;MAPQ20"
    echo -e "normalization\tmean_depth_per_million_eligible_primary_mapped_records"
    echo -e "presence_threshold\tNOT_IMPOSED"
    echo -e "date\t$(date -Iseconds)"
} > "${OUT}/91F_run_metadata.tsv"

echo
echo "======================================================================"
echo "6. VALIDACION"
echo "======================================================================"

for F in \
    "${OUT}/91F_BAM_normalization.tsv" \
    "${OUT}/91F_ARG_coverage_long.tsv" \
    "${OUT}/91F_ARG_annotation_coverage.tsv" \
    "${OUT}/91F_coverage_summary.tsv" \
    "${OUT}/91F_run_metadata.tsv"
do

    [[ -s "$F" ]] || {
        echo "ERROR: salida ausente/vacía: $F" >&2
        exit 6
    }

    echo "OK: $(basename "$F")"

done

echo
echo "======================================================================"
echo "RESUMEN COBERTURA"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/91F_coverage_summary.tsv"

echo
echo "======================================================================"
echo "PASO 91F1B FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "======================================================================"
