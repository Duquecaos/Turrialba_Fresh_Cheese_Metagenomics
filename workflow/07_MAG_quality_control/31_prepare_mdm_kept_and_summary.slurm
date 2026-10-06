#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

#SBATCH --job-name=mdm_prepare_kept
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

MDM="${ROOT}/15_mdmcleaner_rescue/strong_gunc_fail"
STAGE="${ROOT}/15_mdmcleaner_rescue/inputs_writable"

OUT="${ROOT}/16_mdmcleaner_kept"
FASTAS="${OUT}/fastas"

SUMMARY="${OUT}/mdmcleaner_filter_summary.tsv"
MANIFEST="${OUT}/mdmcleaner_kept_manifest.tsv"

mkdir -p "${FASTAS}"

echo "============================================================"
echo "PREPARACIÓN MAGs POST-MDMCLEANER"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job: ${SLURM_JOB_ID}"
echo "============================================================"
echo


srun --exact -n 1 -c 1 \
python3 - \
"${MDM}" \
"${STAGE}" \
"${FASTAS}" \
"${SUMMARY}" \
"${MANIFEST}" <<'PY'
import sys
import gzip
import shutil
import hashlib
from pathlib import Path

mdm = Path(sys.argv[1])
stage = Path(sys.argv[2])
fastas = Path(sys.argv[3])
summary_file = Path(sys.argv[4])
manifest_file = Path(sys.argv[5])

mags = [
    ("L1", "L1_maxbin2.051"),
    ("L1", "concoct_20"),
    ("L1", "concoct_24"),
    ("L2", "L2_maxbin2.004_sub"),
    ("L2", "concoct_55"),
    ("L2", "concoct_72"),
    ("L2", "concoct_91"),
    ("M3", "concoct_37"),
]

categories = {
    "keep": "filtered_kept_contigs.fasta.gz",
    "eval_low": "filtered_need_evaluation_low.fasta.gz",
    "eval_high": "filtered_need_evaluation_high.fasta.gz",
    "removed": "filtered_removed.fasta.gz",
}


def fasta_stats(path):
    opener = gzip.open if str(path).endswith(".gz") else open

    nseq = 0
    nbp = 0

    with opener(path, "rt") as fh:
        seq_len = 0
        seen = False

        for line in fh:
            if line.startswith(">"):
                if seen:
                    nseq += 1
                    nbp += seq_len

                seen = True
                seq_len = 0
            else:
                seq_len += len(line.strip())

        if seen:
            nseq += 1
            nbp += seq_len

    return nseq, nbp


def sha256(path):
    h = hashlib.sha256()

    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)

    return h.hexdigest()


rows = []
manifest_rows = []

for group, mag in mags:

    original = stage / group / f"{mag}.fa"

    if not original.is_file() or original.stat().st_size == 0:
        raise SystemExit(
            f"ERROR: falta MAG original staged: {original}"
        )

    original_n, original_bp = fasta_stats(original)

    cat_stats = {}

    for category, suffix in categories.items():

        path = mdm / mag / f"{mag}_{suffix}"

        if not path.is_file():
            raise SystemExit(
                f"ERROR: falta salida {category}: {path}"
            )

        cat_stats[category] = (*fasta_stats(path), path)

    summed_n = sum(x[0] for x in cat_stats.values())
    summed_bp = sum(x[1] for x in cat_stats.values())

    if summed_n != original_n:
        raise SystemExit(
            f"ERROR {mag}: categorías={summed_n} contigs, "
            f"original={original_n}"
        )

    if summed_bp != original_bp:
        raise SystemExit(
            f"ERROR {mag}: categorías={summed_bp} bp, "
            f"original={original_bp} bp"
        )

    keep_n, keep_bp, keep_file = cat_stats["keep"]

    if keep_n == 0 or keep_bp == 0:
        raise SystemExit(
            f"ERROR: {mag} quedó sin contigs KEEP"
        )

    # Nombre único y trazable para análisis posteriores
    final_name = f"{group}__{mag}.fa"
    final_fasta = fastas / final_name

    with gzip.open(keep_file, "rb") as src, \
         open(final_fasta, "wb") as dst:
        shutil.copyfileobj(src, dst)

    final_n, final_bp = fasta_stats(final_fasta)

    if (final_n, final_bp) != (keep_n, keep_bp):
        raise SystemExit(
            f"ERROR al descomprimir KEEP de {mag}"
        )

    keep_pct_bp = 100 * keep_bp / original_bp
    removed_pct_bp = (
        100 * cat_stats["removed"][1] / original_bp
    )
    eval_low_pct_bp = (
        100 * cat_stats["eval_low"][1] / original_bp
    )
    eval_high_pct_bp = (
        100 * cat_stats["eval_high"][1] / original_bp
    )

    rows.append([
        group,
        mag,
        original_n,
        original_bp,
        cat_stats["keep"][0],
        cat_stats["keep"][1],
        f"{keep_pct_bp:.3f}",
        cat_stats["eval_low"][0],
        cat_stats["eval_low"][1],
        f"{eval_low_pct_bp:.3f}",
        cat_stats["eval_high"][0],
        cat_stats["eval_high"][1],
        f"{eval_high_pct_bp:.3f}",
        cat_stats["removed"][0],
        cat_stats["removed"][1],
        f"{removed_pct_bp:.3f}",
    ])

    manifest_rows.append([
        group,
        mag,
        str(original),
        str(keep_file),
        str(final_fasta),
        sha256(final_fasta),
    ])


header = [
    "group",
    "MAG",
    "original_contigs",
    "original_bp",
    "keep_contigs",
    "keep_bp",
    "keep_bp_pct",
    "eval_low_contigs",
    "eval_low_bp",
    "eval_low_bp_pct",
    "eval_high_contigs",
    "eval_high_bp",
    "eval_high_bp_pct",
    "removed_contigs",
    "removed_bp",
    "removed_bp_pct",
]

with open(summary_file, "w") as out:
    out.write("\t".join(header) + "\n")

    for row in rows:
        out.write("\t".join(map(str, row)) + "\n")


with open(manifest_file, "w") as out:
    out.write(
        "group\tMAG\toriginal_fasta\t"
        "mdm_keep_gz\tpost_mdm_fasta\tsha256\n"
    )

    for row in manifest_rows:
        out.write("\t".join(row) + "\n")


print("MAGs preparados:", len(rows))
print()
print("Resumen:")

for row in rows:
    print(
        f"{row[0]:2s}  {row[1]:22s} "
        f"KEEP={row[4]:4d} contigs "
        f"{row[6]:>7s}% bp | "
        f"LOW={row[7]:3d} | "
        f"HIGH={row[10]:3d} | "
        f"DELETE={row[13]:3d}"
    )
PY


echo
echo "===== TABLA COMPLETA ====="

column -t -s $'\t' "${SUMMARY}"

echo
echo "===== FASTA POST-MDMCLEANER ====="

ls -lh "${FASTAS}"

N_FASTA=$(find "${FASTAS}" \
    -maxdepth 1 \
    -type f \
    -name '*.fa' \
    | wc -l)

echo
echo "Número de MAGs KEEP: ${N_FASTA}"

if [[ "${N_FASTA}" -ne 8 ]]; then
    echo "ERROR: se esperaban 8 MAGs KEEP" >&2
    exit 1
fi

echo
echo "============================================================"
echo "MDMCLEANER KEPT PREPARADOS"
echo "Resumen: ${SUMMARY}"
echo "FASTAs:  ${FASTAS}"
echo "Fin: $(date)"
echo "============================================================"
