#!/bin/bash

#SBATCH --job-name=91f1b_repair
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=91f1b_repair_%j.out
#SBATCH --error=91f1b_repair_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"
OUT="${ROOT}/81b_resistome_ARG_abundance_strict"

LOCI="${ROOT}/80_resistome_community_consolidated/91E_community_AMR_loci.tsv"

mkdir -p "$OUT"

echo "======================================================================"
echo "REPARACION 91F1B - CONSOLIDACION ESTRICTA"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"

[[ -s "$LOCI" ]] || {
    echo "ERROR: falta $LOCI" >&2
    exit 2
}

[[ -s "${OUT}/91F_BAM_normalization.tsv" ]] || {
    echo "ERROR: falta normalización estricta." >&2
    exit 3
}

srun --exact -n 1 -c 1 python3 - <<'PY'
import csv
import os
from pathlib import Path
from collections import Counter, defaultdict

root = Path(
    f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"
)

out = root / "81b_resistome_ARG_abundance_strict"

loci_file = (
    root
    / "80_resistome_community_consolidated"
    / "91E_community_AMR_loci.tsv"
)

coassemblies = ["L1","L2","L3","M1","M2","M3"]

samples = [
    f"{g}_{t}"
    for g in coassemblies
    for t in [1,2,3]
]

# ==========================================================
# 1. Leer catálogo de loci
# ==========================================================

with loci_file.open() as fh:
    r = csv.DictReader(fh, delimiter="\t")
    loci_rows = list(r)

loci = {
    x["community_locus_id"]: x
    for x in loci_rows
}

n_by_group = Counter(
    x["coassembly"]
    for x in loci_rows
)

print("===== LOCI ESPERADOS =====")

for g in coassemblies:
    print(
        g,
        n_by_group[g],
        sep="\t"
    )

if len(loci) != 565:
    raise RuntimeError(
        f"Se esperaban 565 loci; encontrados {len(loci)}"
    )

# ==========================================================
# 2. Validar y reunir los 18 archivos estrictos
# ==========================================================

all_cov = []
header = None

print()
print("===== ARCHIVOS ESTRICTOS =====")

for sample in samples:

    group = sample.split("_")[0]

    f = out / f"{sample}.ARG_coverage.tsv"

    if not f.exists() or f.stat().st_size == 0:
        raise RuntimeError(
            f"Falta archivo estricto: {f}"
        )

    with f.open() as fh:
        r = csv.DictReader(
            fh,
            delimiter="\t"
        )

        rows = list(r)

        if r.fieldnames is None:
            raise RuntimeError(
                f"Sin header: {f}"
            )

        if header is None:
            header = list(r.fieldnames)

        elif list(r.fieldnames) != header:
            raise RuntimeError(
                f"Header incompatible: {f}"
            )

    expected = n_by_group[group] * 2

    if len(rows) != expected:
        raise RuntimeError(
            f"{sample}: esperaba {expected} registros, "
            f"encontró {len(rows)}"
        )

    mode_counts = Counter(
        x["depth_mode"]
        for x in rows
    )

    if mode_counts["ALL_PRIMARY"] != n_by_group[group]:
        raise RuntimeError(
            f"{sample}: ALL_PRIMARY incorrecto"
        )

    if mode_counts["MAPQ20"] != n_by_group[group]:
        raise RuntimeError(
            f"{sample}: MAPQ20 incorrecto"
        )

    keys = [
        (
            x["community_locus_id"],
            x["sample"],
            x["depth_mode"]
        )
        for x in rows
    ]

    if len(keys) != len(set(keys)):
        raise RuntimeError(
            f"{sample}: claves duplicadas"
        )

    all_cov.extend(rows)

    print(
        sample,
        len(rows),
        mode_counts["ALL_PRIMARY"],
        mode_counts["MAPQ20"],
        sep="\t"
    )

expected_total = 565 * 3 * 2

if len(all_cov) != expected_total:
    raise RuntimeError(
        f"Total esperado {expected_total}, "
        f"obtenido {len(all_cov)}"
    )

# ==========================================================
# 3. Sobrescribir tabla long estricta correcta
# ==========================================================

coverage_long = (
    out
    / "91F_ARG_coverage_long.tsv"
)

with coverage_long.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=header,
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(all_cov)

print()
print(
    "coverage_long_rows:",
    len(all_cov)
)

# ==========================================================
# 4. Integrar anotación
# ==========================================================

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

result = []

for row in all_cov:

    lid = row["community_locus_id"]

    if lid not in loci:
        raise RuntimeError(
            f"Locus desconocido: {lid}"
        )

    locus = loci[lid]

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

integrated = (
    out
    / "91F_ARG_annotation_coverage.tsv"
)

fields = (
    header
    + annotation_fields
    + ["structural_flag"]
)

with integrated.open(
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

# ==========================================================
# 5. Resumen
# ==========================================================

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

summary_file = (
    out
    / "91F_coverage_summary.tsv"
)

with summary_file.open("w") as fh:

    fh.write("metric\tvalue\n")

    for k, v in metrics:
        fh.write(
            f"{k}\t{v}\n"
        )

# ==========================================================
# 6. QC adicional
# ==========================================================

key_counts = Counter(
    (
        x["community_locus_id"],
        x["sample"],
        x["depth_mode"]
    )
    for x in result
)

duplicates = sum(
    n - 1
    for n in key_counts.values()
    if n > 1
)

structural_counts = Counter(
    x["structural_flag"]
    for x in result
)

qc_file = (
    out
    / "91F1B_repair_QC.tsv"
)

with qc_file.open("w") as fh:

    fh.write("metric\tvalue\n")
    fh.write(
        f"community_loci\t{len(loci)}\n"
    )
    fh.write(
        f"samples\t{len(samples)}\n"
    )
    fh.write(
        f"integrated_rows\t{len(result)}\n"
    )
    fh.write(
        f"unique_keys\t{len(key_counts)}\n"
    )
    fh.write(
        f"duplicate_keys\t{duplicates}\n"
    )

    for k, v in sorted(
        structural_counts.items()
    ):
        fh.write(
            f"structural_{k}\t{v}\n"
        )

print()
print("===== RESUMEN ESTRICTO =====")

for k, v in metrics:
    print(k, v, sep="\t")

print()
print("integrated_rows:", len(result))
print("unique_keys:", len(key_counts))
print("duplicate_keys:", duplicates)

PY

# ============================================================
# Metadata de reparación
# ============================================================

{
    echo -e "field\tvalue"
    echo -e "analysis\t91F1B_strict_ARG_coverage_repaired_consolidation"
    echo -e "source_per_sample_files\t81b_resistome_ARG_abundance_strict/*.ARG_coverage.tsv"
    echo -e "samtools_depth_overlap_handling\t-s"
    echo -e "samtools_depth_excluded_flags\tUNMAP;SECONDARY;QCFAIL;DUP;SUPPLEMENTARY"
    echo -e "mapping_quality_modes\tALL_PRIMARY;MAPQ20"
    echo -e "depth_recalculated_in_repair\tNO"
    echo -e "repair_reason\toriginal_91F1B_had_hardcoded_81_resistome_ARG_abundance_paths_in_python_consolidation_blocks"
    echo -e "date\t$(date -Iseconds)"
} > "${OUT}/91F1B_repair_metadata.tsv"

echo
echo "======================================================================"
echo "VALIDACION"
echo "======================================================================"

for F in \
    "${OUT}/91F_BAM_normalization.tsv" \
    "${OUT}/91F_ARG_coverage_long.tsv" \
    "${OUT}/91F_ARG_annotation_coverage.tsv" \
    "${OUT}/91F_coverage_summary.tsv" \
    "${OUT}/91F1B_repair_QC.tsv" \
    "${OUT}/91F1B_repair_metadata.tsv"
do

    [[ -s "$F" ]] || {
        echo "ERROR: salida ausente/vacía: $F" >&2
        exit 4
    }

    echo "OK: $(basename "$F")"
done

echo
echo "======================================================================"
echo "RESUMEN FINAL"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/91F_coverage_summary.tsv"

echo
echo "======================================================================"
echo "REPARACION 91F1B FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "======================================================================"
