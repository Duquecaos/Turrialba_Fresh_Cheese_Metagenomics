#!/bin/bash

#SBATCH --job-name=91f2_amrfinal
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=91f2_amrfinal_%j.out
#SBATCH --error=91f2_amrfinal_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

STRICT="${ROOT}/81b_resistome_ARG_abundance_strict"
COMM="${ROOT}/80_resistome_community_consolidated"
FINAL18="${ROOT}/78_resistome_final18_integrated"

OUT="${ROOT}/82_resistome_final_synthesis"

mkdir -p "$OUT"

INPUT="${STRICT}/91F_ARG_annotation_coverage.tsv"
LOCI="${COMM}/91E_community_AMR_loci.tsv"
FINAL18_LOCI="${FINAL18}/91D_integrated_AMR_loci.tsv"

for F in \
    "$INPUT" \
    "$LOCI" \
    "$FINAL18_LOCI"
do
    [[ -s "$F" ]] || {
        echo "ERROR: falta archivo: $F" >&2
        exit 2
    }
done

echo "======================================================================"
echo "PASO 91F2 - SINTESIS FINAL DEL RESISTOMA"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"

srun --exact -n 1 -c 1 python3 - <<'PY'
import csv
import os
import statistics
from pathlib import Path
from collections import Counter, defaultdict

root = Path(
    f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"
)

strict = root / "81b_resistome_ARG_abundance_strict"
comm = root / "80_resistome_community_consolidated"
final18 = root / "78_resistome_final18_integrated"

out = root / "82_resistome_final_synthesis"

src = strict / "91F_ARG_annotation_coverage.tsv"
loci_file = comm / "91E_community_AMR_loci.tsv"
final18_file = final18 / "91D_integrated_AMR_loci.tsv"

# ============================================================
# Helpers
# ============================================================

def clean(x):
    if x is None:
        return ""
    x = str(x).strip()
    if x.upper() == "NA":
        return ""
    return x


def fnum(x):
    x = clean(x)
    return float(x) if x else 0.0


def semijoin(values):
    return ";".join(
        sorted({
            clean(x)
            for x in values
            if clean(x)
        })
    )


def read_support_class(b1, b5):

    if b1 == 0:
        return "no_read_support"

    if b1 >= 0.95 and b5 >= 0.80:
        return "high_depth_broad"

    if b1 >= 0.95:
        return "very_broad_low_depth"

    if b1 >= 0.80:
        return "broad"

    if b1 >= 0.50:
        return "moderate"

    return "limited"


def annotation_group(method):

    method = clean(method)

    if (
        "PARTIAL" in method
        or "INTERNAL_STOP" in method
    ):
        return "partial_or_interrupted"

    if method in {
        "EXACTP",
        "EXACTX",
        "ALLELEP",
        "ALLELEX",
    }:
        return "exact_or_allele"

    if method in {
        "BLASTP",
        "BLASTX",
        "HMM",
    }:
        return "curated_homology"

    return "other"


# ============================================================
# Final18: símbolo -> MAGs / soporte
# ============================================================

final18_by_symbol = defaultdict(
    lambda: {
        "MAGs": set(),
        "taxa": set(),
        "supports": set(),
    }
)

with final18_file.open() as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in r:

        symbols = [
            x.strip()
            for x in row["gene_symbols"].split(";")
            if x.strip()
        ]

        taxon = clean(row.get("species"))

        if not taxon:
            taxon = clean(row.get("genus"))

        for symbol in symbols:

            final18_by_symbol[symbol]["MAGs"].add(
                row["MAG"]
            )

            if taxon:
                final18_by_symbol[symbol]["taxa"].add(
                    taxon
                )

            final18_by_symbol[symbol]["supports"].add(
                row["evidence_support"]
            )


# ============================================================
# Leer cobertura estricta y emparejar modos
# ============================================================

with src.open() as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    raw = list(r)

print("raw_records:", len(raw))

if len(raw) != 3390:
    raise RuntimeError(
        f"Se esperaban 3390 filas; hay {len(raw)}"
    )

paired = defaultdict(dict)

for row in raw:

    key = (
        row["community_locus_id"],
        row["sample"],
    )

    mode = row["depth_mode"]

    if mode in paired[key]:
        raise RuntimeError(
            f"Modo duplicado: {key} {mode}"
        )

    paired[key][mode] = row


# ============================================================
# Tabla locus x muestra usando MAPQ20 como principal
# ============================================================

sample_rows = []

for key, modes in paired.items():

    if set(modes) != {
        "ALL_PRIMARY",
        "MAPQ20",
    }:
        raise RuntimeError(
            f"Faltan modos en {key}: {set(modes)}"
        )

    a = modes["ALL_PRIMARY"]
    m = modes["MAPQ20"]

    all_b1 = fnum(a["breadth_1x"])
    map_b1 = fnum(m["breadth_1x"])

    map_b5 = fnum(m["breadth_5x"])
    map_b10 = fnum(m["breadth_10x"])

    all_depth = fnum(a["mean_depth"])
    map_depth = fnum(m["mean_depth"])

    breadth_loss = all_b1 - map_b1
    depth_loss = all_depth - map_depth

    if all_depth > 0:
        mapq_depth_retention = (
            map_depth / all_depth
        )
    else:
        mapq_depth_retention = ""

    # ------------------------------------------
    # Sensibilidad de mapeo
    # ------------------------------------------

    if map_b1 >= 0.80:

        mapping_flag = (
            "broad_support_MAPQ20"
        )

    elif (
        all_b1 >= 0.80
        and map_b1 < 0.50
    ):

        mapping_flag = (
            "multimapping_sensitive"
        )

    elif all_b1 >= 0.80:

        mapping_flag = (
            "intermediate_MAPQ_loss"
        )

    else:

        mapping_flag = (
            "limited_even_ALL_PRIMARY"
        )

    # ------------------------------------------
    # Evidencia de reads
    # ------------------------------------------

    support = read_support_class(
        map_b1,
        map_b5
    )

    method = m["method"]

    ann_group = annotation_group(
        method
    )

    symbol = m["element_symbol"]

    final18_info = final18_by_symbol.get(
        symbol,
        None
    )

    if final18_info:

        final18_mags = ";".join(
            sorted(final18_info["MAGs"])
        )

        final18_taxa = ";".join(
            sorted(final18_info["taxa"])
        )

        final18_support = ";".join(
            sorted(final18_info["supports"])
        )

    else:

        final18_mags = ""
        final18_taxa = ""
        final18_support = ""

    sample_rows.append({
        "community_locus_id":
            m["community_locus_id"],

        "element_symbol":
            symbol,

        "element_name":
            m["element_name"],

        "class":
            m["class"],

        "subclass":
            m["subclass"],

        "coassembly":
            m["coassembly"],

        "producer":
            m["producer"],

        "biological_unit":
            m["biological_unit"],

        "sample":
            m["sample"],

        "week":
            m["week"],

        "contig":
            m["contig"],

        "start":
            m["start"],

        "stop":
            m["stop"],

        "length_bp":
            m["length_bp"],

        "method":
            method,

        "annotation_evidence_group":
            ann_group,

        "structural_flag":
            m["structural_flag"],

        "scope":
            m["scope"],

        "identity_reference":
            m["identity_reference"],

        "coverage_reference":
            m["coverage_reference"],

        "representation_vs_final18":
            m["representation_vs_final18"],

        "final18_MAGs_same_symbol":
            final18_mags,

        "final18_taxa_same_symbol":
            final18_taxa,

        "final18_evidence_same_symbol":
            final18_support,

        "MAPQ20_breadth_1x":
            map_b1,

        "MAPQ20_breadth_5x":
            map_b5,

        "MAPQ20_breadth_10x":
            map_b10,

        "MAPQ20_mean_depth":
            map_depth,

        "MAPQ20_median_depth":
            fnum(m["median_depth"]),

        "MAPQ20_max_depth":
            fnum(m["max_depth"]),

        "normalized_mean_depth":
            fnum(
                m[
                    "mean_depth_per_million_mapped"
                ]
            ),

        "eligible_primary_mapped_records":
            m[
                "eligible_primary_mapped_records"
            ],

        "ALL_PRIMARY_breadth_1x":
            all_b1,

        "ALL_PRIMARY_mean_depth":
            all_depth,

        "breadth_loss_after_MAPQ20":
            breadth_loss,

        "depth_loss_after_MAPQ20":
            depth_loss,

        "MAPQ20_depth_retention":
            mapq_depth_retention,

        "mapping_support_flag":
            mapping_flag,

        "read_support_class":
            support,
    })


# ============================================================
# Orden
# ============================================================

sample_rows.sort(
    key=lambda x: (
        x["producer"],
        int(x["biological_unit"]),
        int(x["week"]),
        x["element_symbol"],
        x["community_locus_id"],
    )
)


# ============================================================
# Escribir tabla locus x muestra
# ============================================================

sample_file = (
    out
    / "91F2_ARG_sample_evidence.tsv"
)

sample_fields = list(
    sample_rows[0].keys()
)

with sample_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=sample_fields,
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(sample_rows)


# ============================================================
# Resumen por locus
# ============================================================

by_locus = defaultdict(list)

for row in sample_rows:
    by_locus[
        row["community_locus_id"]
    ].append(row)

locus_rows = []

for locus_id, rows in sorted(
    by_locus.items()
):

    if len(rows) != 3:
        raise RuntimeError(
            f"{locus_id}: esperaba 3 tiempos, "
            f"encontró {len(rows)}"
        )

    rows = sorted(
        rows,
        key=lambda x: int(x["week"])
    )

    broad = [
        x
        for x in rows
        if float(
            x["MAPQ20_breadth_1x"]
        ) >= 0.80
    ]

    very_broad = [
        x
        for x in rows
        if float(
            x["MAPQ20_breadth_1x"]
        ) >= 0.95
    ]

    high_depth = [
        x
        for x in rows
        if (
            float(
                x["MAPQ20_breadth_1x"]
            ) >= 0.95
            and
            float(
                x["MAPQ20_breadth_5x"]
            ) >= 0.80
        )
    ]

    n_broad = len(broad)

    if n_broad == 3:
        longitudinal = (
            "broad_support_all_3_weeks"
        )

    elif n_broad == 2:
        longitudinal = (
            "broad_support_2_of_3_weeks"
        )

    elif n_broad == 1:
        longitudinal = (
            "broad_support_1_of_3_weeks"
        )

    else:
        longitudinal = (
            "no_week_with_broad_support"
        )

    r0 = rows[0]

    locus_rows.append({
        "community_locus_id":
            locus_id,

        "element_symbol":
            r0["element_symbol"],

        "element_name":
            r0["element_name"],

        "class":
            r0["class"],

        "subclass":
            r0["subclass"],

        "coassembly":
            r0["coassembly"],

        "producer":
            r0["producer"],

        "biological_unit":
            r0["biological_unit"],

        "method":
            r0["method"],

        "annotation_evidence_group":
            r0["annotation_evidence_group"],

        "structural_flag":
            r0["structural_flag"],

        "representation_vs_final18":
            r0["representation_vs_final18"],

        "final18_MAGs_same_symbol":
            r0["final18_MAGs_same_symbol"],

        "n_weeks_breadth_ge_0.80":
            n_broad,

        "n_weeks_breadth_ge_0.95":
            len(very_broad),

        "n_weeks_high_depth_broad":
            len(high_depth),

        "longitudinal_support":
            longitudinal,

        "week0_breadth":
            rows[0]["MAPQ20_breadth_1x"],

        "week1_breadth":
            rows[1]["MAPQ20_breadth_1x"],

        "week2_breadth":
            rows[2]["MAPQ20_breadth_1x"],

        "week0_norm_depth":
            rows[0]["normalized_mean_depth"],

        "week1_norm_depth":
            rows[1]["normalized_mean_depth"],

        "week2_norm_depth":
            rows[2]["normalized_mean_depth"],

        "max_MAPQ20_breadth":
            max(
                float(
                    x["MAPQ20_breadth_1x"]
                )
                for x in rows
            ),

        "max_normalized_depth":
            max(
                float(
                    x["normalized_mean_depth"]
                )
                for x in rows
            ),

        "median_normalized_depth":
            statistics.median(
                float(
                    x["normalized_mean_depth"]
                )
                for x in rows
            ),

        "n_multimapping_sensitive_weeks":
            sum(
                x["mapping_support_flag"]
                ==
                "multimapping_sensitive"
                for x in rows
            ),
    })


locus_file = (
    out
    / "91F2_ARG_locus_summary.tsv"
)

with locus_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            locus_rows[0].keys()
        ),
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(locus_rows)


# ============================================================
# Resumen por símbolo
# ============================================================

by_symbol = defaultdict(list)

for row in locus_rows:
    by_symbol[
        row["element_symbol"]
    ].append(row)

symbol_rows = []

for symbol, rows in sorted(
    by_symbol.items()
):

    producers = sorted({
        x["producer"]
        for x in rows
    })

    coassemblies = sorted({
        x["coassembly"]
        for x in rows
    })

    persistent = [
        x for x in rows
        if x["n_weeks_breadth_ge_0.80"] == 3
    ]

    recurrent = [
        x for x in rows
        if x["n_weeks_breadth_ge_0.80"] >= 2
    ]

    partial = [
        x for x in rows
        if x["structural_flag"]
        ==
        "partial_or_interrupted_detection"
    ]

    symbol_rows.append({
        "element_symbol":
            symbol,

        "classes":
            semijoin(
                x["class"]
                for x in rows
            ),

        "subclasses":
            semijoin(
                x["subclass"]
                for x in rows
            ),

        "n_community_loci":
            len(rows),

        "n_coassemblies":
            len(coassemblies),

        "coassemblies":
            ";".join(coassemblies),

        "n_producers":
            len(producers),

        "producers":
            ";".join(producers),

        "n_loci_broad_all_3_weeks":
            len(persistent),

        "n_loci_broad_at_least_2_weeks":
            len(recurrent),

        "n_partial_or_interrupted_loci":
            len(partial),

        "annotation_evidence_groups":
            semijoin(
                x["annotation_evidence_group"]
                for x in rows
            ),

        "representation_vs_final18":
            semijoin(
                x["representation_vs_final18"]
                for x in rows
            ),

        "final18_MAGs_same_symbol":
            semijoin(
                x["final18_MAGs_same_symbol"]
                for x in rows
            ),

        "max_normalized_depth":
            max(
                x["max_normalized_depth"]
                for x in rows
            ),
    })


symbol_file = (
    out
    / "91F2_ARG_symbol_summary.tsv"
)

with symbol_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            symbol_rows[0].keys()
        ),
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(symbol_rows)


# ============================================================
# Productor x semana:
# contar símbolos, NO sumar profundidad
# ============================================================

producer_week = defaultdict(list)

for row in sample_rows:

    producer_week[
        (
            row["producer"],
            row["week"]
        )
    ].append(row)

pw_rows = []

for (
    producer,
    week
), rows in sorted(
    producer_week.items()
):

    broad = [
        x for x in rows
        if float(
            x["MAPQ20_breadth_1x"]
        ) >= 0.80
    ]

    high = [
        x for x in rows
        if (
            float(
                x["MAPQ20_breadth_1x"]
            ) >= 0.95
            and
            float(
                x["MAPQ20_breadth_5x"]
            ) >= 0.80
        )
    ]

    pw_rows.append({
        "producer":
            producer,

        "week":
            week,

        "locus_sample_observations":
            len(rows),

        "broad_supported_loci":
            len(broad),

        "broad_supported_unique_symbols":
            len({
                x["element_symbol"]
                for x in broad
            }),

        "high_depth_broad_loci":
            len(high),

        "high_depth_broad_unique_symbols":
            len({
                x["element_symbol"]
                for x in high
            }),
    })


pw_file = (
    out
    / "91F2_producer_week_summary.tsv"
)

with pw_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            pw_rows[0].keys()
        ),
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(pw_rows)


# ============================================================
# Prioritarios
# ============================================================

priority = {
    "mecA",
    "mcr-3",
    "mcr-9.1",
    "erm(B)",
    "sul1",
    "tet(M)",
    "tet(S)",
    "qnrB",
    "qnrC",
    "qnrD1",
    "qnrE",
    "qnrS1",
    "qnrVC5",
    "dfrA14",
    "dfrA31",
    "dfrA51",
    "dfrG",
}

priority_rows = [
    x for x in locus_rows
    if x["element_symbol"] in priority
]

priority_file = (
    out
    / "91F2_priority_ARG_loci.tsv"
)

with priority_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            locus_rows[0].keys()
        ),
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(priority_rows)


# ============================================================
# Resumen global
# ============================================================

support_counts = Counter(
    x["read_support_class"]
    for x in sample_rows
)

mapping_counts = Counter(
    x["mapping_support_flag"]
    for x in sample_rows
)

long_counts = Counter(
    x["longitudinal_support"]
    for x in locus_rows
)

annotation_counts = Counter(
    x["annotation_evidence_group"]
    for x in locus_rows
)

struct_counts = Counter(
    x["structural_flag"]
    for x in locus_rows
)

repr_counts = Counter(
    x["representation_vs_final18"]
    for x in locus_rows
)

summary = [
    (
        "community_ARG_loci",
        len(locus_rows)
    ),
    (
        "unique_ARG_symbols",
        len(symbol_rows)
    ),
    (
        "locus_sample_observations",
        len(sample_rows)
    ),
]

for k, v in sorted(
    support_counts.items()
):
    summary.append(
        (
            f"read_support_{k}",
            v
        )
    )

for k, v in sorted(
    mapping_counts.items()
):
    summary.append(
        (
            f"mapping_{k}",
            v
        )
    )

for k, v in sorted(
    long_counts.items()
):
    summary.append(
        (
            f"longitudinal_{k}",
            v
        )
    )

for k, v in sorted(
    annotation_counts.items()
):
    summary.append(
        (
            f"annotation_{k}",
            v
        )
    )

for k, v in sorted(
    struct_counts.items()
):
    summary.append(
        (
            f"structural_{k}",
            v
        )
    )

for k, v in sorted(
    repr_counts.items()
):
    summary.append(
        (
            f"final18_relation_{k}",
            v
        )
    )


summary_file = (
    out
    / "91F2_resistome_final_summary.tsv"
)

with summary_file.open("w") as fh:

    fh.write("metric\tvalue\n")

    for k, v in summary:
        fh.write(
            f"{k}\t{v}\n"
        )


# ============================================================
# QC
# ============================================================

if len(sample_rows) != 1695:
    raise RuntimeError(
        f"Esperaba 1695 observaciones MAPQ20; "
        f"hay {len(sample_rows)}"
    )

if len(locus_rows) != 565:
    raise RuntimeError(
        f"Esperaba 565 loci; "
        f"hay {len(locus_rows)}"
    )

if len(symbol_rows) != 115:
    raise RuntimeError(
        f"Esperaba 115 símbolos; "
        f"hay {len(symbol_rows)}"
    )

print()
print("=" * 78)
print("RESUMEN FINAL")
print("=" * 78)

for k, v in summary:
    print(k, v, sep="\t")

print()
print("Producer/week:")

for x in pw_rows:
    print(
        x["producer"],
        x["week"],
        x["broad_supported_loci"],
        x["broad_supported_unique_symbols"],
        x["high_depth_broad_loci"],
        x["high_depth_broad_unique_symbols"],
        sep="\t"
    )

PY

# ============================================================
# Metadata
# ============================================================

{
    echo -e "field\tvalue"
    echo -e "analysis\tfinal_resistome_synthesis"
    echo -e "canonical_coverage\t81b_resistome_ARG_abundance_strict"
    echo -e "principal_mapping_filter\tMAPQ20"
    echo -e "mate_overlap_handling\tsamtools_depth_-s"
    echo -e "excluded_alignment_flags\tUNMAP;SECONDARY;QCFAIL;DUP;SUPPLEMENTARY"
    echo -e "presence_absence_call\tNOT_USED"
    echo -e "read_support_high_depth_broad\tbreadth1>=0.95 AND breadth5>=0.80"
    echo -e "read_support_very_broad\tbreadth1>=0.95"
    echo -e "read_support_broad\tbreadth1>=0.80"
    echo -e "read_support_moderate\tbreadth1>=0.50"
    echo -e "read_support_limited\t0<breadth1<0.50"
    echo -e "read_support_none\tbreadth1=0"
    echo -e "threshold_interpretation\tdescriptive_operational_not_biological_presence_threshold"
    echo -e "date\t$(date -Iseconds)"
} > "${OUT}/91F2_run_metadata.tsv"

echo
echo "======================================================================"
echo "VALIDACION"
echo "======================================================================"

for F in \
    "${OUT}/91F2_ARG_sample_evidence.tsv" \
    "${OUT}/91F2_ARG_locus_summary.tsv" \
    "${OUT}/91F2_ARG_symbol_summary.tsv" \
    "${OUT}/91F2_producer_week_summary.tsv" \
    "${OUT}/91F2_priority_ARG_loci.tsv" \
    "${OUT}/91F2_resistome_final_summary.tsv" \
    "${OUT}/91F2_run_metadata.tsv"
do

    [[ -s "$F" ]] || {
        echo "ERROR: salida ausente/vacía: $F" >&2
        exit 3
    }

    echo "OK: $(basename "$F")"

done

echo
echo "======================================================================"
echo "RESUMEN FINAL DEL RESISTOMA"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/91F2_resistome_final_summary.tsv"

echo
echo "======================================================================"
echo "PRODUCTOR x SEMANA"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/91F2_producer_week_summary.tsv"

echo
echo "======================================================================"
echo "PASO 91F2 FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "======================================================================"
