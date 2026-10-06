#!/bin/bash

#SBATCH --job-name=92h1_hosttax
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=92h1_hosttax_%j.out
#SBATCH --error=92h1_hosttax_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

HOSTLINKS="${ROOT}/94_votu_host_taxonomy/audit/92H0b_direct_votu_MAG_links.tsv"
MAGTAX="${ROOT}/43_MAG_abundance_taxonomy/representative_MAG_taxonomy.tsv"

VOTUMETA="${ROOT}/93_votu_coverage/integrated/92G4_votu_metadata.tsv"
RECURRENCE="${ROOT}/93_votu_coverage/integrated/92G4_mapq20_id95_library_recurrence.tsv"
LONG="${ROOT}/93_votu_coverage/integrated/92G4_votu_sample_long.tsv.gz"

OUT="${ROOT}/94_votu_host_taxonomy/direct_host"

mkdir -p "$OUT"

DONE="${OUT}/92H1_COMPLETE.ok"

echo "======================================================================"
echo "92H1 - DIRECT HOST TAXONOMY INTEGRATION"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"

for F in \
    "$HOSTLINKS" \
    "$MAGTAX" \
    "$VOTUMETA" \
    "$RECURRENCE" \
    "$LONG"
do
    [[ -s "$F" ]] || {
        echo "ERROR: falta $F" >&2
        exit 2
    }
done


python3 - \
    "$HOSTLINKS" \
    "$MAGTAX" \
    "$VOTUMETA" \
    "$RECURRENCE" \
    "$LONG" \
    "$OUT" <<'PY'

import csv
import gzip
import sys

from collections import Counter, defaultdict
from pathlib import Path


(
    hostlinks_file,
    magtax_file,
    votumeta_file,
    recurrence_file,
    long_file,
    outdir,
) = map(Path, sys.argv[1:])


# ============================================================
# Helpers
# ============================================================

def read_tsv(path):

    with path.open() as fh:

        return list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


def clean(x):

    return "" if x is None else str(x).strip()


def ffloat(x):

    x = clean(x)

    if x in {
        "",
        "NA",
        "NaN",
        "nan",
    }:
        return None

    return float(x)


def fint(x):

    x = clean(x)

    if x == "":
        return None

    return int(
        float(x)
    )


# ============================================================
# 1. Direct host links
# ============================================================

hostlinks = read_tsv(
    hostlinks_file
)

if len(hostlinks) != 31:

    raise RuntimeError(
        f"Esperaba 31 pares vOTU-MAG; "
        f"encontré {len(hostlinks)}"
    )


pair_keys = {
    (
        x["votu_id"],
        x["final18_MAG"],
    )
    for x in hostlinks
}

if len(pair_keys) != 31:

    raise RuntimeError(
        "Hay pares vOTU-MAG duplicados"
    )


# ============================================================
# 2. MAG taxonomy
# ============================================================

magtax_rows = read_tsv(
    magtax_file
)

if not magtax_rows:

    raise RuntimeError(
        "MAG taxonomy vacío"
    )


required_mag_cols = {
    "MAG",
    "completeness",
    "contamination",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "taxon_label",
    "taxon_label_rank",
    "warnings",
}

missing = (
    required_mag_cols
    - set(
        magtax_rows[0].keys()
    )
)

if missing:

    raise RuntimeError(
        f"Faltan columnas MAG taxonomy: "
        f"{sorted(missing)}"
    )


magtax = {}

for row in magtax_rows:

    mag = row["MAG"]

    if mag in magtax:

        raise RuntimeError(
            f"MAG duplicado en taxonomy: {mag}"
        )

    magtax[mag] = row


linked_mags = {
    x["final18_MAG"]
    for x in hostlinks
}

unmatched_mags = (
    linked_mags
    - set(magtax)
)

if unmatched_mags:

    raise RuntimeError(
        "MAGs sin taxonomía: "
        + ";".join(
            sorted(
                unmatched_mags
            )
        )
    )


# ============================================================
# 3. vOTU metadata
# ============================================================

votu_rows = read_tsv(
    votumeta_file
)

votu_meta = {}

for row in votu_rows:

    votu = row["votu_id"]

    if votu in votu_meta:

        raise RuntimeError(
            f"vOTU duplicado: {votu}"
        )

    votu_meta[votu] = row


linked_votus = {
    x["votu_id"]
    for x in hostlinks
}

missing_votus = (
    linked_votus
    - set(votu_meta)
)

if missing_votus:

    raise RuntimeError(
        "vOTUs ausentes de metadata: "
        + ";".join(
            sorted(
                missing_votus
            )
        )
    )


# ============================================================
# 4. Recurrence MAPQ20_ID95
# ============================================================

rec_rows = read_tsv(
    recurrence_file
)

recurrence = {
    row["votu_id"]: row
    for row in rec_rows
}

if len(recurrence) != 4491:

    raise RuntimeError(
        f"Recurrence rows={len(recurrence)}; "
        "esperaba 4491"
    )


# ============================================================
# 5. MAPQ20_ID95 quantitative metrics for the 31 vOTUs
# ============================================================

quant = defaultdict(dict)

with gzip.open(
    long_file,
    "rt"
) as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in reader:

        if row["mode"] != "MAPQ20_ID95":
            continue

        votu = row["votu_id"]

        if votu not in linked_votus:
            continue

        sample = row["sample"]

        quant[votu][sample] = {
            "mapped_records":
                int(
                    float(
                        row[
                            "mapped_records"
                        ]
                    )
                ),

            "rpm_hostfree":
                float(
                    row[
                        "mapped_records_per_million_hostfree_reads"
                    ]
                ),

            "breadth_1x":
                float(
                    row[
                        "breadth_1x"
                    ]
                ),

            "breadth_5x":
                float(
                    row[
                        "breadth_5x"
                    ]
                ),

            "mean_depth":
                float(
                    row[
                        "mean_depth"
                    ]
                ),

            "mean_depth_per_million_hostfree_reads":
                float(
                    row[
                        "mean_depth_per_million_hostfree_reads"
                    ]
                ),
        }


samples = [
    "L1_1", "L1_2", "L1_3",
    "L2_1", "L2_2", "L2_3",
    "L3_1", "L3_2", "L3_3",
    "M1_1", "M1_2", "M1_3",
    "M2_1", "M2_2", "M2_3",
    "M3_1", "M3_2", "M3_3",
]


for votu in linked_votus:

    if len(
        quant[votu]
    ) != 18:

        raise RuntimeError(
            f"{votu}: solo "
            f"{len(quant[votu])}/18 "
            "bibliotecas en cuantificación"
        )


# ============================================================
# 6. Integrated pair table
# ============================================================

integrated = []


for link in hostlinks:

    votu = link[
        "votu_id"
    ]

    mag = link[
        "final18_MAG"
    ]

    h = magtax[
        mag
    ]

    v = votu_meta[
        votu
    ]

    r = recurrence[
        votu
    ]

    q = quant[
        votu
    ]


    breadth80_samples = [
        s
        for s in samples
        if q[s][
            "breadth_1x"
        ] >= 0.80
    ]

    breadth95_samples = [
        s
        for s in samples
        if q[s][
            "breadth_1x"
        ] >= 0.95
    ]


    max_depth_sample = max(
        samples,
        key=lambda s:
            q[s][
                "mean_depth_per_million_hostfree_reads"
            ]
    )


    integrated.append({

        # ---------------- Host evidence ----------------
        "votu_id":
            votu,

        "final18_MAG":
            mag,

        "host_evidence":
            link[
                "host_evidence"
            ],

        "n_provirus_predictions":
            link[
                "n_provirus_predictions"
            ],

        "n_parent_contigs":
            link[
                "n_parent_contigs"
            ],

        "parent_contigs":
            link[
                "parent_contigs"
            ],

        "coassemblies":
            link[
                "coassemblies"
            ],


        # ---------------- MAG quality ----------------
        "host_MAG_completeness":
            h[
                "completeness"
            ],

        "host_MAG_contamination":
            h[
                "contamination"
            ],

        "host_MAG_warnings":
            h[
                "warnings"
            ],


        # ---------------- Host taxonomy ----------------
        "host_phylum":
            h[
                "phylum"
            ],

        "host_class":
            h[
                "class"
            ],

        "host_order":
            h[
                "order"
            ],

        "host_family":
            h[
                "family"
            ],

        "host_genus":
            h[
                "genus"
            ],

        "host_species":
            h[
                "species"
            ],

        "host_taxon_label":
            h[
                "taxon_label"
            ],

        "host_taxon_label_rank":
            h[
                "taxon_label_rank"
            ],


        # ---------------- Viral quality ----------------
        "votu_seed_length":
            v[
                "seed_length"
            ],

        "votu_n_members":
            v[
                "n_vclean_members"
            ],

        "votu_seed_checkv_quality":
            v[
                "seed_checkv_quality"
            ],

        "votu_seed_completeness":
            v[
                "seed_completeness"
            ],

        "votu_best_checkv_quality":
            v[
                "best_checkv_quality"
            ],

        "votu_best_completeness":
            v[
                "best_completeness"
            ],

        "votu_genomad_topologies":
            v[
                "genomad_topologies"
            ],

        "votu_genomad_terminal_taxa":
            v[
                "genomad_terminal_taxa"
            ],

        "votu_max_genomad_virus_score":
            v[
                "max_genomad_virus_score"
            ],

        "votu_max_genomad_hallmarks":
            v[
                "max_genomad_hallmarks"
            ],


        # ---------------- Recurrence ----------------
        "n_libraries_mapped_gt0":
            r[
                "n_libraries_mapped_gt0"
            ],

        "n_libraries_breadth1_ge_0.80":
            r[
                "n_libraries_breadth1_ge_0.80"
            ],

        "n_libraries_breadth1_ge_0.95":
            r[
                "n_libraries_breadth1_ge_0.95"
            ],

        "n_L_libraries_breadth1_ge_0.80":
            r[
                "n_L_libraries_breadth1_ge_0.80"
            ],

        "n_M_libraries_breadth1_ge_0.80":
            r[
                "n_M_libraries_breadth1_ge_0.80"
            ],

        "n_biological_units_breadth1_ge_0.80":
            r[
                "n_biological_units_breadth1_ge_0.80"
            ],

        "n_producers_breadth1_ge_0.80":
            r[
                "n_producers_breadth1_ge_0.80"
            ],


        # ---------------- Quantitative summary ----------------
        "breadth80_samples":
            ";".join(
                breadth80_samples
            ),

        "breadth95_samples":
            ";".join(
                breadth95_samples
            ),

        "max_normalized_depth_sample":
            max_depth_sample,

        "max_mean_depth_per_million_hostfree_reads":
            q[
                max_depth_sample
            ][
                "mean_depth_per_million_hostfree_reads"
            ],
    })


integrated_file = (
    outdir
    / "92H1_direct_votu_host_taxonomy.tsv"
)

with integrated_file.open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=list(
            integrated[
                0
            ].keys()
        ),
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        integrated
    )


# ============================================================
# 7. Per-library quantitative table
# ============================================================

quant_rows = []


for row in integrated:

    votu = row[
        "votu_id"
    ]

    mag = row[
        "final18_MAG"
    ]

    for sample in samples:

        q = quant[
            votu
        ][sample]

        quant_rows.append({

            "votu_id":
                votu,

            "final18_MAG":
                mag,

            "host_taxon_label":
                row[
                    "host_taxon_label"
                ],

            "host_genus":
                row[
                    "host_genus"
                ],

            "host_species":
                row[
                    "host_species"
                ],

            "sample":
                sample,

            "producer":
                sample[0],

            "biological_unit":
                sample[1],

            "time_index":
                sample[-1],

            "week":
                int(
                    sample[-1]
                ) - 1,

            **q,
        })


quant_file = (
    outdir
    / "92H1_direct_votu_host_quantification.tsv"
)

with quant_file.open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=list(
            quant_rows[
                0
            ].keys()
        ),
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        quant_rows
    )


if len(
    quant_rows
) != 31 * 18:

    raise RuntimeError(
        f"Quant rows={len(quant_rows)}; "
        "esperaba 558"
    )


# ============================================================
# 8. Summary by MAG
# ============================================================

by_mag = defaultdict(
    list
)

for row in integrated:

    by_mag[
        row[
            "final18_MAG"
        ]
    ].append(row)


mag_summary_rows = []


for mag, rows in sorted(
    by_mag.items()
):

    first = rows[0]

    mag_summary_rows.append({

        "final18_MAG":
            mag,

        "host_taxon_label":
            first[
                "host_taxon_label"
            ],

        "host_taxon_label_rank":
            first[
                "host_taxon_label_rank"
            ],

        "host_genus":
            first[
                "host_genus"
            ],

        "host_species":
            first[
                "host_species"
            ],

        "host_MAG_completeness":
            first[
                "host_MAG_completeness"
            ],

        "host_MAG_contamination":
            first[
                "host_MAG_contamination"
            ],

        "n_direct_vOTUs":
            len(rows),

        "n_provirus_predictions":
            sum(
                int(
                    x[
                        "n_provirus_predictions"
                    ]
                )
                for x in rows
            ),

        "n_direct_vOTUs_breadth80_both_producers":
            sum(
                int(
                    x[
                        "n_producers_breadth1_ge_0.80"
                    ]
                ) == 2
                for x in rows
            ),

        "direct_vOTUs":
            ";".join(
                sorted(
                    x[
                        "votu_id"
                    ]
                    for x in rows
                )
            ),
    })


mag_summary_file = (
    outdir
    / "92H1_summary_by_MAG.tsv"
)

with mag_summary_file.open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=list(
            mag_summary_rows[
                0
            ].keys()
        ),
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        mag_summary_rows
    )


# ============================================================
# 9. Summary by host taxon label
# ============================================================

by_taxon = defaultdict(
    list
)

for row in integrated:

    label = clean(
        row[
            "host_taxon_label"
        ]
    )

    if not label:
        label = (
            "UNRESOLVED_HOST_TAXON"
        )

    by_taxon[
        label
    ].append(row)


tax_summary_rows = []


for label, rows in sorted(
    by_taxon.items()
):

    mags = {
        x[
            "final18_MAG"
        ]
        for x in rows
    }

    votus = {
        x[
            "votu_id"
        ]
        for x in rows
    }

    tax_summary_rows.append({

        "host_taxon_label":
            label,

        "n_MAGs":
            len(
                mags
            ),

        "n_direct_vOTUs":
            len(
                votus
            ),

        "n_provirus_predictions":
            sum(
                int(
                    x[
                        "n_provirus_predictions"
                    ]
                )
                for x in rows
            ),

        "MAGs":
            ";".join(
                sorted(
                    mags
                )
            ),

        "vOTUs":
            ";".join(
                sorted(
                    votus
                )
            ),
    })


tax_summary_file = (
    outdir
    / "92H1_summary_by_host_taxon.tsv"
)

with tax_summary_file.open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=list(
            tax_summary_rows[
                0
            ].keys()
        ),
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        tax_summary_rows
    )


# ============================================================
# 10. Global summary
# ============================================================

quality_counts = Counter(
    x[
        "votu_best_checkv_quality"
    ]
    for x in integrated
)

host_rank_counts = Counter(
    x[
        "host_taxon_label_rank"
    ]
    for x in integrated
)


summary = [

    (
        "direct_votu_MAG_pairs",
        len(
            integrated
        ),
    ),

    (
        "unique_direct_vOTUs",
        len(
            {
                x[
                    "votu_id"
                ]
                for x in integrated
            }
        ),
    ),

    (
        "unique_direct_MAGs",
        len(
            {
                x[
                    "final18_MAG"
                ]
                for x in integrated
            }
        ),
    ),

    (
        "unique_host_taxon_labels",
        len(
            {
                x[
                    "host_taxon_label"
                ]
                for x in integrated
            }
        ),
    ),

    (
        "direct_vOTUs_breadth80_both_producers",
        sum(
            int(
                x[
                    "n_producers_breadth1_ge_0.80"
                ]
            ) == 2
            for x in integrated
        ),
    ),

    (
        "direct_vOTUs_breadth95_at_least_one_library",
        sum(
            int(
                x[
                    "n_libraries_breadth1_ge_0.95"
                ]
            ) >= 1
            for x in integrated
        ),
    ),

    (
        "viral_quality_complete",
        quality_counts.get(
            "Complete",
            0
        ),
    ),

    (
        "viral_quality_high_quality",
        quality_counts.get(
            "High-quality",
            0
        ),
    ),

    (
        "viral_quality_medium_quality",
        quality_counts.get(
            "Medium-quality",
            0
        ),
    ),

    (
        "viral_quality_low_quality",
        quality_counts.get(
            "Low-quality",
            0
        ),
    ),

    (
        "viral_quality_not_determined",
        quality_counts.get(
            "Not-determined",
            0
        ),
    ),
]


for rank, count in sorted(
    host_rank_counts.items()
):

    summary.append(
        (
            f"host_label_rank_{rank or 'blank'}",
            count,
        )
    )


summary_file = (
    outdir
    / "92H1_global_summary.tsv"
)

with summary_file.open(
    "w"
) as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for key, value in summary:

        fh.write(
            f"{key}\t{value}\n"
        )


# ============================================================
# 11. Metadata
# ============================================================

metadata = (
    outdir
    / "92H1_run_metadata.tsv"
)

with metadata.open(
    "w"
) as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "analysis\t"
        "direct_provirus_MAG_host_taxonomy_integration\n"
    )

    fh.write(
        "host_evidence\t"
        "exact_parent_contig_membership_in_final18_MAG\n"
    )

    fh.write(
        "host_taxonomy_source\t"
        "43_MAG_abundance_taxonomy/representative_MAG_taxonomy.tsv\n"
    )

    fh.write(
        "viral_quantification_mode\t"
        "MAPQ20_ID95\n"
    )

    fh.write(
        "direct_host_assignment_is_prediction_only\tNO\n"
    )

    fh.write(
        "MAG_taxonomy_uncertainty_retained\tYES\n"
    )

    fh.write(
        "presence_absence_threshold_finalized\tNO\n"
    )


# ============================================================
# Final validation
# ============================================================

required = [
    integrated_file,
    quant_file,
    mag_summary_file,
    tax_summary_file,
    summary_file,
    metadata,
]

for path in required:

    if (
        not path.is_file()
        or path.stat().st_size == 0
    ):

        raise RuntimeError(
            f"Salida vacía/ausente: {path}"
        )


print()
print("=" * 80)
print("92H1 COMPLETED")
print("=" * 80)

for key, value in summary:

    print(
        key,
        value,
        sep="\t"
    )

print()
print("=" * 80)
print("BY MAG")
print("=" * 80)

for x in mag_summary_rows:

    print(
        x[
            "final18_MAG"
        ],
        x[
            "host_taxon_label"
        ],
        x[
            "host_MAG_completeness"
        ],
        x[
            "host_MAG_contamination"
        ],
        x[
            "n_direct_vOTUs"
        ],
        x[
            "n_provirus_predictions"
        ],
        sep="\t"
    )

PY


# ======================================================================
# SHA256
# ======================================================================

find "$OUT" \
    -maxdepth 1 \
    -type f \
    ! -name '92H1_sha256.tsv' \
    ! -name '92H1_COMPLETE.ok' \
    -print0 \
| sort -z \
| xargs -0 sha256sum \
> "${OUT}/92H1_sha256.tsv"


printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"$SLURM_JOB_ID" \
> "$DONE"


echo
echo "======================================================================"
echo "GLOBAL"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92H1_global_summary.tsv"


echo
echo "======================================================================"
echo "BY MAG"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92H1_summary_by_MAG.tsv"


echo
echo "======================================================================"
echo "BY HOST TAXON"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92H1_summary_by_host_taxon.tsv"


echo
echo "======================================================================"
echo "DIRECT ASSOCIATIONS"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92H1_direct_votu_host_taxonomy.tsv"


echo
echo "======================================================================"
echo "92H1 FINALIZO CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "======================================================================"

