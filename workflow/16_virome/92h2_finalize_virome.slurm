#!/bin/bash

#SBATCH --job-name=92h2_virome_final
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=92h2_virome_final_%j.out
#SBATCH --error=92h2_virome_final_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

VOTUMETA="${ROOT}/93_votu_coverage/integrated/92G4_votu_metadata.tsv"
RECURRENCE="${ROOT}/93_votu_coverage/integrated/92G4_mapq20_id95_library_recurrence.tsv"
HOST="${ROOT}/94_votu_host_taxonomy/direct_host/92H1_direct_votu_host_taxonomy.tsv"

OUT="${ROOT}/94_votu_host_taxonomy/final"

mkdir -p "$OUT"

DONE="${OUT}/92H2_COMPLETE.ok"

echo "======================================================================"
echo "92H2 - FINALIZACION DEL VIROMA"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"

for F in \
    "$VOTUMETA" \
    "$RECURRENCE" \
    "$HOST"
do
    [[ -s "$F" ]] || {
        echo "ERROR: falta $F" >&2
        exit 2
    }
done


python3 - \
    "$VOTUMETA" \
    "$RECURRENCE" \
    "$HOST" \
    "$OUT" <<'PY'

import csv
import sys

from collections import Counter, defaultdict
from pathlib import Path


(
    votumeta_file,
    recurrence_file,
    host_file,
    outdir,
) = map(Path, sys.argv[1:])


N_VOTU = 4491


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


def to_float(x):

    x = clean(x)

    if x in {
        "",
        "NA",
        "NaN",
        "nan",
    }:
        return None

    return float(x)


# ======================================================================
# 1. vOTU metadata
# ======================================================================

meta_rows = read_tsv(
    votumeta_file
)

if len(meta_rows) != N_VOTU:

    raise RuntimeError(
        f"vOTU metadata rows={len(meta_rows)}; "
        f"esperaba {N_VOTU}"
    )


meta = {}

for row in meta_rows:

    votu = row[
        "votu_id"
    ]

    if votu in meta:

        raise RuntimeError(
            f"vOTU duplicado: {votu}"
        )

    meta[votu] = row


expected_ids = {
    f"VOTU{i:05d}"
    for i in range(
        1,
        N_VOTU + 1
    )
}

if set(meta) != expected_ids:

    raise RuntimeError(
        "IDs vOTU no son exactamente "
        "VOTU00001..VOTU04491"
    )


# ======================================================================
# 2. Recurrence
# ======================================================================

rec_rows = read_tsv(
    recurrence_file
)

if len(rec_rows) != N_VOTU:

    raise RuntimeError(
        f"recurrence rows={len(rec_rows)}"
    )


rec = {}

for row in rec_rows:

    votu = row[
        "votu_id"
    ]

    if votu in rec:

        raise RuntimeError(
            f"recurrence duplicada: {votu}"
        )

    rec[votu] = row


if set(rec) != expected_ids:

    raise RuntimeError(
        "Recurrence no cubre los 4491 vOTUs"
    )


# ======================================================================
# 3. Direct host
# ======================================================================

host_rows = read_tsv(
    host_file
)

if len(host_rows) != 31:

    raise RuntimeError(
        f"direct host rows={len(host_rows)}; "
        "esperaba 31"
    )


host = {}

for row in host_rows:

    votu = row[
        "votu_id"
    ]

    if votu in host:

        raise RuntimeError(
            f"Más de un host directo para {votu}"
        )

    host[votu] = row


if not set(host).issubset(
    expected_ids
):

    raise RuntimeError(
        "Host table contiene vOTU desconocido"
    )


# ======================================================================
# 4. Master table
# ======================================================================

master = []


for votu in sorted(
    expected_ids
):

    m = meta[
        votu
    ]

    r = rec[
        votu
    ]

    h = host.get(
        votu
    )


    direct = (
        h is not None
    )


    if direct:

        host_status = (
            "DIRECT_PROVIRUS_MAG"
        )

        host_mag = h[
            "final18_MAG"
        ]

        host_taxon = h[
            "host_taxon_label"
        ]

        host_rank = h[
            "host_taxon_label_rank"
        ]

        host_phylum = h[
            "host_phylum"
        ]

        host_genus = h[
            "host_genus"
        ]

        host_species = h[
            "host_species"
        ]

        host_completeness = h[
            "host_MAG_completeness"
        ]

        host_contamination = h[
            "host_MAG_contamination"
        ]

        host_contigs = h[
            "parent_contigs"
        ]

        n_provirus_predictions = h[
            "n_provirus_predictions"
        ]

    else:

        host_status = (
            "UNASSIGNED"
        )

        host_mag = ""
        host_taxon = ""
        host_rank = ""
        host_phylum = ""
        host_genus = ""
        host_species = ""
        host_completeness = ""
        host_contamination = ""
        host_contigs = ""
        n_provirus_predictions = ""


    best_quality = m[
        "best_checkv_quality"
    ]

    length = int(
        m[
            "seed_length"
        ]
    )


    if best_quality in {
        "Complete",
        "High-quality",
    }:

        quality_tier = (
            "GENOME_LEVEL_STRONG"
        )

    elif best_quality == "Medium-quality":

        quality_tier = (
            "GENOME_LEVEL_INTERMEDIATE"
        )

    elif best_quality == "Low-quality":

        quality_tier = (
            "FRAGMENTARY_LOW"
        )

    else:

        quality_tier = (
            "QUALITY_UNDETERMINED"
        )


    # --------------------------------------------------------
    # Descriptive recurrence only.
    # NOT a biological presence/prevalence classification.
    # --------------------------------------------------------

    n_b80 = int(
        r[
            "n_libraries_breadth1_ge_0.80"
        ]
    )

    n_b95 = int(
        r[
            "n_libraries_breadth1_ge_0.95"
        ]
    )

    n_prod80 = int(
        r[
            "n_producers_breadth1_ge_0.80"
        ]
    )


    master.append({

        "votu_id":
            votu,

        # ----------------------------------------------------
        # Sequence / clustering
        # ----------------------------------------------------
        "seed_length_bp":
            length,

        "n_vclean_members":
            m[
                "n_vclean_members"
            ],

        "n_preclean_exact_clusters":
            m[
                "n_preclean_exact_clusters"
            ],

        "n_original_genomad_predictions":
            m[
                "n_original_genomad_predictions"
            ],


        # ----------------------------------------------------
        # Viral quality
        # ----------------------------------------------------
        "seed_checkv_quality":
            m[
                "seed_checkv_quality"
            ],

        "seed_completeness":
            m[
                "seed_completeness"
            ],

        "best_checkv_quality":
            best_quality,

        "best_completeness":
            m[
                "best_completeness"
            ],

        "quality_interpretation_tier":
            quality_tier,


        # ----------------------------------------------------
        # geNomad
        # ----------------------------------------------------
        "genomad_topologies":
            m[
                "genomad_topologies"
            ],

        "genomad_terminal_taxa":
            m[
                "genomad_terminal_taxa"
            ],

        "max_genomad_virus_score":
            m[
                "max_genomad_virus_score"
            ],

        "min_genomad_fdr":
            m[
                "min_genomad_fdr"
            ],

        "max_genomad_hallmarks":
            m[
                "max_genomad_hallmarks"
            ],


        # ----------------------------------------------------
        # Assembly origin
        # ----------------------------------------------------
        "n_source_coassemblies":
            m[
                "n_coassemblies"
            ],

        "source_coassemblies":
            m[
                "coassemblies"
            ],

        "n_source_producers":
            m[
                "n_producers"
            ],

        "source_producers":
            m[
                "producers"
            ],


        # ----------------------------------------------------
        # Read support MAPQ20_ID95
        # ----------------------------------------------------
        "n_libraries_mapped_gt0":
            r[
                "n_libraries_mapped_gt0"
            ],

        "n_libraries_breadth10":
            r[
                "n_libraries_breadth1_ge_0.10"
            ],

        "n_libraries_breadth50":
            r[
                "n_libraries_breadth1_ge_0.50"
            ],

        "n_libraries_breadth80":
            n_b80,

        "n_libraries_breadth95":
            n_b95,

        "n_L_libraries_breadth80":
            r[
                "n_L_libraries_breadth1_ge_0.80"
            ],

        "n_M_libraries_breadth80":
            r[
                "n_M_libraries_breadth1_ge_0.80"
            ],

        "n_L_libraries_breadth95":
            r[
                "n_L_libraries_breadth1_ge_0.95"
            ],

        "n_M_libraries_breadth95":
            r[
                "n_M_libraries_breadth1_ge_0.95"
            ],

        "n_biological_units_breadth80":
            r[
                "n_biological_units_breadth1_ge_0.80"
            ],

        "n_biological_units_breadth95":
            r[
                "n_biological_units_breadth1_ge_0.95"
            ],

        "n_producers_breadth80":
            n_prod80,

        "n_producers_breadth95":
            r[
                "n_producers_breadth1_ge_0.95"
            ],


        # ----------------------------------------------------
        # Direct host
        # ----------------------------------------------------
        "host_assignment_status":
            host_status,

        "direct_host_MAG":
            host_mag,

        "direct_host_taxon":
            host_taxon,

        "direct_host_rank":
            host_rank,

        "direct_host_phylum":
            host_phylum,

        "direct_host_genus":
            host_genus,

        "direct_host_species":
            host_species,

        "direct_host_MAG_completeness":
            host_completeness,

        "direct_host_MAG_contamination":
            host_contamination,

        "direct_host_parent_contigs":
            host_contigs,

        "direct_host_n_provirus_predictions":
            n_provirus_predictions,

        "direct_host_evidence":
            (
                "proviral_sequence_on_exact_contig_of_final18_MAG"
                if direct
                else ""
            ),

        "host_assignment_is_exclusive":
            (
                "NO"
                if direct
                else ""
            ),

        "host_assignment_scope":
            (
                "demonstrated_for_at_least_one_proviral_member"
                if direct
                else ""
            ),
    })


# ======================================================================
# 5. Write master
# ======================================================================

master_file = (
    outdir
    / "92H2_votu_master.tsv"
)

with master_file.open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=list(
            master[0].keys()
        ),
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        master
    )


if len(master) != N_VOTU:

    raise RuntimeError(
        f"master rows={len(master)}"
    )


# ======================================================================
# 6. Direct host subset
# ======================================================================

direct_subset = [
    x
    for x in master
    if x[
        "host_assignment_status"
    ]
    == "DIRECT_PROVIRUS_MAG"
]


direct_file = (
    outdir
    / "92H2_direct_host_votus.tsv"
)

with direct_file.open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=list(
            master[0].keys()
        ),
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        direct_subset
    )


if len(
    direct_subset
) != 31:

    raise RuntimeError(
        f"direct subset={len(direct_subset)}"
    )


# ======================================================================
# 7. Summary
# ======================================================================

quality_counts = Counter(
    x[
        "best_checkv_quality"
    ]
    for x in master
)

direct_quality = Counter(
    x[
        "best_checkv_quality"
    ]
    for x in direct_subset
)

terminal_taxa = Counter()

for x in master:

    taxa = clean(
        x[
            "genomad_terminal_taxa"
        ]
    )

    if not taxa:
        terminal_taxa[
            "UNCLASSIFIED"
        ] += 1
        continue

    # A vOTU may aggregate >1 member taxonomy.
    terminal_taxa[
        taxa
    ] += 1


host_taxa = Counter(
    x[
        "direct_host_taxon"
    ]
    for x in direct_subset
)


summary = [

    (
        "votu_total",
        len(
            master
        ),
    ),

    (
        "direct_host_votus",
        len(
            direct_subset
        ),
    ),

    (
        "unassigned_host_votus",
        len(
            master
        )
        - len(
            direct_subset
        ),
    ),

    (
        "direct_host_fraction",
        len(
            direct_subset
        )
        / len(
            master
        ),
    ),

    (
        "direct_host_unique_MAGs",
        len(
            {
                x[
                    "direct_host_MAG"
                ]
                for x in direct_subset
            }
        ),
    ),

    (
        "direct_host_unique_taxa",
        len(
            {
                x[
                    "direct_host_taxon"
                ]
                for x in direct_subset
            }
        ),
    ),

    (
        "direct_host_species_rank",
        sum(
            x[
                "direct_host_rank"
            ]
            == "species"
            for x in direct_subset
        ),
    ),

    (
        "direct_host_genus_rank",
        sum(
            x[
                "direct_host_rank"
            ]
            == "genus"
            for x in direct_subset
        ),
    ),

    (
        "direct_host_complete_or_HQ",
        sum(
            x[
                "best_checkv_quality"
            ]
            in {
                "Complete",
                "High-quality",
            }
            for x in direct_subset
        ),
    ),

    (
        "direct_host_medium_or_better",
        sum(
            x[
                "best_checkv_quality"
            ]
            in {
                "Complete",
                "High-quality",
                "Medium-quality",
            }
            for x in direct_subset
        ),
    ),

    (
        "direct_host_breadth80_both_producers",
        sum(
            int(
                x[
                    "n_producers_breadth80"
                ]
            ) == 2
            for x in direct_subset
        ),
    ),

    (
        "direct_host_breadth95_any_library",
        sum(
            int(
                x[
                    "n_libraries_breadth95"
                ]
            ) >= 1
            for x in direct_subset
        ),
    ),

    (
        "complete_all_catalog",
        quality_counts.get(
            "Complete",
            0
        ),
    ),

    (
        "high_quality_all_catalog",
        quality_counts.get(
            "High-quality",
            0
        ),
    ),

    (
        "medium_quality_all_catalog",
        quality_counts.get(
            "Medium-quality",
            0
        ),
    ),

    (
        "low_quality_all_catalog",
        quality_counts.get(
            "Low-quality",
            0
        ),
    ),

    (
        "not_determined_all_catalog",
        quality_counts.get(
            "Not-determined",
            0
        ),
    ),
]


summary_file = (
    outdir
    / "92H2_global_summary.tsv"
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


# ======================================================================
# 8. Direct host taxon counts
# ======================================================================

host_counts_file = (
    outdir
    / "92H2_direct_host_taxon_counts.tsv"
)

with host_counts_file.open(
    "w"
) as fh:

    fh.write(
        "host_taxon\tn_direct_vOTUs\n"
    )

    for taxon, count in sorted(
        host_taxa.items(),
        key=lambda x: (
            -x[1],
            x[0],
        )
    ):

        fh.write(
            f"{taxon}\t{count}\n"
        )


# ======================================================================
# 9. Viral terminal taxonomy counts at vOTU level
# ======================================================================

viral_tax_file = (
    outdir
    / "92H2_votu_terminal_taxonomy_counts.tsv"
)

with viral_tax_file.open(
    "w"
) as fh:

    fh.write(
        "genomad_terminal_taxa\t"
        "n_vOTUs\n"
    )

    for taxon, count in sorted(
        terminal_taxa.items(),
        key=lambda x: (
            -x[1],
            x[0],
        )
    ):

        fh.write(
            f"{taxon}\t{count}\n"
        )


# ======================================================================
# 10. Highly recurrent direct-host subset
#
# Descriptive only; not called "core".
# ======================================================================

recurrent_direct = [
    x
    for x in direct_subset
    if int(
        x[
            "n_libraries_breadth80"
        ]
    ) >= 4
]


recurrent_file = (
    outdir
    / "92H2_direct_host_recurrent_breadth80_ge4.tsv"
)

with recurrent_file.open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=list(
            master[0].keys()
        ),
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        recurrent_direct
    )


# ======================================================================
# 11. Metadata
# ======================================================================

metadata_file = (
    outdir
    / "92H2_run_metadata.tsv"
)

with metadata_file.open(
    "w"
) as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "analysis\t"
        "final_integrated_virome_master_table\n"
    )

    fh.write(
        "votu_definition\t"
        "ANI95_AFshorter85_operational_species_level\n"
    )

    fh.write(
        "quantitative_support_mode\t"
        "MAPQ20_ID95\n"
    )

    fh.write(
        "direct_host_definition\t"
        "provirus_parent_contig_exactly_present_in_final18_MAG\n"
    )

    fh.write(
        "direct_host_is_exclusive\tNO\n"
    )

    fh.write(
        "host_unassigned_does_not_mean_no_host\tYES\n"
    )

    fh.write(
        "breadth_recurrence_is_library_recurrence_not_prevalence\tYES\n"
    )

    fh.write(
        "low_quality_vOTUs_retained_in_inclusive_catalog\tYES\n"
    )

    fh.write(
        "presence_absence_threshold_formally_finalized\tNO\n"
    )


# ======================================================================
# 12. Validation
# ======================================================================

required = [
    master_file,
    direct_file,
    summary_file,
    host_counts_file,
    viral_tax_file,
    recurrent_file,
    metadata_file,
]

for path in required:

    if (
        not path.is_file()
        or path.stat().st_size == 0
    ):

        raise RuntimeError(
            f"Salida ausente/vacía: {path}"
        )


print()
print("=" * 80)
print("92H2 FINAL VIROME")
print("=" * 80)

for key, value in summary:

    print(
        key,
        value,
        sep="\t"
    )

print()
print("=" * 80)
print("DIRECT HOST COUNTS")
print("=" * 80)

for taxon, count in sorted(
    host_taxa.items(),
    key=lambda x: (
        -x[1],
        x[0],
    )
):

    print(
        taxon,
        count,
        sep="\t"
    )

print()
print(
    "recurrent_direct_breadth80_ge4",
    len(
        recurrent_direct
    ),
    sep="\t"
)

PY


# ======================================================================
# Checksums + marker
# ======================================================================

find "$OUT" \
    -maxdepth 1 \
    -type f \
    ! -name '92H2_sha256.tsv' \
    ! -name '92H2_COMPLETE.ok' \
    -print0 \
| sort -z \
| xargs -0 sha256sum \
> "${OUT}/92H2_sha256.tsv"


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
    "${OUT}/92H2_global_summary.tsv"


echo
echo "======================================================================"
echo "DIRECT HOST TAXA"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92H2_direct_host_taxon_counts.tsv"


echo
echo "======================================================================"
echo "TOP VIRAL TAXA"
echo "======================================================================"

head -n 25 \
    "${OUT}/92H2_votu_terminal_taxonomy_counts.tsv" \
| column -t -s $'\t'


echo
echo "======================================================================"
echo "RECURRENT DIRECT-HOST vOTUs"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92H2_direct_host_recurrent_breadth80_ge4.tsv"


echo
echo "======================================================================"
echo "92H2 FINALIZO CORRECTAMENTE."
echo "STEP 92 VIROMA CERRADO."
echo "ES SEGURO SALIR."
echo "======================================================================"

