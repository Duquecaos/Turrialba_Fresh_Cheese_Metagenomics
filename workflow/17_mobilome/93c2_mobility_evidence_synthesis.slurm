#!/bin/bash

#SBATCH --job-name=93c2_mobility_synth
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=93c2_mobility_synth_%j.out
#SBATCH --error=93c2_mobility_synth_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

C1="${ROOT}/95_mobilome/93C1_plasmid_cargo"
C0="${ROOT}/95_mobilome/93C0_colocalization_audit"
B="${ROOT}/95_mobilome/93B_plasmid_catalog"

REPL="${ROOT}/25_abricate_consolidated/plasmid_replicon_loci_90_90.tsv"

OUT="${ROOT}/95_mobilome/93C2_mobility_synthesis"

mkdir -p "$OUT"

echo "======================================================================"
echo "93C2 - VALIDACION Y SINTESIS DE EVIDENCIA DE MOVILIDAD"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"


for F in \
    "${C1}/93C1_COMPLETE.ok" \
    "${C1}/93C1_global_summary.tsv" \
    "${C1}/93C1_plasmid_candidate_integrated.tsv" \
    "${C1}/93C1_community_ARG_candidates.tsv" \
    "${C1}/93C1_replicon_supported_candidates.tsv" \
    "${C0}/93C0_traceable_original_mapping.tsv" \
    "${B}/93B_strand_exact_clusters.tsv" \
    "$REPL"
do
    [[ -s "$F" ]] || {
        echo "ERROR: falta $F" >&2
        exit 2
    }
done


python3 - \
    "$C1" \
    "$C0" \
    "$B" \
    "$REPL" \
    "$OUT" <<'PY'

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path


c1 = Path(sys.argv[1])
c0 = Path(sys.argv[2])
b = Path(sys.argv[3])
replicon_file = Path(sys.argv[4])
out = Path(sys.argv[5])


def read_tsv(path):

    with path.open(errors="replace") as fh:
        return list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


def clean(x):
    return "" if x is None else str(x).strip()


def intval(x):
    return int(float(clean(x) or 0))


def yes(row, field):
    return clean(row.get(field)) == "YES"


# =====================================================================
# 1. Validate 93C1 global results
# =====================================================================

global_rows = read_tsv(
    c1 / "93C1_global_summary.tsv"
)

global_values = {
    x["metric"]: x["value"]
    for x in global_rows
}

expected = {
    "traceable_original_sequence_identity_PASS": 2463,
    "plasmid_candidates_total": 13599,
    "plasmidfinder_loci_total": 7,
    "plasmidfinder_loci_on_genomad_candidates": 6,
    "candidates_with_plasmidfinder_replicon": 6,
    "community_AMR_loci_total": 565,
    "community_AMR_loci_on_genomad_candidates": 176,
    "candidates_with_community_AMR": 152,
    "final18_ABRicate_AMR_loci_total": 29,
    "final18_ABRicate_AMR_loci_on_genomad_candidates": 7,
    "VFDB_loci_total": 11,
    "VFDB_loci_on_genomad_candidates": 0,
    "candidates_with_VFDB_virulence": 0,
    "candidates_community_AMR_plus_conjugation": 22,
}

for key, value in expected.items():

    observed = intval(
        global_values.get(key)
    )

    if observed != value:
        raise RuntimeError(
            f"{key}: observed={observed}, expected={value}"
        )


# =====================================================================
# 2. Main integrated catalog
# =====================================================================

integrated = read_tsv(
    c1 / "93C1_plasmid_candidate_integrated.tsv"
)

if len(integrated) != 13599:
    raise RuntimeError(
        f"Integrated rows={len(integrated)}"
    )


candidate = {}

for row in integrated:

    key = (
        row["coassembly"],
        row["original_contig"],
    )

    if key in candidate:
        raise RuntimeError(
            f"Duplicate candidate key: {key}"
        )

    candidate[key] = row


# =====================================================================
# 3. Strand-exact recurrence metadata
# =====================================================================

cluster_rows = read_tsv(
    b / "93B_strand_exact_clusters.tsv"
)

cluster_meta = {
    x["strand_exact_cluster_id"]: x
    for x in cluster_rows
}


# =====================================================================
# 4. Traceable -> original
# =====================================================================

mapping = read_tsv(
    c0 / "93C0_traceable_original_mapping.tsv"
)

trace_map = {
    x["traceable_contig"]: x
    for x in mapping
}

if len(trace_map) != 2463:
    raise RuntimeError(
        f"Trace map rows={len(trace_map)}"
    )


# =====================================================================
# 5. Audit all 7 PlasmidFinder loci
# =====================================================================

replicons = read_tsv(
    replicon_file
)

if len(replicons) != 7:
    raise RuntimeError(
        f"Replicon loci={len(replicons)}"
    )


rep_audit = []


for row in replicons:

    trace = row["contig"]

    if trace not in trace_map:
        raise RuntimeError(
            f"Unmapped traceable contig: {trace}"
        )

    info = trace_map[trace]

    key = (
        info["coassembly"],
        info["original_contig"],
    )

    cand = candidate.get(key)

    rep_audit.append({

        "locus_id":
            row["locus_id"],

        "MAG":
            row["MAG"],

        "coassembly":
            info["coassembly"],

        "traceable_contig":
            trace,

        "original_contig":
            info["original_contig"],

        "replicon":
            row["gene_labels"],

        "identity":
            row["max_identity"],

        "coverage":
            row["max_coverage"],

        "genomad_plasmid_candidate":
            "YES" if cand else "NO",

        "genomad_seq_name":
            cand["seq_name"] if cand else "",

        "genomad_plasmid_score":
            cand["plasmid_score"] if cand else "",

        "genomad_topology":
            cand["topology"] if cand else "",

        "community_AMR_same_contig":
            (
                cand["independent_AMR_colocalized"]
                if cand
                else "NO"
            ),

        "conjugation_annotation":
            (
                cand["has_conjugation_annotation"]
                if cand
                else "NO"
            ),
    })


if sum(
    x["genomad_plasmid_candidate"] == "YES"
    for x in rep_audit
) != 6:

    raise RuntimeError(
        "Expected exactly 6/7 PlasmidFinder loci "
        "inside geNomad candidates"
    )


# =====================================================================
# 6. Build AMR mobility context table
# =====================================================================

amr_rows = [
    x
    for x in integrated
    if yes(
        x,
        "independent_AMR_colocalized"
    )
]

if len(amr_rows) != 152:
    raise RuntimeError(
        f"AMR candidate contexts={len(amr_rows)}"
    )


amr_contexts = []


for row in amr_rows:

    cid = row[
        "strand_exact_cluster_id"
    ]

    cm = cluster_meta.get(
        cid,
        {}
    )


    evidence = [
        "GE_NOMAD_PLASMID_PREDICTION",
        "AMRFINDERPLUS_ARG_SAME_CONTIG",
    ]


    if yes(
        row,
        "has_conjugation_annotation"
    ):
        evidence.append(
            "GENOMAD_CONJUGATION_ANNOTATION"
        )


    if intval(
        row[
            "n_final18_ABRicate_AMR_loci"
        ]
    ) > 0:
        evidence.append(
            "FINAL18_ABRICATE_AMR_SAME_CONTIG"
        )


    if yes(
        row,
        "has_genomad_amr_annotation"
    ):
        evidence.append(
            "GENOMAD_AMR_ANNOTATION"
        )


    if intval(
        row[
            "n_final18_MAGs"
        ]
    ) > 0:
        evidence.append(
            "FINAL18_MAG_CONTIG"
        )


    amr_contexts.append({

        "seq_name":
            row["seq_name"],

        "coassembly":
            row["coassembly"],

        "producer":
            row["producer"],

        "biological_unit":
            row["biological_unit"],

        "original_contig":
            row["original_contig"],

        "length":
            row["length"],

        "topology":
            row["topology"],

        "plasmid_score":
            row["plasmid_score"],

        "fdr":
            row["fdr"],

        "strand_exact_cluster_id":
            cid,

        "exact_cluster_n_members":
            cm.get(
                "n_members",
                ""
            ),

        "exact_cluster_n_coassemblies":
            cm.get(
                "n_coassemblies",
                ""
            ),

        "exact_cluster_coassemblies":
            cm.get(
                "coassemblies",
                ""
            ),

        "exact_cluster_n_producers":
            cm.get(
                "n_producers",
                ""
            ),

        "exact_cluster_producers":
            cm.get(
                "producers",
                ""
            ),

        "n_community_AMR_loci":
            row[
                "n_community_AMR_loci"
            ],

        "community_AMR_symbols":
            row[
                "community_AMR_symbols"
            ],

        "community_AMR_classes":
            row[
                "community_AMR_classes"
            ],

        "conjugation_genes_genomad":
            row[
                "conjugation_genes_genomad"
            ],

        "has_conjugation_annotation":
            row[
                "has_conjugation_annotation"
            ],

        "genomad_amr_genes":
            row[
                "genomad_amr_genes"
            ],

        "has_genomad_amr_annotation":
            row[
                "has_genomad_amr_annotation"
            ],

        "n_final18_ABRicate_AMR_loci":
            row[
                "n_final18_ABRicate_AMR_loci"
            ],

        "final18_ABRicate_AMR_labels":
            row[
                "final18_ABRicate_AMR_labels"
            ],

        "n_final18_MAGs":
            row[
                "n_final18_MAGs"
            ],

        "final18_MAGs":
            row[
                "final18_MAGs"
            ],

        "final18_MAG_taxa":
            row[
                "final18_MAG_taxa"
            ],

        "evidence_features":
            ";".join(
                evidence
            ),
    })


# =====================================================================
# 7. Multi-evidence AMR contexts
#
# No score/ranking:
# requires independent AMR plus at least one additional mobility/context
# feature.
# =====================================================================

multi = [

    x for x in amr_contexts

    if (
        x[
            "has_conjugation_annotation"
        ] == "YES"

        or intval(
            x[
                "n_final18_ABRicate_AMR_loci"
            ]
        ) > 0

        or intval(
            x[
                "n_final18_MAGs"
            ]
        ) > 0
    )
]


amr_conjugation = [

    x for x in amr_contexts

    if x[
        "has_conjugation_annotation"
    ] == "YES"
]


final18_amr_contexts = [

    x for x in amr_contexts

    if intval(
        x[
            "n_final18_ABRicate_AMR_loci"
        ]
    ) > 0
]


# =====================================================================
# 8. AMR symbol summary
# =====================================================================

symbol_data = defaultdict(
    lambda: {
        "candidates": set(),
        "coassemblies": set(),
        "producers": set(),
        "clusters": set(),
        "conjugation_candidates": set(),
        "final18_candidates": set(),
        "MAGs": set(),
        "taxa": set(),
    }
)


for row in amr_contexts:

    symbols = [
        x
        for x in row[
            "community_AMR_symbols"
        ].split(";")
        if x
    ]


    for symbol in symbols:

        x = symbol_data[
            symbol
        ]

        x[
            "candidates"
        ].add(
            row[
                "seq_name"
            ]
        )

        x[
            "coassemblies"
        ].add(
            row[
                "coassembly"
            ]
        )

        x[
            "producers"
        ].add(
            row[
                "producer"
            ]
        )

        x[
            "clusters"
        ].add(
            row[
                "strand_exact_cluster_id"
            ]
        )


        if row[
            "has_conjugation_annotation"
        ] == "YES":

            x[
                "conjugation_candidates"
            ].add(
                row[
                    "seq_name"
                ]
            )


        if intval(
            row[
                "n_final18_MAGs"
            ]
        ) > 0:

            x[
                "final18_candidates"
            ].add(
                row[
                    "seq_name"
                ]
            )


        for mag in row[
            "final18_MAGs"
        ].split(";"):

            if mag:
                x[
                    "MAGs"
                ].add(
                    mag
                )


        for taxon in row[
            "final18_MAG_taxa"
        ].split(";"):

            if taxon:
                x[
                    "taxa"
                ].add(
                    taxon
                )


symbol_summary = []


for symbol, x in symbol_data.items():

    symbol_summary.append({

        "AMR_symbol":
            symbol,

        "n_candidate_contigs":
            len(
                x[
                    "candidates"
                ]
            ),

        "n_exact_clusters":
            len(
                x[
                    "clusters"
                ]
            ),

        "n_coassemblies":
            len(
                x[
                    "coassemblies"
                ]
            ),

        "coassemblies":
            ";".join(
                sorted(
                    x[
                        "coassemblies"
                    ]
                )
            ),

        "n_producers":
            len(
                x[
                    "producers"
                ]
            ),

        "producers":
            ";".join(
                sorted(
                    x[
                        "producers"
                    ]
                )
            ),

        "n_candidates_with_conjugation":
            len(
                x[
                    "conjugation_candidates"
                ]
            ),

        "n_candidates_in_final18_MAGs":
            len(
                x[
                    "final18_candidates"
                ]
            ),

        "final18_MAGs":
            ";".join(
                sorted(
                    x[
                        "MAGs"
                    ]
                )
            ),

        "final18_taxa":
            ";".join(
                sorted(
                    x[
                        "taxa"
                    ]
                )
            ),
    })


symbol_summary.sort(
    key=lambda x: (
        -x[
            "n_candidate_contigs"
        ],
        x[
            "AMR_symbol"
        ],
    )
)


# =====================================================================
# 9. Evidence combination counts across all 13,599 candidates
# =====================================================================

patterns = Counter()


for row in integrated:

    pattern = (
        "replicon="
        + (
            "YES"
            if yes(
                row,
                "replicon_supported"
            )
            else "NO"
        )
        + ";community_AMR="
        + (
            "YES"
            if yes(
                row,
                "independent_AMR_colocalized"
            )
            else "NO"
        )
        + ";conjugation="
        + (
            "YES"
            if yes(
                row,
                "has_conjugation_annotation"
            )
            else "NO"
        )
        + ";genomad_AMR="
        + (
            "YES"
            if yes(
                row,
                "has_genomad_amr_annotation"
            )
            else "NO"
        )
        + ";final18_MAG="
        + (
            "YES"
            if intval(
                row[
                    "n_final18_MAGs"
                ]
            ) > 0
            else "NO"
        )
    )

    patterns[
        pattern
    ] += 1


# =====================================================================
# Writers
# =====================================================================

def write_table(
    path,
    rows,
    fields=None,
):

    with path.open(
        "w",
        newline=""
    ) as fh:

        if fields is None:

            if not rows:
                raise RuntimeError(
                    f"No rows for {path}"
                )

            fields = list(
                rows[0].keys()
            )


        w = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t"
        )

        w.writeheader()

        if rows:
            w.writerows(
                rows
            )


write_table(
    out
    / "93C2_replicon_audit.tsv",
    rep_audit,
)

write_table(
    out
    / "93C2_AMR_plasmid_contexts.tsv",
    amr_contexts,
)

write_table(
    out
    / "93C2_AMR_conjugation_contexts.tsv",
    amr_conjugation,
)

write_table(
    out
    / "93C2_final18_AMR_plasmid_contexts.tsv",
    final18_amr_contexts,
)

write_table(
    out
    / "93C2_multievidence_AMR_contexts.tsv",
    multi,
)

write_table(
    out
    / "93C2_AMR_symbol_summary.tsv",
    symbol_summary,
)


with (
    out
    / "93C2_evidence_pattern_counts.tsv"
).open("w") as fh:

    fh.write(
        "evidence_pattern\t"
        "n_candidates\n"
    )

    for pattern, count in sorted(
        patterns.items(),
        key=lambda x: (
            -x[1],
            x[0],
        )
    ):

        fh.write(
            f"{pattern}\t{count}\n"
        )


# =====================================================================
# 10. Global synthesis
# =====================================================================

exact_recurrent_amr = {

    x[
        "strand_exact_cluster_id"
    ]

    for x in amr_contexts

    if intval(
        x[
            "exact_cluster_n_coassemblies"
        ]
    ) > 1
}


both_producer_exact_amr = {

    x[
        "strand_exact_cluster_id"
    ]

    for x in amr_contexts

    if intval(
        x[
            "exact_cluster_n_producers"
        ]
    ) > 1
}


summary = [

    (
        "source_93C1_scheduler_job",
        "131807",
    ),

    (
        "source_93C1_scheduler_state",
        "FAILED_after_outputs",
    ),

    (
        "source_93C1_complete_marker_present",
        "YES",
    ),

    (
        "source_93C1_outputs_validated",
        "YES",
    ),

    (
        "plasmid_candidates_total",
        len(
            integrated
        ),
    ),

    (
        "PlasmidFinder_loci_total",
        len(
            rep_audit
        ),
    ),

    (
        "PlasmidFinder_on_geNomad_candidates",
        sum(
            x[
                "genomad_plasmid_candidate"
            ] == "YES"
            for x in rep_audit
        ),
    ),

    (
        "PlasmidFinder_not_on_geNomad_candidates",
        sum(
            x[
                "genomad_plasmid_candidate"
            ] == "NO"
            for x in rep_audit
        ),
    ),

    (
        "AMR_candidate_contigs",
        len(
            amr_contexts
        ),
    ),

    (
        "AMR_loci_on_candidate_contigs",
        176,
    ),

    (
        "AMR_candidates_with_conjugation_annotation",
        len(
            amr_conjugation
        ),
    ),

    (
        "AMR_candidates_with_final18_ABRicate_AMR",
        len(
            final18_amr_contexts
        ),
    ),

    (
        "final18_ABRicate_AMR_loci_on_candidates",
        7,
    ),

    (
        "AMR_multievidence_contexts",
        len(
            multi
        ),
    ),

    (
        "AMR_symbols_on_candidate_contigs",
        len(
            symbol_summary
        ),
    ),

    (
        "AMR_exact_clusters_recurrent_across_coassemblies",
        len(
            exact_recurrent_amr
        ),
    ),

    (
        "AMR_exact_clusters_observed_in_both_producers",
        len(
            both_producer_exact_amr
        ),
    ),

    (
        "VFDB_loci_on_geNomad_candidates",
        0,
    ),
]


with (
    out
    / "93C2_global_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for key, value in summary:

        fh.write(
            f"{key}\t{value}\n"
        )


# =====================================================================
# 11. Metadata
# =====================================================================

with (
    out
    / "93C2_run_metadata.tsv"
).open("w") as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "analysis\t"
        "mobilome_AMR_evidence_synthesis_and_93C1_validation\n"
    )

    fh.write(
        "93C1_results_recomputed\tNO\n"
    )

    fh.write(
        "same_contig_colocalization_equals_transfer_proof\tNO\n"
    )

    fh.write(
        "genomad_plasmid_prediction_equals_confirmed_plasmid\tNO\n"
    )

    fh.write(
        "conjugation_gene_annotation_equals_demonstrated_transfer\tNO\n"
    )

    fh.write(
        "cross_coassembly_exact_recurrence_equals_prevalence\tNO\n"
    )

    fh.write(
        "evidence_score_used\tNO\n"
    )


print()
print("=" * 78)
print("93C2 GLOBAL")
print("=" * 78)

for key, value in summary:

    print(
        key,
        value,
        sep="\t"
    )


print()
print("=" * 78)
print("PLASMIDFINDER AUDIT")
print("=" * 78)

for x in rep_audit:

    print(
        x[
            "MAG"
        ],
        x[
            "replicon"
        ],
        x[
            "original_contig"
        ],
        x[
            "genomad_plasmid_candidate"
        ],
        x[
            "genomad_plasmid_score"
        ],
        sep="\t"
    )


print()
print("=" * 78)
print("TOP AMR SYMBOLS")
print("=" * 78)

for x in symbol_summary[:30]:

    print(
        x[
            "AMR_symbol"
        ],
        x[
            "n_candidate_contigs"
        ],
        x[
            "n_coassemblies"
        ],
        x[
            "n_producers"
        ],
        x[
            "n_candidates_with_conjugation"
        ],
        x[
            "n_candidates_in_final18_MAGs"
        ],
        sep="\t"
    )

PY


# ======================================================================
# SHA256 + clean successful validation marker
# ======================================================================

find "$OUT" \
    -maxdepth 1 \
    -type f \
    ! -name '93C2_sha256.tsv' \
    ! -name '93C2_COMPLETE.ok' \
    -print0 \
| sort -z \
| xargs -0 sha256sum \
> "${OUT}/93C2_sha256.tsv"


printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"$SLURM_JOB_ID" \
> "${OUT}/93C2_COMPLETE.ok"


echo
echo "======================================================================"
echo "93C2 FINALIZO CORRECTAMENTE."
echo "93C1 fue validado sin recalcular sus resultados."
echo "ES SEGURO SALIR."
echo "======================================================================"

