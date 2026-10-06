#!/bin/bash

#SBATCH --job-name=93b_plasmid_exact
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=93b_plasmid_exact_%j.out
#SBATCH --error=93b_plasmid_exact_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

GENOMAD="${ROOT}/85_genomad_coassemblies"
OUT="${ROOT}/95_mobilome/93B_plasmid_catalog"

mkdir -p "$OUT"

echo "======================================================================"
echo "93B - CONSOLIDACION EXACTA DE CANDIDATOS PLASMIDICOS"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"

python3 - "$GENOMAD" "$OUT" <<'PY'

import csv
import hashlib
import sys

from collections import Counter, defaultdict
from pathlib import Path


genomad = Path(sys.argv[1])
out = Path(sys.argv[2])

groups = [
    "L1", "L2", "L3",
    "M1", "M2", "M3",
]

EXPECTED = 13599


# ======================================================================
# Helpers
# ======================================================================

RC_TABLE = str.maketrans(
    "ACGTNacgtn",
    "TGCANtgcan"
)


def reverse_complement(seq):

    return seq.translate(
        RC_TABLE
    )[::-1]


def sha256_text(seq):

    return hashlib.sha256(
        seq.encode()
    ).hexdigest()


def canonical_orientation(seq):
    """
    Conservador:
    colapsa una secuencia y su reverse-complement,
    pero NO intenta rotaciones circulares.
    """

    rc = reverse_complement(seq)

    if rc < seq:
        return rc, "REVERSE_COMPLEMENT"

    return seq, "FORWARD"


def read_fasta(path):

    name = None
    chunks = []

    with path.open() as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if name is not None:

                    yield (
                        name,
                        "".join(chunks).upper(),
                    )

                name = (
                    line[1:]
                    .split()[0]
                )

                chunks = []

            else:

                chunks.append(line)


    if name is not None:

        yield (
            name,
            "".join(chunks).upper(),
        )


def empty_annotation(x):

    x = (
        ""
        if x is None
        else str(x).strip()
    )

    return x.lower() in {
        "",
        "na",
        "n/a",
        "none",
        "null",
        "0",
        "-",
    }


# ======================================================================
# Consolidate
# ======================================================================

records = []

sequence_by_id = {}

topology_counts = Counter()

conjugation_raw_counts = Counter()
amr_raw_counts = Counter()


for group in groups:

    summary = (
        genomad
        / group
        / f"{group}_summary"
        / f"{group}_plasmid_summary.tsv"
    )

    fasta = (
        genomad
        / group
        / f"{group}_summary"
        / f"{group}_plasmid.fna"
    )


    for path in [
        summary,
        fasta,
    ]:

        if not path.is_file():

            raise RuntimeError(
                f"Falta {path}"
            )


    # --------------------------------------------------------
    # FASTA
    # --------------------------------------------------------

    seqs = {
        name: seq
        for name, seq
        in read_fasta(fasta)
    }


    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    with summary.open() as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )


        required = {
            "seq_name",
            "length",
            "topology",
            "n_genes",
            "genetic_code",
            "plasmid_score",
            "fdr",
            "n_hallmarks",
            "marker_enrichment",
            "conjugation_genes",
            "amr_genes",
        }


        missing = (
            required
            - set(
                reader.fieldnames or []
            )
        )

        if missing:

            raise RuntimeError(
                f"{summary}: faltan "
                f"{sorted(missing)}"
            )


        seen = set()

        for row in reader:

            seq_id = row[
                "seq_name"
            ]

            if seq_id in seen:

                raise RuntimeError(
                    f"{summary}: duplicado {seq_id}"
                )

            seen.add(seq_id)


            if seq_id not in seqs:

                raise RuntimeError(
                    f"{summary}: {seq_id} "
                    "no está en FASTA"
                )


            seq = seqs[
                seq_id
            ]

            declared_length = int(
                row[
                    "length"
                ]
            )

            if len(seq) != declared_length:

                raise RuntimeError(
                    f"{seq_id}: FASTA len={len(seq)} "
                    f"summary len={declared_length}"
                )


            direct_hash = sha256_text(
                seq
            )

            canon_seq, orientation = (
                canonical_orientation(
                    seq
                )
            )

            canonical_hash = (
                sha256_text(
                    canon_seq
                )
            )


            producer = group[0]

            biological_unit = int(
                group[1]
            )


            topology = (
                row[
                    "topology"
                ].strip()
            )

            conjugation = (
                row[
                    "conjugation_genes"
                ].strip()
            )

            amr = (
                row[
                    "amr_genes"
                ].strip()
            )


            topology_counts[
                topology
            ] += 1


            if not empty_annotation(
                conjugation
            ):

                conjugation_raw_counts[
                    conjugation
                ] += 1


            if not empty_annotation(
                amr
            ):

                amr_raw_counts[
                    amr
                ] += 1


            rec = {

                "seq_name":
                    seq_id,

                "coassembly":
                    group,

                "producer":
                    producer,

                "biological_unit":
                    biological_unit,

                "length":
                    declared_length,

                "topology":
                    topology,

                "n_genes":
                    row[
                        "n_genes"
                    ],

                "genetic_code":
                    row[
                        "genetic_code"
                    ],

                "plasmid_score":
                    row[
                        "plasmid_score"
                    ],

                "fdr":
                    row[
                        "fdr"
                    ],

                "n_hallmarks":
                    row[
                        "n_hallmarks"
                    ],

                "marker_enrichment":
                    row[
                        "marker_enrichment"
                    ],

                "conjugation_genes":
                    conjugation,

                "amr_genes":
                    amr,

                "has_conjugation_annotation":
                    (
                        "NO"
                        if empty_annotation(
                            conjugation
                        )
                        else "YES"
                    ),

                "has_genomad_amr_annotation":
                    (
                        "NO"
                        if empty_annotation(
                            amr
                        )
                        else "YES"
                    ),

                "sequence_sha256":
                    direct_hash,

                "strand_invariant_sha256":
                    canonical_hash,

                "canonical_orientation":
                    orientation,
            }


            records.append(
                rec
            )

            sequence_by_id[
                seq_id
            ] = seq


        if set(seqs) != seen:

            extra = (
                set(seqs)
                - seen
            )

            raise RuntimeError(
                f"{group}: FASTA contiene "
                f"{len(extra)} secuencias "
                "sin fila summary"
            )


# ======================================================================
# Global consistency
# ======================================================================

if len(records) != EXPECTED:

    raise RuntimeError(
        f"records={len(records)}; "
        f"esperaba {EXPECTED}"
    )


ids = [
    x[
        "seq_name"
    ]
    for x in records
]

if len(ids) != len(
    set(ids)
):

    raise RuntimeError(
        "seq_name global duplicado"
    )


# ======================================================================
# Exact direct-sequence clusters
# ======================================================================

by_direct_hash = defaultdict(
    list
)

for x in records:

    by_direct_hash[
        x[
            "sequence_sha256"
        ]
    ].append(x)


direct_clusters = []

direct_members = []


for i, (
    seq_hash,
    members
) in enumerate(
    sorted(
        by_direct_hash.items()
    ),
    start=1
):

    cid = (
        f"PEXACT{i:06d}"
    )

    representative = sorted(
        x[
            "seq_name"
        ]
        for x in members
    )[0]


    coassemblies = sorted(
        {
            x[
                "coassembly"
            ]
            for x in members
        }
    )

    producers = sorted(
        {
            x[
                "producer"
            ]
            for x in members
        }
    )

    topologies = sorted(
        {
            x[
                "topology"
            ]
            for x in members
        }
    )


    direct_clusters.append({

        "exact_cluster_id":
            cid,

        "sequence_sha256":
            seq_hash,

        "length":
            members[
                0
            ][
                "length"
            ],

        "n_members":
            len(
                members
            ),

        "representative_seq_name":
            representative,

        "n_coassemblies":
            len(
                coassemblies
            ),

        "coassemblies":
            ";".join(
                coassemblies
            ),

        "n_producers":
            len(
                producers
            ),

        "producers":
            ";".join(
                producers
            ),

        "topologies":
            ";".join(
                topologies
            ),
    })


    for x in members:

        direct_members.append({

            "exact_cluster_id":
                cid,

            "seq_name":
                x[
                    "seq_name"
                ],

            "coassembly":
                x[
                    "coassembly"
                ],

            "producer":
                x[
                    "producer"
                ],

            "length":
                x[
                    "length"
                ],

            "is_exact_representative":
                (
                    "YES"
                    if x[
                        "seq_name"
                    ]
                    == representative
                    else "NO"
                ),
        })


# ======================================================================
# Strand-invariant exact clusters
#
# Same complete sequence in either orientation.
# Does NOT collapse circular rotations.
# ======================================================================

by_strand_hash = defaultdict(
    list
)

for x in records:

    by_strand_hash[
        x[
            "strand_invariant_sha256"
        ]
    ].append(x)


strand_clusters = []

strand_members = []


for i, (
    strand_hash,
    members
) in enumerate(
    sorted(
        by_strand_hash.items()
    ),
    start=1
):

    cid = (
        f"PSIEXACT{i:06d}"
    )


    # Longitudes must be identical.
    lengths = {
        x[
            "length"
        ]
        for x in members
    }

    if len(lengths) != 1:

        raise RuntimeError(
            f"{cid}: longitudes distintas"
        )


    representative = sorted(
        x[
            "seq_name"
        ]
        for x in members
    )[0]


    coassemblies = sorted(
        {
            x[
                "coassembly"
            ]
            for x in members
        }
    )

    producers = sorted(
        {
            x[
                "producer"
            ]
            for x in members
        }
    )


    n_forward = sum(
        x[
            "canonical_orientation"
        ]
        == "FORWARD"
        for x in members
    )

    n_rc = sum(
        x[
            "canonical_orientation"
        ]
        == "REVERSE_COMPLEMENT"
        for x in members
    )


    strand_clusters.append({

        "strand_exact_cluster_id":
            cid,

        "strand_invariant_sha256":
            strand_hash,

        "length":
            next(
                iter(
                    lengths
                )
            ),

        "n_members":
            len(
                members
            ),

        "representative_seq_name":
            representative,

        "n_coassemblies":
            len(
                coassemblies
            ),

        "coassemblies":
            ";".join(
                coassemblies
            ),

        "n_producers":
            len(
                producers
            ),

        "producers":
            ";".join(
                producers
            ),

        "n_forward_canonical":
            n_forward,

        "n_reverse_complement_canonical":
            n_rc,
    })


    for x in members:

        strand_members.append({

            "strand_exact_cluster_id":
                cid,

            "seq_name":
                x[
                    "seq_name"
                ],

            "coassembly":
                x[
                    "coassembly"
                ],

            "producer":
                x[
                    "producer"
                ],

            "length":
                x[
                    "length"
                ],

            "canonical_orientation":
                x[
                    "canonical_orientation"
                ],

            "is_strand_exact_representative":
                (
                    "YES"
                    if x[
                        "seq_name"
                    ]
                    == representative
                    else "NO"
                ),
        })


# ======================================================================
# Map seq_name -> cluster IDs
# ======================================================================

exact_for_seq = {}

for x in direct_members:

    exact_for_seq[
        x[
            "seq_name"
        ]
    ] = x[
        "exact_cluster_id"
    ]


strand_for_seq = {}

for x in strand_members:

    strand_for_seq[
        x[
            "seq_name"
        ]
    ] = x[
        "strand_exact_cluster_id"
    ]


for x in records:

    x[
        "exact_cluster_id"
    ] = exact_for_seq[
        x[
            "seq_name"
        ]
    ]

    x[
        "strand_exact_cluster_id"
    ] = strand_for_seq[
        x[
            "seq_name"
        ]
    ]


# ======================================================================
# Write consolidated catalog
# ======================================================================

catalog_file = (
    out
    / "93B_plasmid_catalog_all.tsv"
)

with catalog_file.open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=list(
            records[
                0
            ].keys()
        ),
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        records
    )


# ======================================================================
# Write exact tables
# ======================================================================

def write_table(
    path,
    rows
):

    with path.open(
        "w",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=list(
                rows[
                    0
                ].keys()
            ),
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(
            rows
        )


write_table(
    out
    / "93B_exact_sequence_clusters.tsv",
    direct_clusters,
)

write_table(
    out
    / "93B_exact_sequence_members.tsv",
    direct_members,
)

write_table(
    out
    / "93B_strand_exact_clusters.tsv",
    strand_clusters,
)

write_table(
    out
    / "93B_strand_exact_members.tsv",
    strand_members,
)


# ======================================================================
# Representative FASTA:
# one representative per strand-invariant exact cluster.
# ======================================================================

rep_fasta = (
    out
    / "93B_strand_exact_representatives.fna"
)

with rep_fasta.open(
    "w"
) as fh:

    for cluster in (
        strand_clusters
    ):

        cid = cluster[
            "strand_exact_cluster_id"
        ]

        seq_id = cluster[
            "representative_seq_name"
        ]

        seq = sequence_by_id[
            seq_id
        ]

        fh.write(
            f">{cid} "
            f"original={seq_id} "
            f"length={len(seq)}\n"
        )

        for i in range(
            0,
            len(seq),
            80
        ):

            fh.write(
                seq[
                    i:i+80
                ]
                + "\n"
            )


# ======================================================================
# Per-coassembly summary
# ======================================================================

coassembly_summary = []


for group in groups:

    rows = [
        x
        for x in records
        if x[
            "coassembly"
        ] == group
    ]

    lengths = sorted(
        x[
            "length"
        ]
        for x in rows
    )


    coassembly_summary.append({

        "coassembly":
            group,

        "predictions":
            len(
                rows
            ),

        "total_bp":
            sum(
                lengths
            ),

        "min_length":
            min(
                lengths
            ),

        "median_length":
            lengths[
                len(
                    lengths
                ) // 2
            ],

        "max_length":
            max(
                lengths
            ),

        "length_ge_5000":
            sum(
                x >= 5000
                for x in lengths
            ),

        "length_ge_10000":
            sum(
                x >= 10000
                for x in lengths
            ),

        "length_ge_20000":
            sum(
                x >= 20000
                for x in lengths
            ),

        "length_ge_50000":
            sum(
                x >= 50000
                for x in lengths
            ),

        "length_ge_100000":
            sum(
                x >= 100000
                for x in lengths
            ),

        "with_conjugation_annotation":
            sum(
                x[
                    "has_conjugation_annotation"
                ] == "YES"
                for x in rows
            ),

        "with_genomad_amr_annotation":
            sum(
                x[
                    "has_genomad_amr_annotation"
                ] == "YES"
                for x in rows
            ),
    })


write_table(
    out
    / "93B_summary_by_coassembly.tsv",
    coassembly_summary,
)


# ======================================================================
# Topology
# ======================================================================

with (
    out
    / "93B_topology_counts.tsv"
).open("w") as fh:

    fh.write(
        "topology\tcount\n"
    )

    for key, value in sorted(
        topology_counts.items(),
        key=lambda x: (
            -x[1],
            x[0],
        )
    ):

        fh.write(
            f"{key}\t{value}\n"
        )


# ======================================================================
# Raw conjugation / AMR annotation signatures
#
# These are geNomad annotations, not independent validation.
# ======================================================================

with (
    out
    / "93B_conjugation_annotation_counts.tsv"
).open("w") as fh:

    fh.write(
        "conjugation_genes_raw\tcount\n"
    )

    for key, value in sorted(
        conjugation_raw_counts.items(),
        key=lambda x: (
            -x[1],
            x[0],
        )
    ):

        fh.write(
            f"{key}\t{value}\n"
        )


with (
    out
    / "93B_genomad_amr_annotation_counts.tsv"
).open("w") as fh:

    fh.write(
        "amr_genes_raw\tcount\n"
    )

    for key, value in sorted(
        amr_raw_counts.items(),
        key=lambda x: (
            -x[1],
            x[0],
        )
    ):

        fh.write(
            f"{key}\t{value}\n"
        )


# ======================================================================
# Global summary
# ======================================================================

total_bp = sum(
    x[
        "length"
    ]
    for x in records
)

direct_multimember = sum(
    x[
        "n_members"
    ] > 1
    for x in direct_clusters
)

strand_multimember = sum(
    x[
        "n_members"
    ] > 1
    for x in strand_clusters
)

direct_cross_coassembly = sum(
    x[
        "n_coassemblies"
    ] > 1
    for x in direct_clusters
)

strand_cross_coassembly = sum(
    x[
        "n_coassemblies"
    ] > 1
    for x in strand_clusters
)

direct_both_producers = sum(
    x[
        "n_producers"
    ] > 1
    for x in direct_clusters
)

strand_both_producers = sum(
    x[
        "n_producers"
    ] > 1
    for x in strand_clusters
)


summary = [

    (
        "plasmid_predictions",
        len(
            records
        ),
    ),

    (
        "total_bp",
        total_bp,
    ),

    (
        "exact_sequence_clusters",
        len(
            direct_clusters
        ),
    ),

    (
        "exact_sequence_multimember_clusters",
        direct_multimember,
    ),

    (
        "exact_sequence_cross_coassembly_clusters",
        direct_cross_coassembly,
    ),

    (
        "exact_sequence_both_producers_clusters",
        direct_both_producers,
    ),

    (
        "strand_exact_clusters",
        len(
            strand_clusters
        ),
    ),

    (
        "strand_exact_multimember_clusters",
        strand_multimember,
    ),

    (
        "strand_exact_cross_coassembly_clusters",
        strand_cross_coassembly,
    ),

    (
        "strand_exact_both_producers_clusters",
        strand_both_producers,
    ),

    (
        "sequences_collapsed_by_exact",
        len(
            records
        )
        - len(
            direct_clusters
        ),
    ),

    (
        "sequences_collapsed_by_strand_exact",
        len(
            records
        )
        - len(
            strand_clusters
        ),
    ),

    (
        "additional_collapses_due_to_reverse_complement",
        len(
            direct_clusters
        )
        - len(
            strand_clusters
        ),
    ),

    (
        "with_conjugation_annotation",
        sum(
            x[
                "has_conjugation_annotation"
            ] == "YES"
            for x in records
        ),
    ),

    (
        "with_genomad_amr_annotation",
        sum(
            x[
                "has_genomad_amr_annotation"
            ] == "YES"
            for x in records
        ),
    ),
]


with (
    out
    / "93B_global_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for key, value in summary:

        fh.write(
            f"{key}\t{value}\n"
        )


# ======================================================================
# Metadata
# ======================================================================

with (
    out
    / "93B_run_metadata.tsv"
).open("w") as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "analysis\t"
        "conservative_exact_plasmid_candidate_consolidation\n"
    )

    fh.write(
        "input_predictions\t13599\n"
    )

    fh.write(
        "prediction_source\t"
        "geNomad_1.11.0\n"
    )

    fh.write(
        "primary_exact_definition\t"
        "identical_nucleotide_sequence\n"
    )

    fh.write(
        "strand_exact_definition\t"
        "identical_sequence_allowing_reverse_complement\n"
    )

    fh.write(
        "circular_rotation_equivalence_tested\tNO\n"
    )

    fh.write(
        "ANI_clustering_performed\tNO\n"
    )

    fh.write(
        "plasmid_OTUs_defined\tNO\n"
    )

    fh.write(
        "cross_coassembly_exact_recurrence_is_prevalence\tNO\n"
    )

    fh.write(
        "genomad_AMR_annotation_is_final_AMR_validation\tNO\n"
    )

    fh.write(
        "genomad_conjugation_annotation_is_transfer_proof\tNO\n"
    )


# ======================================================================
# Final validation
# ======================================================================

required = [
    out
    / "93B_plasmid_catalog_all.tsv",

    out
    / "93B_exact_sequence_clusters.tsv",

    out
    / "93B_exact_sequence_members.tsv",

    out
    / "93B_strand_exact_clusters.tsv",

    out
    / "93B_strand_exact_members.tsv",

    out
    / "93B_strand_exact_representatives.fna",

    out
    / "93B_summary_by_coassembly.tsv",

    out
    / "93B_topology_counts.tsv",

    out
    / "93B_conjugation_annotation_counts.tsv",

    out
    / "93B_genomad_amr_annotation_counts.tsv",

    out
    / "93B_global_summary.tsv",

    out
    / "93B_run_metadata.tsv",
]


for path in required:

    if (
        not path.is_file()
        or path.stat().st_size == 0
    ):

        raise RuntimeError(
            f"Salida ausente/vacía: {path}"
        )


# Representative FASTA count
nrep = 0

with rep_fasta.open() as fh:

    for line in fh:

        if line.startswith(">"):
            nrep += 1


if nrep != len(
    strand_clusters
):

    raise RuntimeError(
        f"FASTA reps={nrep}; "
        f"clusters={len(strand_clusters)}"
    )


print()
print("=" * 80)
print("93B GLOBAL")
print("=" * 80)

for key, value in summary:

    print(
        key,
        value,
        sep="\t"
    )


print()
print("=" * 80)
print("BY COASSEMBLY")
print("=" * 80)

for row in coassembly_summary:

    print(
        row[
            "coassembly"
        ],
        row[
            "predictions"
        ],
        row[
            "total_bp"
        ],
        row[
            "median_length"
        ],
        row[
            "max_length"
        ],
        row[
            "with_conjugation_annotation"
        ],
        row[
            "with_genomad_amr_annotation"
        ],
        sep="\t"
    )

PY


# ======================================================================
# SHA256 + marker
# ======================================================================

find "$OUT" \
    -maxdepth 1 \
    -type f \
    ! -name '93B_sha256.tsv' \
    ! -name '93B_COMPLETE.ok' \
    -print0 \
| sort -z \
| xargs -0 sha256sum \
> "${OUT}/93B_sha256.tsv"


printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"$SLURM_JOB_ID" \
> "${OUT}/93B_COMPLETE.ok"


echo
echo "======================================================================"
echo "GLOBAL SUMMARY"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93B_global_summary.tsv"


echo
echo "======================================================================"
echo "BY COASSEMBLY"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93B_summary_by_coassembly.tsv"


echo
echo "======================================================================"
echo "TOPOLOGIES"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93B_topology_counts.tsv"


echo
echo "======================================================================"
echo "TOP CONJUGATION ANNOTATIONS"
echo "======================================================================"

head -n 25 \
    "${OUT}/93B_conjugation_annotation_counts.tsv" \
| column -t -s $'\t'


echo
echo "======================================================================"
echo "TOP geNomad AMR ANNOTATIONS"
echo "======================================================================"

head -n 25 \
    "${OUT}/93B_genomad_amr_annotation_counts.tsv" \
| column -t -s $'\t'


echo
echo "======================================================================"
echo "93B FINALIZO CORRECTAMENTE."
echo "No se definieron plasmid OTUs."
echo "ES SEGURO SALIR."
echo "======================================================================"

