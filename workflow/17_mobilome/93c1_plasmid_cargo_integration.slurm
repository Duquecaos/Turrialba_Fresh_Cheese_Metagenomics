#!/bin/bash

#SBATCH --job-name=93c1_plasmid_cargo
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=01:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=93c1_plasmid_cargo_%j.out
#SBATCH --error=93c1_plasmid_cargo_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

CATALOG="${ROOT}/95_mobilome/93B_plasmid_catalog/93B_plasmid_catalog_all.tsv"

MAP="${ROOT}/95_mobilome/93C0_colocalization_audit/93C0_traceable_original_mapping.tsv"

TRACE="${ROOT}/22_final_representative_mags/representative18_traceable.fasta"
MAGDIR="${ROOT}/22_final_representative_mags/fastas"

REPL="${ROOT}/25_abricate_consolidated/plasmid_replicon_loci_90_90.tsv"
AMR18="${ROOT}/25_abricate_consolidated/amr_loci_90_90.tsv"
VIR18="${ROOT}/25_abricate_consolidated/virulence_loci_90_90.tsv"

COMMAMR="${ROOT}/80_resistome_community_consolidated/91E_community_AMR_loci.tsv"

MAGTAX="${ROOT}/43_MAG_abundance_taxonomy/representative_MAG_taxonomy.tsv"

OUT="${ROOT}/95_mobilome/93C1_plasmid_cargo"

mkdir -p "$OUT"


echo "======================================================================"
echo "93C1 - INTEGRACION PLASMIDIO / REPLICON / ARG / VIRULENCIA"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"


for F in \
    "$CATALOG" \
    "$MAP" \
    "$TRACE" \
    "$REPL" \
    "$AMR18" \
    "$VIR18" \
    "$COMMAMR" \
    "$MAGTAX"
do

    [[ -s "$F" ]] || {
        echo "ERROR: falta $F" >&2
        exit 2
    }

done


[[ -d "$MAGDIR" ]] || {
    echo "ERROR: falta $MAGDIR" >&2
    exit 3
}


python3 - \
    "$CATALOG" \
    "$MAP" \
    "$TRACE" \
    "$MAGDIR" \
    "$REPL" \
    "$AMR18" \
    "$VIR18" \
    "$COMMAMR" \
    "$MAGTAX" \
    "$OUT" <<'PY'

import csv
import sys

from collections import Counter, defaultdict
from pathlib import Path


(
    catalog_file,
    map_file,
    trace_fasta,
    magdir,
    replicon_file,
    amr18_file,
    vir18_file,
    community_amr_file,
    magtax_file,
    outdir,
) = map(Path, sys.argv[1:])


EXPECTED_CANDIDATES = 13599
EXPECTED_COMM_AMR = 565
EXPECTED_FINAL18_AMR = 29
EXPECTED_VIR = 11
EXPECTED_REPLICON = 7
EXPECTED_MAG_CONTIGS = 2463


# ======================================================================
# Helpers
# ======================================================================

def read_tsv(path):

    with path.open(
        errors="replace"
    ) as fh:

        return list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


def clean(value):

    if value is None:
        return ""

    return str(value).strip()


def read_fasta(path):

    result = {}

    name = None
    chunks = []

    with path.open() as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if name is not None:

                    if name in result:
                        raise RuntimeError(
                            f"FASTA duplicado: {name}"
                        )

                    result[name] = (
                        "".join(chunks)
                        .upper()
                    )


                name = (
                    line[1:]
                    .split()[0]
                )

                chunks = []

            else:

                chunks.append(
                    line
                )


    if name is not None:

        if name in result:
            raise RuntimeError(
                f"FASTA duplicado: {name}"
            )

        result[name] = (
            "".join(chunks)
            .upper()
        )


    return result


def interval_distance(
    a_start,
    a_end,
    b_start,
    b_end,
):

    a_start = int(a_start)
    a_end = int(a_end)

    b_start = int(b_start)
    b_end = int(b_end)


    if a_end < b_start:

        return (
            b_start
            - a_end
            - 1
        )


    if b_end < a_start:

        return (
            a_start
            - b_end
            - 1
        )


    return 0


def semicolon(values):

    vals = sorted(
        {
            clean(x)
            for x in values
            if clean(x)
        }
    )

    return ";".join(vals)


# ======================================================================
# 1. Validate traceable ↔ original nucleotide identity
# ======================================================================

mapping_rows = read_tsv(
    map_file
)


if len(mapping_rows) != EXPECTED_MAG_CONTIGS:

    raise RuntimeError(
        f"mapping rows={len(mapping_rows)}; "
        f"esperaba {EXPECTED_MAG_CONTIGS}"
    )


trace_sequences = read_fasta(
    trace_fasta
)


mag_files = sorted(
    list(
        magdir.glob("*.fa")
    )
    + list(
        magdir.glob("*.fna")
    )
    + list(
        magdir.glob("*.fasta")
    )
)


if len(mag_files) != 18:

    raise RuntimeError(
        f"MAG FASTA={len(mag_files)}; esperaba 18"
    )


original_sequences = {}


for fasta in mag_files:

    mag = fasta.stem

    seqs = read_fasta(
        fasta
    )

    for contig, sequence in seqs.items():

        key = (
            mag,
            contig,
        )

        if key in original_sequences:

            raise RuntimeError(
                f"Original duplicado: {key}"
            )

        original_sequences[
            key
        ] = sequence


trace_to_original = {}

original_to_mags = defaultdict(
    set
)

sequence_identity_fail = []


for row in mapping_rows:

    mag = row[
        "MAG"
    ]

    coassembly = row[
        "coassembly"
    ]

    trace = row[
        "traceable_contig"
    ]

    original = row[
        "original_contig"
    ]


    if trace in trace_to_original:

        raise RuntimeError(
            f"Traceable duplicado: {trace}"
        )


    t_seq = trace_sequences.get(
        trace
    )

    o_seq = original_sequences.get(
        (
            mag,
            original,
        )
    )


    if t_seq is None:

        raise RuntimeError(
            f"Falta trace sequence: {trace}"
        )


    if o_seq is None:

        raise RuntimeError(
            f"Falta original: {mag}/{original}"
        )


    if t_seq != o_seq:

        sequence_identity_fail.append(
            (
                mag,
                trace,
                original,
                len(t_seq),
                len(o_seq),
            )
        )


    trace_to_original[
        trace
    ] = {

        "MAG":
            mag,

        "coassembly":
            coassembly,

        "original_contig":
            original,

        "length":
            len(o_seq),
    }


    original_to_mags[
        (
            coassembly,
            original,
        )
    ].add(
        mag
    )


if sequence_identity_fail:

    raise RuntimeError(
        "Traceable/original sequence mismatch: "
        + repr(
            sequence_identity_fail[:5]
        )
    )


if len(
    trace_to_original
) != EXPECTED_MAG_CONTIGS:

    raise RuntimeError(
        "Trace mapping incompleto"
    )


# ======================================================================
# 2. MAG taxonomy
# ======================================================================

mag_tax_rows = read_tsv(
    magtax_file
)

mag_tax = {}


for row in mag_tax_rows:

    mag = row[
        "MAG"
    ]

    if mag in mag_tax:

        raise RuntimeError(
            f"MAG taxonomy duplicado: {mag}"
        )

    mag_tax[
        mag
    ] = row


# ======================================================================
# 3. geNomad plasmid candidate catalog
# ======================================================================

catalog_rows = read_tsv(
    catalog_file
)


if len(
    catalog_rows
) != EXPECTED_CANDIDATES:

    raise RuntimeError(
        f"catalog rows={len(catalog_rows)}; "
        f"esperaba {EXPECTED_CANDIDATES}"
    )


candidate_by_key = {}


for row in catalog_rows:

    seq_name = row[
        "seq_name"
    ]

    coassembly = row[
        "coassembly"
    ]


    expected_prefix = (
        coassembly
        + "__"
    )


    if not seq_name.startswith(
        expected_prefix
    ):

        raise RuntimeError(
            f"{seq_name}: prefijo no coincide "
            f"con {coassembly}"
        )


    original_contig = (
        seq_name[
            len(
                expected_prefix
            ):
        ]
    )


    key = (
        coassembly,
        original_contig,
    )


    if key in candidate_by_key:

        raise RuntimeError(
            f"Candidato duplicado por contig: {key}"
        )


    row[
        "original_contig"
    ] = original_contig

    candidate_by_key[
        key
    ] = row


# ======================================================================
# 4. Map final18 ABRicate loci back to original coassembly contigs
# ======================================================================

def map_final18_loci(
    path,
    expected_rows,
    label,
):

    rows = read_tsv(
        path
    )


    if len(rows) != expected_rows:

        raise RuntimeError(
            f"{label}: rows={len(rows)}; "
            f"esperaba {expected_rows}"
        )


    mapped = []


    for row in rows:

        trace = row[
            "contig"
        ]


        if trace not in trace_to_original:

            raise RuntimeError(
                f"{label}: traceable no mapeado "
                f"{trace}"
            )


        info = trace_to_original[
            trace
        ]


        if row[
            "MAG"
        ] != info[
            "MAG"
        ]:

            raise RuntimeError(
                f"{label}: MAG mismatch "
                f"{trace}"
            )


        new = dict(
            row
        )

        new[
            "traceable_contig"
        ] = trace

        new[
            "coassembly"
        ] = info[
            "coassembly"
        ]

        new[
            "original_contig"
        ] = info[
            "original_contig"
        ]

        mapped.append(
            new
        )


    return mapped


replicon_rows = map_final18_loci(
    replicon_file,
    EXPECTED_REPLICON,
    "PlasmidFinder",
)

amr18_rows = map_final18_loci(
    amr18_file,
    EXPECTED_FINAL18_AMR,
    "ABRicate-AMR",
)

vir_rows = map_final18_loci(
    vir18_file,
    EXPECTED_VIR,
    "VFDB",
)


# ======================================================================
# 5. Community AMR
# ======================================================================

community_amr_rows = read_tsv(
    community_amr_file
)


if len(
    community_amr_rows
) != EXPECTED_COMM_AMR:

    raise RuntimeError(
        f"community AMR rows="
        f"{len(community_amr_rows)}; "
        f"esperaba {EXPECTED_COMM_AMR}"
    )


# ======================================================================
# 6. Group all evidence by original assembly contig
# ======================================================================

replicons_by_key = defaultdict(
    list
)

amr18_by_key = defaultdict(
    list
)

vir_by_key = defaultdict(
    list
)

community_amr_by_key = defaultdict(
    list
)


for row in replicon_rows:

    key = (
        row[
            "coassembly"
        ],
        row[
            "original_contig"
        ],
    )

    replicons_by_key[
        key
    ].append(
        row
    )


for row in amr18_rows:

    key = (
        row[
            "coassembly"
        ],
        row[
            "original_contig"
        ],
    )

    amr18_by_key[
        key
    ].append(
        row
    )


for row in vir_rows:

    key = (
        row[
            "coassembly"
        ],
        row[
            "original_contig"
        ],
    )

    vir_by_key[
        key
    ].append(
        row
    )


for row in community_amr_rows:

    key = (
        row[
            "coassembly"
        ],
        row[
            "contig"
        ],
    )

    community_amr_by_key[
        key
    ].append(
        row
    )


# ======================================================================
# 7. Integrated candidate table
# ======================================================================

integrated = []

cargo_links = []

replicon_cargo_pairs = []


for key, candidate in sorted(
    candidate_by_key.items()
):

    (
        coassembly,
        original_contig,
    ) = key


    reps = replicons_by_key.get(
        key,
        []
    )

    comm_amr = community_amr_by_key.get(
        key,
        []
    )

    final_amr = amr18_by_key.get(
        key,
        []
    )

    virs = vir_by_key.get(
        key,
        []
    )


    mags = sorted(
        original_to_mags.get(
            key,
            set(),
        )
    )


    host_taxa = []

    for mag in mags:

        if mag in mag_tax:

            host_taxa.append(
                mag_tax[
                    mag
                ].get(
                    "taxon_label",
                    ""
                )
            )


    evidence = [
        "GE_NOMAD_PLASMID_PREDICTION"
    ]


    if reps:

        evidence.append(
            "PLASMIDFINDER_REPLICON_SAME_CONTIG"
        )


    if comm_amr:

        evidence.append(
            "AMRFINDERPLUS_ARG_SAME_CONTIG"
        )


    if final_amr:

        evidence.append(
            "ABRICATE_AMR_SAME_CONTIG"
        )


    if virs:

        evidence.append(
            "VFDB_VIRULENCE_SAME_CONTIG"
        )


    if (
        candidate[
            "has_conjugation_annotation"
        ]
        == "YES"
    ):

        evidence.append(
            "GENOMAD_CONJUGATION_ANNOTATION"
        )


    if (
        candidate[
            "has_genomad_amr_annotation"
        ]
        == "YES"
    ):

        evidence.append(
            "GENOMAD_AMR_ANNOTATION"
        )


    if mags:

        evidence.append(
            "FINAL18_MAG_CONTIG"
        )


    integrated.append({

        "seq_name":
            candidate[
                "seq_name"
            ],

        "coassembly":
            coassembly,

        "producer":
            candidate[
                "producer"
            ],

        "biological_unit":
            candidate[
                "biological_unit"
            ],

        "original_contig":
            original_contig,

        "length":
            candidate[
                "length"
            ],

        "topology":
            candidate[
                "topology"
            ],

        "plasmid_score":
            candidate[
                "plasmid_score"
            ],

        "fdr":
            candidate[
                "fdr"
            ],

        "n_hallmarks":
            candidate[
                "n_hallmarks"
            ],

        "strand_exact_cluster_id":
            candidate[
                "strand_exact_cluster_id"
            ],

        "conjugation_genes_genomad":
            candidate[
                "conjugation_genes"
            ],

        "genomad_amr_genes":
            candidate[
                "amr_genes"
            ],

        "has_conjugation_annotation":
            candidate[
                "has_conjugation_annotation"
            ],

        "has_genomad_amr_annotation":
            candidate[
                "has_genomad_amr_annotation"
            ],

        "n_plasmidfinder_replicons":
            len(
                reps
            ),

        "plasmidfinder_replicons":
            semicolon(
                x[
                    "gene_labels"
                ]
                for x in reps
            ),

        "n_community_AMR_loci":
            len(
                comm_amr
            ),

        "community_AMR_symbols":
            semicolon(
                x[
                    "element_symbol"
                ]
                for x in comm_amr
            ),

        "community_AMR_classes":
            semicolon(
                x[
                    "class"
                ]
                for x in comm_amr
            ),

        "n_final18_ABRicate_AMR_loci":
            len(
                final_amr
            ),

        "final18_ABRicate_AMR_labels":
            semicolon(
                x[
                    "gene_labels"
                ]
                for x in final_amr
            ),

        "n_VFDB_virulence_loci":
            len(
                virs
            ),

        "VFDB_virulence_labels":
            semicolon(
                x[
                    "gene_labels"
                ]
                for x in virs
            ),

        "n_final18_MAGs":
            len(
                mags
            ),

        "final18_MAGs":
            ";".join(
                mags
            ),

        "final18_MAG_taxa":
            semicolon(
                host_taxa
            ),

        "evidence_features":
            ";".join(
                evidence
            ),

        "replicon_supported":
            (
                "YES"
                if reps
                else "NO"
            ),

        "independent_AMR_colocalized":
            (
                "YES"
                if comm_amr
                else "NO"
            ),

        "virulence_colocalized":
            (
                "YES"
                if virs
                else "NO"
            ),
    })


    # --------------------------------------------------------
    # Long-form cargo evidence
    # --------------------------------------------------------

    for row in reps:

        cargo_links.append({

            "seq_name":
                candidate[
                    "seq_name"
                ],

            "coassembly":
                coassembly,

            "original_contig":
                original_contig,

            "evidence_source":
                "PlasmidFinder",

            "evidence_type":
                "plasmid_replicon",

            "locus_id":
                row[
                    "locus_id"
                ],

            "symbol":
                row[
                    "gene_labels"
                ],

            "start":
                row[
                    "start"
                ],

            "end":
                row[
                    "end"
                ],

            "strand":
                row[
                    "strand"
                ],

            "method_or_database":
                row[
                    "databases"
                ],

            "MAG":
                row[
                    "MAG"
                ],
        })


    for row in comm_amr:

        cargo_links.append({

            "seq_name":
                candidate[
                    "seq_name"
                ],

            "coassembly":
                coassembly,

            "original_contig":
                original_contig,

            "evidence_source":
                "AMRFinderPlus",

            "evidence_type":
                "AMR",

            "locus_id":
                row[
                    "community_locus_id"
                ],

            "symbol":
                row[
                    "element_symbol"
                ],

            "start":
                row[
                    "start"
                ],

            "end":
                row[
                    "stop"
                ],

            "strand":
                row[
                    "strand"
                ],

            "method_or_database":
                row[
                    "method"
                ],

            "MAG":
                "",
        })


    for row in final_amr:

        cargo_links.append({

            "seq_name":
                candidate[
                    "seq_name"
                ],

            "coassembly":
                coassembly,

            "original_contig":
                original_contig,

            "evidence_source":
                "ABRicate",

            "evidence_type":
                "AMR",

            "locus_id":
                row[
                    "locus_id"
                ],

            "symbol":
                row[
                    "gene_labels"
                ],

            "start":
                row[
                    "start"
                ],

            "end":
                row[
                    "end"
                ],

            "strand":
                row[
                    "strand"
                ],

            "method_or_database":
                row[
                    "databases"
                ],

            "MAG":
                row[
                    "MAG"
                ],
        })


    for row in virs:

        cargo_links.append({

            "seq_name":
                candidate[
                    "seq_name"
                ],

            "coassembly":
                coassembly,

            "original_contig":
                original_contig,

            "evidence_source":
                "VFDB_ABRicate",

            "evidence_type":
                "virulence",

            "locus_id":
                row[
                    "locus_id"
                ],

            "symbol":
                row[
                    "gene_labels"
                ],

            "start":
                row[
                    "start"
                ],

            "end":
                row[
                    "end"
                ],

            "strand":
                row[
                    "strand"
                ],

            "method_or_database":
                row[
                    "databases"
                ],

            "MAG":
                row[
                    "MAG"
                ],
        })


    # --------------------------------------------------------
    # Replicon ↔ ARG / virulence pairs on exact same contig
    # --------------------------------------------------------

    for rep in reps:

        for cargo in comm_amr:

            replicon_cargo_pairs.append({

                "seq_name":
                    candidate[
                        "seq_name"
                    ],

                "coassembly":
                    coassembly,

                "original_contig":
                    original_contig,

                "pair_type":
                    "REPLICON_COMMUNITY_AMR",

                "replicon_locus_id":
                    rep[
                        "locus_id"
                    ],

                "replicon_label":
                    rep[
                        "gene_labels"
                    ],

                "replicon_start":
                    rep[
                        "start"
                    ],

                "replicon_end":
                    rep[
                        "end"
                    ],

                "cargo_locus_id":
                    cargo[
                        "community_locus_id"
                    ],

                "cargo_symbol":
                    cargo[
                        "element_symbol"
                    ],

                "cargo_start":
                    cargo[
                        "start"
                    ],

                "cargo_end":
                    cargo[
                        "stop"
                    ],

                "distance_bp":
                    interval_distance(
                        rep[
                            "start"
                        ],
                        rep[
                            "end"
                        ],
                        cargo[
                            "start"
                        ],
                        cargo[
                            "stop"
                        ],
                    ),

                "MAG":
                    rep[
                        "MAG"
                    ],
            })


        for cargo in virs:

            replicon_cargo_pairs.append({

                "seq_name":
                    candidate[
                        "seq_name"
                    ],

                "coassembly":
                    coassembly,

                "original_contig":
                    original_contig,

                "pair_type":
                    "REPLICON_VIRULENCE",

                "replicon_locus_id":
                    rep[
                        "locus_id"
                    ],

                "replicon_label":
                    rep[
                        "gene_labels"
                    ],

                "replicon_start":
                    rep[
                        "start"
                    ],

                "replicon_end":
                    rep[
                        "end"
                    ],

                "cargo_locus_id":
                    cargo[
                        "locus_id"
                    ],

                "cargo_symbol":
                    cargo[
                        "gene_labels"
                    ],

                "cargo_start":
                    cargo[
                        "start"
                    ],

                "cargo_end":
                    cargo[
                        "end"
                    ],

                "distance_bp":
                    interval_distance(
                        rep[
                            "start"
                        ],
                        rep[
                            "end"
                        ],
                        cargo[
                            "start"
                        ],
                        cargo[
                            "end"
                        ],
                    ),

                "MAG":
                    rep[
                        "MAG"
                    ],
            })


# ======================================================================
# 8. Write main table
# ======================================================================

integrated_file = (
    outdir
    / "93C1_plasmid_candidate_integrated.tsv"
)


with integrated_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            integrated[
                0
            ].keys()
        ),
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(
        integrated
    )


# ======================================================================
# 9. Long cargo/evidence table
# ======================================================================

cargo_file = (
    outdir
    / "93C1_colocalized_loci_long.tsv"
)


if cargo_links:

    with cargo_file.open(
        "w",
        newline=""
    ) as fh:

        w = csv.DictWriter(
            fh,
            fieldnames=list(
                cargo_links[
                    0
                ].keys()
            ),
            delimiter="\t"
        )

        w.writeheader()
        w.writerows(
            cargo_links
        )

else:

    cargo_file.write_text(
        "seq_name\tcoassembly\toriginal_contig\t"
        "evidence_source\tevidence_type\tlocus_id\t"
        "symbol\tstart\tend\tstrand\t"
        "method_or_database\tMAG\n"
    )


# ======================================================================
# 10. Replicon/cargo pairs
# ======================================================================

pair_file = (
    outdir
    / "93C1_replicon_cargo_pairs.tsv"
)


pair_fields = [
    "seq_name",
    "coassembly",
    "original_contig",
    "pair_type",
    "replicon_locus_id",
    "replicon_label",
    "replicon_start",
    "replicon_end",
    "cargo_locus_id",
    "cargo_symbol",
    "cargo_start",
    "cargo_end",
    "distance_bp",
    "MAG",
]


with pair_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=pair_fields,
        delimiter="\t"
    )

    w.writeheader()

    if replicon_cargo_pairs:

        w.writerows(
            replicon_cargo_pairs
        )


# ======================================================================
# 11. Useful subsets
# ======================================================================

rep_supported = [
    x
    for x in integrated
    if x[
        "replicon_supported"
    ]
    == "YES"
]


arg_candidates = [
    x
    for x in integrated
    if x[
        "independent_AMR_colocalized"
    ]
    == "YES"
]


vir_candidates = [
    x
    for x in integrated
    if x[
        "virulence_colocalized"
    ]
    == "YES"
]


independent_support_subset = [
    x
    for x in integrated
    if (
        x[
            "replicon_supported"
        ] == "YES"
        or x[
            "independent_AMR_colocalized"
        ] == "YES"
        or x[
            "virulence_colocalized"
        ] == "YES"
    )
]


def write_subset(
    name,
    rows,
):

    path = (
        outdir
        / name
    )


    with path.open(
        "w",
        newline=""
    ) as fh:

        w = csv.DictWriter(
            fh,
            fieldnames=list(
                integrated[
                    0
                ].keys()
            ),
            delimiter="\t"
        )

        w.writeheader()

        if rows:

            w.writerows(
                rows
            )


write_subset(
    "93C1_replicon_supported_candidates.tsv",
    rep_supported,
)

write_subset(
    "93C1_community_ARG_candidates.tsv",
    arg_candidates,
)

write_subset(
    "93C1_virulence_candidates.tsv",
    vir_candidates,
)

write_subset(
    "93C1_candidates_with_independent_support_or_cargo.tsv",
    independent_support_subset,
)


# ======================================================================
# 12. Summary by coassembly
# ======================================================================

summary_by_group = []


for group in [
    "L1", "L2", "L3",
    "M1", "M2", "M3",
]:

    rows = [
        x
        for x in integrated
        if x[
            "coassembly"
        ] == group
    ]


    summary_by_group.append({

        "coassembly":
            group,

        "plasmid_candidates":
            len(
                rows
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

        "with_plasmidfinder_replicon":
            sum(
                x[
                    "replicon_supported"
                ] == "YES"
                for x in rows
            ),

        "with_community_AMR":
            sum(
                x[
                    "independent_AMR_colocalized"
                ] == "YES"
                for x in rows
            ),

        "with_VFDB_virulence":
            sum(
                x[
                    "virulence_colocalized"
                ] == "YES"
                for x in rows
            ),

        "within_final18_MAG_contig":
            sum(
                int(
                    x[
                        "n_final18_MAGs"
                    ]
                ) > 0
                for x in rows
            ),

        "replicon_and_community_AMR":
            sum(
                (
                    x[
                        "replicon_supported"
                    ] == "YES"
                )
                and (
                    x[
                        "independent_AMR_colocalized"
                    ] == "YES"
                )
                for x in rows
            ),

        "community_AMR_and_conjugation":
            sum(
                (
                    x[
                        "independent_AMR_colocalized"
                    ] == "YES"
                )
                and (
                    x[
                        "has_conjugation_annotation"
                    ] == "YES"
                )
                for x in rows
            ),
    })


group_file = (
    outdir
    / "93C1_summary_by_coassembly.tsv"
)


with group_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            summary_by_group[
                0
            ].keys()
        ),
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(
        summary_by_group
    )


# ======================================================================
# 13. Global summary
# ======================================================================

replicon_loci_on_candidates = sum(
    len(
        replicons_by_key.get(
            key,
            []
        )
    )
    for key in candidate_by_key
)

community_amr_loci_on_candidates = sum(
    len(
        community_amr_by_key.get(
            key,
            []
        )
    )
    for key in candidate_by_key
)

final18_amr_loci_on_candidates = sum(
    len(
        amr18_by_key.get(
            key,
            []
        )
    )
    for key in candidate_by_key
)

vir_loci_on_candidates = sum(
    len(
        vir_by_key.get(
            key,
            []
        )
    )
    for key in candidate_by_key
)


global_summary = [

    (
        "traceable_original_sequence_identity_PASS",
        EXPECTED_MAG_CONTIGS,
    ),

    (
        "plasmid_candidates_total",
        len(
            integrated
        ),
    ),

    (
        "candidates_with_conjugation_annotation",
        sum(
            x[
                "has_conjugation_annotation"
            ] == "YES"
            for x in integrated
        ),
    ),

    (
        "candidates_with_genomad_AMR_annotation",
        sum(
            x[
                "has_genomad_amr_annotation"
            ] == "YES"
            for x in integrated
        ),
    ),

    (
        "plasmidfinder_loci_total",
        len(
            replicon_rows
        ),
    ),

    (
        "plasmidfinder_loci_on_genomad_candidates",
        replicon_loci_on_candidates,
    ),

    (
        "candidates_with_plasmidfinder_replicon",
        len(
            rep_supported
        ),
    ),

    (
        "community_AMR_loci_total",
        len(
            community_amr_rows
        ),
    ),

    (
        "community_AMR_loci_on_genomad_candidates",
        community_amr_loci_on_candidates,
    ),

    (
        "candidates_with_community_AMR",
        len(
            arg_candidates
        ),
    ),

    (
        "final18_ABRicate_AMR_loci_total",
        len(
            amr18_rows
        ),
    ),

    (
        "final18_ABRicate_AMR_loci_on_genomad_candidates",
        final18_amr_loci_on_candidates,
    ),

    (
        "VFDB_loci_total",
        len(
            vir_rows
        ),
    ),

    (
        "VFDB_loci_on_genomad_candidates",
        vir_loci_on_candidates,
    ),

    (
        "candidates_with_VFDB_virulence",
        len(
            vir_candidates
        ),
    ),

    (
        "candidates_on_final18_MAG_contigs",
        sum(
            int(
                x[
                    "n_final18_MAGs"
                ]
            ) > 0
            for x in integrated
        ),
    ),

    (
        "candidates_replicon_plus_community_AMR",
        sum(
            (
                x[
                    "replicon_supported"
                ] == "YES"
            )
            and (
                x[
                    "independent_AMR_colocalized"
                ] == "YES"
            )
            for x in integrated
        ),
    ),

    (
        "candidates_replicon_plus_VFDB_virulence",
        sum(
            (
                x[
                    "replicon_supported"
                ] == "YES"
            )
            and (
                x[
                    "virulence_colocalized"
                ] == "YES"
            )
            for x in integrated
        ),
    ),

    (
        "candidates_replicon_plus_conjugation",
        sum(
            (
                x[
                    "replicon_supported"
                ] == "YES"
            )
            and (
                x[
                    "has_conjugation_annotation"
                ] == "YES"
            )
            for x in integrated
        ),
    ),

    (
        "candidates_community_AMR_plus_conjugation",
        sum(
            (
                x[
                    "independent_AMR_colocalized"
                ] == "YES"
            )
            and (
                x[
                    "has_conjugation_annotation"
                ] == "YES"
            )
            for x in integrated
        ),
    ),

    (
        "replicon_cargo_pairs",
        len(
            replicon_cargo_pairs
        ),
    ),

    (
        "candidates_with_independent_support_or_cargo",
        len(
            independent_support_subset
        ),
    ),
]


global_file = (
    outdir
    / "93C1_global_summary.tsv"
)


with global_file.open(
    "w"
) as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for key, value in (
        global_summary
    ):

        fh.write(
            f"{key}\t{value}\n"
        )


# ======================================================================
# 14. Metadata
# ======================================================================

metadata_file = (
    outdir
    / "93C1_run_metadata.tsv"
)


with metadata_file.open(
    "w"
) as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "analysis\t"
        "plasmid_candidate_independent_evidence_and_cargo_colocalization\n"
    )

    fh.write(
        "primary_candidate_source\t"
        "geNomad_1.11.0\n"
    )

    fh.write(
        "plasmid_replicon_source\t"
        "PlasmidFinder_via_ABRicate_identity90_coverage90\n"
    )

    fh.write(
        "community_AMR_source\t"
        "AMRFinderPlus_community_coassembly_loci\n"
    )

    fh.write(
        "virulence_source\t"
        "VFDB_via_ABRicate_final18_identity90_coverage90\n"
    )

    fh.write(
        "colocalization_definition\t"
        "same_original_coassembly_contig\n"
    )

    fh.write(
        "traceable_original_sequences_verified_identical\tYES\n"
    )

    fh.write(
        "genomad_plasmid_candidate_equals_confirmed_plasmid\tNO\n"
    )

    fh.write(
        "same_contig_ARG_equals_demonstrated_transfer\tNO\n"
    )

    fh.write(
        "plasmidfinder_negative_equals_nonplasmid\tNO\n"
    )

    fh.write(
        "genomad_conjugation_annotation_equals_demonstrated_conjugation\tNO\n"
    )


# ======================================================================
# 15. Validation
# ======================================================================

if len(
    integrated
) != EXPECTED_CANDIDATES:

    raise RuntimeError(
        "Integrated candidate count mismatch"
    )


required = [
    integrated_file,
    cargo_file,
    pair_file,
    outdir
    / "93C1_replicon_supported_candidates.tsv",
    outdir
    / "93C1_community_ARG_candidates.tsv",
    outdir
    / "93C1_virulence_candidates.tsv",
    outdir
    / "93C1_candidates_with_independent_support_or_cargo.tsv",
    group_file,
    global_file,
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
print("93C1 GLOBAL")
print("=" * 80)

for key, value in global_summary:

    print(
        key,
        value,
        sep="\t"
    )


print()
print("=" * 80)
print("BY COASSEMBLY")
print("=" * 80)

for row in summary_by_group:

    print(
        row[
            "coassembly"
        ],
        row[
            "plasmid_candidates"
        ],
        row[
            "with_plasmidfinder_replicon"
        ],
        row[
            "with_community_AMR"
        ],
        row[
            "with_VFDB_virulence"
        ],
        row[
            "with_conjugation_annotation"
        ],
        row[
            "replicon_and_community_AMR"
        ],
        row[
            "community_AMR_and_conjugation"
        ],
        sep="\t"
    )

PY


# ======================================================================
# Checksums
# ======================================================================

find "$OUT" \
    -maxdepth 1 \
    -type f \
    ! -name '93C1_sha256.tsv' \
    ! -name '93C1_COMPLETE.ok' \
    -print0 \
| sort -z \
| xargs -0 sha256sum \
> "${OUT}/93C1_sha256.tsv"


printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"$SLURM_JOB_ID" \
> "${OUT}/93C1_COMPLETE.ok"


echo
echo "======================================================================"
echo "GLOBAL"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93C1_global_summary.tsv"


echo
echo "======================================================================"
echo "BY COASSEMBLY"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93C1_summary_by_coassembly.tsv"


echo
echo "======================================================================"
echo "REPLICON-SUPPORTED CANDIDATES"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93C1_replicon_supported_candidates.tsv"


echo
echo "======================================================================"
echo "REPLICON-CARGO PAIRS"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93C1_replicon_cargo_pairs.tsv"


echo
echo "======================================================================"
echo "ARG-CARRYING CANDIDATES - PREVIEW"
echo "======================================================================"

{
    head -n 1 \
        "${OUT}/93C1_community_ARG_candidates.tsv"

    tail -n +2 \
        "${OUT}/93C1_community_ARG_candidates.tsv" \
    | head -n 30

} | column -t -s $'\t'


echo
echo "======================================================================"
echo "93C1 FINALIZO CORRECTAMENTE."
echo "Co-localizacion = mismo contig; no implica transferencia demostrada."
echo "ES SEGURO SALIR."
echo "======================================================================"

