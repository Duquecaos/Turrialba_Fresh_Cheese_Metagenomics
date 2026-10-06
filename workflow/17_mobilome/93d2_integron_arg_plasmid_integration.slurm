#!/bin/bash

#SBATCH --job-name=93d2_int_arg_plasmid
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=01:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=93d2_int_arg_plasmid_%j.out
#SBATCH --error=93d2_int_arg_plasmid_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

D1="${ROOT}/95_mobilome/93D1_integronfinder"

AMR="${ROOT}/80_resistome_community_consolidated/91E_community_AMR_loci.tsv"

PLASMIDS="${ROOT}/95_mobilome/93C1_plasmid_cargo/93C1_plasmid_candidate_integrated.tsv"

TRACE_MAP="${ROOT}/95_mobilome/93C0_colocalization_audit/93C0_traceable_original_mapping.tsv"

MAGTAX="${ROOT}/43_MAG_abundance_taxonomy/representative_MAG_taxonomy.tsv"

OUT="${ROOT}/95_mobilome/93D2_integron_ARG_plasmid"

mkdir -p "$OUT"

echo "======================================================================"
echo "93D2 - INTEGRON x ARG x PLASMID x MAG"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"


for F in \
    "$AMR" \
    "$PLASMIDS" \
    "$TRACE_MAP" \
    "$MAGTAX"
do
    [[ -s "$F" ]] || {
        echo "ERROR: falta $F" >&2
        exit 2
    }
done


for G in L1 L2 L3 M1 M2 M3
do
    [[ -s "${D1}/${G}/93D1_COMPLETE.ok" ]] || {
        echo "ERROR: 93D1 no está cerrado para $G" >&2
        exit 3
    }
done


python3 - \
    "$D1" \
    "$AMR" \
    "$PLASMIDS" \
    "$TRACE_MAP" \
    "$MAGTAX" \
    "$OUT" <<'PY'

import csv
import math
import sys

from collections import Counter, defaultdict
from pathlib import Path


D1 = Path(sys.argv[1])
AMR_FILE = Path(sys.argv[2])
PLASMID_FILE = Path(sys.argv[3])
TRACE_FILE = Path(sys.argv[4])
MAGTAX_FILE = Path(sys.argv[5])
OUT = Path(sys.argv[6])

GROUPS = [
    "L1", "L2", "L3",
    "M1", "M2", "M3",
]

EXPECTED_STRUCTURES = 375
EXPECTED_COMPLETE = 7
EXPECTED_IN0 = 10
EXPECTED_CALIN = 358
EXPECTED_ELEMENT_ROWS = 2088
EXPECTED_AMR = 565
EXPECTED_PLASMIDS = 13599


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


def clean(x):

    if x is None:
        return ""

    return str(x).strip()


def as_int(x):

    return int(
        float(
            clean(x)
        )
    )


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
        return b_start - a_end - 1

    if b_end < a_start:
        return a_start - b_end - 1

    return 0


def join_unique(values):

    return ";".join(
        sorted(
            {
                clean(x)
                for x in values
                if clean(x)
            }
        )
    )


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
                    f"No rows and no fields: {path}"
                )

            fields = list(
                rows[0].keys()
            )

        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t"
        )

        writer.writeheader()

        if rows:
            writer.writerows(rows)


# ======================================================================
# 1. Read IntegronFinder outputs
# ======================================================================

element_rows = []

for group in GROUPS:

    candidates = sorted(
        D1.glob(
            f"{group}/**/*.integrons"
        )
    )

    candidates = [
        x
        for x in candidates
        if x.is_file()
    ]

    if len(candidates) != 1:

        raise RuntimeError(
            f"{group}: esperaba 1 .integrons; "
            f"encontré {len(candidates)}"
        )

    path = candidates[0]

    with path.open(
        errors="replace"
    ) as fh:

        reader = csv.DictReader(
            (
                line
                for line in fh
                if (
                    line.strip()
                    and not line.startswith("#")
                )
            ),
            delimiter="\t"
        )

        required = {
            "ID_integron",
            "ID_replicon",
            "element",
            "pos_beg",
            "pos_end",
            "strand",
            "evalue",
            "type_elt",
            "annotation",
            "model",
            "type",
            "default",
            "distance_2attC",
            "considered_topology",
        }

        missing = (
            required
            - set(
                reader.fieldnames or []
            )
        )

        if missing:

            raise RuntimeError(
                f"{group}: faltan columnas "
                f"{sorted(missing)}"
            )

        for row in reader:

            x = dict(row)

            x[
                "coassembly"
            ] = group

            x[
                "producer"
            ] = group[0]

            x[
                "biological_unit"
            ] = group[1]

            x[
                "pos_beg"
            ] = as_int(
                x["pos_beg"]
            )

            x[
                "pos_end"
            ] = as_int(
                x["pos_end"]
            )

            element_rows.append(x)


if len(element_rows) != EXPECTED_ELEMENT_ROWS:

    raise RuntimeError(
        f"IntegronFinder element rows="
        f"{len(element_rows)}; "
        f"esperaba {EXPECTED_ELEMENT_ROWS}"
    )


# ======================================================================
# 2. Build structures
# ======================================================================

structure_elements = defaultdict(
    list
)

for row in element_rows:

    key = (
        row["coassembly"],
        row["ID_replicon"],
        row["ID_integron"],
        row["type"],
    )

    structure_elements[
        key
    ].append(row)


if len(structure_elements) != EXPECTED_STRUCTURES:

    raise RuntimeError(
        f"structures={len(structure_elements)}; "
        f"esperaba {EXPECTED_STRUCTURES}"
    )


type_counts = Counter(
    key[3]
    for key in structure_elements
)


expected_types = {
    "complete": EXPECTED_COMPLETE,
    "In0": EXPECTED_IN0,
    "CALIN": EXPECTED_CALIN,
}

for typ, expected in expected_types.items():

    observed = type_counts.get(
        typ,
        0
    )

    if observed != expected:

        raise RuntimeError(
            f"{typ}: observed={observed}; "
            f"expected={expected}"
        )


# ======================================================================
# 3. Read plasmid candidate catalog
# ======================================================================

plasmids = read_tsv(
    PLASMID_FILE
)

if len(plasmids) != EXPECTED_PLASMIDS:

    raise RuntimeError(
        f"plasmids={len(plasmids)}"
    )


plasmid_by_key = {}

for row in plasmids:

    key = (
        row["coassembly"],
        row["original_contig"],
    )

    if key in plasmid_by_key:

        raise RuntimeError(
            f"Duplicate plasmid key {key}"
        )

    plasmid_by_key[
        key
    ] = row


# ======================================================================
# 4. Final18 MAG mapping
# ======================================================================

trace_rows = read_tsv(
    TRACE_FILE
)

mags_by_contig = defaultdict(
    set
)

for row in trace_rows:

    key = (
        row["coassembly"],
        row["original_contig"],
    )

    mags_by_contig[
        key
    ].add(
        row["MAG"]
    )


mag_tax_rows = read_tsv(
    MAGTAX_FILE
)

tax_by_mag = {
    row["MAG"]:
        row.get(
            "taxon_label",
            ""
        )
    for row in mag_tax_rows
}


# ======================================================================
# 5. Read community AMR
# ======================================================================

amr_rows = read_tsv(
    AMR_FILE
)

if len(amr_rows) != EXPECTED_AMR:

    raise RuntimeError(
        f"AMR rows={len(amr_rows)}; "
        f"esperaba {EXPECTED_AMR}"
    )


amr_by_contig = defaultdict(
    list
)

amr_by_protein = defaultdict(
    list
)

for row in amr_rows:

    key = (
        row["coassembly"],
        row["contig"],
    )

    amr_by_contig[
        key
    ].append(row)


    protein = clean(
        row.get(
            "protein_id"
        )
    )

    if protein:

        amr_by_protein[
            (
                row["coassembly"],
                protein,
            )
        ].append(row)


# ======================================================================
# 6. Construct integrated structure table
# ======================================================================

structures = []

structure_raw = {}

counter = 0


for key in sorted(
    structure_elements
):

    (
        group,
        contig,
        iid,
        typ,
    ) = key

    elements = structure_elements[
        key
    ]

    counter += 1

    sid = (
        f"INTSTRUCT{counter:06d}"
    )


    start = min(
        x["pos_beg"]
        for x in elements
    )

    end = max(
        x["pos_end"]
        for x in elements
    )


    attcs = [
        x
        for x in elements
        if clean(
            x["type_elt"]
        ) == "attC"
    ]

    proteins = [
        x
        for x in elements
        if clean(
            x["type_elt"]
        ) == "protein"
    ]

    promoters = [
        x
        for x in elements
        if clean(
            x["type_elt"]
        ) == "Promoter"
    ]


    protein_ids = {
        clean(
            x["element"]
        )
        for x in proteins
        if clean(
            x["element"]
        )
    }


    ckey = (
        group,
        contig,
    )

    plasmid = plasmid_by_key.get(
        ckey
    )

    mags = sorted(
        mags_by_contig.get(
            ckey,
            set(),
        )
    )

    taxa = sorted(
        {
            tax_by_mag.get(
                mag,
                ""
            )
            for mag in mags
            if tax_by_mag.get(
                mag,
                ""
            )
        }
    )


    amrs_same_contig = amr_by_contig.get(
        ckey,
        []
    )


    overlapping = []

    exact_protein = []


    for amr in amrs_same_contig:

        astart = as_int(
            amr["start"]
        )

        aend = as_int(
            amr["stop"]
        )

        if interval_distance(
            astart,
            aend,
            start,
            end,
        ) == 0:

            overlapping.append(
                amr
            )


        protein_id = clean(
            amr.get(
                "protein_id"
            )
        )

        if (
            protein_id
            and protein_id in protein_ids
        ):

            exact_protein.append(
                amr
            )


    combined_direct = {}

    for amr in (
        overlapping
        + exact_protein
    ):

        combined_direct[
            amr[
                "community_locus_id"
            ]
        ] = amr


    direct_amrs = list(
        combined_direct.values()
    )


    row = {

        "structure_id":
            sid,

        "coassembly":
            group,

        "producer":
            group[0],

        "biological_unit":
            group[1],

        "contig":
            contig,

        "ID_integron":
            iid,

        "integron_type":
            typ,

        "structure_start":
            start,

        "structure_end":
            end,

        "structure_length_bp":
            end - start + 1,

        "n_element_rows":
            len(
                elements
            ),

        "n_attC":
            len(
                attcs
            ),

        "n_proteins":
            len(
                proteins
            ),

        "n_promoters":
            len(
                promoters
            ),

        "protein_ids":
            ";".join(
                sorted(
                    protein_ids
                )
            ),

        "plasmid_candidate":
            (
                "YES"
                if plasmid
                else "NO"
            ),

        "plasmid_seq_name":
            (
                plasmid[
                    "seq_name"
                ]
                if plasmid
                else ""
            ),

        "plasmid_score":
            (
                plasmid[
                    "plasmid_score"
                ]
                if plasmid
                else ""
            ),

        "plasmid_topology":
            (
                plasmid[
                    "topology"
                ]
                if plasmid
                else ""
            ),

        "plasmid_conjugation_annotation":
            (
                plasmid[
                    "has_conjugation_annotation"
                ]
                if plasmid
                else "NO"
            ),

        "plasmid_conjugation_genes":
            (
                plasmid[
                    "conjugation_genes_genomad"
                ]
                if plasmid
                else ""
            ),

        "n_AMR_same_contig":
            len(
                amrs_same_contig
            ),

        "n_AMR_direct_structure":
            len(
                direct_amrs
            ),

        "n_AMR_exact_protein_ID":
            len(
                {
                    x[
                        "community_locus_id"
                    ]
                    for x in exact_protein
                }
            ),

        "AMR_symbols_direct":
            join_unique(
                x[
                    "element_symbol"
                ]
                for x in direct_amrs
            ),

        "AMR_classes_direct":
            join_unique(
                x[
                    "class"
                ]
                for x in direct_amrs
            ),

        "n_final18_MAGs":
            len(
                mags
            ),

        "final18_MAGs":
            ";".join(
                mags
            ),

        "final18_taxa":
            ";".join(
                taxa
            ),
    }


    structures.append(row)

    structure_raw[
        sid
    ] = {
        "row":
            row,

        "elements":
            elements,

        "attcs":
            attcs,

        "proteins":
            proteins,

        "protein_ids":
            protein_ids,
    }


# ======================================================================
# 7. ARG-to-structure links
# ======================================================================

direct_links = []

nearest_rows = []


structures_by_contig = defaultdict(
    list
)

for row in structures:

    structures_by_contig[
        (
            row[
                "coassembly"
            ],
            row[
                "contig"
            ],
        )
    ].append(row)


for amr in amr_rows:

    group = amr[
        "coassembly"
    ]

    contig = amr[
        "contig"
    ]

    astart = as_int(
        amr[
            "start"
        ]
    )

    aend = as_int(
        amr[
            "stop"
        ]
    )

    protein_id = clean(
        amr.get(
            "protein_id"
        )
    )


    candidates = structures_by_contig.get(
        (
            group,
            contig,
        ),
        []
    )


    nearest = None


    for structure in candidates:

        sid = structure[
            "structure_id"
        ]

        raw = structure_raw[
            sid
        ]


        dist = interval_distance(
            astart,
            aend,
            structure[
                "structure_start"
            ],
            structure[
                "structure_end"
            ],
        )


        exact_id = (
            protein_id
            and protein_id
            in raw[
                "protein_ids"
            ]
        )


        coord_overlap = (
            dist == 0
        )


        attc_distances = [
            interval_distance(
                astart,
                aend,
                x[
                    "pos_beg"
                ],
                x[
                    "pos_end"
                ],
            )
            for x in raw[
                "attcs"
            ]
        ]


        nearest_attc = (
            min(
                attc_distances
            )
            if attc_distances
            else None
        )


        if exact_id:

            evidence = (
                "EXACT_PROTEIN_ID_IN_STRUCTURE"
            )

        elif coord_overlap:

            evidence = (
                "COORDINATE_OVERLAP_STRUCTURE"
            )

        else:

            evidence = ""


        if evidence:

            direct_links.append({

                "community_locus_id":
                    amr[
                        "community_locus_id"
                    ],

                "coassembly":
                    group,

                "producer":
                    amr[
                        "producer"
                    ],

                "biological_unit":
                    amr[
                        "biological_unit"
                    ],

                "contig":
                    contig,

                "AMR_start":
                    astart,

                "AMR_end":
                    aend,

                "protein_id":
                    protein_id,

                "AMR_symbol":
                    amr[
                        "element_symbol"
                    ],

                "AMR_name":
                    amr[
                        "element_name"
                    ],

                "AMR_class":
                    amr[
                        "class"
                    ],

                "AMR_method":
                    amr[
                        "method"
                    ],

                "structure_id":
                    sid,

                "integron_type":
                    structure[
                        "integron_type"
                    ],

                "structure_start":
                    structure[
                        "structure_start"
                    ],

                "structure_end":
                    structure[
                        "structure_end"
                    ],

                "link_evidence":
                    evidence,

                "distance_to_structure_bp":
                    dist,

                "distance_to_nearest_attC_bp":
                    (
                        nearest_attc
                        if nearest_attc is not None
                        else ""
                    ),

                "n_attC":
                    structure[
                        "n_attC"
                    ],

                "plasmid_candidate":
                    structure[
                        "plasmid_candidate"
                    ],

                "plasmid_seq_name":
                    structure[
                        "plasmid_seq_name"
                    ],

                "plasmid_score":
                    structure[
                        "plasmid_score"
                    ],

                "plasmid_conjugation_annotation":
                    structure[
                        "plasmid_conjugation_annotation"
                    ],

                "plasmid_conjugation_genes":
                    structure[
                        "plasmid_conjugation_genes"
                    ],

                "final18_MAGs":
                    structure[
                        "final18_MAGs"
                    ],

                "final18_taxa":
                    structure[
                        "final18_taxa"
                    ],
            })


        near_key = (
            dist,
            0 if exact_id else 1,
            sid,
        )


        if (
            nearest is None
            or near_key
            < nearest[
                0
            ]
        ):

            nearest = (
                near_key,
                structure,
                nearest_attc,
                exact_id,
            )


    if nearest is None:

        nearest_rows.append({

            "community_locus_id":
                amr[
                    "community_locus_id"
                ],

            "coassembly":
                group,

            "contig":
                contig,

            "protein_id":
                protein_id,

            "AMR_symbol":
                amr[
                    "element_symbol"
                ],

            "AMR_class":
                amr[
                    "class"
                ],

            "integron_structure_same_contig":
                "NO",

            "nearest_structure_id":
                "",

            "nearest_integron_type":
                "",

            "distance_to_structure_bp":
                "",

            "distance_to_nearest_attC_bp":
                "",

            "exact_protein_ID_in_structure":
                "NO",
        })

    else:

        (
            near_key,
            structure,
            nearest_attc,
            exact_id,
        ) = nearest

        dist = near_key[
            0
        ]


        nearest_rows.append({

            "community_locus_id":
                amr[
                    "community_locus_id"
                ],

            "coassembly":
                group,

            "contig":
                contig,

            "protein_id":
                protein_id,

            "AMR_symbol":
                amr[
                    "element_symbol"
                ],

            "AMR_class":
                amr[
                    "class"
                ],

            "integron_structure_same_contig":
                "YES",

            "nearest_structure_id":
                structure[
                    "structure_id"
                ],

            "nearest_integron_type":
                structure[
                    "integron_type"
                ],

            "distance_to_structure_bp":
                dist,

            "distance_to_nearest_attC_bp":
                (
                    nearest_attc
                    if nearest_attc is not None
                    else ""
                ),

            "exact_protein_ID_in_structure":
                (
                    "YES"
                    if exact_id
                    else "NO"
                ),
        })


# ======================================================================
# 8. Deduplicate direct locus-structure pair
# ======================================================================

pair_seen = {}
dedup_links = []


for row in direct_links:

    key = (
        row[
            "community_locus_id"
        ],
        row[
            "structure_id"
        ],
    )


    # Prefer exact protein-ID evidence if duplicate
    if key not in pair_seen:

        pair_seen[
            key
        ] = row

    elif (
        row[
            "link_evidence"
        ]
        == "EXACT_PROTEIN_ID_IN_STRUCTURE"
    ):

        pair_seen[
            key
        ] = row


dedup_links = list(
    pair_seen.values()
)

dedup_links.sort(
    key=lambda x: (
        x[
            "coassembly"
        ],
        x[
            "contig"
        ],
        x[
            "structure_id"
        ],
        x[
            "community_locus_id"
        ],
    )
)


# ======================================================================
# 9. Important subsets
# ======================================================================

complete_structures = [
    x
    for x in structures
    if x[
        "integron_type"
    ] == "complete"
]


structures_with_arg = [
    x
    for x in structures
    if int(
        x[
            "n_AMR_direct_structure"
        ]
    ) > 0
]


plasmid_structures = [
    x
    for x in structures
    if x[
        "plasmid_candidate"
    ] == "YES"
]


triple_contexts = [
    x
    for x in structures
    if (
        x[
            "plasmid_candidate"
        ] == "YES"
        and int(
            x[
                "n_AMR_direct_structure"
            ]
        ) > 0
    )
]


complete_arg = [
    x
    for x in structures
    if (
        x[
            "integron_type"
        ] == "complete"
        and int(
            x[
                "n_AMR_direct_structure"
            ]
        ) > 0
    )
]


complete_arg_plasmid = [
    x
    for x in structures
    if (
        x[
            "integron_type"
        ] == "complete"
        and x[
            "plasmid_candidate"
        ] == "YES"
        and int(
            x[
                "n_AMR_direct_structure"
            ]
        ) > 0
    )
]


# ======================================================================
# 10. Summary by coassembly/type
# ======================================================================

summary_group_type = []


for group in GROUPS:

    for typ in [
        "complete",
        "In0",
        "CALIN",
    ]:

        rows = [
            x
            for x in structures
            if (
                x[
                    "coassembly"
                ] == group
                and x[
                    "integron_type"
                ] == typ
            )
        ]


        summary_group_type.append({

            "coassembly":
                group,

            "integron_type":
                typ,

            "structures":
                len(
                    rows
                ),

            "with_direct_AMR":
                sum(
                    int(
                        x[
                            "n_AMR_direct_structure"
                        ]
                    ) > 0
                    for x in rows
                ),

            "with_plasmid_candidate":
                sum(
                    x[
                        "plasmid_candidate"
                    ] == "YES"
                    for x in rows
                ),

            "with_AMR_and_plasmid":
                sum(
                    (
                        int(
                            x[
                                "n_AMR_direct_structure"
                            ]
                        ) > 0
                    )
                    and (
                        x[
                            "plasmid_candidate"
                        ] == "YES"
                    )
                    for x in rows
                ),

            "within_final18_MAG":
                sum(
                    int(
                        x[
                            "n_final18_MAGs"
                        ]
                    ) > 0
                    for x in rows
                ),
        })


# ======================================================================
# 11. AMR symbol summary
# ======================================================================

symbol_stats = defaultdict(
    lambda: {
        "loci": set(),
        "structures": set(),
        "types": set(),
        "groups": set(),
        "producers": set(),
        "plasmid_structures": set(),
        "exact_protein_loci": set(),
        "MAGs": set(),
        "taxa": set(),
    }
)


for row in dedup_links:

    symbol = row[
        "AMR_symbol"
    ]

    x = symbol_stats[
        symbol
    ]

    x[
        "loci"
    ].add(
        row[
            "community_locus_id"
        ]
    )

    x[
        "structures"
    ].add(
        row[
            "structure_id"
        ]
    )

    x[
        "types"
    ].add(
        row[
            "integron_type"
        ]
    )

    x[
        "groups"
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


    if row[
        "plasmid_candidate"
    ] == "YES":

        x[
            "plasmid_structures"
        ].add(
            row[
                "structure_id"
            ]
        )


    if (
        row[
            "link_evidence"
        ]
        == "EXACT_PROTEIN_ID_IN_STRUCTURE"
    ):

        x[
            "exact_protein_loci"
        ].add(
            row[
                "community_locus_id"
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
        "final18_taxa"
    ].split(";"):

        if taxon:
            x[
                "taxa"
            ].add(
                taxon
            )


symbol_summary = []


for symbol, x in (
    symbol_stats.items()
):

    symbol_summary.append({

        "AMR_symbol":
            symbol,

        "n_AMR_loci":
            len(
                x[
                    "loci"
                ]
            ),

        "n_integron_structures":
            len(
                x[
                    "structures"
                ]
            ),

        "integron_types":
            ";".join(
                sorted(
                    x[
                        "types"
                    ]
                )
            ),

        "n_coassemblies":
            len(
                x[
                    "groups"
                ]
            ),

        "coassemblies":
            ";".join(
                sorted(
                    x[
                        "groups"
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

        "n_loci_exact_protein_ID":
            len(
                x[
                    "exact_protein_loci"
                ]
            ),

        "n_plasmid_integron_structures":
            len(
                x[
                    "plasmid_structures"
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
            "n_AMR_loci"
        ],
        x[
            "AMR_symbol"
        ],
    )
)


# ======================================================================
# 12. Global metrics
# ======================================================================

unique_direct_AMR = {
    x[
        "community_locus_id"
    ]
    for x in dedup_links
}

unique_exact_protein_AMR = {
    x[
        "community_locus_id"
    ]
    for x in dedup_links
    if (
        x[
            "link_evidence"
        ]
        == "EXACT_PROTEIN_ID_IN_STRUCTURE"
    )
}

unique_coordinate_only_AMR = {
    x[
        "community_locus_id"
    ]
    for x in dedup_links
    if (
        x[
            "link_evidence"
        ]
        == "COORDINATE_OVERLAP_STRUCTURE"
    )
}


same_contig_AMR = {
    x[
        "community_locus_id"
    ]
    for x in nearest_rows
    if (
        x[
            "integron_structure_same_contig"
        ]
        == "YES"
    )
}


within_1000_not_direct = {
    x[
        "community_locus_id"
    ]
    for x in nearest_rows
    if (
        x[
            "integron_structure_same_contig"
        ]
        == "YES"
        and clean(
            x[
                "distance_to_structure_bp"
            ]
        )
        and int(
            x[
                "distance_to_structure_bp"
            ]
        ) > 0
        and int(
            x[
                "distance_to_structure_bp"
            ]
        ) <= 1000
    )
}


global_summary = [

    (
        "integron_element_rows",
        len(
            element_rows
        ),
    ),

    (
        "integron_structures_total",
        len(
            structures
        ),
    ),

    (
        "complete_integrons",
        type_counts[
            "complete"
        ],
    ),

    (
        "In0_structures",
        type_counts[
            "In0"
        ],
    ),

    (
        "CALIN_structures",
        type_counts[
            "CALIN"
        ],
    ),

    (
        "structures_with_plasmid_candidate",
        len(
            plasmid_structures
        ),
    ),

    (
        "structures_with_direct_AMR",
        len(
            structures_with_arg
        ),
    ),

    (
        "unique_AMR_loci_directly_in_integron_structure",
        len(
            unique_direct_AMR
        ),
    ),

    (
        "unique_AMR_loci_exact_protein_ID_in_structure",
        len(
            unique_exact_protein_AMR
        ),
    ),

    (
        "unique_AMR_loci_coordinate_overlap_only",
        len(
            unique_coordinate_only_AMR
        ),
    ),

    (
        "AMR_loci_on_contigs_with_integron_structure",
        len(
            same_contig_AMR
        ),
    ),

    (
        "AMR_loci_within_1kb_but_not_directly_in_structure",
        len(
            within_1000_not_direct
        ),
    ),

    (
        "structures_with_AMR_and_plasmid_candidate",
        len(
            triple_contexts
        ),
    ),

    (
        "complete_integrons_with_direct_AMR",
        len(
            complete_arg
        ),
    ),

    (
        "complete_integrons_with_AMR_and_plasmid_candidate",
        len(
            complete_arg_plasmid
        ),
    ),

    (
        "AMR_symbols_directly_associated_with_structures",
        len(
            symbol_summary
        ),
    ),
]


# ======================================================================
# 13. Write
# ======================================================================

write_table(
    OUT
    / "93D2_integron_structures_integrated.tsv",
    structures,
)

write_table(
    OUT
    / "93D2_AMR_integron_direct_links.tsv",
    dedup_links,
    fields=[
        "community_locus_id",
        "coassembly",
        "producer",
        "biological_unit",
        "contig",
        "AMR_start",
        "AMR_end",
        "protein_id",
        "AMR_symbol",
        "AMR_name",
        "AMR_class",
        "AMR_method",
        "structure_id",
        "integron_type",
        "structure_start",
        "structure_end",
        "link_evidence",
        "distance_to_structure_bp",
        "distance_to_nearest_attC_bp",
        "n_attC",
        "plasmid_candidate",
        "plasmid_seq_name",
        "plasmid_score",
        "plasmid_conjugation_annotation",
        "plasmid_conjugation_genes",
        "final18_MAGs",
        "final18_taxa",
    ],
)

write_table(
    OUT
    / "93D2_AMR_nearest_integron_structure.tsv",
    nearest_rows,
)

write_table(
    OUT
    / "93D2_complete_integrons.tsv",
    complete_structures,
)

write_table(
    OUT
    / "93D2_integron_structures_with_AMR.tsv",
    structures_with_arg,
)

write_table(
    OUT
    / "93D2_plasmid_integron_ARG_contexts.tsv",
    triple_contexts,
)

write_table(
    OUT
    / "93D2_complete_integron_ARG_contexts.tsv",
    complete_arg,
    fields=list(
        structures[
            0
        ].keys()
    ),
)

write_table(
    OUT
    / "93D2_complete_integron_plasmid_ARG_contexts.tsv",
    complete_arg_plasmid,
    fields=list(
        structures[
            0
        ].keys()
    ),
)

write_table(
    OUT
    / "93D2_summary_by_coassembly_type.tsv",
    summary_group_type,
)

write_table(
    OUT
    / "93D2_AMR_symbol_integron_summary.tsv",
    symbol_summary,
    fields=[
        "AMR_symbol",
        "n_AMR_loci",
        "n_integron_structures",
        "integron_types",
        "n_coassemblies",
        "coassemblies",
        "n_producers",
        "producers",
        "n_loci_exact_protein_ID",
        "n_plasmid_integron_structures",
        "final18_MAGs",
        "final18_taxa",
    ],
)


with (
    OUT
    / "93D2_global_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for key, value in global_summary:

        fh.write(
            f"{key}\t{value}\n"
        )


# ======================================================================
# 14. Metadata
# ======================================================================

with (
    OUT
    / "93D2_run_metadata.tsv"
).open("w") as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "integron_source\t"
        "IntegronFinder_2.0.6_local_max\n"
    )

    fh.write(
        "AMR_source\t"
        "AMRFinderPlus_community_565_loci\n"
    )

    fh.write(
        "primary_direct_link\t"
        "exact_AMRFinder_protein_id_matches_"
        "IntegronFinder_protein_element\n"
    )

    fh.write(
        "secondary_direct_link\t"
        "AMR_coordinates_overlap_operational_"
        "IntegronFinder_structure_span\n"
    )

    fh.write(
        "same_contig_without_overlap_is_integron_membership\tNO\n"
    )

    fh.write(
        "CALIN_equals_complete_integron\tNO\n"
    )

    fh.write(
        "In0_equals_complete_integron\tNO\n"
    )

    fh.write(
        "plasmid_source\t"
        "geNomad_plasmid_candidate\n"
    )

    fh.write(
        "plasmid_candidate_equals_confirmed_plasmid\tNO\n"
    )

    fh.write(
        "integron_ARG_colocalization_equals_transfer_proof\tNO\n"
    )

    fh.write(
        "cross_coassembly_structure_count_equals_prevalence\tNO\n"
    )


# ======================================================================
# 15. Validation
# ======================================================================

required = [
    OUT
    / "93D2_integron_structures_integrated.tsv",

    OUT
    / "93D2_AMR_integron_direct_links.tsv",

    OUT
    / "93D2_AMR_nearest_integron_structure.tsv",

    OUT
    / "93D2_complete_integrons.tsv",

    OUT
    / "93D2_integron_structures_with_AMR.tsv",

    OUT
    / "93D2_plasmid_integron_ARG_contexts.tsv",

    OUT
    / "93D2_complete_integron_ARG_contexts.tsv",

    OUT
    / "93D2_complete_integron_plasmid_ARG_contexts.tsv",

    OUT
    / "93D2_summary_by_coassembly_type.tsv",

    OUT
    / "93D2_AMR_symbol_integron_summary.tsv",

    OUT
    / "93D2_global_summary.tsv",

    OUT
    / "93D2_run_metadata.tsv",
]


for path in required:

    if (
        not path.is_file()
        or path.stat().st_size == 0
    ):

        raise RuntimeError(
            f"Salida faltante/vacía: {path}"
        )


print()
print("=" * 78)
print("93D2 GLOBAL")
print("=" * 78)

for key, value in global_summary:

    print(
        key,
        value,
        sep="\t"
    )


print()
print("=" * 78)
print("COMPLETE INTEGRONS")
print("=" * 78)

for row in complete_structures:

    print(
        row[
            "structure_id"
        ],
        row[
            "coassembly"
        ],
        row[
            "contig"
        ],
        row[
            "n_attC"
        ],
        row[
            "n_AMR_direct_structure"
        ],
        row[
            "AMR_symbols_direct"
        ],
        row[
            "plasmid_candidate"
        ],
        row[
            "final18_taxa"
        ],
        sep="\t"
    )


print()
print("=" * 78)
print("DIRECT ARG LINKS")
print("=" * 78)

for row in dedup_links:

    print(
        row[
            "AMR_symbol"
        ],
        row[
            "coassembly"
        ],
        row[
            "contig"
        ],
        row[
            "integron_type"
        ],
        row[
            "link_evidence"
        ],
        row[
            "distance_to_nearest_attC_bp"
        ],
        row[
            "plasmid_candidate"
        ],
        sep="\t"
    )

PY


# ======================================================================
# Checksums + marker
# ======================================================================

find "$OUT" \
    -maxdepth 1 \
    -type f \
    ! -name '93D2_sha256.tsv' \
    ! -name '93D2_COMPLETE.ok' \
    -print0 \
| sort -z \
| xargs -0 sha256sum \
> "${OUT}/93D2_sha256.tsv"


printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"$SLURM_JOB_ID" \
> "${OUT}/93D2_COMPLETE.ok"


echo
echo "======================================================================"
echo "GLOBAL"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93D2_global_summary.tsv"


echo
echo "======================================================================"
echo "93D2 FINALIZO CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "======================================================================"

