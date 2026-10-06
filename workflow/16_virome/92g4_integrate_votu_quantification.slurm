#!/bin/bash

#SBATCH --job-name=92g4_votu_integrate
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=92g4_votu_integrate_%j.out
#SBATCH --error=92g4_votu_integrate_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

META="${ROOT}/91_votu_catalog/92F_votu_summary.tsv"

COVROOT="${ROOT}/93_votu_coverage"
IDROOT="${COVROOT}/identity95"

OUT="${COVROOT}/integrated"
MATRIX="${OUT}/matrices"

mkdir -p \
    "$OUT" \
    "$MATRIX"

DONE="${OUT}/92G4_COMPLETE.ok"

echo "======================================================================"
echo "92G4 - INTEGRACION FINAL DE CUANTIFICACION vOTU"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"

[[ -s "$META" ]] || {
    echo "ERROR: falta $META" >&2
    exit 2
}


python3 - <<'PY'

import csv
import gzip
import os
import re
from collections import Counter, defaultdict
from pathlib import Path


# ======================================================================
# Paths
# ======================================================================

root = Path(
    f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"
)

meta_file = (
    root
    / "91_votu_catalog"
    / "92F_votu_summary.tsv"
)

covroot = (
    root
    / "93_votu_coverage"
)

idroot = (
    covroot
    / "identity95"
)

out = (
    covroot
    / "integrated"
)

matrix_dir = (
    out
    / "matrices"
)

out.mkdir(
    parents=True,
    exist_ok=True
)

matrix_dir.mkdir(
    parents=True,
    exist_ok=True
)


# ======================================================================
# Experimental design
# ======================================================================

samples = [
    "L1_1", "L1_2", "L1_3",
    "L2_1", "L2_2", "L2_3",
    "L3_1", "L3_2", "L3_3",
    "M1_1", "M1_2", "M1_3",
    "M2_1", "M2_2", "M2_3",
    "M3_1", "M3_2", "M3_3",
]

modes = [
    "ALL_PRIMARY",
    "MAPQ20",
    "MAPQ20_ID95",
]

N_VOTU = 4491
EXPECTED_ROWS = (
    N_VOTU
    * len(samples)
    * len(modes)
)


def sample_info(sample):

    m = re.fullmatch(
        r"([LM])([123])_([123])",
        sample
    )

    if not m:
        raise RuntimeError(
            f"ID de muestra inesperado: {sample}"
        )

    producer = m.group(1)
    biological_unit = int(
        m.group(2)
    )

    time_index = int(
        m.group(3)
    )

    # Codificación experimental:
    # _1 = semana 0
    # _2 = semana 1
    # _3 = semana 2
    week = (
        time_index
        - 1
    )

    return {
        "sample":
            sample,

        "producer":
            producer,

        "biological_unit":
            biological_unit,

        "time_index":
            time_index,

        "week":
            week,
    }


sample_meta = {
    s: sample_info(s)
    for s in samples
}


# ======================================================================
# Helpers
# ======================================================================

def clean(x):

    if x is None:
        return ""

    return str(x).strip()


def to_int(x):

    x = clean(x)

    if x == "":
        return None

    return int(
        float(x)
    )


def to_float(x):

    x = clean(x)

    if x == "":
        return None

    return float(x)


def read_tsv(path):

    with path.open() as fh:

        return list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


# ======================================================================
# vOTU metadata
# ======================================================================

meta_rows = read_tsv(
    meta_file
)

if len(meta_rows) != N_VOTU:

    raise RuntimeError(
        f"92F metadata rows={len(meta_rows)}; "
        f"esperaba {N_VOTU}"
    )


votu_meta = {}

for row in meta_rows:

    votu = row[
        "votu_id"
    ]

    if votu in votu_meta:

        raise RuntimeError(
            f"vOTU duplicado en 92F: {votu}"
        )

    votu_meta[votu] = row


expected_ids = {
    f"VOTU{i:05d}"
    for i in range(
        1,
        N_VOTU + 1
    )
}

if set(votu_meta) != expected_ids:

    raise RuntimeError(
        "Los IDs de 92F no corresponden "
        "exactamente a los 4491 VOTU."
    )


# ======================================================================
# Read quantitative data
#
# data[mode][sample][votu] = metrics
# ======================================================================

data = {
    mode: {}
    for mode in modes
}


def normalize_row(
    row,
    sample,
    mode,
):

    votu = row[
        "votu_id"
    ]

    length = to_int(
        row[
            "length_bp"
        ]
    )

    mapped = to_int(
        row[
            "mapped_records"
        ]
    )

    cov1 = to_int(
        row[
            "covered_bp_1x"
        ]
    )

    cov5 = to_int(
        row[
            "covered_bp_5x"
        ]
    )

    depth_sum = to_int(
        row[
            "depth_sum"
        ]
    )

    max_depth = to_int(
        row[
            "max_depth"
        ]
    )

    mean_depth = to_float(
        row[
            "mean_depth"
        ]
    )

    b1 = to_float(
        row[
            "breadth_1x"
        ]
    )

    b5 = to_float(
        row[
            "breadth_5x"
        ]
    )

    rpm_input = to_float(
        row[
            "mapped_records_per_million_input_read_records"
        ]
    )

    norm_depth_input = to_float(
        row[
            "mean_depth_per_million_input_read_records"
        ]
    )

    input_reads = to_int(
        row[
            "input_read_records"
        ]
    )

    mode_total = to_int(
        row[
            "mode_total_mapped_records"
        ]
    )


    required_numeric = {
        "length": length,
        "mapped": mapped,
        "cov1": cov1,
        "cov5": cov5,
        "depth_sum": depth_sum,
        "max_depth": max_depth,
        "mean_depth": mean_depth,
        "breadth1": b1,
        "breadth5": b5,
        "rpm_input": rpm_input,
        "norm_depth_input": norm_depth_input,
        "input_reads": input_reads,
        "mode_total": mode_total,
    }

    for key, value in (
        required_numeric.items()
    ):

        if value is None:

            raise RuntimeError(
                f"{sample}/{mode}/{votu}: "
                f"{key} ausente"
            )


    if not (
        0
        <= cov5
        <= cov1
        <= length
    ):

        raise RuntimeError(
            f"{sample}/{mode}/{votu}: "
            "coverage inconsistente"
        )


    return {
        "sample":
            sample,

        "mode":
            mode,

        "votu_id":
            votu,

        "length_bp":
            length,

        "mapped_records":
            mapped,

        "mapped_records_per_million_hostfree_reads":
            rpm_input,

        "covered_bp_1x":
            cov1,

        "breadth_1x":
            b1,

        "covered_bp_5x":
            cov5,

        "breadth_5x":
            b5,

        "depth_sum":
            depth_sum,

        "mean_depth":
            mean_depth,

        "max_depth":
            max_depth,

        "mean_depth_per_million_hostfree_reads":
            norm_depth_input,

        "input_hostfree_read_records":
            input_reads,

        "mode_total_mapped_records":
            mode_total,
    }


for sample in samples:

    previous_file = (
        covroot
        / sample
        / f"{sample}.votu_metrics.tsv"
    )

    id95_file = (
        idroot
        / sample
        / f"{sample}.MAPQ20_ID95.votu_metrics.tsv"
    )


    for path in [
        previous_file,
        id95_file,
    ]:

        if not path.is_file():

            raise RuntimeError(
                f"Falta {path}"
            )


    previous_rows = read_tsv(
        previous_file
    )

    if len(previous_rows) != (
        N_VOTU * 2
    ):

        raise RuntimeError(
            f"{sample}: 92G2 rows="
            f"{len(previous_rows)}"
        )


    for mode in [
        "ALL_PRIMARY",
        "MAPQ20",
    ]:

        selected = [
            r
            for r in previous_rows
            if r[
                "mode"
            ] == mode
        ]

        if len(selected) != N_VOTU:

            raise RuntimeError(
                f"{sample}/{mode}: "
                f"rows={len(selected)}"
            )


        d = {}

        for row in selected:

            x = normalize_row(
                row,
                sample,
                mode,
            )

            votu = x[
                "votu_id"
            ]

            if votu in d:

                raise RuntimeError(
                    f"{sample}/{mode}: "
                    f"duplicado {votu}"
                )

            d[votu] = x


        if set(d) != expected_ids:

            raise RuntimeError(
                f"{sample}/{mode}: "
                "IDs incompletos"
            )

        data[
            mode
        ][sample] = d


    id95_rows = read_tsv(
        id95_file
    )

    if len(id95_rows) != N_VOTU:

        raise RuntimeError(
            f"{sample}/ID95 rows="
            f"{len(id95_rows)}"
        )


    d = {}

    for row in id95_rows:

        x = normalize_row(
            row,
            sample,
            "MAPQ20_ID95",
        )

        votu = x[
            "votu_id"
        ]

        if votu in d:

            raise RuntimeError(
                f"{sample}/ID95 duplicate "
                f"{votu}"
            )

        d[votu] = x


    if set(d) != expected_ids:

        raise RuntimeError(
            f"{sample}/ID95 IDs incompletos"
        )

    data[
        "MAPQ20_ID95"
    ][sample] = d


# ======================================================================
# Validate length against canonical vOTU seed
# ======================================================================

for sample in samples:

    for mode in modes:

        for votu in expected_ids:

            observed = (
                data[
                    mode
                ][sample][votu][
                    "length_bp"
                ]
            )

            expected = int(
                votu_meta[
                    votu
                ][
                    "seed_length"
                ]
            )

            if observed != expected:

                raise RuntimeError(
                    f"{sample}/{mode}/{votu}: "
                    f"length {observed} != "
                    f"92F {expected}"
                )


# ======================================================================
# Validate nested filtering
#
# MAPQ20_ID95 ⊆ MAPQ20 ⊆ ALL_PRIMARY
# ======================================================================

nested_tests = 0

for sample in samples:

    for votu in expected_ids:

        allp = (
            data[
                "ALL_PRIMARY"
            ][sample][votu]
        )

        q20 = (
            data[
                "MAPQ20"
            ][sample][votu]
        )

        i95 = (
            data[
                "MAPQ20_ID95"
            ][sample][votu]
        )


        for field in [
            "mapped_records",
            "covered_bp_1x",
            "covered_bp_5x",
            "depth_sum",
        ]:

            if not (
                i95[field]
                <= q20[field]
                <= allp[field]
            ):

                raise RuntimeError(
                    f"Nesting FAIL "
                    f"{sample}/{votu}/{field}: "
                    f"{i95[field]} > "
                    f"{q20[field]} > "
                    f"{allp[field]}"
                )

            nested_tests += 1


# ======================================================================
# Sample metadata
# ======================================================================

sample_meta_file = (
    out
    / "92G4_sample_metadata.tsv"
)

with sample_meta_file.open(
    "w",
    newline=""
) as fh:

    fields = [
        "sample",
        "producer",
        "biological_unit",
        "time_index",
        "week",
    ]

    w = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    w.writeheader()

    for sample in samples:
        w.writerow(
            sample_meta[
                sample
            ]
        )


# ======================================================================
# vOTU metadata
# ======================================================================

votu_meta_file = (
    out
    / "92G4_votu_metadata.tsv"
)


votu_meta_fields = [
    "votu_id",
    "aniclust_seed",
    "seed_length",
    "n_vclean_members",
    "best_quality_member",
    "seed_checkv_quality",
    "seed_completeness",
    "best_checkv_quality",
    "best_completeness",
    "n_preclean_exact_clusters",
    "n_original_genomad_predictions",
    "n_coassemblies",
    "coassemblies",
    "n_producers",
    "producers",
    "genomad_topologies",
    "genomad_terminal_taxa",
    "max_genomad_virus_score",
    "min_genomad_fdr",
    "max_genomad_hallmarks",
    "best_quality_HQ_or_complete",
    "length_ge_2500",
    "length_ge_5000",
    "length_ge_10000",
    "length_ge_20000",
    "length_ge_50000",
]


with votu_meta_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=votu_meta_fields,
        delimiter="\t"
    )

    w.writeheader()


    for votu in sorted(
        expected_ids
    ):

        m = votu_meta[
            votu
        ]

        length = int(
            m[
                "seed_length"
            ]
        )

        bestq = m[
            "best_checkv_quality"
        ]


        row = {
            key:
                m.get(
                    key,
                    ""
                )
            for key in votu_meta_fields
        }


        row[
            "best_quality_HQ_or_complete"
        ] = (
            "YES"
            if bestq in {
                "Complete",
                "High-quality",
            }
            else "NO"
        )

        row[
            "length_ge_2500"
        ] = (
            "YES"
            if length >= 2500
            else "NO"
        )

        row[
            "length_ge_5000"
        ] = (
            "YES"
            if length >= 5000
            else "NO"
        )

        row[
            "length_ge_10000"
        ] = (
            "YES"
            if length >= 10000
            else "NO"
        )

        row[
            "length_ge_20000"
        ] = (
            "YES"
            if length >= 20000
            else "NO"
        )

        row[
            "length_ge_50000"
        ] = (
            "YES"
            if length >= 50000
            else "NO"
        )

        w.writerow(row)


# ======================================================================
# Long table
# ======================================================================

long_file = (
    out
    / "92G4_votu_sample_long.tsv.gz"
)


long_fields = [
    "sample",
    "producer",
    "biological_unit",
    "time_index",
    "week",
    "mode",
    "votu_id",
    "length_bp",
    "mapped_records",
    "mapped_records_per_million_hostfree_reads",
    "covered_bp_1x",
    "breadth_1x",
    "covered_bp_5x",
    "breadth_5x",
    "depth_sum",
    "mean_depth",
    "max_depth",
    "mean_depth_per_million_hostfree_reads",
    "input_hostfree_read_records",
    "mode_total_mapped_records",
    "best_checkv_quality",
    "best_completeness",
    "assembly_coassemblies",
    "assembly_producers",
    "genomad_topologies",
    "genomad_terminal_taxa",
]


n_long = 0


with gzip.open(
    long_file,
    "wt",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=long_fields,
        delimiter="\t"
    )

    w.writeheader()


    for sample in samples:

        sinfo = sample_meta[
            sample
        ]


        for mode in modes:

            for votu in sorted(
                expected_ids
            ):

                q = (
                    data[
                        mode
                    ][sample][votu]
                )

                m = (
                    votu_meta[
                        votu
                    ]
                )


                row = {
                    "sample":
                        sample,

                    "producer":
                        sinfo[
                            "producer"
                        ],

                    "biological_unit":
                        sinfo[
                            "biological_unit"
                        ],

                    "time_index":
                        sinfo[
                            "time_index"
                        ],

                    "week":
                        sinfo[
                            "week"
                        ],

                    "mode":
                        mode,

                    "votu_id":
                        votu,

                    "length_bp":
                        q[
                            "length_bp"
                        ],

                    "mapped_records":
                        q[
                            "mapped_records"
                        ],

                    "mapped_records_per_million_hostfree_reads":
                        q[
                            "mapped_records_per_million_hostfree_reads"
                        ],

                    "covered_bp_1x":
                        q[
                            "covered_bp_1x"
                        ],

                    "breadth_1x":
                        q[
                            "breadth_1x"
                        ],

                    "covered_bp_5x":
                        q[
                            "covered_bp_5x"
                        ],

                    "breadth_5x":
                        q[
                            "breadth_5x"
                        ],

                    "depth_sum":
                        q[
                            "depth_sum"
                        ],

                    "mean_depth":
                        q[
                            "mean_depth"
                        ],

                    "max_depth":
                        q[
                            "max_depth"
                        ],

                    "mean_depth_per_million_hostfree_reads":
                        q[
                            "mean_depth_per_million_hostfree_reads"
                        ],

                    "input_hostfree_read_records":
                        q[
                            "input_hostfree_read_records"
                        ],

                    "mode_total_mapped_records":
                        q[
                            "mode_total_mapped_records"
                        ],

                    "best_checkv_quality":
                        m[
                            "best_checkv_quality"
                        ],

                    "best_completeness":
                        m[
                            "best_completeness"
                        ],

                    "assembly_coassemblies":
                        m[
                            "coassemblies"
                        ],

                    "assembly_producers":
                        m[
                            "producers"
                        ],

                    "genomad_topologies":
                        m[
                            "genomad_topologies"
                        ],

                    "genomad_terminal_taxa":
                        m[
                            "genomad_terminal_taxa"
                        ],
                }

                w.writerow(row)

                n_long += 1


if n_long != EXPECTED_ROWS:

    raise RuntimeError(
        f"Long rows={n_long}; "
        f"esperaba {EXPECTED_ROWS}"
    )


# ======================================================================
# Sample × mode summary
# ======================================================================

sample_mode_summary = []


for sample in samples:

    all_total = (
        data[
            "ALL_PRIMARY"
        ][sample][
            "VOTU00001"
        ][
            "mode_total_mapped_records"
        ]
    )

    q20_total = (
        data[
            "MAPQ20"
        ][sample][
            "VOTU00001"
        ][
            "mode_total_mapped_records"
        ]
    )

    id95_total = (
        data[
            "MAPQ20_ID95"
        ][sample][
            "VOTU00001"
        ][
            "mode_total_mapped_records"
        ]
    )


    for mode in modes:

        rows = [
            data[
                mode
            ][sample][v]
            for v in sorted(
                expected_ids
            )
        ]


        total = rows[
            0
        ][
            "mode_total_mapped_records"
        ]


        if any(
            x[
                "mode_total_mapped_records"
            ] != total
            for x in rows
        ):

            raise RuntimeError(
                f"{sample}/{mode}: "
                "mode total varía entre vOTUs"
            )


        sample_mode_summary.append({

            "sample":
                sample,

            "producer":
                sample_meta[
                    sample
                ][
                    "producer"
                ],

            "biological_unit":
                sample_meta[
                    sample
                ][
                    "biological_unit"
                ],

            "week":
                sample_meta[
                    sample
                ][
                    "week"
                ],

            "mode":
                mode,

            "mapped_records":
                total,

            "retention_vs_ALL_PRIMARY":
                (
                    total / all_total
                    if all_total
                    else 0
                ),

            "retention_vs_MAPQ20":
                (
                    total / q20_total
                    if q20_total
                    else 0
                ),

            "votu_mapped_gt0":
                sum(
                    x[
                        "mapped_records"
                    ] > 0
                    for x in rows
                ),

            "breadth1_ge_0.10":
                sum(
                    x[
                        "breadth_1x"
                    ] >= 0.10
                    for x in rows
                ),

            "breadth1_ge_0.50":
                sum(
                    x[
                        "breadth_1x"
                    ] >= 0.50
                    for x in rows
                ),

            "breadth1_ge_0.80":
                sum(
                    x[
                        "breadth_1x"
                    ] >= 0.80
                    for x in rows
                ),

            "breadth1_ge_0.95":
                sum(
                    x[
                        "breadth_1x"
                    ] >= 0.95
                    for x in rows
                ),

            "breadth5_ge_0.10":
                sum(
                    x[
                        "breadth_5x"
                    ] >= 0.10
                    for x in rows
                ),

            "breadth5_ge_0.50":
                sum(
                    x[
                        "breadth_5x"
                    ] >= 0.50
                    for x in rows
                ),

            "breadth5_ge_0.80":
                sum(
                    x[
                        "breadth_5x"
                    ] >= 0.80
                    for x in rows
                ),

            "HQ_or_complete_breadth1_ge_0.80":
                sum(
                    (
                        x[
                            "breadth_1x"
                        ] >= 0.80
                    )
                    and (
                        votu_meta[
                            x[
                                "votu_id"
                            ]
                        ][
                            "best_checkv_quality"
                        ]
                        in {
                            "Complete",
                            "High-quality",
                        }
                    )
                    for x in rows
                ),
        })


summary_file = (
    out
    / "92G4_sample_mode_summary.tsv"
)

with summary_file.open(
    "w",
    newline=""
) as fh:

    fields = list(
        sample_mode_summary[
            0
        ].keys()
    )

    w = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(
        sample_mode_summary
    )


# ======================================================================
# Matrices
# ======================================================================

def matrix_filename(
    mode,
    metric,
):

    mode_slug = (
        mode.lower()
    )

    return (
        matrix_dir
        / f"92G4_{mode_slug}_{metric}.tsv.gz"
    )


def write_matrix(
    mode,
    metric,
):

    path = matrix_filename(
        mode,
        metric,
    )

    with gzip.open(
        path,
        "wt",
        newline=""
    ) as fh:

        w = csv.writer(
            fh,
            delimiter="\t"
        )

        w.writerow(
            [
                "votu_id",
                *samples,
            ]
        )


        for votu in sorted(
            expected_ids
        ):

            w.writerow(
                [
                    votu,
                    *[
                        data[
                            mode
                        ][sample][votu][
                            metric
                        ]
                        for sample in samples
                    ],
                ]
            )


matrix_metrics = [
    "mapped_records",
    "mapped_records_per_million_hostfree_reads",
    "mean_depth",
    "mean_depth_per_million_hostfree_reads",
    "breadth_1x",
    "breadth_5x",
]


for mode in modes:

    for metric in matrix_metrics:

        write_matrix(
            mode,
            metric,
        )


# ======================================================================
# Breadth threshold indicator matrices
#
# Descriptive sensitivity matrices.
# NOT biological presence calls.
# ======================================================================

for threshold in [
    0.10,
    0.50,
    0.80,
    0.95,
]:

    label = str(
        int(
            threshold * 100
        )
    )

    path = (
        matrix_dir
        / (
            "92G4_mapq20_id95_"
            f"breadth1_ge_{label}_indicator.tsv.gz"
        )
    )

    with gzip.open(
        path,
        "wt",
        newline=""
    ) as fh:

        w = csv.writer(
            fh,
            delimiter="\t"
        )

        w.writerow(
            [
                "votu_id",
                *samples,
            ]
        )


        for votu in sorted(
            expected_ids
        ):

            w.writerow(
                [
                    votu,
                    *[
                        int(
                            data[
                                "MAPQ20_ID95"
                            ][sample][votu][
                                "breadth_1x"
                            ]
                            >= threshold
                        )
                        for sample in samples
                    ],
                ]
            )


# ======================================================================
# Recurrence across libraries
#
# IMPORTANT:
# counts libraries, not independent biological prevalence.
# ======================================================================

recurrence_file = (
    out
    / "92G4_mapq20_id95_library_recurrence.tsv"
)


recurrence_fields = [
    "votu_id",

    "n_libraries_mapped_gt0",

    "n_libraries_breadth1_ge_0.10",
    "n_libraries_breadth1_ge_0.50",
    "n_libraries_breadth1_ge_0.80",
    "n_libraries_breadth1_ge_0.95",

    "n_L_libraries_breadth1_ge_0.80",
    "n_M_libraries_breadth1_ge_0.80",

    "n_L_libraries_breadth1_ge_0.95",
    "n_M_libraries_breadth1_ge_0.95",

    "n_biological_units_breadth1_ge_0.80",
    "n_biological_units_breadth1_ge_0.95",

    "n_producers_breadth1_ge_0.80",
    "n_producers_breadth1_ge_0.95",
]


with recurrence_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=recurrence_fields,
        delimiter="\t"
    )

    w.writeheader()


    for votu in sorted(
        expected_ids
    ):

        values = {
            sample:
                data[
                    "MAPQ20_ID95"
                ][sample][votu]
            for sample in samples
        }


        def passing(
            threshold
        ):

            return [
                sample
                for sample in samples
                if values[
                    sample
                ][
                    "breadth_1x"
                ] >= threshold
            ]


        p80 = passing(
            0.80
        )

        p95 = passing(
            0.95
        )


        bio80 = {
            (
                sample_meta[
                    s
                ][
                    "producer"
                ],
                sample_meta[
                    s
                ][
                    "biological_unit"
                ],
            )
            for s in p80
        }

        bio95 = {
            (
                sample_meta[
                    s
                ][
                    "producer"
                ],
                sample_meta[
                    s
                ][
                    "biological_unit"
                ],
            )
            for s in p95
        }


        prod80 = {
            sample_meta[
                s
            ][
                "producer"
            ]
            for s in p80
        }

        prod95 = {
            sample_meta[
                s
            ][
                "producer"
            ]
            for s in p95
        }


        w.writerow({

            "votu_id":
                votu,

            "n_libraries_mapped_gt0":
                sum(
                    values[
                        s
                    ][
                        "mapped_records"
                    ] > 0
                    for s in samples
                ),

            "n_libraries_breadth1_ge_0.10":
                len(
                    passing(
                        0.10
                    )
                ),

            "n_libraries_breadth1_ge_0.50":
                len(
                    passing(
                        0.50
                    )
                ),

            "n_libraries_breadth1_ge_0.80":
                len(
                    p80
                ),

            "n_libraries_breadth1_ge_0.95":
                len(
                    p95
                ),

            "n_L_libraries_breadth1_ge_0.80":
                sum(
                    s.startswith(
                        "L"
                    )
                    for s in p80
                ),

            "n_M_libraries_breadth1_ge_0.80":
                sum(
                    s.startswith(
                        "M"
                    )
                    for s in p80
                ),

            "n_L_libraries_breadth1_ge_0.95":
                sum(
                    s.startswith(
                        "L"
                    )
                    for s in p95
                ),

            "n_M_libraries_breadth1_ge_0.95":
                sum(
                    s.startswith(
                        "M"
                    )
                    for s in p95
                ),

            "n_biological_units_breadth1_ge_0.80":
                len(
                    bio80
                ),

            "n_biological_units_breadth1_ge_0.95":
                len(
                    bio95
                ),

            "n_producers_breadth1_ge_0.80":
                len(
                    prod80
                ),

            "n_producers_breadth1_ge_0.95":
                len(
                    prod95
                ),
        })


# ======================================================================
# Global vOTU summary
# ======================================================================

best_quality_counts = Counter(
    votu_meta[
        v
    ][
        "best_checkv_quality"
    ]
    for v in expected_ids
)


id95_retentions = []

for sample in samples:

    q20 = (
        data[
            "MAPQ20"
        ][sample][
            "VOTU00001"
        ][
            "mode_total_mapped_records"
        ]
    )

    i95 = (
        data[
            "MAPQ20_ID95"
        ][sample][
            "VOTU00001"
        ][
            "mode_total_mapped_records"
        ]
    )

    id95_retentions.append(
        i95 / q20
    )


global_summary = [
    (
        "votu_total",
        N_VOTU,
    ),

    (
        "libraries",
        len(samples),
    ),

    (
        "quantification_modes",
        len(modes),
    ),

    (
        "long_table_rows",
        n_long,
    ),

    (
        "nested_filter_tests",
        nested_tests,
    ),

    (
        "MAPQ20_ID95_min_retention_vs_MAPQ20",
        min(
            id95_retentions
        ),
    ),

    (
        "MAPQ20_ID95_max_retention_vs_MAPQ20",
        max(
            id95_retentions
        ),
    ),

    (
        "best_quality_complete",
        best_quality_counts.get(
            "Complete",
            0,
        ),
    ),

    (
        "best_quality_high_quality",
        best_quality_counts.get(
            "High-quality",
            0,
        ),
    ),

    (
        "best_quality_medium_quality",
        best_quality_counts.get(
            "Medium-quality",
            0,
        ),
    ),

    (
        "best_quality_low_quality",
        best_quality_counts.get(
            "Low-quality",
            0,
        ),
    ),

    (
        "best_quality_not_determined",
        best_quality_counts.get(
            "Not-determined",
            0,
        ),
    ),
]


global_file = (
    out
    / "92G4_global_summary.tsv"
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
# Metadata
# ======================================================================

metadata_file = (
    out
    / "92G4_run_metadata.tsv"
)

with metadata_file.open(
    "w"
) as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "analysis\t"
        "integrated_vOTU_quantification\n"
    )

    fh.write(
        "votu_reference\t"
        "4491_species_rank_operational_vOTUs\n"
    )

    fh.write(
        "modes\t"
        "ALL_PRIMARY;MAPQ20;MAPQ20_ID95\n"
    )

    fh.write(
        "MAPQ20_threshold\t20\n"
    )

    fh.write(
        "ID95_definition\t"
        "read_edit_identity_ge_0.95\n"
    )

    fh.write(
        "abundance_normalization\t"
        "per_million_hostfree_read_records\n"
    )

    fh.write(
        "breadth_threshold_indicators_are_presence_calls\tNO\n"
    )

    fh.write(
        "recurrence_unit\t"
        "sequencing_library_not_independent_prevalence\n"
    )

    fh.write(
        "taxonomy_field_scope\t"
        "geNomad_terminal_taxa_aggregated_across_vOTU_members\n"
    )

    fh.write(
        "quality_field_for_flags\t"
        "best_CheckV_member_within_vOTU\n"
    )

    fh.write(
        "week_mapping\t"
        "_1=week0;_2=week1;_3=week2\n"
    )


# ======================================================================
# Final validation
# ======================================================================

required = [
    sample_meta_file,
    votu_meta_file,
    long_file,
    summary_file,
    recurrence_file,
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


matrix_files = list(
    matrix_dir.glob(
        "*.tsv.gz"
    )
)

# 3 modes × 6 metrics = 18
# + 4 indicator matrices = 22
if len(matrix_files) != 22:

    raise RuntimeError(
        f"Matrix files={len(matrix_files)}; "
        "esperaba 22"
    )


print()
print("=" * 80)
print("92G4 COMPLETED")
print("=" * 80)

for key, value in global_summary:

    print(
        key,
        value,
        sep="\t"
    )

print()
print(
    "matrix_files",
    len(matrix_files),
    sep="\t"
)

print()
print("=" * 80)
print("SAMPLE × MODE SUMMARY")
print("=" * 80)

print(
    "sample\tmode\tmapped\t"
    "ret_all\tret_q20\t"
    "gt0\tb1_50\tb1_80\tb1_95\t"
    "HQ_b1_80"
)

for x in sample_mode_summary:

    print(
        x[
            "sample"
        ],
        x[
            "mode"
        ],
        x[
            "mapped_records"
        ],
        f"{x['retention_vs_ALL_PRIMARY']:.6f}",
        f"{x['retention_vs_MAPQ20']:.6f}",
        x[
            "votu_mapped_gt0"
        ],
        x[
            "breadth1_ge_0.50"
        ],
        x[
            "breadth1_ge_0.80"
        ],
        x[
            "breadth1_ge_0.95"
        ],
        x[
            "HQ_or_complete_breadth1_ge_0.80"
        ],
        sep="\t"
    )

PY


# ======================================================================
# Checksums
# ======================================================================

find "$OUT" \
    -maxdepth 2 \
    -type f \
    ! -name '92G4_sha256.tsv' \
    ! -name '92G4_COMPLETE.ok' \
    -print0 \
| sort -z \
| xargs -0 sha256sum \
> "${OUT}/92G4_sha256.tsv"


printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"$SLURM_JOB_ID" \
> "$DONE"


# ======================================================================
# Shell validation
# ======================================================================

echo
echo "======================================================================"
echo "GLOBAL SUMMARY"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92G4_global_summary.tsv"


echo
echo "======================================================================"
echo "SAMPLE MODE SUMMARY"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92G4_sample_mode_summary.tsv"


echo
echo "======================================================================"
echo "FILES"
echo "======================================================================"

echo "Matrices:"
find "$MATRIX" \
    -maxdepth 1 \
    -type f \
    -name '*.tsv.gz' \
    | sort \
    | wc -l

echo
echo "Integrated outputs:"
find "$OUT" \
    -maxdepth 1 \
    -type f \
    -printf '%f\n' \
    | sort


echo
echo "======================================================================"
echo "92G4 FINALIZO CORRECTAMENTE."
echo "Viroma cuantitativo integrado."
echo "NO se definió presencia/ausencia biológica."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "======================================================================"

