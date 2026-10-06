#!/bin/bash

#SBATCH --job-name=92g2b_rebuild
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=92g2b_rebuild_%j.out
#SBATCH --error=92g2b_rebuild_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

MAPROOT="${ROOT}/92_votu_abundance/mapping"
OUTROOT="${ROOT}/93_votu_coverage"

mkdir -p "$OUTROOT"

echo "======================================================================"
echo "92G2b - RECONSTRUCCION ROBUSTA DE METRICAS"
echo "NO repite Bowtie2"
echo "NO repite samtools depth"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"

python3 - <<'PY'

import csv
import os
import re
from pathlib import Path


root = Path(
    f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"
)

maproot = (
    root
    / "92_votu_abundance"
    / "mapping"
)

outroot = (
    root
    / "93_votu_coverage"
)

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
]

N_VOTU = 4491
TOTAL_BP = 26784056


# ============================================================
# Helpers
# ============================================================

def mapping_metric(path, key):

    with path.open() as fh:

        r = csv.reader(
            fh,
            delimiter="\t"
        )

        next(r)

        for row in r:

            if (
                len(row) >= 2
                and row[0] == key
            ):
                return row[1]

    raise RuntimeError(
        f"No encontré {key} en {path}"
    )


def input_pairs_from_bt2(path):

    with path.open() as fh:

        for line in fh:

            m = re.match(
                r"^\s*([0-9]+)\s+reads;",
                line
            )

            if m:
                return int(
                    m.group(1)
                )

    raise RuntimeError(
        f"No pude recuperar reads de {path}"
    )


def read_counts(path):

    result = {}

    with path.open() as fh:

        r = csv.DictReader(
            fh,
            delimiter="\t"
        )

        if set(
            r.fieldnames or []
        ) != {
            "votu_id",
            "mapped_records",
        }:

            raise RuntimeError(
                f"Schema inesperado: {path} "
                f"{r.fieldnames}"
            )

        for row in r:

            votu = row[
                "votu_id"
            ]

            value = row[
                "mapped_records"
            ].strip()

            if not value:

                raise RuntimeError(
                    f"{path}: mapped_records vacío "
                    f"para {votu}"
                )

            if votu in result:

                raise RuntimeError(
                    f"{path}: ID duplicado {votu}"
                )

            result[votu] = int(
                value
            )

    return result


def parse_int_required(
    value,
    field,
    path,
    votu,
):

    value = (
        ""
        if value is None
        else value.strip()
    )

    if value == "":

        raise RuntimeError(
            f"{path}: {field} vacío "
            f"para {votu}"
        )

    return int(
        float(value)
    )


def parse_depth(path, sample, mode):

    rows = []

    patches = {
        "covered_bp_1x": 0,
        "covered_bp_5x": 0,
        "max_depth": 0,
    }

    max_b1_rounding_diff = 0.0
    max_b5_rounding_diff = 0.0
    max_mean_rounding_diff = 0.0


    with path.open() as fh:

        r = csv.DictReader(
            fh,
            delimiter="\t"
        )

        required = {
            "sample",
            "mode",
            "votu_id",
            "length_bp",
            "depth_sum",
            "mean_depth",
            "covered_bp_1x",
            "breadth_1x",
            "breadth_1x_percent",
            "covered_bp_5x",
            "breadth_5x",
            "breadth_5x_percent",
            "max_depth",
        }

        missing = (
            required
            - set(
                r.fieldnames or []
            )
        )

        if missing:

            raise RuntimeError(
                f"{path}: faltan columnas "
                f"{sorted(missing)}"
            )


        for index, row in enumerate(r):

            votu = row[
                "votu_id"
            ]

            if row["sample"] != sample:

                raise RuntimeError(
                    f"{path}: sample "
                    f"{row['sample']} != {sample}"
                )

            if row["mode"] != mode:

                raise RuntimeError(
                    f"{path}: mode "
                    f"{row['mode']} != {mode}"
                )


            length = parse_int_required(
                row["length_bp"],
                "length_bp",
                path,
                votu,
            )

            depth_sum = parse_int_required(
                row["depth_sum"],
                "depth_sum",
                path,
                votu,
            )


            # ------------------------------------------------
            # BUG conocido del AWK:
            # un contador cero no inicializado puede aparecer
            # como cadena vacía SOLO en la primera referencia.
            # ------------------------------------------------

            primitive = {}

            for field in [
                "covered_bp_1x",
                "covered_bp_5x",
                "max_depth",
            ]:

                value = (
                    ""
                    if row[field] is None
                    else row[field].strip()
                )

                if value == "":

                    if index != 0:

                        raise RuntimeError(
                            f"{path}: {field} vacío "
                            f"fuera de primera fila: "
                            f"{votu}"
                        )

                    primitive[
                        field
                    ] = 0

                    patches[
                        field
                    ] += 1

                else:

                    primitive[
                        field
                    ] = int(
                        float(value)
                    )


            cov1 = primitive[
                "covered_bp_1x"
            ]

            cov5 = primitive[
                "covered_bp_5x"
            ]

            maxdepth = primitive[
                "max_depth"
            ]


            # ------------------------------------------------
            # Validaciones de las métricas enteras.
            # ------------------------------------------------

            if length <= 0:

                raise RuntimeError(
                    f"{path}: longitud <=0 "
                    f"para {votu}"
                )

            if not (
                0
                <= cov5
                <= cov1
                <= length
            ):

                raise RuntimeError(
                    f"{sample}/{mode}/{votu}: "
                    f"cov5={cov5}, cov1={cov1}, "
                    f"length={length}"
                )

            if depth_sum < cov1:

                raise RuntimeError(
                    f"{sample}/{mode}/{votu}: "
                    f"depth_sum={depth_sum} "
                    f"< covered_bp_1x={cov1}"
                )

            if cov1 == 0:

                if depth_sum != 0:

                    raise RuntimeError(
                        f"{sample}/{mode}/{votu}: "
                        "cov1=0 pero depth_sum>0"
                    )

                if maxdepth != 0:

                    raise RuntimeError(
                        f"{sample}/{mode}/{votu}: "
                        "cov1=0 pero maxdepth>0"
                    )

            else:

                if maxdepth < 1:

                    raise RuntimeError(
                        f"{sample}/{mode}/{votu}: "
                        "cov1>0 pero maxdepth<1"
                    )


            # ------------------------------------------------
            # RECONSTRUIR LOS DECIMALES DESDE ENTEROS.
            #
            # No usamos los breadth/mean originales para
            # el resultado final porque fueron impresos por
            # AWK con precisión limitada.
            # ------------------------------------------------

            mean_depth = (
                depth_sum
                / length
            )

            breadth1 = (
                cov1
                / length
            )

            breadth5 = (
                cov5
                / length
            )


            # ------------------------------------------------
            # Auditoría de cuánto redondeó AWK.
            # No se usa como criterio de rechazo.
            # ------------------------------------------------

            original_mean = (
                row[
                    "mean_depth"
                ].strip()
            )

            original_b1 = (
                row[
                    "breadth_1x"
                ].strip()
            )

            original_b5 = (
                row[
                    "breadth_5x"
                ].strip()
            )


            if original_mean:

                diff = abs(
                    float(
                        original_mean
                    )
                    - mean_depth
                )

                max_mean_rounding_diff = max(
                    max_mean_rounding_diff,
                    diff,
                )


            if original_b1:

                diff = abs(
                    float(
                        original_b1
                    )
                    - breadth1
                )

                max_b1_rounding_diff = max(
                    max_b1_rounding_diff,
                    diff,
                )


            if original_b5:

                diff = abs(
                    float(
                        original_b5
                    )
                    - breadth5
                )

                max_b5_rounding_diff = max(
                    max_b5_rounding_diff,
                    diff,
                )


            rows.append({
                "sample":
                    sample,

                "mode":
                    mode,

                "votu_id":
                    votu,

                "length_bp":
                    length,

                "depth_sum":
                    depth_sum,

                "mean_depth":
                    mean_depth,

                "covered_bp_1x":
                    cov1,

                "breadth_1x":
                    breadth1,

                "breadth_1x_percent":
                    breadth1 * 100,

                "covered_bp_5x":
                    cov5,

                "breadth_5x":
                    breadth5,

                "breadth_5x_percent":
                    breadth5 * 100,

                "max_depth":
                    maxdepth,
            })


    # --------------------------------------------------------
    # QC global de la referencia.
    # --------------------------------------------------------

    if len(rows) != N_VOTU:

        raise RuntimeError(
            f"{path}: filas={len(rows)}, "
            f"esperaba {N_VOTU}"
        )


    ids = [
        x["votu_id"]
        for x in rows
    ]

    if len(ids) != len(
        set(ids)
    ):

        raise RuntimeError(
            f"{path}: vOTUs duplicados"
        )


    if sum(
        x["length_bp"]
        for x in rows
    ) != TOTAL_BP:

        raise RuntimeError(
            f"{path}: longitud de referencia "
            "no coincide"
        )


    audit = {
        "patched_blank_covered_bp_1x":
            patches[
                "covered_bp_1x"
            ],

        "patched_blank_covered_bp_5x":
            patches[
                "covered_bp_5x"
            ],

        "patched_blank_max_depth":
            patches[
                "max_depth"
            ],

        "max_abs_rounding_diff_mean_depth":
            max_mean_rounding_diff,

        "max_abs_rounding_diff_breadth1":
            max_b1_rounding_diff,

        "max_abs_rounding_diff_breadth5":
            max_b5_rounding_diff,
    }


    return rows, audit


# ============================================================
# Main
# ============================================================

repair_audit = []
global_summary = []


for sample in samples:

    print()
    print(
        "=" * 72
    )
    print(
        "SAMPLE",
        sample
    )
    print(
        "=" * 72
    )

    sout = (
        outroot
        / sample
    )

    sout.mkdir(
        parents=True,
        exist_ok=True
    )


    bt2log = (
        maproot
        / sample
        / f"{sample}.bowtie2.log"
    )

    oldmetrics = (
        maproot
        / sample
        / f"{sample}.mapping_metrics.tsv"
    )


    if not bt2log.is_file():

        raise RuntimeError(
            f"Falta {bt2log}"
        )

    if not oldmetrics.is_file():

        raise RuntimeError(
            f"Falta {oldmetrics}"
        )


    input_pairs = (
        input_pairs_from_bt2(
            bt2log
        )
    )

    input_read_records = (
        input_pairs
        * 2
    )


    expected_mapped = {

        "ALL_PRIMARY":
            int(
                mapping_metric(
                    oldmetrics,
                    "primary_mapped_alignment_records",
                )
            ),

        "MAPQ20":
            int(
                mapping_metric(
                    oldmetrics,
                    "mapq20_primary_mapped_alignment_records",
                )
            ),
    }


    parsed = {}
    counts_by_mode = {}


    # ========================================================
    # Leer los dos modos.
    # ========================================================

    for mode in modes:

        counts_file = (
            sout
            / f"{sample}.{mode}.counts.tsv"
        )

        depth_file = (
            sout
            / f"{sample}.{mode}.depth.tsv"
        )


        for path in [
            counts_file,
            depth_file,
        ]:

            if not path.is_file():

                raise RuntimeError(
                    f"Falta intermedio {path}"
                )


        counts = read_counts(
            counts_file
        )

        rows, audit = parse_depth(
            depth_file,
            sample,
            mode,
        )


        mapped_sum = sum(
            counts.values()
        )


        if mapped_sum != expected_mapped[
            mode
        ]:

            raise RuntimeError(
                f"{sample}/{mode}: "
                f"counts sum={mapped_sum}, "
                f"92G1={expected_mapped[mode]}"
            )


        row_ids = {
            x["votu_id"]
            for x in rows
        }

        if (
            set(counts)
            - row_ids
        ):

            raise RuntimeError(
                f"{sample}/{mode}: "
                "counts contiene IDs "
                "fuera de referencia"
            )


        parsed[
            mode
        ] = {
            x["votu_id"]: x
            for x in rows
        }

        counts_by_mode[
            mode
        ] = counts


        repair_audit.append({

            "sample":
                sample,

            "mode":
                mode,

            "depth_rows":
                len(rows),

            "reference_total_bp":
                sum(
                    x[
                        "length_bp"
                    ]
                    for x in rows
                ),

            "mapped_records_counts":
                mapped_sum,

            "mapped_records_expected_92G1":
                expected_mapped[
                    mode
                ],

            **audit,

            "validation":
                "PASS",
        })


    # ========================================================
    # Sensibilidad MAPQ20:
    # debe ser subconjunto de ALL_PRIMARY.
    # ========================================================

    for votu in parsed[
        "ALL_PRIMARY"
    ]:

        a = parsed[
            "ALL_PRIMARY"
        ][votu]

        q = parsed[
            "MAPQ20"
        ][votu]


        if (
            q[
                "covered_bp_1x"
            ]
            > a[
                "covered_bp_1x"
            ]
        ):

            raise RuntimeError(
                f"{sample}/{votu}: "
                "MAPQ20 cov1 > ALL_PRIMARY"
            )


        if (
            q[
                "covered_bp_5x"
            ]
            > a[
                "covered_bp_5x"
            ]
        ):

            raise RuntimeError(
                f"{sample}/{votu}: "
                "MAPQ20 cov5 > ALL_PRIMARY"
            )


        if (
            q[
                "depth_sum"
            ]
            > a[
                "depth_sum"
            ]
        ):

            raise RuntimeError(
                f"{sample}/{votu}: "
                "MAPQ20 depth_sum > ALL_PRIMARY"
            )


        if (
            counts_by_mode[
                "MAPQ20"
            ].get(
                votu,
                0
            )
            >
            counts_by_mode[
                "ALL_PRIMARY"
            ].get(
                votu,
                0
            )
        ):

            raise RuntimeError(
                f"{sample}/{votu}: "
                "MAPQ20 records > ALL_PRIMARY"
            )


    # ========================================================
    # Construir salida final.
    # ========================================================

    final_rows = []
    summary_rows = []


    for mode in modes:

        mapped_total = (
            expected_mapped[
                mode
            ]
        )

        rows = [
            parsed[
                mode
            ][v]
            for v in sorted(
                parsed[
                    mode
                ]
            )
        ]


        for row in rows:

            votu = row[
                "votu_id"
            ]

            mapped = (
                counts_by_mode[
                    mode
                ].get(
                    votu,
                    0
                )
            )

            final_rows.append({

                "sample":
                    sample,

                "mode":
                    mode,

                "votu_id":
                    votu,

                "length_bp":
                    row[
                        "length_bp"
                    ],

                "mapped_records":
                    mapped,

                "mapped_records_per_million_input_read_records":
                    (
                        mapped
                        * 1_000_000
                        / input_read_records
                    ),

                "mapped_records_per_million_mode_mapped_records":
                    (
                        mapped
                        * 1_000_000
                        / mapped_total
                        if mapped_total
                        else 0
                    ),

                "covered_bp_1x":
                    row[
                        "covered_bp_1x"
                    ],

                "breadth_1x":
                    row[
                        "breadth_1x"
                    ],

                "breadth_1x_percent":
                    row[
                        "breadth_1x_percent"
                    ],

                "covered_bp_5x":
                    row[
                        "covered_bp_5x"
                    ],

                "breadth_5x":
                    row[
                        "breadth_5x"
                    ],

                "breadth_5x_percent":
                    row[
                        "breadth_5x_percent"
                    ],

                "depth_sum":
                    row[
                        "depth_sum"
                    ],

                "mean_depth":
                    row[
                        "mean_depth"
                    ],

                "max_depth":
                    row[
                        "max_depth"
                    ],

                "mean_depth_per_million_input_read_records":
                    (
                        row[
                            "mean_depth"
                        ]
                        * 1_000_000
                        / input_read_records
                    ),

                "mean_depth_per_million_mode_mapped_records":
                    (
                        row[
                            "mean_depth"
                        ]
                        * 1_000_000
                        / mapped_total
                        if mapped_total
                        else 0
                    ),

                "input_read_pairs":
                    input_pairs,

                "input_read_records":
                    input_read_records,

                "mode_total_mapped_records":
                    mapped_total,
            })


        metrics = {

            "mode_total_mapped_records":
                mapped_total,

            "votu_mapped_records_gt0":
                sum(
                    counts_by_mode[
                        mode
                    ].get(
                        x[
                            "votu_id"
                        ],
                        0
                    ) > 0
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
        }


        global_summary.append({
            "sample":
                sample,

            "mode":
                mode,

            **metrics,
        })


        for metric, value in (
            metrics.items()
        ):

            summary_rows.append({
                "sample":
                    sample,

                "mode":
                    mode,

                "metric":
                    metric,

                "value":
                    value,
            })


    if len(final_rows) != 8982:

        raise RuntimeError(
            f"{sample}: final_rows="
            f"{len(final_rows)}"
        )


    # ========================================================
    # Write sample outputs.
    # ========================================================

    final_file = (
        sout
        / f"{sample}.votu_metrics.tsv"
    )

    with final_file.open(
        "w",
        newline=""
    ) as fh:

        w = csv.DictWriter(
            fh,
            fieldnames=list(
                final_rows[
                    0
                ].keys()
            ),
            delimiter="\t"
        )

        w.writeheader()
        w.writerows(
            final_rows
        )


    summary_file = (
        sout
        / f"{sample}.92G2_summary.tsv"
    )

    with summary_file.open(
        "w",
        newline=""
    ) as fh:

        w = csv.DictWriter(
            fh,
            fieldnames=[
                "sample",
                "mode",
                "metric",
                "value",
            ],
            delimiter="\t"
        )

        w.writeheader()
        w.writerows(
            summary_rows
        )


    metadata_file = (
        sout
        / f"{sample}.92G2_metadata.tsv"
    )

    with metadata_file.open(
        "w"
    ) as fh:

        fh.write(
            "field\tvalue\n"
        )

        fh.write(
            "analysis\t"
            "vOTU_breadth_depth_rebuilt\n"
        )

        fh.write(
            f"sample\t{sample}\n"
        )

        fh.write(
            "source_depth_job\t131710\n"
        )

        fh.write(
            "mapping_recomputed\tNO\n"
        )

        fh.write(
            "samtools_depth_recomputed\tNO\n"
        )

        fh.write(
            "derived_float_metrics\t"
            "recalculated_from_integer_primitives\n"
        )

        fh.write(
            "mean_depth_formula\t"
            "depth_sum/length_bp\n"
        )

        fh.write(
            "breadth_1x_formula\t"
            "covered_bp_1x/length_bp\n"
        )

        fh.write(
            "breadth_5x_formula\t"
            "covered_bp_5x/length_bp\n"
        )

        fh.write(
            "modes\t"
            "ALL_PRIMARY;MAPQ20\n"
        )

        fh.write(
            "presence_threshold\t"
            "NOT_DEFINED\n"
        )


    marker = (
        sout
        / f"{sample}.92G2_COMPLETE.ok"
    )

    marker.write_text(
        "COMPLETED_REBUILT\t"
        "source_depth_job=131710\t"
        f"sample={sample}\n"
    )


    print(
        sample,
        "PASS",
        sep="\t"
    )


# ============================================================
# Global files
# ============================================================

audit_file = (
    outroot
    / "92G2_repair_audit.tsv"
)

with audit_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            repair_audit[
                0
            ].keys()
        ),
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(
        repair_audit
    )


summary_file = (
    outroot
    / "92G2_global_sample_summary.tsv"
)

with summary_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            global_summary[
                0
            ].keys()
        ),
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(
        global_summary
    )


if len(
    global_summary
) != 36:

    raise RuntimeError(
        f"Resumen global tiene "
        f"{len(global_summary)} filas; "
        "esperaba 36"
    )


# ============================================================
# Final summary
# ============================================================

print()
print("=" * 78)
print("92G2b GLOBAL")
print("=" * 78)

print(
    "sample\tmode\tmapped\tgt0\t"
    "b1_10\tb1_50\tb1_80\tb1_95\t"
    "b5_10\tb5_50\tb5_80"
)

for x in global_summary:

    print(
        x["sample"],
        x["mode"],
        x[
            "mode_total_mapped_records"
        ],
        x[
            "votu_mapped_records_gt0"
        ],
        x[
            "breadth1_ge_0.10"
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
            "breadth5_ge_0.10"
        ],
        x[
            "breadth5_ge_0.50"
        ],
        x[
            "breadth5_ge_0.80"
        ],
        sep="\t"
    )

PY


# ======================================================================
# VALIDACION FINAL
# ======================================================================

echo
echo "======================================================================"
echo "VALIDACION FINAL"
echo "======================================================================"

SAMPLES=(
    L1_1 L1_2 L1_3
    L2_1 L2_2 L2_3
    L3_1 L3_2 L3_3
    M1_1 M1_2 M1_3
    M2_1 M2_2 M2_3
    M3_1 M3_2 M3_3
)

for SAMPLE in "${SAMPLES[@]}"
do

    F="${OUTROOT}/${SAMPLE}/${SAMPLE}.votu_metrics.tsv"
    S="${OUTROOT}/${SAMPLE}/${SAMPLE}.92G2_summary.tsv"
    M="${OUTROOT}/${SAMPLE}/${SAMPLE}.92G2_COMPLETE.ok"

    [[ -s "$F" ]] || {
        echo "ERROR: falta $F" >&2
        exit 3
    }

    [[ -s "$S" ]] || {
        echo "ERROR: falta $S" >&2
        exit 4
    }

    [[ -s "$M" ]] || {
        echo "ERROR: falta $M" >&2
        exit 5
    }

    N="$(
        awk 'END{print NR-1}' "$F"
    )"

    [[ "$N" -eq 8982 ]] || {
        echo "ERROR: ${SAMPLE} rows=${N}" >&2
        exit 6
    }

    echo "OK  ${SAMPLE}  rows=${N}"

done


echo
echo "======================================================================"
echo "AUDITORIA"
echo "======================================================================"

column -t -s $'\t' \
    "${OUTROOT}/92G2_repair_audit.tsv"


echo
echo "======================================================================"
echo "RESUMEN GLOBAL"
echo "======================================================================"

column -t -s $'\t' \
    "${OUTROOT}/92G2_global_sample_summary.tsv"


echo
echo "======================================================================"
echo "92G2b FINALIZO CORRECTAMENTE."
echo "Los decimales se reconstruyeron desde métricas enteras."
echo "No se repitió Bowtie2."
echo "No se repitió samtools depth."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "======================================================================"

