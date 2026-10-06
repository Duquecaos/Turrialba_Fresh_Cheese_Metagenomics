#!/usr/bin/env python3

import csv
import statistics
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(sys.argv[1])

D0 = (
    ROOT
    / "97_16S_shotgun"
    / "95D0_mapping_audit"
)

D1 = (
    ROOT
    / "97_16S_shotgun"
    / "95D1_locus_quantification"
)

TARGETS_FILE = (
    D0
    / "95D0_target_loci_manifest.tsv"
)

BAM_MANIFEST = (
    D0
    / "95D0_bam_manifest.tsv"
)

MAPPING_FILE = (
    D0
    / "95D0_mapping_strategy.tsv"
)


def read_tsv(path):

    if not path.is_file():
        raise RuntimeError(
            f"Archivo faltante: {path}"
        )

    with path.open(
        newline="",
        errors="replace"
    ) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        return list(reader)


def write_tsv(path, rows, fields):

    with path.open(
        "w",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)


def pct(n, d):

    if d == 0:
        return 0.0

    return 100.0 * n / d


def ratio(a, b):

    if b == 0:
        return 0.0

    return a / b


def depth_summary(values):

    n = len(values)

    if n == 0:
        raise RuntimeError(
            "Vector de profundidad vacío"
        )

    total = sum(values)

    return {
        "sum_depth":
            total,

        "mean_depth":
            total / n,

        "median_depth":
            statistics.median(
                values
            ),

        "max_depth":
            max(values),

        "breadth_1x_pct":
            pct(
                sum(x >= 1 for x in values),
                n
            ),

        "breadth_5x_pct":
            pct(
                sum(x >= 5 for x in values),
                n
            ),

        "breadth_10x_pct":
            pct(
                sum(x >= 10 for x in values),
                n
            ),
    }


# ============================================================
# Load inputs
# ============================================================

targets = read_tsv(
    TARGETS_FILE
)

bam_manifest = read_tsv(
    BAM_MANIFEST
)

mapping_rows = read_tsv(
    MAPPING_FILE
)


if len(targets) != 121:
    raise RuntimeError(
        f"Targets={len(targets)} expected=121"
    )

if len(bam_manifest) != 18:
    raise RuntimeError(
        f"BAM manifest={len(bam_manifest)} expected=18"
    )

if len(mapping_rows) != 18:
    raise RuntimeError(
        f"Mapping rows={len(mapping_rows)} expected=18"
    )


targets_by_group = defaultdict(list)

for r in targets:
    targets_by_group[
        r["coassembly"]
    ].append(r)


mapping_by_sample = {
    r["sample"]: r
    for r in mapping_rows
}


# ============================================================
# Parse one depth file into locus-specific vectors
# ============================================================

def parse_depth(
    path,
    group_targets
):

    if not path.is_file():
        raise RuntimeError(
            f"Depth faltante: {path}"
        )

    # Empty depth is technically possible but unexpected here.
    if path.stat().st_size == 0:
        raise RuntimeError(
            f"Depth vacío: {path}"
        )


    vectors = {}

    intervals = defaultdict(list)


    for r in group_targets:

        locus = r["locus_id"]

        start = int(
            r["start"]
        )

        end = int(
            r["end"]
        )

        length = int(
            r["length_bp"]
        )


        if end - start + 1 != length:
            raise RuntimeError(
                f"{locus}: coordenadas incompatibles "
                f"con length_bp"
            )


        vectors[locus] = [
            0
        ] * length


        intervals[
            r["contig"]
        ].append(
            (
                start,
                end,
                locus,
            )
        )


    for contig in intervals:
        intervals[
            contig
        ].sort()


    depth_rows = 0
    assigned_rows = 0
    unassigned_rows = 0
    multiassigned_positions = 0


    with path.open() as fh:

        for line in fh:

            if not line.strip():
                continue

            parts = (
                line.rstrip("\n")
                .split("\t")
            )

            if len(parts) < 3:
                continue


            contig = parts[0]
            pos = int(
                parts[1]
            )

            depth = int(
                parts[2]
            )

            depth_rows += 1


            matches = []

            for start, end, locus in (
                intervals.get(
                    contig,
                    []
                )
            ):

                if start <= pos <= end:
                    matches.append(
                        (
                            start,
                            locus,
                        )
                    )


            if not matches:

                unassigned_rows += 1
                continue


            if len(matches) > 1:
                multiassigned_positions += 1


            for start, locus in matches:

                idx = pos - start

                vectors[
                    locus
                ][idx] = depth

                assigned_rows += 1


    return (
        vectors,
        {
            "depth_rows":
                depth_rows,

            "assigned_rows":
                assigned_rows,

            "unassigned_rows":
                unassigned_rows,

            "multiassigned_positions":
                multiassigned_positions,
        }
    )


# ============================================================
# Recover all 18 samples
# ============================================================

global_metrics = []
sample_summaries = []


for bamrow in bam_manifest:

    sample = bamrow[
        "sample"
    ]

    group = bamrow[
        "coassembly"
    ]

    if sample not in mapping_by_sample:
        raise RuntimeError(
            f"Mapping metadata missing for {sample}"
        )


    mapping = mapping_by_sample[
        sample
    ]

    if mapping[
        "coassembly"
    ] != group:
        raise RuntimeError(
            f"{sample}: group mismatch"
        )


    mapped_alignments = int(
        mapping[
            "idxstats_mapped_alignments"
        ]
    )

    if mapped_alignments <= 0:
        raise RuntimeError(
            f"{sample}: mapped alignments <=0"
        )


    group_targets = (
        targets_by_group[
            group
        ]
    )


    outdir = (
        D1
        / "results"
        / sample
    )

    if not outdir.is_dir():
        raise RuntimeError(
            f"Directorio faltante: {outdir}"
        )


    depth_all = (
        outdir
        / f"{sample}.depth.all.tsv"
    )

    depth_q10 = (
        outdir
        / f"{sample}.depth.mapq10.tsv"
    )


    vectors_all, audit_all = (
        parse_depth(
            depth_all,
            group_targets
        )
    )

    vectors_q10, audit_q10 = (
        parse_depth(
            depth_q10,
            group_targets
        )
    )


    # samtools depth -b BED should only emit positions
    # belonging to the BED targets.
    if (
        audit_all[
            "unassigned_rows"
        ] != 0
    ):
        raise RuntimeError(
            f"{sample}: ALL contains "
            f"{audit_all['unassigned_rows']} "
            f"unassigned positions"
        )


    if (
        audit_q10[
            "unassigned_rows"
        ] != 0
    ):
        raise RuntimeError(
            f"{sample}: Q10 contains "
            f"{audit_q10['unassigned_rows']} "
            f"unassigned positions"
        )


    million_mapped = (
        mapped_alignments
        / 1_000_000.0
    )


    sample_metrics = []


    for r in sorted(
        group_targets,
        key=lambda x:
            x["locus_id"]
    ):

        locus = r[
            "locus_id"
        ]

        a = depth_summary(
            vectors_all[
                locus
            ]
        )

        q = depth_summary(
            vectors_q10[
                locus
            ]
        )


        rec = {
            "sample":
                sample,

            "coassembly":
                group,

            "locus_id":
                locus,

            "exact_cluster_id":
                r[
                    "taxonomy_exact_cluster_id"
                ],

            "contig":
                r["contig"],

            "start":
                r["start"],

            "end":
                r["end"],

            "strand":
                r["strand"],

            "length_bp":
                r["length_bp"],

            "deepest_informative_LCA_taxon":
                r[
                    "deepest_informative_LCA_taxon"
                ],

            "SILVA_LCA_terminal":
                r[
                    "SILVA_LCA_terminal"
                ],

            "resolution_tier":
                r[
                    "revised_resolution_tier"
                ],

            "MAG_associated":
                r[
                    "exact_final18_contig_membership"
                ],

            "final18_MAGs":
                r[
                    "final18_MAGs_exact"
                ],

            "mapped_alignments_library":
                mapped_alignments,

            "all_sum_depth":
                a[
                    "sum_depth"
                ],

            "all_mean_depth":
                f"{a['mean_depth']:.8f}",

            "all_median_depth":
                f"{a['median_depth']:.8f}",

            "all_max_depth":
                a[
                    "max_depth"
                ],

            "all_breadth_1x_pct":
                f"{a['breadth_1x_pct']:.6f}",

            "all_breadth_5x_pct":
                f"{a['breadth_5x_pct']:.6f}",

            "all_breadth_10x_pct":
                f"{a['breadth_10x_pct']:.6f}",

            "all_mean_depth_per_million_mapped":
                f"{a['mean_depth'] / million_mapped:.10f}",

            "all_sum_depth_per_million_mapped":
                f"{a['sum_depth'] / million_mapped:.10f}",

            "q10_sum_depth":
                q[
                    "sum_depth"
                ],

            "q10_mean_depth":
                f"{q['mean_depth']:.8f}",

            "q10_median_depth":
                f"{q['median_depth']:.8f}",

            "q10_max_depth":
                q[
                    "max_depth"
                ],

            "q10_breadth_1x_pct":
                f"{q['breadth_1x_pct']:.6f}",

            "q10_breadth_5x_pct":
                f"{q['breadth_5x_pct']:.6f}",

            "q10_breadth_10x_pct":
                f"{q['breadth_10x_pct']:.6f}",

            "q10_mean_depth_per_million_mapped":
                f"{q['mean_depth'] / million_mapped:.10f}",

            "q10_sum_depth_per_million_mapped":
                f"{q['sum_depth'] / million_mapped:.10f}",

            "q10_to_all_mean_depth_ratio":
                f"{ratio(q['mean_depth'], a['mean_depth']):.8f}",

            "q10_to_all_sum_depth_ratio":
                f"{ratio(q['sum_depth'], a['sum_depth']):.8f}",

            "interpretation":
                "mapping_derived_16S_locus_signal_not_cell_abundance",
        }


        sample_metrics.append(
            rec
        )

        global_metrics.append(
            rec
        )


    fields = list(
        sample_metrics[0].keys()
    )


    write_tsv(
        outdir
        / f"{sample}.95D1_16S_locus_metrics.tsv",
        sample_metrics,
        fields
    )


    summary = {
        "sample":
            sample,

        "coassembly":
            group,

        "target_loci":
            len(
                sample_metrics
            ),

        "mapped_alignments_library":
            mapped_alignments,

        "loci_with_any_all_coverage":
            sum(
                float(
                    r[
                        "all_breadth_1x_pct"
                    ]
                ) > 0
                for r in sample_metrics
            ),

        "loci_all_breadth_ge50pct":
            sum(
                float(
                    r[
                        "all_breadth_1x_pct"
                    ]
                ) >= 50
                for r in sample_metrics
            ),

        "loci_all_breadth_ge90pct":
            sum(
                float(
                    r[
                        "all_breadth_1x_pct"
                    ]
                ) >= 90
                for r in sample_metrics
            ),

        "loci_q10_breadth_ge50pct":
            sum(
                float(
                    r[
                        "q10_breadth_1x_pct"
                    ]
                ) >= 50
                for r in sample_metrics
            ),

        "loci_q10_breadth_ge90pct":
            sum(
                float(
                    r[
                        "q10_breadth_1x_pct"
                    ]
                ) >= 90
                for r in sample_metrics
            ),

        "depth_rows_all":
            audit_all[
                "depth_rows"
            ],

        "depth_rows_q10":
            audit_q10[
                "depth_rows"
            ],

        "multiassigned_depth_positions_all":
            audit_all[
                "multiassigned_positions"
            ],

        "multiassigned_depth_positions_q10":
            audit_q10[
                "multiassigned_positions"
            ],

        "interpretation":
            "breadth_counts_descriptive_not_presence_calls",
    }


    write_tsv(
        outdir
        / f"{sample}.95D1_summary.tsv",
        [summary],
        list(
            summary.keys()
        )
    )


    sample_summaries.append(
        summary
    )


    with (
        outdir
        / "95D1R_RECOVERED.ok"
    ).open("w") as fh:

        fh.write(
            f"RECOVERED\t{sample}\t"
            f"from_existing_depth_files\n"
        )


# ============================================================
# Consolidated outputs
# ============================================================

write_tsv(
    D1
    / "95D1_all_samples_locus_metrics.tsv",
    global_metrics,
    list(
        global_metrics[0].keys()
    )
)


write_tsv(
    D1
    / "95D1_all_samples_summary.tsv",
    sample_summaries,
    list(
        sample_summaries[0].keys()
    )
)


expected_rows = sum(
    len(
        targets_by_group[
            r["coassembly"]
        ]
    )
    for r in bam_manifest
)


if len(global_metrics) != expected_rows:
    raise RuntimeError(
        f"Global metric rows={len(global_metrics)} "
        f"expected={expected_rows}"
    )


with (
    D1
    / "95D1R_global_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    fh.write(
        f"samples_recovered\t{len(sample_summaries)}\n"
    )

    fh.write(
        f"locus_sample_rows\t{len(global_metrics)}\n"
    )

    fh.write(
        f"expected_locus_sample_rows\t{expected_rows}\n"
    )

    fh.write(
        "samtools_depth_rerun\tNO\n"
    )

    fh.write(
        "reused_existing_depth_files\tYES\n"
    )

    fh.write(
        "original_array_job\t132799\n"
    )

    fh.write(
        "original_array_status\t"
        "FAILED_after_depth_generation_before_metrics_validation\n"
    )


print(
    f"SAMPLES={len(sample_summaries)}"
)

print(
    f"LOCUS_SAMPLE_ROWS={len(global_metrics)}"
)

print(
    f"EXPECTED_ROWS={expected_rows}"
)

print(
    "SAMTOOLS_DEPTH_RERUN=NO"
)

print(
    "95D1R_RECOVERY=PASS"
)

