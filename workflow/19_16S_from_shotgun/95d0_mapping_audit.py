#!/usr/bin/env python3

import csv
import re
import shlex
import sys
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(sys.argv[1])

CATALOG = (
    ROOT
    / "97_16S_shotgun"
    / "95C2C_final_curated_bacterial_catalog"
)

PRIMARY = (
    CATALOG
    / "95C2C_primary_158_loci_final_QC.tsv"
)

OUT = (
    ROOT
    / "97_16S_shotgun"
    / "95D0_mapping_audit"
)

HEADERS = (
    OUT
    / "headers"
)

IDXSTATS = (
    OUT
    / "idxstats"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)

HEADERS.mkdir(
    parents=True,
    exist_ok=True
)

IDXSTATS.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# Helpers
# ============================================================

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

        return (
            list(reader),
            reader.fieldnames or []
        )


def write_tsv(
    path,
    rows,
    fields
):

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


def yes(value):

    return (
        str(value)
        .strip()
        .lower()
        in {
            "yes",
            "true",
            "1",
            "y",
        }
    )


def classify_bowtie2_command(cl):

    if not cl:
        return (
            "bowtie2_command_not_recovered"
        )

    try:
        tokens = shlex.split(
            cl
        )
    except Exception:
        tokens = cl.split()

    # --all or -a
    if (
        "-a" in tokens
        or "--all" in tokens
    ):
        return (
            "all_alignments_explicit"
        )

    # -k N or -kN
    for i, token in enumerate(tokens):

        if token == "-k":
            return (
                "multiple_alignments_k_explicit"
            )

        if re.fullmatch(
            r"-k[0-9]+",
            token
        ):
            return (
                "multiple_alignments_k_explicit"
            )

    return (
        "no_explicit_a_or_k_detected"
    )


def extract_pg_fields(line):

    fields = {}

    for field in line.rstrip(
        "\n"
    ).split("\t")[1:]:

        if ":" not in field:
            continue

        key, value = field.split(
            ":",
            1
        )

        fields[key] = value

    return fields


# ============================================================
# 1. Final 121 primary loci
# ============================================================

primary, primary_fields = read_tsv(
    PRIMARY
)

if len(primary) != 158:
    raise RuntimeError(
        f"Primary rows={len(primary)} expected=158"
    )


targets = [
    r
    for r in primary
    if yes(
        r[
            "final_bacterial_quantification_target_v2"
        ]
    )
]


if len(targets) != 121:
    raise RuntimeError(
        f"Final bacterial loci={len(targets)} expected=121"
    )


target_fields = [
    "locus_id",
    "coassembly",
    "contig",
    "start",
    "end",
    "strand",
    "length_bp",
    "taxonomy_exact_cluster_id",
    "SILVA_LCA_lineage",
    "SILVA_LCA_terminal",
    "deepest_informative_LCA_taxon",
    "revised_resolution_tier",
    "exact_final18_contig_membership",
    "final18_MAGs_exact",
]


write_tsv(
    OUT
    / "95D0_target_loci_manifest.tsv",
    targets,
    target_fields
)


# ============================================================
# 2. BED files by coassembly
# ============================================================

targets_by_group = defaultdict(list)

for row in targets:

    targets_by_group[
        row["coassembly"]
    ].append(row)


BEDDIR = (
    OUT
    / "beds"
)

BEDDIR.mkdir(
    parents=True,
    exist_ok=True
)


for group, rows in sorted(
    targets_by_group.items()
):

    bed = (
        BEDDIR
        / f"{group}.16S_targets.bed"
    )

    with bed.open("w") as fh:

        for row in sorted(
            rows,
            key=lambda x: (
                x["contig"],
                int(x["start"]),
                int(x["end"]),
            )
        ):

            # BED = 0-based, half-open
            start0 = (
                int(row["start"])
                - 1
            )

            end = int(
                row["end"]
            )

            fh.write(
                f"{row['contig']}\t"
                f"{start0}\t"
                f"{end}\t"
                f"{row['locus_id']}\n"
            )


# ============================================================
# 3. BAM manifest
# ============================================================

bam_root = (
    ROOT
    / "06_mapping"
)

bam_rows = []


for bam in sorted(
    bam_root.glob(
        "*/*/*.sorted.bam"
    )
):

    sample = bam.parent.name
    group = bam.parent.parent.name

    bai1 = Path(
        str(bam)
        + ".bai"
    )

    bai2 = bam.with_suffix(
        ".bai"
    )

    indexed = (
        bai1.is_file()
        or bai2.is_file()
    )

    bam_rows.append({
        "coassembly":
            group,

        "sample":
            sample,

        "bam":
            str(bam),

        "bam_size_bytes":
            bam.stat().st_size,

        "index_exists":
            "YES"
            if indexed
            else "NO",

        "header_file":
            str(
                HEADERS
                / f"{sample}.header.sam"
            ),

        "idxstats_file":
            str(
                IDXSTATS
                / f"{sample}.idxstats.tsv"
            ),
    })


if len(bam_rows) != 18:
    raise RuntimeError(
        f"BAMs={len(bam_rows)} expected=18"
    )


write_tsv(
    OUT
    / "95D0_bam_manifest.tsv",
    bam_rows,
    [
        "coassembly",
        "sample",
        "bam",
        "bam_size_bytes",
        "index_exists",
        "header_file",
        "idxstats_file",
    ]
)


# ============================================================
# 4. Parse headers / Bowtie2 commands
# ============================================================

mapping_rows = []


for bamrow in bam_rows:

    sample = bamrow[
        "sample"
    ]

    header = Path(
        bamrow[
            "header_file"
        ]
    )

    idxstats = Path(
        bamrow[
            "idxstats_file"
        ]
    )

    if not header.is_file():
        raise RuntimeError(
            f"Header faltante: {header}"
        )

    if not idxstats.is_file():
        raise RuntimeError(
            f"idxstats faltante: {idxstats}"
        )


    pg_records = []

    bowtie_records = []

    with header.open(
        errors="replace"
    ) as fh:

        for line in fh:

            if not line.startswith(
                "@PG\t"
            ):
                continue

            pg = extract_pg_fields(
                line
            )

            pg_records.append(
                pg
            )

            searchable = (
                " ".join(
                    [
                        pg.get(
                            "ID",
                            ""
                        ),
                        pg.get(
                            "PN",
                            ""
                        ),
                        pg.get(
                            "CL",
                            ""
                        ),
                    ]
                )
            ).lower()

            if "bowtie2" in searchable:
                bowtie_records.append(
                    pg
                )


    bowtie_cl = ""

    if bowtie_records:

        # Prefer record explicitly PN=bowtie2.
        explicit = [
            x
            for x in bowtie_records
            if x.get(
                "PN",
                ""
            ).lower()
            == "bowtie2"
        ]

        chosen = (
            explicit[0]
            if explicit
            else bowtie_records[0]
        )

        bowtie_cl = chosen.get(
            "CL",
            ""
        )


    multimap_mode = (
        classify_bowtie2_command(
            bowtie_cl
        )
    )


    reference_names = set()

    mapped_idxstats = 0
    unmapped_idxstats = 0

    with idxstats.open() as fh:

        for line in fh:

            parts = (
                line.rstrip("\n")
                .split("\t")
            )

            if len(parts) < 4:
                continue

            ref = parts[0]

            if ref != "*":
                reference_names.add(
                    ref
                )

            mapped_idxstats += int(
                parts[2]
            )

            unmapped_idxstats += int(
                parts[3]
            )


    group_targets = (
        targets_by_group[
            bamrow["coassembly"]
        ]
    )

    target_contigs = {
        x["contig"]
        for x in group_targets
    }

    missing_contigs = sorted(
        target_contigs
        - reference_names
    )

    loci_on_present_refs = sum(
        r["contig"]
        in reference_names
        for r in group_targets
    )


    mapping_rows.append({
        "coassembly":
            bamrow[
                "coassembly"
            ],

        "sample":
            sample,

        "bam":
            bamrow[
                "bam"
            ],

        "index_exists":
            bamrow[
                "index_exists"
            ],

        "PG_records":
            len(pg_records),

        "bowtie2_PG_detected":
            "YES"
            if bowtie_records
            else "NO",

        "bowtie2_command":
            bowtie_cl,

        "multimap_reporting_mode":
            multimap_mode,

        "bam_reference_sequences":
            len(reference_names),

        "target_16S_loci_expected":
            len(group_targets),

        "target_16S_loci_on_present_refs":
            loci_on_present_refs,

        "target_contigs_missing_from_BAM_header":
            len(
                missing_contigs
            ),

        "missing_target_contigs":
            ";".join(
                missing_contigs
            ),

        "idxstats_mapped_alignments":
            mapped_idxstats,

        "idxstats_unmapped_records":
            unmapped_idxstats,
    })


mapping_fields = list(
    mapping_rows[0].keys()
)


write_tsv(
    OUT
    / "95D0_mapping_strategy.tsv",
    mapping_rows,
    mapping_fields
)


# ============================================================
# 5. Target summary by coassembly
# ============================================================

target_summary = []


for group in [
    "L1",
    "L2",
    "L3",
    "M1",
    "M2",
    "M3",
]:

    rr = targets_by_group[
        group
    ]

    target_summary.append({
        "coassembly":
            group,

        "final_bacterial_16S_loci":
            len(rr),

        "unique_exact_clusters":
            len({
                r[
                    "taxonomy_exact_cluster_id"
                ]
                for r in rr
            }),

        "target_contigs":
            len({
                r[
                    "contig"
                ]
                for r in rr
            }),

        "MAG_associated_loci":
            sum(
                r[
                    "exact_final18_contig_membership"
                ]
                == "YES"
                for r in rr
            ),
    })


write_tsv(
    OUT
    / "95D0_target_summary_by_coassembly.tsv",
    target_summary,
    list(
        target_summary[0].keys()
    )
)


# ============================================================
# 6. Global summary
# ============================================================

mode_counts = Counter(
    r[
        "multimap_reporting_mode"
    ]
    for r in mapping_rows
)


metrics = [
    (
        "final_primary_bacterial_16S_loci",
        len(targets)
    ),
    (
        "coassemblies_with_targets",
        len(
            targets_by_group
        )
    ),
    (
        "BAMs_total",
        len(
            bam_rows
        )
    ),
    (
        "BAMs_indexed",
        sum(
            r[
                "index_exists"
            ]
            == "YES"
            for r in bam_rows
        )
    ),
    (
        "BAMs_with_Bowtie2_PG",
        sum(
            r[
                "bowtie2_PG_detected"
            ]
            == "YES"
            for r in mapping_rows
        )
    ),
    (
        "BAMs_all_target_refs_present",
        sum(
            int(
                r[
                    "target_contigs_missing_from_BAM_header"
                ]
            )
            == 0
            for r in mapping_rows
        )
    ),
    (
        "BAMs_explicit_all_alignments",
        mode_counts.get(
            "all_alignments_explicit",
            0
        )
    ),
    (
        "BAMs_explicit_k_multiple",
        mode_counts.get(
            "multiple_alignments_k_explicit",
            0
        )
    ),
    (
        "BAMs_no_explicit_a_or_k",
        mode_counts.get(
            "no_explicit_a_or_k_detected",
            0
        )
    ),
    (
        "BAMs_Bowtie2_command_not_recovered",
        mode_counts.get(
            "bowtie2_command_not_recovered",
            0
        )
    ),
]


with (
    OUT
    / "95D0_global_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for key, value in metrics:
        fh.write(
            f"{key}\t{value}\n"
        )


# ============================================================
# 7. Methodological scope
# ============================================================

with (
    OUT
    / "95D0_methodological_scope.tsv"
).open("w") as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "quantification_source\texisting_whole_coassembly_BAMs\n"
    )

    fh.write(
        "target_unit\t121_curated_primary_bacterial_16S_loci\n"
    )

    fh.write(
        "existing_BAMs_assumed_valid_before_audit\tNO\n"
    )

    fh.write(
        "Bowtie2_multimapping_to_be_interpreted_from_header\tYES\n"
    )

    fh.write(
        "locus_depth_equals_taxon_abundance\tNO\n"
    )

    fh.write(
        "multiple_16S_operons_possible\tYES\n"
    )

    fh.write(
        "multimapping_between_similar_16S_loci_expected\tYES\n"
    )

    fh.write(
        "planned_primary_metrics\tmean_depth_median_depth_breadth_1x_5x_10x\n"
    )

    fh.write(
        "planned_sensitivity_metric\tMAPQ_filtered_depth\n"
    )

    fh.write(
        "library_size_normalization_planned\tYES\n"
    )


print(
    f"TARGET_LOCI={len(targets)}"
)

print(
    f"BAMS={len(bam_rows)}"
)

print(
    "95D0_PARSE=PASS"
)

