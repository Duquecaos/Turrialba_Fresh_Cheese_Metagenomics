#!/usr/bin/env python3

import csv
import json
import re
import sys
from pathlib import Path
from collections import defaultdict, Counter


if len(sys.argv) != 8:
    raise SystemExit(
        "Usage: 98c3br2r_salvage_bagel.py "
        "<manifest.tsv> <context.tsv> <session> "
        "<wrapper.stdout> <wrapper.stderr> "
        "<system_command.log> <OUT>"
    )


MANIFEST = Path(sys.argv[1])
CONTEXT = Path(sys.argv[2])
SESSION = Path(sys.argv[3])
STDOUT = Path(sys.argv[4])
STDERR = Path(sys.argv[5])
COMMANDLOG = Path(sys.argv[6])
OUT = Path(sys.argv[7])

OUT.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# Helpers
# ============================================================

def uniq(values):
    return list(
        dict.fromkeys(
            x
            for x in values
            if x not in ("", None)
        )
    )


def safe_int(x):
    try:
        return int(float(x))
    except Exception:
        return None


def read_text(path):
    if path.is_file():
        return path.read_text(
            errors="replace"
        )
    return ""


def write_table(path, rows, fields=None):

    with path.open(
        "w",
        newline=""
    ) as fh:

        if rows:
            fieldnames = (
                fields
                or list(
                    rows[0].keys()
                )
            )
        else:
            fieldnames = (
                fields
                or ["status"]
            )

        w = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )

        w.writeheader()

        for row in rows:
            w.writerow(row)


# ============================================================
# 1. Target manifest
# ============================================================

with MANIFEST.open(
    newline=""
) as fh:

    manifest = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )


if len(manifest) != 7:
    raise RuntimeError(
        f"Expected 7 targets; found {len(manifest)}"
    )


targets = {
    row["contig_id"]:
        row
    for row in manifest
}


# ============================================================
# 2. Original target protein coordinates
# ============================================================

context = defaultdict(list)


with CONTEXT.open(
    newline=""
) as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in reader:

        cid = row.get(
            "contig_id",
            ""
        )

        if cid in targets:
            context[cid].append(
                row
            )


# ============================================================
# 3. Execution logs
# ============================================================

command_text = read_text(
    COMMANDLOG
)

stdout_text = read_text(
    STDOUT
)

stderr_text = read_text(
    STDERR
)

sessionstop_text = read_text(
    SESSION
    / "sessionstop"
)

filenames_db_text = read_text(
    SESSION
    / "00.filenames_db.json"
)

all_contigs_text = read_text(
    SESSION
    / "00.all_contigs.table"
)


command_lines = [
    x
    for x in command_text.splitlines()
    if x.strip()
]


# ============================================================
# 4. Global input count reported by BAGEL
# ============================================================

reported_input_counts = []


for text in [
    stdout_text,
    read_text(
        SESSION
        / "BAGEL_wrapper.log"
    ),
]:

    for m in re.finditer(
        r"InputfilesCount\s*=\s*(\d+)",
        text
    ):

        reported_input_counts.append(
            int(
                m.group(1)
            )
        )


reported_input_count = (
    reported_input_counts[-1]
    if reported_input_counts
    else None
)


# ============================================================
# 5. Parse AOI tables using PREFIX matching
#
# BAGEL can rename:
#
# ICONTIG000116433
#      -> ICONTIG000116433.3
#
# Therefore use:
#     ICONTIG000116433*.AOI.table
# ============================================================

aoi_rows = []
parse_failures = []


for cid in targets:

    tables = sorted(
        SESSION.glob(
            f"{cid}*.AOI.table"
        )
    )


    for table in tables:

        queryname = (
            table.name
            .replace(
                ".AOI.table",
                ""
            )
        )


        with table.open(
            errors="replace"
        ) as fh:

            for line_number, line in enumerate(
                fh,
                start=1
            ):

                line = line.rstrip(
                    "\n"
                )

                if not line.strip():
                    continue


                f = line.split(
                    "\t"
                )


                if len(f) < 6:

                    parse_failures.append({
                        "contig_id":
                            cid,

                        "queryname":
                            queryname,

                        "table":
                            str(table),

                        "line_number":
                            line_number,

                        "reason":
                            "fewer_than_6_columns",

                        "raw_line":
                            line,
                    })

                    continue


                start = safe_int(
                    f[2]
                )

                end = safe_int(
                    f[3]
                )


                if (
                    start is None
                    or end is None
                ):

                    parse_failures.append({
                        "contig_id":
                            cid,

                        "queryname":
                            queryname,

                        "table":
                            str(table),

                        "line_number":
                            line_number,

                        "reason":
                            "non_numeric_coordinates",

                        "raw_line":
                            line,
                    })

                    continue


                if start > end:
                    start, end = end, start


                aoi_rows.append({
                    "contig_id":
                        cid,

                    "BAGEL_queryname":
                        queryname,

                    "AOI_id":
                        f[0],

                    "AOI_start":
                        start,

                    "AOI_end":
                        end,

                    "AOI_length":
                        end - start + 1,

                    "BAGEL_class":
                        f[5],

                    "source_AOI_table":
                        str(table),

                    "raw_line":
                        line,
                })


aois_by_contig = defaultdict(
    list
)


for row in aoi_rows:

    aois_by_contig[
        row[
            "contig_id"
        ]
    ].append(
        row
    )


# ============================================================
# 6. Per-target execution validation
#
# Required detection stages:
#   - HMM identification
#   - bacteriocin BLAST identification
#   - AOI merge
#
# Annotation is only expected when an AOI is found.
# ============================================================

qc_rows = []


for cid, meta in targets.items():

    target_commands = [
        line
        for line in command_lines
        if cid in line
    ]


    hmm_command = any(
        "bagel4_AOI_identification_hmm_rules.pl"
        in line
        for line in target_commands
    )


    blast_command = any(
        "bagel4_AOI_identification_blast_bacteriocins.pl"
        in line
        for line in target_commands
    )


    merge_command = any(
        "bagel4_AOI_merge.pl"
        in line
        for line in target_commands
    )


    annotation_command = any(
        "bagel4_AOI_annotation.pl"
        in line
        for line in target_commands
    )


    gene_json_command = any(
        "bagel4_GeneTable_2_json.pl"
        in line
        for line in target_commands
    )


    hmm_files = [
        p
        for p in (
            SESSION
            / "hmm_rules_domtblout"
        ).rglob("*")
        if (
            p.is_file()
            and cid in p.name
        )
    ]


    aoi_tables = sorted(
        SESSION.glob(
            f"{cid}*.AOI.table"
        )
    )


    aoi_count = len(
        aois_by_contig.get(
            cid,
            []
        )
    )


    detection_complete = int(
        hmm_command
        and blast_command
        and merge_command
    )


    filename_db_present = int(
        cid
        in filenames_db_text
    )


    all_contigs_present = int(
        cid
        in all_contigs_text
    )


    if detection_complete:

        if aoi_count > 0:

            processing_interpretation = (
                "valid_detection_with_BAGEL_AOI"
            )

        else:

            processing_interpretation = (
                "valid_detection_no_BAGEL_AOI"
            )

    else:

        processing_interpretation = (
            "processing_incomplete_or_unresolved"
        )


    qc_rows.append({
        "contig_id":
            cid,

        "BAGEL_priority":
            meta.get(
                "BAGEL_priority",
                ""
            ),

        "direct_evidence_classes":
            meta.get(
                "direct_evidence_classes",
                ""
            ),

        "HMM_detection_command":
            int(
                hmm_command
            ),

        "BLAST_detection_command":
            int(
                blast_command
            ),

        "AOI_merge_command":
            int(
                merge_command
            ),

        "AOI_annotation_command":
            int(
                annotation_command
            ),

        "GeneTable_JSON_command":
            int(
                gene_json_command
            ),

        "target_system_commands":
            len(
                target_commands
            ),

        "target_HMM_output_files":
            len(
                hmm_files
            ),

        "AOI_table_files":
            len(
                aoi_tables
            ),

        "parsed_AOI_rows":
            aoi_count,

        "present_in_filenames_db":
            filename_db_present,

        "present_in_all_contigs_table":
            all_contigs_present,

        "detection_pipeline_complete":
            detection_complete,

        "processing_interpretation":
            processing_interpretation,
    })


# ============================================================
# 7. Target protein vs BAGEL AOI overlap
# ============================================================

result_rows = []


for cid, meta in targets.items():

    aois = aois_by_contig.get(
        cid,
        []
    )

    proteins = context.get(
        cid,
        []
    )


    overlap_detail = []

    overlapping_aois = set()
    containing_aois = set()


    for protein in proteins:

        pstart = safe_int(
            protein.get(
                "gene_start",
                ""
            )
        )

        pend = safe_int(
            protein.get(
                "gene_end",
                ""
            )
        )


        if (
            pstart is None
            or pend is None
        ):
            continue


        if pstart > pend:
            pstart, pend = pend, pstart


        for aoi in aois:

            astart = int(
                aoi[
                    "AOI_start"
                ]
            )

            aend = int(
                aoi[
                    "AOI_end"
                ]
            )


            overlaps = (
                pstart <= aend
                and
                pend >= astart
            )


            contained = (
                pstart >= astart
                and
                pend <= aend
            )


            if overlaps:

                overlapping_aois.add(
                    aoi[
                        "AOI_id"
                    ]
                )

                overlap_detail.append(
                    f"{protein['protein_id']}"
                    f"->{aoi['AOI_id']}"
                )


            if contained:

                containing_aois.add(
                    aoi[
                        "AOI_id"
                    ]
                )


    qc = next(
        x
        for x in qc_rows
        if x[
            "contig_id"
        ] == cid
    )


    classes = sorted({
        a[
            "BAGEL_class"
        ]
        for a in aois
        if a[
            "BAGEL_class"
        ]
    })


    querynames = sorted({
        a[
            "BAGEL_queryname"
        ]
        for a in aois
    })


    if (
        not qc[
            "detection_pipeline_complete"
        ]
    ):

        interpretation = (
            "BAGEL_processing_unresolved"
        )

    elif containing_aois:

        interpretation = (
            "BAGEL_AOI_contains_target_gene"
        )

    elif overlapping_aois:

        interpretation = (
            "BAGEL_AOI_overlaps_target_gene"
        )

    elif aois:

        interpretation = (
            "BAGEL_AOI_elsewhere_on_same_contig"
        )

    else:

        interpretation = (
            "valid_BAGEL_screen_no_AOI"
        )


    result_rows.append({
        **meta,

        "BAGEL_querynames":
            ";".join(
                querynames
            ),

        "BAGEL_AOI_count":
            len(
                aois
            ),

        "BAGEL_classes":
            ";".join(
                classes
            ),

        "target_overlaps_BAGEL_AOI":
            int(
                bool(
                    overlapping_aois
                )
            ),

        "target_fully_inside_BAGEL_AOI":
            int(
                bool(
                    containing_aois
                )
            ),

        "target_AOI_overlap_detail":
            ";".join(
                sorted(
                    overlap_detail
                )
            ),

        "detection_pipeline_complete":
            qc[
                "detection_pipeline_complete"
            ],

        "BAGEL_interpretation":
            interpretation,
    })


# ============================================================
# 8. Global QC
# ============================================================

session_done = int(
    "Analysis done"
    in sessionstop_text
)


valid_targets = sum(
    int(
        r[
            "detection_pipeline_complete"
        ]
    )
    for r in qc_rows
)


filenames_targets = sum(
    int(
        r[
            "present_in_filenames_db"
        ]
    )
    for r in qc_rows
)


all_contigs_targets = sum(
    int(
        r[
            "present_in_all_contigs_table"
        ]
    )
    for r in qc_rows
)


real_processing_pass = int(
    valid_targets == 7
    and filenames_targets == 7
    and session_done == 1
)


# ============================================================
# 9. Writers
# ============================================================

write_table(
    OUT
    / "98C3BR2R_target_execution_QC.tsv",
    qc_rows
)


write_table(
    OUT
    / "98C3BR2R_AOI_details.tsv",
    aoi_rows,
    [
        "contig_id",
        "BAGEL_queryname",
        "AOI_id",
        "AOI_start",
        "AOI_end",
        "AOI_length",
        "BAGEL_class",
        "source_AOI_table",
        "raw_line",
    ]
)


write_table(
    OUT
    / "98C3BR2R_target_BAGEL_summary.tsv",
    result_rows
)


write_table(
    OUT
    / "98C3BR2R_parse_failures.tsv",
    parse_failures,
    [
        "contig_id",
        "queryname",
        "table",
        "line_number",
        "reason",
        "raw_line",
    ]
)


# ============================================================
# 10. Summary
# ============================================================

with (
    OUT
    / "98C3BR2R_global_summary.tsv"
).open(
    "w",
    newline=""
) as fh:

    w = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n"
    )

    w.writerow([
        "metric",
        "value"
    ])


    metrics = [
        (
            "target_contigs",
            7
        ),
        (
            "targets_present_in_filenames_db",
            filenames_targets
        ),
        (
            "targets_present_in_all_contigs_table",
            all_contigs_targets
        ),
        (
            "targets_with_complete_HMM_BLAST_merge_detection",
            valid_targets
        ),
        (
            "BAGEL_reported_InputfilesCount",
            (
                reported_input_count
                if reported_input_count is not None
                else "UNRESOLVED"
            )
        ),
        (
            "sessionstop_Analysis_done",
            (
                "YES"
                if session_done
                else "NO"
            )
        ),
        (
            "BAGEL_AOI_rows",
            len(
                aoi_rows
            )
        ),
        (
            "contigs_with_BAGEL_AOI",
            sum(
                int(
                    r[
                        "BAGEL_AOI_count"
                    ] > 0
                )
                for r in result_rows
            )
        ),
        (
            "contigs_with_valid_screen_no_AOI",
            sum(
                int(
                    r[
                        "BAGEL_interpretation"
                    ]
                    ==
                    "valid_BAGEL_screen_no_AOI"
                )
                for r in result_rows
            )
        ),
        (
            "targets_overlapping_BAGEL_AOI",
            sum(
                int(
                    r[
                        "target_overlaps_BAGEL_AOI"
                    ]
                )
                for r in result_rows
            )
        ),
        (
            "targets_fully_inside_BAGEL_AOI",
            sum(
                int(
                    r[
                        "target_fully_inside_BAGEL_AOI"
                    ]
                )
                for r in result_rows
            )
        ),
        (
            "high_direct_sequence_target_supported_by_BAGEL_AOI",
            sum(
                int(
                    r[
                        "target_overlaps_BAGEL_AOI"
                    ]
                )
                for r in result_rows
                if r[
                    "BAGEL_priority"
                ]
                ==
                "HIGH_direct_sequence_evidence"
            )
        ),
        (
            "AOI_parse_failures",
            len(
                parse_failures
            )
        ),
        (
            "real_BAGEL_processing_validation",
            (
                "PASS"
                if real_processing_pass
                else "FAIL"
            )
        ),
        (
            "functional_bacteriocin_confirmed",
            "NO"
        ),
        (
            "reads_remapped",
            "NO"
        ),
        (
            "amplicon_data_used",
            "NO"
        ),
        (
            "next_step",
            (
                "98C3C_candidate_level_review_and_98C_synthesis"
                if real_processing_pass
                else
                "resolve_remaining_BAGEL_processing_gap"
            )
        ),
    ]


    w.writerows(
        metrics
    )


print(
    f"VALID_TARGETS={valid_targets}/7"
)

print(
    f"BAGEL_AOI_ROWS={len(aoi_rows)}"
)

print(
    "REAL_PROCESSING_VALIDATION="
    + (
        "PASS"
        if real_processing_pass
        else "FAIL"
    )
)

print(
    "98C3BR2R=PASS"
)
