#!/usr/bin/env python3

import csv
import sys
from pathlib import Path
from collections import defaultdict, Counter


if len(sys.argv) < 2:
    raise SystemExit(
        "Usage:\n"
        "  prepare <targets.fna> <manifest.tsv> <query_dir>\n"
        "  parse <manifest.tsv> <context.tsv> <session_dir> <OUT>"
    )


MODE = sys.argv[1]


# ============================================================
# FASTA helper
# ============================================================

def read_fasta(path):

    records = {}

    current = None
    seq = []

    def save():

        if current is not None:

            records[
                current
            ] = "".join(
                seq
            )

    with path.open(
        errors="replace"
    ) as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                save()

                current = (
                    line[1:]
                    .split()[0]
                )

                seq = []

            else:

                seq.append(
                    line
                )

    save()

    return records


# ============================================================
# PREPARE
# ============================================================

if MODE == "prepare":

    if len(sys.argv) != 5:
        raise SystemExit(
            "prepare <targets.fna> <manifest.tsv> <query_dir>"
        )

    FASTA = Path(sys.argv[2])
    MANIFEST = Path(sys.argv[3])
    QUERY = Path(sys.argv[4])

    QUERY.mkdir(
        parents=True,
        exist_ok=True
    )


    with MANIFEST.open(
        newline=""
    ) as fh:

        targets = list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


    wanted = [
        r[
            "contig_id"
        ]
        for r in targets
    ]


    if len(wanted) != 7:

        raise RuntimeError(
            f"Expected 7 targets; found {len(wanted)}"
        )


    records = read_fasta(
        FASTA
    )


    missing = [
        cid
        for cid in wanted
        if cid not in records
    ]


    if missing:

        raise RuntimeError(
            "Missing FASTA records: "
            + ",".join(
                missing
            )
        )


    for row in targets:

        cid = row[
            "contig_id"
        ]

        seq = records[
            cid
        ]


        expected = int(
            row[
                "length_bp"
            ]
        )


        if len(seq) != expected:

            raise RuntimeError(
                f"Length mismatch {cid}: "
                f"{len(seq)} != {expected}"
            )


        out = (
            QUERY
            / f"{cid}.fna"
        )


        with out.open(
            "w"
        ) as fh:

            fh.write(
                f">{cid}\n"
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


    print(
        f"QUERY_FASTA_FILES={len(wanted)}"
    )

    print(
        "98C3B_PREPARE=PASS"
    )


# ============================================================
# PARSE
# ============================================================

elif MODE == "parse":

    if len(sys.argv) != 6:
        raise SystemExit(
            "parse <manifest.tsv> <context.tsv> "
            "<session_dir> <OUT>"
        )

    MANIFEST = Path(sys.argv[2])
    CONTEXT = Path(sys.argv[3])
    SESSION = Path(sys.argv[4])
    OUT = Path(sys.argv[5])

    OUT.mkdir(
        parents=True,
        exist_ok=True
    )


    # --------------------------------------------------------
    # Target manifest
    # --------------------------------------------------------

    with MANIFEST.open(
        newline=""
    ) as fh:

        manifest = list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


    targets = {
        r[
            "contig_id"
        ]:
        r
        for r in manifest
    }


    # --------------------------------------------------------
    # Target protein coordinates
    # --------------------------------------------------------

    context_by_contig = defaultdict(
        list
    )


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

                context_by_contig[
                    cid
                ].append(
                    row
                )


    # --------------------------------------------------------
    # Parse BAGEL AOI tables
    #
    # BAGEL wrapper uses:
    #   item[0] = AOI/query identifier
    #   item[2] = start
    #   item[3] = end
    #   item[5] = class
    #
    # according to the actual patched wrapper.
    # --------------------------------------------------------

    aoi_rows = []

    parse_failures = []


    for cid in targets:

        table = (
            SESSION
            / f"{cid}.AOI.table"
        )


        if (
            not table.is_file()
            or table.stat().st_size == 0
        ):
            continue


        with table.open(
            errors="replace"
        ) as fh:

            for line_number, line in enumerate(
                fh,
                start=1
            ):

                line = line.rstrip("\n")

                if not line.strip():
                    continue

                fields = line.split(
                    "\t"
                )


                if len(fields) < 6:

                    parse_failures.append({
                        "contig_id":
                            cid,

                        "line_number":
                            line_number,

                        "reason":
                            "fewer_than_6_columns",

                        "raw_line":
                            line,
                    })

                    continue


                try:

                    start = int(
                        float(
                            fields[2]
                        )
                    )

                    end = int(
                        float(
                            fields[3]
                        )
                    )

                except ValueError:

                    parse_failures.append({
                        "contig_id":
                            cid,

                        "line_number":
                            line_number,

                        "reason":
                            "non_numeric_AOI_coordinates",

                        "raw_line":
                            line,
                    })

                    continue


                if start > end:
                    start, end = end, start


                aoi_rows.append({
                    "contig_id":
                        cid,

                    "AOI_id":
                        fields[0],

                    "AOI_start":
                        start,

                    "AOI_end":
                        end,

                    "AOI_length":
                        end - start + 1,

                    "BAGEL_class":
                        fields[5],

                    "raw_field_1":
                        fields[1],

                    "raw_field_4":
                        fields[4],

                    "source_table":
                        str(
                            table
                        ),

                    "raw_line":
                        line,
                })


    # --------------------------------------------------------
    # Associate target protein coordinates with AOIs
    # --------------------------------------------------------

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


    target_rows = []


    for cid, manifest_row in targets.items():

        aois = aois_by_contig.get(
            cid,
            []
        )

        proteins = context_by_contig.get(
            cid,
            []
        )


        overlapping_AOIs = set()
        containing_AOIs = set()

        protein_overlap_detail = []


        for protein in proteins:

            try:

                pstart = int(
                    float(
                        protein[
                            "gene_start"
                        ]
                    )
                )

                pend = int(
                    float(
                        protein[
                            "gene_end"
                        ]
                    )
                )

            except Exception:
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


                overlap = (
                    pstart <= aend
                    and
                    pend >= astart
                )


                contained = (
                    pstart >= astart
                    and
                    pend <= aend
                )


                if overlap:

                    overlapping_AOIs.add(
                        aoi[
                            "AOI_id"
                        ]
                    )


                if contained:

                    containing_AOIs.add(
                        aoi[
                            "AOI_id"
                        ]
                    )


                if overlap:

                    protein_overlap_detail.append(
                        f"{protein['protein_id']}"
                        f"->{aoi['AOI_id']}"
                    )


        classes = sorted({
            x[
                "BAGEL_class"
            ]
            for x in aois
            if x[
                "BAGEL_class"
            ]
        })


        coordinates = sorted({
            f"{x['AOI_start']}-{x['AOI_end']}"
            for x in aois
        })


        target_rows.append({
            **manifest_row,

            "BAGEL_AOI_count":
                len(
                    aois
                ),

            "BAGEL_AOI_detected":
                int(
                    len(
                        aois
                    ) > 0
                ),

            "BAGEL_classes":
                ";".join(
                    classes
                ),

            "BAGEL_AOI_coordinates":
                ";".join(
                    coordinates
                ),

            "target_protein_count":
                len(
                    proteins
                ),

            "target_overlaps_BAGEL_AOI":
                int(
                    len(
                        overlapping_AOIs
                    ) > 0
                ),

            "target_fully_inside_BAGEL_AOI":
                int(
                    len(
                        containing_AOIs
                    ) > 0
                ),

            "target_AOI_overlap_detail":
                ";".join(
                    sorted(
                        protein_overlap_detail
                    )
                ),

            "BAGEL_interpretation":
                (
                    "BAGEL_AOI_contains_target"
                    if containing_AOIs
                    else
                    "BAGEL_AOI_overlaps_target"
                    if overlapping_AOIs
                    else
                    "BAGEL_AOI_elsewhere_on_contig"
                    if aois
                    else
                    "no_BAGEL_AOI_detected"
                ),
        })


    # --------------------------------------------------------
    # BAGEL evidence artifact inventory
    # --------------------------------------------------------

    artifact_patterns = {
        "bacteriocin_hmmsearch":
            "*.bacteriocin_hmmsearch",

        "bacteriocin_blast":
            "*.bacteriocin.blast_*",

        "AOI_predict":
            "*.AOI_*.predict",

        "AOI_faa":
            "*.AOI_*.faa",

        "AOI_gbk":
            "*.AOI_*.gbk",
    }


    artifact_rows = []


    for label, pattern in artifact_patterns.items():

        files = list(
            SESSION.glob(
                pattern
            )
        )

        artifact_rows.append({
            "artifact_type":
                label,

            "file_count":
                len(
                    files
                ),

            "nonempty_file_count":
                sum(
                    int(
                        f.stat().st_size > 0
                    )
                    for f in files
                ),
        })


    # --------------------------------------------------------
    # Writers
    # --------------------------------------------------------

    def write_table(
        path,
        rows,
        fields=None
    ):

        with path.open(
            "w",
            newline=""
        ) as fh:

            if rows:
                fieldnames = (
                    fields
                    or
                    list(
                        rows[0].keys()
                    )
                )

            else:
                fieldnames = (
                    fields
                    or
                    ["status"]
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
                w.writerow(
                    row
                )


    write_table(
        OUT
        / "98C3B_target_BAGEL_summary.tsv",
        target_rows
    )


    write_table(
        OUT
        / "98C3B_AOI_details.tsv",
        aoi_rows,
        [
            "contig_id",
            "AOI_id",
            "AOI_start",
            "AOI_end",
            "AOI_length",
            "BAGEL_class",
            "raw_field_1",
            "raw_field_4",
            "source_table",
            "raw_line",
        ]
    )


    write_table(
        OUT
        / "98C3B_parse_failures.tsv",
        parse_failures,
        [
            "contig_id",
            "line_number",
            "reason",
            "raw_line",
        ]
    )


    write_table(
        OUT
        / "98C3B_BAGEL_artifacts.tsv",
        artifact_rows
    )


    # --------------------------------------------------------
    # Global summary
    # --------------------------------------------------------

    class_counts = Counter()

    for aoi in aoi_rows:

        class_counts[
            aoi[
                "BAGEL_class"
            ]
        ] += 1


    with (
        OUT
        / "98C3B_global_summary.tsv"
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
                len(
                    target_rows
                )
            ),
            (
                "BAGEL_AOIs_total",
                len(
                    aoi_rows
                )
            ),
            (
                "contigs_with_any_BAGEL_AOI",
                sum(
                    int(
                        r[
                            "BAGEL_AOI_detected"
                        ]
                    )
                    for r in target_rows
                )
            ),
            (
                "contigs_without_BAGEL_AOI",
                sum(
                    int(
                        r[
                            "BAGEL_AOI_detected"
                        ]
                    ) == 0
                    for r in target_rows
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
                    for r in target_rows
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
                    for r in target_rows
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
                    for r in target_rows
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
        ]


        for klass, n in sorted(
            class_counts.items()
        ):

            metrics.append(
                (
                    f"BAGEL_class__{klass}",
                    n
                )
            )


        metrics += [
            (
                "BAGEL_AOI_support_equivalent_to_candidate_level_core_detection",
                "NO"
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
                "98C3C_exact_candidate_level_BAGEL_review_and_98C_synthesis"
            ),
        ]


        w.writerows(
            metrics
        )


    print(
        f"BAGEL_AOIS={len(aoi_rows)}"
    )

    print(
        "TARGETS_WITH_AOI="
        + str(
            sum(
                int(
                    r[
                        "BAGEL_AOI_detected"
                    ]
                )
                for r in target_rows
            )
        )
    )

    print(
        "TARGETS_OVERLAPPING_AOI="
        + str(
            sum(
                int(
                    r[
                        "target_overlaps_BAGEL_AOI"
                    ]
                )
                for r in target_rows
            )
        )
    )

    print(
        "98C3B_PARSE=PASS"
    )


else:

    raise SystemExit(
        f"Unknown mode: {MODE}"
    )
