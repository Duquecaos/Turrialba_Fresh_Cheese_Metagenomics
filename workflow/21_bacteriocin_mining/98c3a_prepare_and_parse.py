#!/usr/bin/env python3

import csv
import re
import sys
from pathlib import Path


if len(sys.argv) != 5:
    raise SystemExit(
        "Usage:\n"
        "  prepare: 98c3a_prepare_and_parse.py prepare <ROOT> <PRIMARY.tsv> <OUT>\n"
        "  parse:   98c3a_prepare_and_parse.py parse   <PRIMARY.tsv> <ANTISMASH_OUT> <OUT>"
    )


MODE = sys.argv[1]


def read_fasta(path):
    records = {}

    header = None
    seq = []

    with path.open(errors="replace") as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if header is not None:
                    records[
                        header.split()[0]
                    ] = "".join(seq)

                header = line[1:]
                seq = []

            else:
                seq.append(line)

    if header is not None:
        records[
            header.split()[0]
        ] = "".join(seq)

    return records


def read_fasta_selected(path, wanted):

    found = {}

    current_id = None
    seq = []
    retain = False

    def save():

        if retain and current_id:
            found[
                current_id
            ] = "".join(seq)

    with path.open(errors="replace") as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                save()

                current_id = (
                    line[1:]
                    .split()[0]
                )

                retain = (
                    current_id in wanted
                )

                seq = []

            elif retain:

                seq.append(line)

    save()

    return found


if MODE == "prepare":

    ROOT = Path(sys.argv[2])
    PRIMARY = Path(sys.argv[3])
    OUT = Path(sys.argv[4])

    OUT.mkdir(
        parents=True,
        exist_ok=True
    )


    with PRIMARY.open(
        newline=""
    ) as fh:

        rows = list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


    targets = [
        r["contig_id"]
        for r in rows
    ]

    targets = list(
        dict.fromkeys(
            targets
        )
    )

    wanted = set(
        targets
    )


    if len(targets) != 40:

        raise RuntimeError(
            f"Expected 40 primary contigs; found {len(targets)}"
        )


    # ========================================================
    # Discover nucleotide FASTA candidates
    # ========================================================

    candidates = []

    for subdir in [
        ROOT / "69_incomplete_unique_contigs",
        ROOT / "70_incomplete_gene_catalog",
    ]:

        if not subdir.is_dir():
            continue

        for path in subdir.rglob("*"):

            if (
                path.is_file()
                and path.suffix.lower()
                in {
                    ".fa",
                    ".fna",
                    ".fasta",
                    ".fas",
                }
            ):

                candidates.append(
                    path
                )


    if not candidates:

        raise RuntimeError(
            "No nucleotide FASTA candidates found in steps 69/70"
        )


    audit = []

    best_path = None
    best_found = {}


    for path in sorted(candidates):

        found = read_fasta_selected(
            path,
            wanted
        )

        audit.append({
            "path":
                str(path),

            "targets_found":
                len(found),

            "targets_expected":
                len(wanted),
        })


        if len(found) > len(best_found):

            best_found = found
            best_path = path


    with (
        OUT
        / "98C3A_fasta_discovery.tsv"
    ).open(
        "w",
        newline=""
    ) as fh:

        w = csv.DictWriter(
            fh,
            fieldnames=[
                "path",
                "targets_found",
                "targets_expected",
            ],
            delimiter="\t",
            lineterminator="\n"
        )

        w.writeheader()
        w.writerows(
            audit
        )


    if (
        best_path is None
        or len(best_found) != len(wanted)
    ):

        missing = sorted(
            wanted
            - set(best_found)
        )

        (
            OUT
            / "98C3A_missing_contigs.txt"
        ).write_text(
            "\n".join(
                missing
            )
            + "\n"
        )

        raise RuntimeError(
            "No single FASTA contains all 40 target contigs. "
            f"Best={best_path}; found={len(best_found)}/40"
        )


    # ========================================================
    # Validate lengths against 98C2
    # ========================================================

    expected_lengths = {
        r[
            "contig_id"
        ]:
        int(
            r[
                "length_bp"
            ]
        )
        for r in rows
    }


    manifest = []

    for cid in targets:

        sequence = (
            best_found[
                cid
            ]
            .replace(" ", "")
            .upper()
        )

        observed = len(
            sequence
        )

        expected = expected_lengths[
            cid
        ]


        if observed != expected:

            raise RuntimeError(
                f"Length mismatch {cid}: "
                f"FASTA={observed}, table={expected}"
            )


        manifest.append({
            "contig_id":
                cid,

            "length_bp":
                observed,

            "source_fasta":
                str(
                    best_path
                ),
        })


    # ========================================================
    # Output target FASTA
    # ========================================================

    target_fasta = (
        OUT
        / "98C3A_primary40_contigs.fna"
    )


    with target_fasta.open(
        "w"
    ) as fh:

        for cid in targets:

            seq = best_found[
                cid
            ].upper()

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


    with (
        OUT
        / "98C3A_primary40_manifest.tsv"
    ).open(
        "w",
        newline=""
    ) as fh:

        w = csv.DictWriter(
            fh,
            fieldnames=[
                "contig_id",
                "length_bp",
                "source_fasta",
            ],
            delimiter="\t",
            lineterminator="\n"
        )

        w.writeheader()
        w.writerows(
            manifest
        )


    (
        OUT
        / "98C3A_selected_source_fasta.txt"
    ).write_text(
        str(
            best_path
        )
        + "\n"
    )


    print(
        f"TARGETS={len(targets)}"
    )

    print(
        f"SOURCE_FASTA={best_path}"
    )

    print(
        "98C3A_PREPARE=PASS"
    )


elif MODE == "parse":

    PRIMARY = Path(sys.argv[2])
    ASOUT = Path(sys.argv[3])
    OUT = Path(sys.argv[4])

    OUT.mkdir(
        parents=True,
        exist_ok=True
    )


    with PRIMARY.open(
        newline=""
    ) as fh:

        primary = list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


    primary_by_contig = {
        r["contig_id"]:
            r
        for r in primary
    }


    RIPP_RE = re.compile(
        r"bacterioc|"
        r"RiPP|"
        r"lanthipept|"
        r"lantipept|"
        r"lassopept|"
        r"thiopept|"
        r"microviridin|"
        r"ranthipept|"
        r"sactipept|"
        r"glycocin|"
        r"RaS-RiPP|"
        r"RRE-containing",
        re.I
    )


    # ========================================================
    # Parse antiSMASH region GBKs
    # ========================================================

    regions = []


    for path in sorted(
        ASOUT.rglob(
            "*.region*.gbk"
        )
    ):

        basename = path.name

        contig = basename.split(
            ".region",
            1
        )[0]


        text = path.read_text(
            errors="replace"
        )

        lines = text.splitlines()

        region_location = ""
        products = []
        contig_edge = ""


        in_region = False

        i = 0

        while i < len(lines):

            line = lines[i]

            feature_key = (
                line[5:21].strip()
                if len(line) >= 21
                else ""
            )


            if feature_key:

                if feature_key == "region":

                    if in_region:
                        break

                    in_region = True

                    region_location = (
                        line[21:]
                        .strip()
                    )

                elif in_region:
                    break


            elif in_region:

                stripped = (
                    line.strip()
                )


                if stripped.startswith(
                    "/product="
                ):

                    value = (
                        stripped
                        .split(
                            "=",
                            1
                        )[1]
                        .strip()
                    )


                    if value.startswith('"'):

                        value = value[1:]

                        while (
                            '"' not in value
                            and i + 1 < len(lines)
                        ):

                            i += 1

                            value += (
                                lines[i]
                                .strip()
                            )


                        value = (
                            value.split(
                                '"',
                                1
                            )[0]
                        )

                    else:

                        value = (
                            value.strip('"')
                        )


                    products.append(
                        value
                    )


                elif stripped.startswith(
                    "/contig_edge="
                ):

                    contig_edge = (
                        stripped
                        .split(
                            "=",
                            1
                        )[1]
                        .strip()
                        .strip('"')
                    )


            i += 1


        products = list(
            dict.fromkeys(
                products
            )
        )


        searchable = " ".join(
            products
        )


        regions.append({
            "contig_id":
                contig,

            "region_file":
                str(path),

            "region_location":
                region_location,

            "products":
                ";".join(
                    products
                ),

            "bacteriocin_RiPP_relevant":
                int(
                    bool(
                        RIPP_RE.search(
                            searchable
                        )
                    )
                ),

            "contig_edge":
                contig_edge,
        })


    # ========================================================
    # Summarize per target contig
    # ========================================================

    by_contig = {}

    for cid in primary_by_contig:

        rr = [
            x
            for x in regions
            if x[
                "contig_id"
            ] == cid
        ]


        ripp = [
            x
            for x in rr
            if x[
                "bacteriocin_RiPP_relevant"
            ] == 1
        ]


        source = primary_by_contig[
            cid
        ]


        by_contig[
            cid
        ] = {
            "contig_id":
                cid,

            "primary_positive_protein_count":
                source[
                    "primary_positive_protein_count"
                ],

            "primary_protein_ids":
                source[
                    "primary_protein_ids"
                ],

            "direct_evidence_classes":
                source[
                    "direct_evidence_classes"
                ],

            "best_query_candidates":
                source[
                    "best_query_candidates"
                ],

            "GA_HMM_models":
                source[
                    "GA_HMM_models"
                ],

            "length_bp":
                source[
                    "length_bp"
                ],

            "producers":
                source[
                    "producers"
                ],

            "strongest_incomplete_context":
                source[
                    "strongest_incomplete_context"
                ],

            "n_antismash_regions":
                len(
                    rr
                ),

            "n_RiPP_relevant_regions":
                len(
                    ripp
                ),

            "antismash_products":
                ";".join(
                    sorted({
                        x[
                            "products"
                        ]
                        for x in rr
                        if x[
                            "products"
                        ]
                    })
                ),

            "RiPP_products":
                ";".join(
                    sorted({
                        x[
                            "products"
                        ]
                        for x in ripp
                        if x[
                            "products"
                        ]
                    })
                ),

            "antismash_same_contig_support":
                int(
                    len(
                        ripp
                    ) > 0
                ),

            "interpretation":
                (
                    "orthogonal_antismash_RiPP_support"
                    if ripp
                    else
                    "no_antismash_RiPP_region_detected"
                ),
        }


    # ========================================================
    # Write tables
    # ========================================================

    region_out = (
        OUT
        / "98C3A_antismash_regions.tsv"
    )


    region_fields = [
        "contig_id",
        "region_file",
        "region_location",
        "products",
        "bacteriocin_RiPP_relevant",
        "contig_edge",
    ]


    with region_out.open(
        "w",
        newline=""
    ) as fh:

        w = csv.DictWriter(
            fh,
            fieldnames=region_fields,
            delimiter="\t",
            lineterminator="\n"
        )

        w.writeheader()

        w.writerows(
            regions
        )


    contig_rows = [
        by_contig[
            cid
        ]
        for cid in sorted(
            by_contig
        )
    ]


    contig_out = (
        OUT
        / "98C3A_primary_contig_antismash_summary.tsv"
    )


    with contig_out.open(
        "w",
        newline=""
    ) as fh:

        w = csv.DictWriter(
            fh,
            fieldnames=list(
                contig_rows[0]
                .keys()
            ),
            delimiter="\t",
            lineterminator="\n"
        )

        w.writeheader()

        w.writerows(
            contig_rows
        )


    total_regions = len(
        regions
    )

    ripp_regions = sum(
        int(
            x[
                "bacteriocin_RiPP_relevant"
            ]
        )
        for x in regions
    )

    supported_contigs = sum(
        int(
            x[
                "antismash_same_contig_support"
            ]
        )
        for x in contig_rows
    )


    with (
        OUT
        / "98C3A_global_summary.tsv"
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
                "primary_contigs_submitted",
                len(
                    contig_rows
                )
            ),
            (
                "antismash_regions_total",
                total_regions
            ),
            (
                "RiPP_relevant_regions",
                ripp_regions
            ),
            (
                "primary_contigs_with_RiPP_antismash_support",
                supported_contigs
            ),
            (
                "primary_contigs_without_RiPP_antismash_support",
                len(
                    contig_rows
                )
                -
                supported_contigs
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
                "98C3B_targeted_BAGEL_on_priority_antismash_negative_contexts"
            ),
        ]

        w.writerows(
            metrics
        )


    print(
        f"REGIONS={total_regions}"
    )

    print(
        f"RIPP_REGIONS={ripp_regions}"
    )

    print(
        f"SUPPORTED_CONTIGS={supported_contigs}"
    )

    print(
        "98C3A_PARSE=PASS"
    )


else:

    raise SystemExit(
        f"Unknown mode: {MODE}"
    )
