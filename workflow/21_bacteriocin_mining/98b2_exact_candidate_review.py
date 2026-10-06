#!/usr/bin/env python3

import csv
import re
import sys
from pathlib import Path
from collections import defaultdict


if len(sys.argv) != 5:
    raise SystemExit(
        "Usage: 98b2_exact_candidate_review.py "
        "<ROOT> <ATTRLOC.tsv> <ANTISMASH.tsv> <OUT>"
    )


ROOT = Path(sys.argv[1])
ATTRLOC = Path(sys.argv[2])
ANTI_TSV = Path(sys.argv[3])
OUT = Path(sys.argv[4])

OUT.mkdir(parents=True, exist_ok=True)


BAGEL35 = ROOT / "35_bagel4_test_lactis"
BAGEL36 = ROOT / "36_bagel4_BAL4"

COMP = ROOT / "29_comparippson_lactococcin"

MASTER_FASTA = (
    ROOT
    / "38_bacteriocin_master"
    / "structural_candidates_protein.faa"
)

MASTER_META = (
    ROOT
    / "38_bacteriocin_master"
    / "structural_candidates_metadata.tsv"
)


# ============================================================
# Helpers
# ============================================================

def norm(s):
    return re.sub(
        r"[^A-Za-z0-9]+",
        "",
        str(s)
    ).lower()


def parse_fasta(path):
    records = []

    if not path.is_file():
        return records

    header = None
    seq = []

    with path.open(
        errors="replace"
    ) as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if header is not None:
                    records.append(
                        (
                            header,
                            "".join(seq)
                            .replace("*", "")
                            .upper()
                        )
                    )

                header = line[1:]
                seq = []

            else:
                seq.append(line)

    if header is not None:

        records.append(
            (
                header,
                "".join(seq)
                .replace("*", "")
                .upper()
            )
        )

    return records


def sequence_relationship(a, b):

    a = a.replace("*", "").upper()
    b = b.replace("*", "").upper()

    if not a or not b:
        return "none"

    if a == b:
        return "exact"

    # Require meaningful length for containment.
    if min(len(a), len(b)) >= 20:

        if a in b or b in a:
            return "contained"

    return "none"


def candidate_local_token(candidate, contig):

    pos = candidate.find(contig)

    if pos >= 0:

        tail = candidate[
            pos + len(contig):
        ]

        tail = tail.strip("_- ")

        if tail:
            return tail

    parts = candidate.split("_")

    for x in reversed(parts):

        low = x.lower()

        if (
            low.startswith("orf")
            or low.startswith("sorf")
        ):
            return x

    return candidate


def read_small_text(path):

    try:

        if (
            not path.is_file()
            or path.stat().st_size > 20_000_000
        ):
            return ""

        return path.read_text(
            errors="replace"
        )

    except Exception:
        return ""


def genbank_translations(path):

    text = read_small_text(path)

    if not text:
        return []

    translations = []

    for match in re.finditer(
        r'/translation="([^"]+)"',
        text,
        flags=re.S
    ):

        seq = re.sub(
            r"\s+",
            "",
            match.group(1)
        ).replace("*", "").upper()

        if seq:
            translations.append(seq)

    return translations


# ============================================================
# Input tables
# ============================================================

with ATTRLOC.open(
    newline=""
) as fh:

    loci = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )


anti = {}

with ANTI_TSV.open(
    newline=""
) as fh:

    for row in csv.DictReader(
        fh,
        delimiter="\t"
    ):

        anti[
            row["locus_id"]
        ] = row


# ============================================================
# Master candidate sequences
# ============================================================

master_records = parse_fasta(
    MASTER_FASTA
)

master_meta_text = read_small_text(
    MASTER_META
)


def get_master_sequence(candidate):

    exact_header = []
    normalized_header = []

    cn = norm(candidate)

    for header, seq in master_records:

        if candidate in header:
            exact_header.append(
                (
                    header,
                    seq
                )
            )

        elif cn and cn in norm(header):
            normalized_header.append(
                (
                    header,
                    seq
                )
            )

    hits = (
        exact_header
        if exact_header
        else normalized_header
    )

    if hits:

        return (
            hits[0][1],
            hits[0][0],
            len(hits)
        )

    return (
        "",
        "",
        0
    )


# ============================================================
# Identify BAGEL run belonging to the actual MAG
# ============================================================

def bagel_roots_for_mag(mag):

    roots = []

    # BAGEL36 contains per-MAG directories.
    p = BAGEL36 / mag

    if p.is_dir():
        roots.append(
            (
                "BAGEL36",
                p
            )
        )

    # BAGEL35 was a specific run. Only count it as same-MAG
    # if its queryfolder actually contains this MAG.
    q = (
        BAGEL35
        / "session"
        / "queryfolder"
        / f"{mag}.fna"
    )

    if q.is_file():
        roots.append(
            (
                "BAGEL35",
                BAGEL35
            )
        )

    return roots


# ============================================================
# Comparippson inventory
# ============================================================

comp_files = []

if COMP.is_dir():

    for p in COMP.rglob("*"):

        if p.is_file():
            comp_files.append(p)


comp_fastas = []

for p in comp_files:

    if p.suffix.lower() in {
        ".faa",
        ".fa",
        ".fasta",
        ".fas"
    }:

        for h, s in parse_fasta(p):

            comp_fastas.append(
                (
                    str(p),
                    h,
                    s
                )
            )


# ============================================================
# Per-candidate review
# ============================================================

candidate_rows = []


for locus in loci:

    locus_id = locus["locus_id"]
    mag = locus["MAG"]
    contig = locus["contig"]

    candidates = [
        x.strip()
        for x in locus[
            "candidates"
        ].split(";")
        if x.strip()
    ]

    anti_row = anti.get(
        locus_id,
        {}
    )

    anti_context = int(
        anti_row.get(
            "antismash_bacteriocin_RiPP_region_same_contig",
            "0"
        )
        == "1"
    )

    anti_products = anti_row.get(
        "antismash_region_products",
        ""
    )

    anti_region_paths = []

    raw_paths = anti_row.get(
        "antismash_region_files",
        ""
    )

    for x in re.split(
        r"[;,]",
        raw_paths
    ):

        x = x.strip()

        if x:
            p = Path(x)

            if p.is_file():
                anti_region_paths.append(p)


    anti_translations = []

    for p in anti_region_paths:

        for seq in genbank_translations(p):

            anti_translations.append(
                (
                    str(p),
                    seq
                )
            )


    actual_bagel_roots = bagel_roots_for_mag(
        mag
    )


    for candidate in candidates:

        local_token = candidate_local_token(
            candidate,
            contig
        )

        candidate_seq, candidate_header, master_hit_count = (
            get_master_sequence(
                candidate
            )
        )


        # ----------------------------------------------------
        # antiSMASH exact candidate sequence
        # ----------------------------------------------------

        anti_exact = 0
        anti_contained = 0
        anti_match_file = ""

        if candidate_seq:

            for path, seq in anti_translations:

                rel = sequence_relationship(
                    candidate_seq,
                    seq
                )

                if rel == "exact":

                    anti_exact = 1
                    anti_match_file = path
                    break

                elif rel == "contained":

                    anti_contained = 1

                    if not anti_match_file:
                        anti_match_file = path


        # ----------------------------------------------------
        # BAGEL: restrict to SAME MAG
        # ----------------------------------------------------

        bagel_same_contig_files = 0
        bagel_text_hit = 0
        bagel_exact_seq = 0
        bagel_contained_seq = 0

        bagel_text_files = []
        bagel_seq_files = []

        candidate_norm = norm(
            candidate
        )

        token_norm = norm(
            local_token
        )


        for bagel_label, bagel_root in actual_bagel_roots:

            for p in bagel_root.rglob("*"):

                if not p.is_file():
                    continue

                path_str = str(p)
                path_norm = norm(
                    path_str
                )

                same_contig_path = (
                    contig in path_str
                )

                if same_contig_path:
                    bagel_same_contig_files += 1


                # --------------------------------------------
                # Text / filename candidate evidence.
                # Must belong to actual MAG and preserve
                # contig + candidate relationship.
                # --------------------------------------------

                hit_here = False

                if (
                    same_contig_path
                    and token_norm
                    and token_norm in path_norm
                ):
                    hit_here = True


                text = ""

                if p.suffix.lower() in {
                    ".tsv",
                    ".txt",
                    ".table",
                    ".json",
                    ".gbk",
                    ".predict",
                    ".detail",
                    ".faa",
                    ".fna",
                    ".log"
                }:

                    text = read_small_text(p)

                    text_norm = norm(
                        text
                    )

                    # Full candidate identifier is safest.
                    if (
                        candidate_norm
                        and candidate_norm in text_norm
                    ):
                        hit_here = True

                    # Or both exact contig and local ORF token
                    # must occur in the same file.
                    elif (
                        contig in text
                        and token_norm
                        and token_norm in text_norm
                    ):
                        hit_here = True


                if hit_here:

                    bagel_text_hit = 1

                    if len(
                        bagel_text_files
                    ) < 10:

                        bagel_text_files.append(
                            path_str
                        )


                # --------------------------------------------
                # Protein sequence comparison
                # --------------------------------------------

                if (
                    candidate_seq
                    and p.suffix.lower()
                    in {
                        ".faa",
                        ".fa",
                        ".fasta"
                    }
                ):

                    for h, seq in parse_fasta(p):

                        rel = sequence_relationship(
                            candidate_seq,
                            seq
                        )

                        if rel == "exact":

                            bagel_exact_seq = 1

                            if len(
                                bagel_seq_files
                            ) < 10:

                                bagel_seq_files.append(
                                    path_str
                                )

                        elif rel == "contained":

                            bagel_contained_seq = 1

                            if len(
                                bagel_seq_files
                            ) < 10:

                                bagel_seq_files.append(
                                    path_str
                                )


        # ----------------------------------------------------
        # Comparippson candidate-level evidence
        # ----------------------------------------------------

        comp_text_hit = 0
        comp_exact_seq = 0
        comp_contained_seq = 0

        comp_text_files = []
        comp_seq_files = []


        for p in comp_files:

            if p.stat().st_size > 20_000_000:
                continue

            text = read_small_text(p)

            if (
                candidate_norm
                and
                candidate_norm in norm(text)
            ):

                comp_text_hit = 1

                if len(comp_text_files) < 10:
                    comp_text_files.append(
                        str(p)
                    )


        if candidate_seq:

            for path, header, seq in comp_fastas:

                rel = sequence_relationship(
                    candidate_seq,
                    seq
                )

                if rel == "exact":

                    comp_exact_seq = 1

                    if len(comp_seq_files) < 10:
                        comp_seq_files.append(
                            path
                        )

                elif rel == "contained":

                    comp_contained_seq = 1

                    if len(comp_seq_files) < 10:
                        comp_seq_files.append(
                            path
                        )


        anti_candidate = int(
            anti_exact
            or anti_contained
        )

        bagel_candidate = int(
            bagel_text_hit
            or bagel_exact_seq
            or bagel_contained_seq
        )

        comp_candidate = int(
            comp_text_hit
            or comp_exact_seq
            or comp_contained_seq
        )


        exact_method_families = (
            anti_candidate
            + bagel_candidate
            + comp_candidate
        )


        if exact_method_families >= 2:

            pattern = (
                "candidate_supported_by_multiple_specialized_methods"
            )

        elif (
            exact_method_families == 1
            and anti_context
            and not anti_candidate
        ):

            pattern = (
                "single_candidate_method_plus_antismash_context"
            )

        elif exact_method_families == 1:

            pattern = (
                "single_candidate_level_method"
            )

        elif anti_context:

            pattern = (
                "antismash_context_only_no_exact_candidate_match"
            )

        else:

            pattern = (
                "no_exact_candidate_level_support_detected"
            )


        candidate_rows.append({
            "locus_id":
                locus_id,

            "MAG":
                mag,

            "contig":
                contig,

            "candidate_id":
                candidate,

            "local_candidate_token":
                local_token,

            "candidate_sequence_available":
                int(
                    bool(
                        candidate_seq
                    )
                ),

            "candidate_aa_length":
                len(
                    candidate_seq
                )
                if candidate_seq
                else 0,

            "candidate_master_header":
                candidate_header,

            "master_sequence_matches":
                master_hit_count,

            "antismash_RiPP_context_same_contig":
                anti_context,

            "antismash_products":
                anti_products,

            "antismash_candidate_exact_translation":
                anti_exact,

            "antismash_candidate_contained_translation":
                anti_contained,

            "antismash_candidate_match_file":
                anti_match_file,

            "BAGEL_same_MAG_run_available":
                int(
                    bool(
                        actual_bagel_roots
                    )
                ),

            "BAGEL_same_MAG_runs":
                ";".join(
                    x[0]
                    for x in actual_bagel_roots
                ),

            "BAGEL_same_contig_files":
                bagel_same_contig_files,

            "BAGEL_candidate_text_hit":
                bagel_text_hit,

            "BAGEL_candidate_exact_sequence":
                bagel_exact_seq,

            "BAGEL_candidate_contained_sequence":
                bagel_contained_seq,

            "BAGEL_candidate_files":
                ";".join(
                    dict.fromkeys(
                        bagel_text_files
                        + bagel_seq_files
                    )
                ),

            "Comparippson_candidate_text_hit":
                comp_text_hit,

            "Comparippson_candidate_exact_sequence":
                comp_exact_seq,

            "Comparippson_candidate_contained_sequence":
                comp_contained_seq,

            "Comparippson_candidate_files":
                ";".join(
                    dict.fromkeys(
                        comp_text_files
                        + comp_seq_files
                    )
                ),

            "candidate_level_method_families":
                exact_method_families,

            "evidence_pattern":
                pattern,

            "robust_context_wide":
                locus.get(
                    "robust_context_wide",
                    ""
                ),

            "context_evidence_strict":
                locus.get(
                    "context_evidence_strict",
                    ""
                ),

            "left_edge_truncated":
                locus.get(
                    "left_edge_truncated",
                    ""
                ),

            "right_edge_truncated":
                locus.get(
                    "right_edge_truncated",
                    ""
                ),

            "contains_interrupted_candidate":
                locus.get(
                    "contains_interrupted_candidate",
                    ""
                ),
        })


# ============================================================
# Locus-level summary
# ============================================================

by_locus = defaultdict(list)

for row in candidate_rows:

    by_locus[
        row["locus_id"]
    ].append(row)


locus_rows = []


for locus in loci:

    lid = locus["locus_id"]
    rows = by_locus[lid]

    max_methods = max(
        int(
            r[
                "candidate_level_method_families"
            ]
        )
        for r in rows
    )

    n_multi = sum(
        int(
            r[
                "candidate_level_method_families"
            ]
        ) >= 2
        for r in rows
    )

    n_single = sum(
        int(
            r[
                "candidate_level_method_families"
            ]
        ) == 1
        for r in rows
    )

    n_zero = sum(
        int(
            r[
                "candidate_level_method_families"
            ]
        ) == 0
        for r in rows
    )

    anti_context = max(
        int(
            r[
                "antismash_RiPP_context_same_contig"
            ]
        )
        for r in rows
    )

    bagel_same_mag = max(
        int(
            r[
                "BAGEL_same_MAG_run_available"
            ]
        )
        for r in rows
    )


    if max_methods >= 2:

        locus_pattern = (
            "at_least_one_candidate_with_multimethod_support"
        )

    elif max_methods == 1 and anti_context:

        locus_pattern = (
            "candidate_level_single_method_with_antismash_context"
        )

    elif max_methods == 1:

        locus_pattern = (
            "candidate_level_single_method_only"
        )

    elif anti_context:

        locus_pattern = (
            "antismash_context_without_exact_candidate_match"
        )

    else:

        locus_pattern = (
            "no_exact_candidate_level_support_in_reviewed_sources"
        )


    locus_rows.append({
        "locus_id":
            lid,

        "MAG":
            locus["MAG"],

        "contig":
            locus["contig"],

        "candidate_count":
            len(rows),

        "candidate_sequences_available":
            sum(
                int(
                    r[
                        "candidate_sequence_available"
                    ]
                )
                for r in rows
            ),

        "candidates_multimethod":
            n_multi,

        "candidates_single_method":
            n_single,

        "candidates_no_exact_support":
            n_zero,

        "max_candidate_method_families":
            max_methods,

        "antismash_RiPP_context_same_contig":
            anti_context,

        "BAGEL_same_MAG_run_available":
            bagel_same_mag,

        "robust_context_wide":
            locus.get(
                "robust_context_wide",
                ""
            ),

        "context_evidence_strict":
            locus.get(
                "context_evidence_strict",
                ""
            ),

        "left_edge_truncated":
            locus.get(
                "left_edge_truncated",
                ""
            ),

        "right_edge_truncated":
            locus.get(
                "right_edge_truncated",
                ""
            ),

        "contains_interrupted_candidate":
            locus.get(
                "contains_interrupted_candidate",
                ""
            ),

        "locus_evidence_pattern":
            locus_pattern,
    })


# ============================================================
# Write tables
# ============================================================

candidate_fields = [
    "locus_id",
    "MAG",
    "contig",
    "candidate_id",
    "local_candidate_token",
    "candidate_sequence_available",
    "candidate_aa_length",
    "candidate_master_header",
    "master_sequence_matches",
    "antismash_RiPP_context_same_contig",
    "antismash_products",
    "antismash_candidate_exact_translation",
    "antismash_candidate_contained_translation",
    "antismash_candidate_match_file",
    "BAGEL_same_MAG_run_available",
    "BAGEL_same_MAG_runs",
    "BAGEL_same_contig_files",
    "BAGEL_candidate_text_hit",
    "BAGEL_candidate_exact_sequence",
    "BAGEL_candidate_contained_sequence",
    "BAGEL_candidate_files",
    "Comparippson_candidate_text_hit",
    "Comparippson_candidate_exact_sequence",
    "Comparippson_candidate_contained_sequence",
    "Comparippson_candidate_files",
    "candidate_level_method_families",
    "evidence_pattern",
    "robust_context_wide",
    "context_evidence_strict",
    "left_edge_truncated",
    "right_edge_truncated",
    "contains_interrupted_candidate",
]


with (
    OUT
    / "98B2_candidate_exact_evidence.tsv"
).open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=candidate_fields,
        delimiter="\t",
        lineterminator="\n"
    )

    w.writeheader()
    w.writerows(candidate_rows)


locus_fields = list(
    locus_rows[0].keys()
)


with (
    OUT
    / "98B2_locus_exact_summary.tsv"
).open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=locus_fields,
        delimiter="\t",
        lineterminator="\n"
    )

    w.writeheader()
    w.writerows(locus_rows)


# ============================================================
# Global summary
# ============================================================

with (
    OUT
    / "98B2_global_summary.tsv"
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

    w.writerow([
        "ATTRLOC_total",
        len(locus_rows)
    ])

    w.writerow([
        "candidate_hypotheses_total",
        len(candidate_rows)
    ])

    w.writerow([
        "candidate_sequences_available",
        sum(
            int(
                r[
                    "candidate_sequence_available"
                ]
            )
            for r in candidate_rows
        )
    ])

    w.writerow([
        "candidate_multimethod_support",
        sum(
            int(
                r[
                    "candidate_level_method_families"
                ]
            ) >= 2
            for r in candidate_rows
        )
    ])

    w.writerow([
        "candidate_single_method_support",
        sum(
            int(
                r[
                    "candidate_level_method_families"
                ]
            ) == 1
            for r in candidate_rows
        )
    ])

    w.writerow([
        "candidate_no_exact_support",
        sum(
            int(
                r[
                    "candidate_level_method_families"
                ]
            ) == 0
            for r in candidate_rows
        )
    ])

    w.writerow([
        "loci_with_multimethod_candidate_support",
        sum(
            int(
                r[
                    "max_candidate_method_families"
                ]
            ) >= 2
            for r in locus_rows
        )
    ])

    w.writerow([
        "loci_with_actual_BAGEL_same_MAG_run",
        sum(
            int(
                r[
                    "BAGEL_same_MAG_run_available"
                ]
            )
            for r in locus_rows
        )
    ])

    w.writerow([
        "new_prediction_performed",
        "NO"
    ])

    w.writerow([
        "reads_remapped",
        "NO"
    ])

    w.writerow([
        "amplicon_data_used",
        "NO"
    ])

    w.writerow([
        "next_step",
        "98C_incomplete_bin_positive_evidence_expansion"
    ])


print(
    f"candidate_hypotheses={len(candidate_rows)}"
)

print(
    "loci="
    + str(
        len(locus_rows)
    )
)

print(
    "98B2=PASS"
)
