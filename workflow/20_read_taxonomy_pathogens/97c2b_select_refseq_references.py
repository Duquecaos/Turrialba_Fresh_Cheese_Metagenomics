#!/usr/bin/env python3

import csv
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path


if len(sys.argv) != 6:
    raise SystemExit(
        "Usage: 97c2b_select_refseq_references.py "
        "<assembly_summary> "
        "<requested_species.tsv> "
        "<selected.tsv> "
        "<missing.tsv> "
        "<audit.tsv>"
    )


ASSEMBLY = Path(sys.argv[1])
REQUESTS = Path(sys.argv[2])
SELECTED = Path(sys.argv[3])
MISSING = Path(sys.argv[4])
AUDIT = Path(sys.argv[5])


# ============================================================
# Requests
# ============================================================

with REQUESTS.open(
    newline=""
) as fh:

    requests = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )


if not requests:
    raise RuntimeError(
        "Requested-species table is empty"
    )


required_request_fields = {
    "panel",
    "role",
    "species",
}


if not required_request_fields.issubset(
    requests[0].keys()
):
    raise RuntimeError(
        "Requested-species TSV lacks required columns"
    )


requests_by_species = defaultdict(list)

for r in requests:

    species = r["species"].strip()

    requests_by_species[
        species
    ].append(r)


requested_species = set(
    requests_by_species
)


# ============================================================
# Read RefSeq assembly summary
# ============================================================

header = None
idx = None

required_columns = [
    "assembly_accession",
    "refseq_category",
    "taxid",
    "species_taxid",
    "organism_name",
    "version_status",
    "assembly_level",
    "genome_rep",
    "ftp_path",
]


candidates = defaultdict(list)

field_counts = Counter()

data_rows = 0
rows_shorter_than_required = 0
latest_complete_full_rows = 0
matching_rows = 0


with ASSEMBLY.open(
    errors="replace"
) as fh:

    for line in fh:

        # ----------------------------------------------------
        # Header
        # ----------------------------------------------------

        if line.startswith("#"):

            h = line.lstrip(
                "# "
            ).rstrip(
                "\r\n"
            )

            if h.startswith(
                "assembly_accession\t"
            ):

                header = h.split(
                    "\t"
                )

                idx = {
                    name: i
                    for i, name
                    in enumerate(header)
                }

                missing_columns = [
                    x
                    for x in required_columns
                    if x not in idx
                ]

                if missing_columns:

                    raise RuntimeError(
                        "Missing assembly-summary columns: "
                        + ", ".join(
                            missing_columns
                        )
                    )

            continue


        if not line.strip():
            continue


        if header is None:
            continue


        parts = line.rstrip(
            "\r\n"
        ).split(
            "\t"
        )


        data_rows += 1

        field_counts[
            len(parts)
        ] += 1


        max_required_index = max(
            idx[x]
            for x in required_columns
        )


        if len(parts) <= max_required_index:

            rows_shorter_than_required += 1
            continue


        def get(name):

            i = idx.get(name)

            if i is None:
                return ""

            if i >= len(parts):
                return ""

            return parts[i].strip()


        # ----------------------------------------------------
        # Only useful RefSeq records
        # ----------------------------------------------------

        if get(
            "version_status"
        ) != "latest":

            continue


        if get(
            "assembly_level"
        ) != "Complete Genome":

            continue


        if get(
            "genome_rep"
        ) != "Full":

            continue


        ftp_path = get(
            "ftp_path"
        )


        if ftp_path in (
            "",
            "na",
        ):

            continue


        latest_complete_full_rows += 1


        organism = get(
            "organism_name"
        )


        # Requested entries are species binomials.
        words = organism.split()

        if len(words) < 2:
            continue


        species_key = (
            words[0]
            + " "
            + words[1]
        )


        if species_key not in requested_species:
            continue


        # Extra safeguard
        if not (
            organism == species_key
            or organism.startswith(
                species_key + " "
            )
        ):
            continue


        matching_rows += 1


        seq_rel_date = ""

        if "seq_rel_date" in idx:

            seq_rel_date = get(
                "seq_rel_date"
            )

        elif "release_date" in idx:

            seq_rel_date = get(
                "release_date"
            )


        candidates[
            species_key
        ].append({
            "assembly_accession":
                get(
                    "assembly_accession"
                ),

            "refseq_category":
                get(
                    "refseq_category"
                ),

            "organism_name":
                organism,

            "taxid":
                get(
                    "taxid"
                ),

            "species_taxid":
                get(
                    "species_taxid"
                ),

            "assembly_level":
                get(
                    "assembly_level"
                ),

            "release_date":
                seq_rel_date,

            "ftp_path":
                ftp_path,
        })


if header is None:

    raise RuntimeError(
        "Could not locate '# assembly_accession' header"
    )


# ============================================================
# Selection functions
# ============================================================

def category_rank(value):

    value = value.strip().lower()

    if value == "reference genome":
        return 0

    if value == "representative genome":
        return 1

    return 2


def parse_date(value):

    for fmt in (
        "%Y/%m/%d",
        "%Y-%m-%d",
    ):

        try:

            return datetime.strptime(
                value,
                fmt
            ).timestamp()

        except Exception:
            pass

    return 0.0


def normalize_ncbi_url(ftp_path):

    x = ftp_path.rstrip("/")

    if x.startswith(
        "ftp://ftp.ncbi.nlm.nih.gov/"
    ):

        x = (
            "https://ftp.ncbi.nlm.nih.gov/"
            + x.split(
                "ftp://ftp.ncbi.nlm.nih.gov/",
                1
            )[1]
        )

    return x


# ============================================================
# Select one reference per requested species
# ============================================================

selected_rows = []
missing_rows = []
candidate_count_rows = []


for req in requests:

    panel = req[
        "panel"
    ]

    role = req[
        "role"
    ]

    species = req[
        "species"
    ].strip()


    cc = candidates.get(
        species,
        []
    )


    candidate_count_rows.append({
        "panel":
            panel,

        "role":
            role,

        "species":
            species,

        "complete_refseq_candidates":
            len(cc),
    })


    if not cc:

        missing_rows.append({
            "panel":
                panel,

            "role":
                role,

            "species":
                species,

            "status":
                "NO_COMPLETE_REFSEQ_ASSEMBLY_FOUND",
        })

        continue


    cc = sorted(
        cc,
        key=lambda r: (
            category_rank(
                r[
                    "refseq_category"
                ]
            ),
            -parse_date(
                r[
                    "release_date"
                ]
            ),
            r[
                "assembly_accession"
            ],
        )
    )


    best = cc[0]


    ftp_path = normalize_ncbi_url(
        best[
            "ftp_path"
        ]
    )


    basename = ftp_path.split(
        "/"
    )[-1]


    download_url = (
        ftp_path
        + "/"
        + basename
        + "_genomic.fna.gz"
    )


    # IMPORTANT:
    # Keep EXACTLY the original 12-column schema.
    # The downstream Bash read command expects 12 fields.
    selected_rows.append({
        "panel":
            panel,

        "role":
            role,

        "requested_species":
            species,

        "assembly_accession":
            best[
                "assembly_accession"
            ],

        "refseq_category":
            best[
                "refseq_category"
            ],

        "organism_name":
            best[
                "organism_name"
            ],

        "taxid":
            best[
                "taxid"
            ],

        "species_taxid":
            best[
                "species_taxid"
            ],

        "assembly_level":
            best[
                "assembly_level"
            ],

        "release_date":
            best[
                "release_date"
            ],

        "ftp_path":
            ftp_path,

        "download_url":
            download_url,
    })


# ============================================================
# Validate all nine target species
# ============================================================

target_requested = {
    (
        r["panel"],
        r["species"]
    )
    for r in requests
    if r["role"] == "target"
}


target_found = {
    (
        r["panel"],
        r["requested_species"]
    )
    for r in selected_rows
    if r["role"] == "target"
}


missing_targets = sorted(
    target_requested
    - target_found
)


# ============================================================
# Write selected references
# ============================================================

selected_fields = [
    "panel",
    "role",
    "requested_species",
    "assembly_accession",
    "refseq_category",
    "organism_name",
    "taxid",
    "species_taxid",
    "assembly_level",
    "release_date",
    "ftp_path",
    "download_url",
]


with SELECTED.open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=selected_fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(
        selected_rows
    )


# ============================================================
# Write missing
# ============================================================

missing_fields = [
    "panel",
    "role",
    "species",
    "status",
]


with MISSING.open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=missing_fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(
        missing_rows
    )


# ============================================================
# Write audit
# ============================================================

with AUDIT.open(
    "w",
    newline=""
) as fh:

    writer = csv.writer(
        fh,
        delimiter="\t",
        lineterminator="\n"
    )


    writer.writerow([
        "metric",
        "value"
    ])


    audit_rows = [
        (
            "assembly_summary_header_columns",
            len(header)
        ),
        (
            "assembly_summary_data_rows",
            data_rows
        ),
        (
            "rows_shorter_than_required",
            rows_shorter_than_required
        ),
        (
            "latest_complete_full_rows",
            latest_complete_full_rows
        ),
        (
            "matching_requested_species_rows",
            matching_rows
        ),
        (
            "requested_species_rows",
            len(requests)
        ),
        (
            "selected_reference_rows",
            len(selected_rows)
        ),
        (
            "missing_reference_rows",
            len(missing_rows)
        ),
        (
            "target_species_requested",
            len(target_requested)
        ),
        (
            "target_species_found",
            len(target_found)
        ),
        (
            "missing_target_species",
            len(missing_targets)
        ),
    ]


    writer.writerows(
        audit_rows
    )


    for n_fields, n_rows in sorted(
        field_counts.items()
    ):

        writer.writerow([
            f"rows_with_{n_fields}_fields",
            n_rows
        ])


# Separate useful candidate-count table
candidate_audit = (
    AUDIT.parent
    / "97C2B_candidate_counts.tsv"
)


with candidate_audit.open(
    "w",
    newline=""
) as fh:

    fields = [
        "panel",
        "role",
        "species",
        "complete_refseq_candidates",
    ]

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(
        candidate_count_rows
    )


# ============================================================
# Final validation
# ============================================================

if missing_targets:

    raise RuntimeError(
        "Required target species lack a complete RefSeq genome: "
        + "; ".join(
            f"{panel}:{species}"
            for panel, species
            in missing_targets
        )
    )


if len(
    selected_rows
) == 0:

    raise RuntimeError(
        "No references selected"
    )


print(
    f"HEADER_COLUMNS={len(header)}"
)

print(
    f"DATA_ROWS={data_rows}"
)

print(
    f"LATEST_COMPLETE_FULL={latest_complete_full_rows}"
)

print(
    f"REQUESTED={len(requests)}"
)

print(
    f"SELECTED={len(selected_rows)}"
)

print(
    f"MISSING_OPTIONAL={len(missing_rows)}"
)

print(
    f"TARGETS_FOUND={len(target_found)}"
)

print(
    "97C2B_REFSEQ_SELECTION=PASS"
)
