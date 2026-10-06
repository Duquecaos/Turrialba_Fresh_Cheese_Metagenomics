#!/usr/bin/env python3

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(sys.argv[1])

COMM = ROOT / "96_virulence/94B1_VFDB_community"
FINAL18 = ROOT / "22_final_representative_mags/fastas"
AUDIT94A = ROOT / "96_virulence/94A_audit"

OUT = ROOT / "96_virulence/94B2_VFDB_consolidated"

OUT.mkdir(
    parents=True,
    exist_ok=True
)

GROUPS = [
    "L1",
    "L2",
    "L3",
    "M1",
    "M2",
    "M3",
]


# ======================================================================
# Helpers
# ======================================================================

def read_tsv(path):

    if not path.is_file():
        raise RuntimeError(
            f"Archivo faltante: {path}"
        )

    with path.open(errors="replace") as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        rows = list(reader)
        fields = reader.fieldnames or []

    return rows, fields


def write_tsv(path, rows, fields):

    with path.open(
        "w",
        newline=""
    ) as fh:

        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)


def number(value):

    try:
        return float(value)
    except Exception:
        return None


def integer(value):

    return int(float(value))


def unique_join(values):

    return ";".join(
        sorted(
            {
                str(x).strip()
                for x in values
                if str(x).strip()
            }
        )
    )


def parse_product(product):

    """
    VFDB PRODUCT suele tener:
      (gene) descripción [factor (VFxxxx)] [organismo referencia]

    No usamos el organismo referencia como taxonomía de la muestra.
    """

    groups = re.findall(
        r"\[([^\]]+)\]",
        product or ""
    )

    if len(groups) >= 2:
        factor = groups[-2]
        reference_organism = groups[-1]

    elif len(groups) == 1:
        factor = groups[0]
        reference_organism = ""

    else:
        factor = ""
        reference_organism = ""

    return factor, reference_organism


# ======================================================================
# 1. Consolidar 94B1
# ======================================================================

all_rows = []
base_fields = None

for group in GROUPS:

    path = (
        COMM
        / group
        / f"{group}_VFDB_screen80_80.tsv"
    )

    rows, fields = read_tsv(path)

    if base_fields is None:
        base_fields = fields

    required = {
        "community_vfdb_locus_id",
        "coassembly",
        "SEQUENCE",
        "START",
        "END",
        "STRAND",
        "GENE",
        "%COVERAGE",
        "%IDENTITY",
        "ACCESSION",
        "PRODUCT",
    }

    missing = required - set(fields)

    if missing:
        raise RuntimeError(
            f"{group}: faltan columnas "
            + ",".join(sorted(missing))
        )

    for row in rows:

        cov = number(
            row["%COVERAGE"]
        )

        ident = number(
            row["%IDENTITY"]
        )

        if cov is None or ident is None:
            raise RuntimeError(
                f"{row['community_vfdb_locus_id']}: "
                "cov/id no numérica"
            )

        if cov < 80 or ident < 80:
            raise RuntimeError(
                f"{row['community_vfdb_locus_id']}: "
                "hit <80/80"
            )

        factor, ref_org = parse_product(
            row.get(
                "PRODUCT",
                ""
            )
        )

        rec = dict(row)

        rec["producer"] = (
            group[0]
            if group
            else ""
        )

        rec["screen_class_recomputed"] = (
            "high_90_90"
            if cov >= 90 and ident >= 90
            else "screen_80_80"
        )

        rec["vf_factor_annotation"] = factor

        rec[
            "VFDB_reference_organism"
        ] = ref_org

        rec[
            "VFDB_reference_organism_is_sample_taxonomy"
        ] = "NO"

        all_rows.append(rec)


if len(all_rows) != 1075:

    raise RuntimeError(
        "Esperaba 1075 hits comunitarios 80/80; "
        f"observé {len(all_rows)}"
    )


# ======================================================================
# 2. Mapa exacto coassembly + original contig -> final18 MAG
# ======================================================================

mag_contig_map = defaultdict(list)

fasta_files = []

for pattern in [
    "*.fa",
    "*.fna",
    "*.fasta",
]:

    fasta_files.extend(
        FINAL18.glob(pattern)
    )

fasta_files = sorted(
    set(fasta_files)
)

if len(fasta_files) != 18:

    raise RuntimeError(
        f"Esperaba 18 FASTA final18; "
        f"observé {len(fasta_files)}"
    )


for fasta in fasta_files:

    mag = fasta.stem

    if "__" not in mag:

        raise RuntimeError(
            f"No puedo inferir coassembly: {mag}"
        )

    coassembly = mag.split(
        "__",
        1
    )[0]

    with fasta.open() as fh:

        for line in fh:

            if not line.startswith(">"):
                continue

            contig = (
                line[1:]
                .strip()
                .split()[0]
            )

            mag_contig_map[
                (
                    coassembly,
                    contig,
                )
            ].append(mag)


# ======================================================================
# 3. Añadir asignación exacta de contig a MAG
# ======================================================================

for row in all_rows:

    key = (
        row["coassembly"],
        row["SEQUENCE"],
    )

    mags = sorted(
        set(
            mag_contig_map.get(
                key,
                []
            )
        )
    )

    row["n_final18_MAGs_exact"] = len(mags)

    row["final18_MAGs_exact"] = (
        ";".join(mags)
    )

    row["assignment_level"] = (
        "exact_final18_contig"
        if mags
        else "community_contig_only"
    )


# ======================================================================
# 4. Exact-coordinate deduplication INSIDE each coassembly
#
# No collapse across L1/L2/etc.
# Cross-coassembly recurrence is NOT prevalence.
# ======================================================================

dedup_groups = defaultdict(list)

for row in all_rows:

    key = (
        row["coassembly"],
        row["SEQUENCE"],
        row["START"],
        row["END"],
        row["STRAND"],
        row["GENE"],
        row["ACCESSION"],
    )

    dedup_groups[key].append(row)


deduplicated = []
duplicate_groups = []

for key, rows in dedup_groups.items():

    rows = sorted(
        rows,
        key=lambda x: x[
            "community_vfdb_locus_id"
        ]
    )

    rec = dict(
        rows[0]
    )

    rec["exact_coordinate_duplicate_count"] = (
        len(rows)
    )

    rec["source_locus_ids"] = unique_join(
        x["community_vfdb_locus_id"]
        for x in rows
    )

    deduplicated.append(rec)

    if len(rows) > 1:

        duplicate_groups.append({
            "coassembly": key[0],
            "SEQUENCE": key[1],
            "START": key[2],
            "END": key[3],
            "STRAND": key[4],
            "GENE": key[5],
            "ACCESSION": key[6],
            "duplicate_count": len(rows),
            "source_locus_ids": unique_join(
                x["community_vfdb_locus_id"]
                for x in rows
            ),
        })


deduplicated.sort(
    key=lambda x: (
        x["coassembly"],
        x["SEQUENCE"],
        integer(x["START"]),
        integer(x["END"]),
        x["GENE"],
    )
)


# ======================================================================
# 5. High 90/90
# ======================================================================

high = [
    x
    for x in deduplicated
    if x[
        "screen_class_recomputed"
    ] == "high_90_90"
]


if len([
    x for x in all_rows
    if x["screen_class_recomputed"]
    == "high_90_90"
]) != 51:

    raise RuntimeError(
        "El conjunto RAW high90/90 no suma 51"
    )


high_assigned = [
    x
    for x in high
    if x["assignment_level"]
    == "exact_final18_contig"
]

high_community_only = [
    x
    for x in high
    if x["assignment_level"]
    == "community_contig_only"
]


# ======================================================================
# 6. Gene/reference summary high90/90
# ======================================================================

gene_groups = defaultdict(list)

for row in high:

    key = (
        row["GENE"],
        row["ACCESSION"],
        row[
            "vf_factor_annotation"
        ],
    )

    gene_groups[key].append(row)


gene_summary = []

for key, rows in gene_groups.items():

    identities = [
        number(x["%IDENTITY"])
        for x in rows
    ]

    coverages = [
        number(x["%COVERAGE"])
        for x in rows
    ]

    coassemblies = sorted(
        {
            x["coassembly"]
            for x in rows
        }
    )

    producers = sorted(
        {
            x["producer"]
            for x in rows
        }
    )

    contigs = {
        (
            x["coassembly"],
            x["SEQUENCE"],
        )
        for x in rows
    }

    exact_mags = sorted(
        {
            mag
            for x in rows
            for mag in x[
                "final18_MAGs_exact"
            ].split(";")
            if mag
        }
    )

    gene_summary.append({
        "GENE": key[0],
        "ACCESSION": key[1],
        "vf_factor_annotation": key[2],
        "n_hits_high90_90": len(rows),
        "n_contigs": len(contigs),
        "n_coassemblies": len(coassemblies),
        "coassemblies": ";".join(
            coassemblies
        ),
        "n_producers": len(producers),
        "producers": ";".join(
            producers
        ),
        "both_producers": (
            "YES"
            if len(producers) == 2
            else "NO"
        ),
        "min_percent_identity": min(
            identities
        ),
        "max_percent_identity": max(
            identities
        ),
        "min_percent_coverage": min(
            coverages
        ),
        "max_percent_coverage": max(
            coverages
        ),
        "n_exact_final18_MAGs": len(
            exact_mags
        ),
        "exact_final18_MAGs": ";".join(
            exact_mags
        ),
        "interpretation_guardrail":
            "VFDB_homology_not_pathogenic_phenotype",
    })


gene_summary.sort(
    key=lambda x: (
        -x["n_hits_high90_90"],
        x["GENE"],
        x["ACCESSION"],
    )
)


# ======================================================================
# 7. Resumen por coassembly
# ======================================================================

coassembly_summary = []

for group in GROUPS:

    raw_group = [
        x
        for x in all_rows
        if x["coassembly"] == group
    ]

    dedup_group = [
        x
        for x in deduplicated
        if x["coassembly"] == group
    ]

    high_group = [
        x
        for x in high
        if x["coassembly"] == group
    ]

    coassembly_summary.append({
        "coassembly": group,
        "producer": group[0],
        "raw_hits_80_80": len(
            raw_group
        ),
        "deduplicated_hits_80_80": len(
            dedup_group
        ),
        "high_hits_90_90": len(
            high_group
        ),
        "high_unique_gene_labels": len(
            {
                x["GENE"]
                for x in high_group
            }
        ),
        "high_contigs": len(
            {
                x["SEQUENCE"]
                for x in high_group
            }
        ),
        "high_hits_exact_final18_contig": sum(
            x["assignment_level"]
            == "exact_final18_contig"
            for x in high_group
        ),
        "high_hits_community_only": sum(
            x["assignment_level"]
            == "community_contig_only"
            for x in high_group
        ),
    })


# ======================================================================
# 8. Recalcular los 18 hits VFDB de final18
# ======================================================================

final18_file = (
    AUDIT94A
    / "94A_existing_final18_vfdb.tsv"
)

final18_rows, final18_fields = read_tsv(
    final18_file
)

if len(final18_rows) != 18:

    raise RuntimeError(
        f"VFDB final18 esperado=18 observado={len(final18_rows)}"
    )


final18_reclassified = []

for row in final18_rows:

    cov = number(
        row["%COVERAGE"]
    )

    ident = number(
        row["%IDENTITY"]
    )

    if cov is None or ident is None:
        raise RuntimeError(
            "final18 cov/id no numérica"
        )

    recomputed = (
        "high_90_90"
        if cov >= 90 and ident >= 90
        else "screen_80_80"
    )

    factor, ref_org = parse_product(
        row.get(
            "PRODUCT",
            ""
        )
    )

    rec = dict(row)

    rec[
        "screen_class_recomputed"
    ] = recomputed

    rec[
        "screen_class_matches_recomputed"
    ] = (
        "YES"
        if row.get(
            "screen_class",
            ""
        ) == recomputed
        else "NO"
    )

    rec[
        "vf_factor_annotation"
    ] = factor

    rec[
        "VFDB_reference_organism"
    ] = ref_org

    rec[
        "VFDB_reference_organism_is_sample_taxonomy"
    ] = "NO"

    final18_reclassified.append(
        rec
    )


final18_high = [
    x
    for x in final18_reclassified
    if x[
        "screen_class_recomputed"
    ] == "high_90_90"
]

final18_disagreement = sum(
    x[
        "screen_class_matches_recomputed"
    ] == "NO"
    for x in final18_reclassified
)


# ======================================================================
# 9. Write tables
# ======================================================================

derived_fields = [
    "producer",
    "screen_class_recomputed",
    "vf_factor_annotation",
    "VFDB_reference_organism",
    "VFDB_reference_organism_is_sample_taxonomy",
    "n_final18_MAGs_exact",
    "final18_MAGs_exact",
    "assignment_level",
]

all_fields = list(
    base_fields
) + derived_fields


write_tsv(
    OUT
    / "94B2_community_all80_80.tsv",
    all_rows,
    all_fields
)


dedup_fields = (
    all_fields
    + [
        "exact_coordinate_duplicate_count",
        "source_locus_ids",
    ]
)

write_tsv(
    OUT
    / "94B2_community_deduplicated80_80.tsv",
    deduplicated,
    dedup_fields
)

write_tsv(
    OUT
    / "94B2_community_high90_90.tsv",
    high,
    dedup_fields
)

write_tsv(
    OUT
    / "94B2_high90_90_exact_final18_contig.tsv",
    high_assigned,
    dedup_fields
)

write_tsv(
    OUT
    / "94B2_high90_90_community_only.tsv",
    high_community_only,
    dedup_fields
)


duplicate_fields = [
    "coassembly",
    "SEQUENCE",
    "START",
    "END",
    "STRAND",
    "GENE",
    "ACCESSION",
    "duplicate_count",
    "source_locus_ids",
]

write_tsv(
    OUT
    / "94B2_exact_coordinate_duplicate_groups.tsv",
    duplicate_groups,
    duplicate_fields
)


gene_fields = [
    "GENE",
    "ACCESSION",
    "vf_factor_annotation",
    "n_hits_high90_90",
    "n_contigs",
    "n_coassemblies",
    "coassemblies",
    "n_producers",
    "producers",
    "both_producers",
    "min_percent_identity",
    "max_percent_identity",
    "min_percent_coverage",
    "max_percent_coverage",
    "n_exact_final18_MAGs",
    "exact_final18_MAGs",
    "interpretation_guardrail",
]

write_tsv(
    OUT
    / "94B2_high90_90_gene_summary.tsv",
    gene_summary,
    gene_fields
)


coassembly_fields = list(
    coassembly_summary[0].keys()
)

write_tsv(
    OUT
    / "94B2_summary_by_coassembly.tsv",
    coassembly_summary,
    coassembly_fields
)


final18_extra = [
    "screen_class_recomputed",
    "screen_class_matches_recomputed",
    "vf_factor_annotation",
    "VFDB_reference_organism",
    "VFDB_reference_organism_is_sample_taxonomy",
]

write_tsv(
    OUT
    / "94B2_final18_reclassified.tsv",
    final18_reclassified,
    final18_fields
    + final18_extra
)


# ======================================================================
# 10. Summaries
# ======================================================================

high_raw_count = sum(
    x["screen_class_recomputed"]
    == "high_90_90"
    for x in all_rows
)

global_metrics = [
    (
        "community_raw_hits_80_80",
        len(all_rows)
    ),
    (
        "community_deduplicated_hits_80_80",
        len(deduplicated)
    ),
    (
        "exact_coordinate_duplicate_groups",
        len(duplicate_groups)
    ),
    (
        "community_raw_hits_high_90_90",
        high_raw_count
    ),
    (
        "community_deduplicated_hits_high_90_90",
        len(high)
    ),
    (
        "fraction_raw_80_80_retained_high_90_90",
        high_raw_count
        / len(all_rows)
    ),
    (
        "unique_gene_labels_high_90_90",
        len(
            {
                x["GENE"]
                for x in high
            }
        )
    ),
    (
        "unique_VFDB_accessions_high_90_90",
        len(
            {
                x["ACCESSION"]
                for x in high
            }
        )
    ),
    (
        "high90_90_hits_exact_final18_contig",
        len(high_assigned)
    ),
    (
        "high90_90_hits_community_only",
        len(high_community_only)
    ),
    (
        "high90_90_gene_reference_signatures_both_producers",
        sum(
            x["both_producers"]
            == "YES"
            for x in gene_summary
        )
    ),
    (
        "final18_VFDB_hits_total",
        len(final18_reclassified)
    ),
    (
        "final18_VFDB_hits_high_90_90",
        len(final18_high)
    ),
    (
        "final18_VFDB_hits_screen80_not90",
        len(final18_reclassified)
        - len(final18_high)
    ),
    (
        "final18_screen_class_disagreements",
        final18_disagreement
    ),
]


with (
    OUT
    / "94B2_global_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for key, value in global_metrics:
        fh.write(
            f"{key}\t{value}\n"
        )


with (
    OUT
    / "94B2_final18_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    fh.write(
        f"VFDB_hits_total\t"
        f"{len(final18_reclassified)}\n"
    )

    fh.write(
        f"VFDB_hits_high_90_90\t"
        f"{len(final18_high)}\n"
    )

    fh.write(
        f"VFDB_hits_80_80_but_not_90_90\t"
        f"{len(final18_reclassified)-len(final18_high)}\n"
    )

    fh.write(
        f"existing_vs_recomputed_class_disagreements\t"
        f"{final18_disagreement}\n"
    )


# ======================================================================
# 11. Methodological scope
# ======================================================================

with (
    OUT
    / "94B2_methodological_scope.tsv"
).open("w") as fh:

    fh.write(
        "field\tvalue\n"
    )

    fh.write(
        "database\tVFDB_2024-Dec-15\n"
    )

    fh.write(
        "broad_screen\tidentity>=80_and_coverage>=80\n"
    )

    fh.write(
        "high_similarity_screen\tidentity>=90_and_coverage>=90\n"
    )

    fh.write(
        "deduplication_scope\twithin_same_coassembly_only\n"
    )

    fh.write(
        "exact_duplicate_definition\t"
        "same_coassembly_contig_start_end_strand_gene_accession\n"
    )

    fh.write(
        "cross_coassembly_hits_collapsed\tNO\n"
    )

    fh.write(
        "cross_coassembly_recurrence_equals_prevalence\tNO\n"
    )

    fh.write(
        "MAG_assignment_method\t"
        "exact_original_contig_membership_in_final18_FASTA\n"
    )

    fh.write(
        "VFDB_reference_organism_used_as_sample_taxonomy\tNO\n"
    )

    fh.write(
        "VFDB_hit_equals_pathogenic_phenotype\tNO\n"
    )

    fh.write(
        "VFDB_hit_equals_pathogenic_organism\tNO\n"
    )


# ======================================================================
# 12. Validation
# ======================================================================

required_outputs = [
    "94B2_community_all80_80.tsv",
    "94B2_community_deduplicated80_80.tsv",
    "94B2_community_high90_90.tsv",
    "94B2_high90_90_exact_final18_contig.tsv",
    "94B2_high90_90_community_only.tsv",
    "94B2_exact_coordinate_duplicate_groups.tsv",
    "94B2_high90_90_gene_summary.tsv",
    "94B2_summary_by_coassembly.tsv",
    "94B2_final18_reclassified.tsv",
    "94B2_final18_summary.tsv",
    "94B2_global_summary.tsv",
    "94B2_methodological_scope.tsv",
]

for name in required_outputs:

    path = OUT / name

    if not path.is_file():

        raise RuntimeError(
            f"Salida faltante: {path}"
        )


print(
    f"COMMUNITY_80_80={len(all_rows)}"
)

print(
    f"COMMUNITY_HIGH_90_90={high_raw_count}"
)

print(
    f"FINAL18_TOTAL={len(final18_reclassified)}"
)

print(
    f"FINAL18_HIGH_90_90={len(final18_high)}"
)

print(
    "94B2_CONSOLIDATION=PASS"
)

