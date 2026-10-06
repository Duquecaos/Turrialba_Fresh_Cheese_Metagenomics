#!/bin/bash

#SBATCH --job-name=92d_viralcat
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=01:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=92d_viralcat_%j.out
#SBATCH --error=92d_viralcat_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

SRC="${ROOT}/85_genomad_coassemblies"
OUT="${ROOT}/86_viral_catalog_preQC"

mkdir -p "$OUT"

echo "======================================================================"
echo "PASO 92D - CONSOLIDACION Y PRE-QC DEL CATALOGO VIRAL"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "======================================================================"

srun --exact -n 1 -c 1 python3 - <<'PY'

import csv
import hashlib
import os
import statistics
from pathlib import Path
from collections import Counter, defaultdict

root = Path(
    f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"
)

src = root / "85_genomad_coassemblies"
out = root / "86_viral_catalog_preQC"

groups = [
    "L1",
    "L2",
    "L3",
    "M1",
    "M2",
    "M3",
]


# ============================================================
# Helpers
# ============================================================

def clean(x):
    if x is None:
        return ""
    x = str(x).strip()

    if x.upper() in {"NA", "NAN", "NONE"}:
        return ""

    return x


def fnum(x):
    x = clean(x)

    if not x:
        return None

    try:
        return float(x)
    except ValueError:
        return None


def inum(x):
    x = clean(x)

    if not x:
        return None

    return int(float(x))


def read_fasta(path):
    seqs = {}

    name = None
    chunks = []

    with path.open() as fh:

        for line in fh:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if name is not None:
                    if name in seqs:
                        raise RuntimeError(
                            f"FASTA ID duplicado: {name}"
                        )

                    seqs[name] = "".join(
                        chunks
                    ).upper()

                name = line[1:].split()[0]
                chunks = []

            else:
                chunks.append(line)

        if name is not None:

            if name in seqs:
                raise RuntimeError(
                    f"FASTA ID duplicado: {name}"
                )

            seqs[name] = "".join(
                chunks
            ).upper()

    return seqs


def write_fasta_record(fh, name, seq, width=80):

    fh.write(f">{name}\n")

    for i in range(
        0,
        len(seq),
        width
    ):
        fh.write(
            seq[i:i+width] + "\n"
        )


def sequence_sha256(seq):

    return hashlib.sha256(
        seq.encode("ascii")
    ).hexdigest()


def terminal_taxon(taxonomy):

    taxonomy = clean(taxonomy)

    if not taxonomy:
        return ""

    levels = [
        x.strip()
        for x in taxonomy.split(";")
        if x.strip()
    ]

    if not levels:
        return ""

    return levels[-1]


# ============================================================
# 1. Leer resultados por coensamblaje
# ============================================================

catalog = []
sequence_by_id = {}

expected_total = 0

for group in groups:

    summary_file = (
        src
        / group
        / f"{group}_summary"
        / f"{group}_virus_summary.tsv"
    )

    fasta_file = (
        src
        / group
        / f"{group}_summary"
        / f"{group}_virus.fna"
    )

    contig_map_file = (
        src
        / "contig_maps"
        / f"{group}.contig_map.tsv"
    )

    for f in [
        summary_file,
        fasta_file,
        contig_map_file,
    ]:

        if not f.exists():
            raise RuntimeError(
                f"Falta archivo requerido: {f}"
            )

    # --------------------------------------------------------
    # Mapa prefixed -> original
    # --------------------------------------------------------

    pref_to_original = {}

    with contig_map_file.open() as fh:

        r = csv.DictReader(
            fh,
            delimiter="\t"
        )

        for row in r:

            pref_to_original[
                row["prefixed_contig"]
            ] = row["original_contig"]

    # --------------------------------------------------------
    # FASTA viral
    # --------------------------------------------------------

    seqs = read_fasta(
        fasta_file
    )

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    with summary_file.open() as fh:

        r = csv.DictReader(
            fh,
            delimiter="\t"
        )

        rows = list(r)

    expected_total += len(rows)

    if len(rows) != len(seqs):

        raise RuntimeError(
            f"{group}: summary={len(rows)} "
            f"pero FASTA={len(seqs)}"
        )

    summary_ids = {
        row["seq_name"]
        for row in rows
    }

    fasta_ids = set(
        seqs.keys()
    )

    if summary_ids != fasta_ids:

        missing_fasta = (
            summary_ids - fasta_ids
        )

        missing_summary = (
            fasta_ids - summary_ids
        )

        raise RuntimeError(
            f"{group}: IDs summary/FASTA no coinciden. "
            f"missing_fasta={len(missing_fasta)} "
            f"missing_summary={len(missing_summary)}"
        )

    for row in rows:

        seq_name = row["seq_name"]

        # geNomad:
        # parent|provirus_START_END
        parent_prefixed = (
            seq_name.split(
                "|provirus_",
                1
            )[0]
        )

        if parent_prefixed not in pref_to_original:

            raise RuntimeError(
                f"{group}: no se encontró contig padre "
                f"{parent_prefixed} en mapa"
            )

        seq = seqs[seq_name]

        reported_length = inum(
            row.get("length")
        )

        if (
            reported_length is not None
            and reported_length != len(seq)
        ):
            raise RuntimeError(
                f"{seq_name}: geNomad length="
                f"{reported_length}, FASTA={len(seq)}"
            )

        taxonomy = clean(
            row.get("taxonomy")
        )

        score = fnum(
            row.get("virus_score")
        )

        fdr = fnum(
            row.get("fdr")
        )

        hallmarks = inum(
            row.get("n_hallmarks")
        )

        topology = clean(
            row.get("topology")
        )

        catalog.append({
            "seq_name":
                seq_name,

            "coassembly":
                group,

            "producer":
                group[0],

            "biological_unit":
                group[1:],

            "parent_contig_prefixed":
                parent_prefixed,

            "parent_contig_original":
                pref_to_original[
                    parent_prefixed
                ],

            "topology":
                topology,

            "coordinates":
                clean(
                    row.get("coordinates")
                ),

            "length":
                len(seq),

            "n_genes":
                clean(
                    row.get("n_genes")
                ),

            "genetic_code":
                clean(
                    row.get("genetic_code")
                ),

            "virus_score":
                (
                    score
                    if score is not None
                    else ""
                ),

            "fdr":
                (
                    fdr
                    if fdr is not None
                    else ""
                ),

            "n_hallmarks":
                (
                    hallmarks
                    if hallmarks is not None
                    else ""
                ),

            "marker_enrichment":
                clean(
                    row.get(
                        "marker_enrichment"
                    )
                ),

            "taxonomy":
                taxonomy,

            "terminal_taxon":
                terminal_taxon(
                    taxonomy
                ),

            "taxonomy_beyond_Viruses":
                (
                    "YES"
                    if terminal_taxon(
                        taxonomy
                    )
                    not in {"", "Viruses"}
                    else "NO"
                ),

            "caudoviricetes_in_taxonomy":
                (
                    "YES"
                    if "Caudoviricetes"
                    in taxonomy
                    else "NO"
                ),

            "viral_context":
                (
                    "PROVIRUS"
                    if topology
                    == "Provirus"
                    else
                    "NON_INTEGRATED_PREDICTION"
                ),

            "sequence_sha256":
                sequence_sha256(
                    seq
                ),
        })

        if seq_name in sequence_by_id:
            raise RuntimeError(
                f"ID viral global duplicado: "
                f"{seq_name}"
            )

        sequence_by_id[
            seq_name
        ] = seq


# ============================================================
# 2. QC global esperado
# ============================================================

if expected_total != 6888:
    raise RuntimeError(
        f"Se esperaban 6888 predicciones; "
        f"se encontraron {expected_total}"
    )

if len(catalog) != 6888:
    raise RuntimeError(
        f"Catálogo inesperado: {len(catalog)}"
    )

if len(sequence_by_id) != 6888:
    raise RuntimeError(
        f"FASTA global inesperado: "
        f"{len(sequence_by_id)}"
    )


# ============================================================
# 3. Escribir catálogo crudo completo
# ============================================================

catalog.sort(
    key=lambda x: (
        x["coassembly"],
        x["seq_name"],
    )
)

catalog_file = (
    out
    / "92D_viral_catalog_all.tsv"
)

catalog_fields = list(
    catalog[0].keys()
)

with catalog_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=catalog_fields,
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(catalog)


# ============================================================
# 4. FASTA completo
# ============================================================

all_fasta = (
    out
    / "92D_viral_sequences_all.fna"
)

with all_fasta.open("w") as fh:

    for row in catalog:

        write_fasta_record(
            fh,
            row["seq_name"],
            sequence_by_id[
                row["seq_name"]
            ]
        )


# ============================================================
# 5. Clustering SOLO por secuencia exacta
# ============================================================

by_hash = defaultdict(list)

for row in catalog:

    by_hash[
        row["sequence_sha256"]
    ].append(row)


hashes = sorted(
    by_hash.keys()
)

cluster_rows = []
member_rows = []

rep_fasta = (
    out
    / "92D_exact_representatives.fna"
)

with rep_fasta.open("w") as repfh:

    for idx, h in enumerate(
        hashes,
        start=1
    ):

        cluster_id = (
            f"VEXACT{idx:06d}"
        )

        members = sorted(
            by_hash[h],
            key=lambda x: (
                x["coassembly"],
                x["seq_name"],
            )
        )

        # Hash idéntico debe implicar
        # secuencia idéntica.
        first_seq = sequence_by_id[
            members[0]["seq_name"]
        ]

        for m in members[1:]:

            if (
                sequence_by_id[
                    m["seq_name"]
                ]
                != first_seq
            ):
                raise RuntimeError(
                    "Colisión SHA256 detectada"
                )

        representative = (
            members[0]["seq_name"]
        )

        write_fasta_record(
            repfh,
            cluster_id,
            first_seq
        )

        coassemblies = sorted({
            m["coassembly"]
            for m in members
        })

        producers = sorted({
            m["producer"]
            for m in members
        })

        topologies = sorted({
            m["topology"]
            for m in members
        })

        terminal_taxa = sorted({
            m["terminal_taxon"]
            for m in members
            if m["terminal_taxon"]
        })

        scores = [
            float(m["virus_score"])
            for m in members
            if m["virus_score"] != ""
        ]

        fdrs = [
            float(m["fdr"])
            for m in members
            if m["fdr"] != ""
        ]

        cluster_rows.append({
            "exact_cluster_id":
                cluster_id,

            "sequence_sha256":
                h,

            "length":
                len(first_seq),

            "n_members":
                len(members),

            "representative_original_id":
                representative,

            "n_coassemblies":
                len(coassemblies),

            "coassemblies":
                ";".join(
                    coassemblies
                ),

            "n_producers":
                len(producers),

            "producers":
                ";".join(
                    producers
                ),

            "topologies":
                ";".join(
                    topologies
                ),

            "terminal_taxa":
                ";".join(
                    terminal_taxa
                ),

            "max_virus_score":
                (
                    max(scores)
                    if scores
                    else ""
                ),

            "min_fdr":
                (
                    min(fdrs)
                    if fdrs
                    else ""
                ),

            "max_n_hallmarks":
                max(
                    [
                        int(m["n_hallmarks"])
                        for m in members
                        if m["n_hallmarks"]
                        != ""
                    ]
                    or [0]
                ),
        })

        for m in members:

            member_rows.append({
                "exact_cluster_id":
                    cluster_id,

                "seq_name":
                    m["seq_name"],

                "coassembly":
                    m["coassembly"],

                "producer":
                    m["producer"],

                "biological_unit":
                    m[
                        "biological_unit"
                    ],

                "parent_contig_original":
                    m[
                        "parent_contig_original"
                    ],

                "topology":
                    m["topology"],

                "coordinates":
                    m["coordinates"],

                "length":
                    m["length"],

                "virus_score":
                    m["virus_score"],

                "fdr":
                    m["fdr"],

                "n_hallmarks":
                    m["n_hallmarks"],

                "taxonomy":
                    m["taxonomy"],

                "is_exact_representative":
                    (
                        "YES"
                        if m["seq_name"]
                        == representative
                        else "NO"
                    ),
            })


cluster_file = (
    out
    / "92D_exact_sequence_clusters.tsv"
)

with cluster_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            cluster_rows[0].keys()
        ),
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(cluster_rows)


members_file = (
    out
    / "92D_exact_sequence_members.tsv"
)

with members_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            member_rows[0].keys()
        ),
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(member_rows)


# ============================================================
# 6. Resumen por coensamblaje
# ============================================================

by_group = defaultdict(list)

for row in catalog:

    by_group[
        row["coassembly"]
    ].append(row)


group_rows = []

for group in groups:

    rows = by_group[group]

    scores = [
        float(x["virus_score"])
        for x in rows
        if x["virus_score"] != ""
    ]

    fdrs = [
        float(x["fdr"])
        for x in rows
        if x["fdr"] != ""
    ]

    group_rows.append({
        "coassembly":
            group,

        "viral_predictions":
            len(rows),

        "proviruses":
            sum(
                x["topology"]
                == "Provirus"
                for x in rows
            ),

        "DTR":
            sum(
                x["topology"]
                == "DTR"
                for x in rows
            ),

        "ITR":
            sum(
                x["topology"]
                == "ITR"
                for x in rows
            ),

        "no_terminal_repeats":
            sum(
                x["topology"]
                == "No terminal repeats"
                for x in rows
            ),

        "total_viral_bp":
            sum(
                x["length"]
                for x in rows
            ),

        "median_length":
            statistics.median(
                x["length"]
                for x in rows
            ),

        "max_length":
            max(
                x["length"]
                for x in rows
            ),

        "score_ge_0.90":
            sum(
                float(x["virus_score"])
                >= 0.90
                for x in rows
                if x["virus_score"] != ""
            ),

        "score_ge_0.95":
            sum(
                float(x["virus_score"])
                >= 0.95
                for x in rows
                if x["virus_score"] != ""
            ),

        "score_ge_0.99":
            sum(
                float(x["virus_score"])
                >= 0.99
                for x in rows
                if x["virus_score"] != ""
            ),

        "fdr_nonNA":
            len(fdrs),

        "fdr_le_0.05":
            sum(
                x <= 0.05
                for x in fdrs
            ),

        "fdr_le_0.01":
            sum(
                x <= 0.01
                for x in fdrs
            ),

        "hallmark_ge_1":
            sum(
                (
                    x["n_hallmarks"]
                    != ""
                    and int(
                        x["n_hallmarks"]
                    ) >= 1
                )
                for x in rows
            ),

        "hallmark_ge_2":
            sum(
                (
                    x["n_hallmarks"]
                    != ""
                    and int(
                        x["n_hallmarks"]
                    ) >= 2
                )
                for x in rows
            ),

        "taxonomy_beyond_Viruses":
            sum(
                x[
                    "taxonomy_beyond_Viruses"
                ]
                == "YES"
                for x in rows
            ),

        "caudoviricetes_taxonomy":
            sum(
                x[
                    "caudoviricetes_in_taxonomy"
                ]
                == "YES"
                for x in rows
            ),
    })


group_file = (
    out
    / "92D_summary_by_coassembly.tsv"
)

with group_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            group_rows[0].keys()
        ),
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(group_rows)


# ============================================================
# 7. Distribución de longitudes
# ============================================================

length_bins = Counter()

for row in catalog:

    L = row["length"]

    if L < 2500:
        key = "lt_2.5kb"

    elif L < 5000:
        key = "2.5_to_5kb"

    elif L < 10000:
        key = "5_to_10kb"

    elif L < 20000:
        key = "10_to_20kb"

    elif L < 50000:
        key = "20_to_50kb"

    else:
        key = "ge_50kb"

    length_bins[key] += 1


# ============================================================
# 8. Terminal taxonomy
# ============================================================

terminal_counts = Counter(
    row["terminal_taxon"]
    if row["terminal_taxon"]
    else "UNASSIGNED"
    for row in catalog
)

taxonomy_file = (
    out
    / "92D_terminal_taxonomy_counts.tsv"
)

with taxonomy_file.open("w") as fh:

    fh.write(
        "terminal_taxon\tcount\n"
    )

    for taxon, n in sorted(
        terminal_counts.items(),
        key=lambda x: (
            -x[1],
            x[0],
        )
    ):
        fh.write(
            f"{taxon}\t{n}\n"
        )


# ============================================================
# 9. Resumen global
# ============================================================

n_clusters = len(
    cluster_rows
)

singleton_clusters = sum(
    x["n_members"] == 1
    for x in cluster_rows
)

multi_clusters = sum(
    x["n_members"] > 1
    for x in cluster_rows
)

multi_coassembly = sum(
    x["n_coassemblies"] > 1
    for x in cluster_rows
)

cross_producer = sum(
    x["n_producers"] > 1
    for x in cluster_rows
)

global_fdrs = [
    float(x["fdr"])
    for x in catalog
    if x["fdr"] != ""
]

global_scores = [
    float(x["virus_score"])
    for x in catalog
    if x["virus_score"] != ""
]

summary = [
    (
        "viral_predictions_total",
        len(catalog)
    ),

    (
        "viral_total_bp",
        sum(
            x["length"]
            for x in catalog
        )
    ),

    (
        "provirus_predictions",
        sum(
            x["topology"]
            == "Provirus"
            for x in catalog
        )
    ),

    (
        "non_integrated_predictions",
        sum(
            x["topology"]
            != "Provirus"
            for x in catalog
        )
    ),

    (
        "DTR_predictions",
        sum(
            x["topology"]
            == "DTR"
            for x in catalog
        )
    ),

    (
        "ITR_predictions",
        sum(
            x["topology"]
            == "ITR"
            for x in catalog
        )
    ),

    (
        "score_ge_0.90",
        sum(
            x >= 0.90
            for x in global_scores
        )
    ),

    (
        "score_ge_0.95",
        sum(
            x >= 0.95
            for x in global_scores
        )
    ),

    (
        "score_ge_0.99",
        sum(
            x >= 0.99
            for x in global_scores
        )
    ),

    (
        "fdr_nonNA",
        len(global_fdrs)
    ),

    (
        "fdr_le_0.05",
        sum(
            x <= 0.05
            for x in global_fdrs
        )
    ),

    (
        "fdr_le_0.01",
        sum(
            x <= 0.01
            for x in global_fdrs
        )
    ),

    (
        "hallmark_ge_1",
        sum(
            (
                x["n_hallmarks"]
                != ""
                and int(
                    x["n_hallmarks"]
                ) >= 1
            )
            for x in catalog
        )
    ),

    (
        "hallmark_ge_2",
        sum(
            (
                x["n_hallmarks"]
                != ""
                and int(
                    x["n_hallmarks"]
                ) >= 2
            )
            for x in catalog
        )
    ),

    (
        "taxonomy_beyond_Viruses",
        sum(
            x[
                "taxonomy_beyond_Viruses"
            ]
            == "YES"
            for x in catalog
        )
    ),

    (
        "caudoviricetes_taxonomy",
        sum(
            x[
                "caudoviricetes_in_taxonomy"
            ]
            == "YES"
            for x in catalog
        )
    ),

    (
        "exact_sequence_clusters",
        n_clusters
    ),

    (
        "exact_singleton_clusters",
        singleton_clusters
    ),

    (
        "exact_multimember_clusters",
        multi_clusters
    ),

    (
        "exact_clusters_multiple_coassemblies",
        multi_coassembly
    ),

    (
        "exact_clusters_both_producers",
        cross_producer
    ),
]

for key in [
    "lt_2.5kb",
    "2.5_to_5kb",
    "5_to_10kb",
    "10_to_20kb",
    "20_to_50kb",
    "ge_50kb",
]:

    summary.append(
        (
            f"length_{key}",
            length_bins[key]
        )
    )


summary_file = (
    out
    / "92D_global_summary.tsv"
)

with summary_file.open("w") as fh:

    fh.write("metric\tvalue\n")

    for k, v in summary:
        fh.write(
            f"{k}\t{v}\n"
        )


# ============================================================
# 10. Metadata
# ============================================================

metadata_file = (
    out
    / "92D_run_metadata.tsv"
)

with metadata_file.open("w") as fh:

    fh.write("field\tvalue\n")
    fh.write(
        "analysis\tviral_catalog_consolidation_preQC\n"
    )
    fh.write(
        "source\tgeNomad_1.11.0_default_summary_filters\n"
    )
    fh.write(
        "score_calibration\tENABLED\n"
    )
    fh.write(
        "sequence_dereplication\texact_SHA256_only\n"
    )
    fh.write(
        "vOTU_clustering\tNOT_PERFORMED\n"
    )
    fh.write(
        "CheckV\tNOT_PERFORMED\n"
    )
    fh.write(
        "biological_presence_threshold\tNOT_DEFINED\n"
    )


# ============================================================
# 11. Print
# ============================================================

print()
print("=" * 78)
print("GLOBAL")
print("=" * 78)

for k, v in summary:
    print(k, v, sep="\t")

print()
print("=" * 78)
print("BY COASSEMBLY")
print("=" * 78)

for row in group_rows:

    print(
        row["coassembly"],
        row["viral_predictions"],
        row["proviruses"],
        row["DTR"],
        row["ITR"],
        row["score_ge_0.95"],
        row["fdr_le_0.05"],
        row["hallmark_ge_1"],
        row["taxonomy_beyond_Viruses"],
        row["caudoviricetes_taxonomy"],
        sep="\t"
    )

PY


echo
echo "======================================================================"
echo "VALIDACION"
echo "======================================================================"

for F in \
    "${OUT}/92D_viral_catalog_all.tsv" \
    "${OUT}/92D_viral_sequences_all.fna" \
    "${OUT}/92D_exact_sequence_clusters.tsv" \
    "${OUT}/92D_exact_sequence_members.tsv" \
    "${OUT}/92D_exact_representatives.fna" \
    "${OUT}/92D_summary_by_coassembly.tsv" \
    "${OUT}/92D_terminal_taxonomy_counts.tsv" \
    "${OUT}/92D_global_summary.tsv" \
    "${OUT}/92D_run_metadata.tsv"
do

    [[ -s "$F" ]] || {
        echo "ERROR: salida ausente/vacía: $F" >&2
        exit 2
    }

    echo "OK: $(basename "$F")"

done


echo
echo "======================================================================"
echo "RESUMEN GLOBAL"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92D_global_summary.tsv"

echo
echo "======================================================================"
echo "POR COENSAMBLAJE"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92D_summary_by_coassembly.tsv"

echo
echo "======================================================================"
echo "TOP 30 TAXONOMIAS TERMINALES"
echo "======================================================================"

head -n 31 \
    "${OUT}/92D_terminal_taxonomy_counts.tsv" \
    | column -t -s $'\t'

echo
echo "======================================================================"
echo "PASO 92D FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "======================================================================"

