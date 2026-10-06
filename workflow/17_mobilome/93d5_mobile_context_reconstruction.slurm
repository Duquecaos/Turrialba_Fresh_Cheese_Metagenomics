#!/bin/bash

#SBATCH --job-name=93d5_mobile_ctx
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=93d5_mobile_ctx_%j.out
#SBATCH --error=93d5_mobile_ctx_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"
MINIMAP2_SIF="${MINIMAP2_SIF:-/opt/ohpc/pub/containers/BIO/minimap2-2.28.sif}"

PAIR_FILE="${ROOT}/95_mobilome/93D4_IS_ARG_integron_plasmid/93D4_AMR_IS_pairs_within5kb.tsv"
IS_FILE="${ROOT}/95_mobilome/93D4_IS_ARG_integron_plasmid/93D4_IS_integrated.tsv"
AMR_FILE="${ROOT}/80_resistome_community_consolidated/91E_community_AMR_loci.tsv"
INT_FILE="${ROOT}/95_mobilome/93D2_integron_ARG_plasmid/93D2_integron_structures_integrated.tsv"
DIRECT_INT_AMR="${ROOT}/95_mobilome/93D2_integron_ARG_plasmid/93D2_AMR_integron_direct_links.tsv"

OUT="${ROOT}/95_mobilome/93D5_mobile_contexts"

rm -rf "$OUT"
mkdir -p "$OUT"

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

echo "======================================================================"
echo "93D5 - RECONSTRUCCION Y COMPARACION DE CONTEXTOS MOVILES"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"

for F in \
    "$PAIR_FILE" \
    "$IS_FILE" \
    "$AMR_FILE" \
    "$INT_FILE" \
    "$DIRECT_INT_AMR"
do
    [[ -s "$F" ]] || {
        echo "ERROR: falta entrada $F" >&2
        exit 2
    }
done

[[ -x "$APPTAINER" ]] || {
    echo "ERROR: Apptainer no disponible" >&2
    exit 3
}

[[ -s "$MINIMAP2_SIF" ]] || {
    echo "ERROR: falta $MINIMAP2_SIF" >&2
    exit 4
}


# ======================================================================
# FASE A
# Construcción de regiones locales
# ======================================================================

python3 - \
    "$ROOT" \
    "$PAIR_FILE" \
    "$IS_FILE" \
    "$AMR_FILE" \
    "$INT_FILE" \
    "$DIRECT_INT_AMR" \
    "$OUT" <<'PY'

import csv
import hashlib
import sys
from collections import defaultdict
from pathlib import Path


ROOT = Path(sys.argv[1])
PAIR_FILE = Path(sys.argv[2])
IS_FILE = Path(sys.argv[3])
AMR_FILE = Path(sys.argv[4])
INT_FILE = Path(sys.argv[5])
DIRECT_FILE = Path(sys.argv[6])
OUT = Path(sys.argv[7])

GROUPS = ["L1", "L2", "L3", "M1", "M2", "M3"]

FLANK_BP = 2000


# =====================================================================
# Utilidades
# =====================================================================

def read_tsv(path):
    with path.open(errors="replace") as fh:
        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )
        rows = list(reader)
        fields = reader.fieldnames or []
    return rows, fields


def write_tsv(path, rows, fields):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def to_int(x):
    x = str(x).strip()

    if x in {"", "NA", "None", "nan"}:
        return None

    return int(float(x))


def revcomp(seq):
    table = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan"
    )
    return seq.translate(table)[::-1]


def parse_needed_fasta(path, needed):
    found = {}

    current = None
    chunks = []

    def save():
        nonlocal current, chunks

        if current is not None and current in needed:
            found[current] = "".join(chunks).upper()

    with path.open() as fh:

        for line in fh:

            if line.startswith(">"):

                save()

                current = (
                    line[1:]
                    .strip()
                    .split()[0]
                )

                chunks = []

            else:

                if current in needed:
                    chunks.append(line.strip())

        save()

    return found


def overlaps(a1, a2, b1, b2):
    return max(a1, b1) <= min(a2, b2)


# =====================================================================
# 1. Semillas ARG–IS de 93D4
# =====================================================================

pairs, pair_fields = read_tsv(PAIR_FILE)

required_pair = {
    "community_locus_id",
    "coassembly",
    "contig",
    "AMR_start",
    "AMR_end",
    "AMR_symbol",
    "IS_id",
    "IS_family",
    "IS_type",
    "IS_start",
    "IS_end",
    "distance_bp",
    "AMR_IS_relation",
    "plasmid_candidate",
}

missing = required_pair - set(pair_fields)

if missing:
    raise RuntimeError(
        f"Faltan columnas en 93D4 pairs: {sorted(missing)}"
    )


if len(pairs) != 37:
    raise RuntimeError(
        f"Esperaba 37 pares ARG-IS <=5 kb; observé {len(pairs)}"
    )


unique_seed_arg = {
    r["community_locus_id"]
    for r in pairs
}

if len(unique_seed_arg) != 33:
    raise RuntimeError(
        f"Esperaba 33 ARG únicos; observé {len(unique_seed_arg)}"
    )


for i, row in enumerate(pairs, 1):

    row["_pair_id"] = f"PAIR{i:04d}"

    row["_amr_start"] = min(
        to_int(row["AMR_start"]),
        to_int(row["AMR_end"])
    )

    row["_amr_end"] = max(
        to_int(row["AMR_start"]),
        to_int(row["AMR_end"])
    )

    row["_is_start"] = min(
        to_int(row["IS_start"]),
        to_int(row["IS_end"])
    )

    row["_is_end"] = max(
        to_int(row["IS_start"]),
        to_int(row["IS_end"])
    )


needed_contigs = defaultdict(set)

for row in pairs:
    needed_contigs[row["coassembly"]].add(
        row["contig"]
    )


# =====================================================================
# 2. Recuperar solamente contigs necesarios
# =====================================================================

sequences = {}

for group in GROUPS:

    fasta = (
        ROOT
        / "05_assembly"
        / group
        / "final.contigs.fa"
    )

    if not fasta.is_file():
        raise RuntimeError(
            f"Falta assembly {fasta}"
        )

    found = parse_needed_fasta(
        fasta,
        needed_contigs[group]
    )

    missing_contigs = (
        needed_contigs[group]
        - set(found)
    )

    if missing_contigs:
        raise RuntimeError(
            f"{group}: contigs no encontrados: "
            f"{sorted(missing_contigs)}"
        )

    for contig, seq in found.items():
        sequences[(group, contig)] = seq


# =====================================================================
# 3. Intervalo inicial por pareja
# =====================================================================

pair_intervals = defaultdict(list)

for row in pairs:

    key = (
        row["coassembly"],
        row["contig"]
    )

    contig_len = len(sequences[key])

    start = max(
        1,
        min(
            row["_amr_start"],
            row["_is_start"]
        ) - FLANK_BP
    )

    end = min(
        contig_len,
        max(
            row["_amr_end"],
            row["_is_end"]
        ) + FLANK_BP
    )

    pair_intervals[key].append({
        "start": start,
        "end": end,
        "pair": row,
    })


# =====================================================================
# 4. Fusionar contextos solapados en el mismo contig
# =====================================================================

contexts_raw = []

for key, intervals in pair_intervals.items():

    group, contig = key

    intervals = sorted(
        intervals,
        key=lambda x: (
            x["start"],
            x["end"]
        )
    )

    current = None

    for item in intervals:

        if current is None:

            current = {
                "coassembly": group,
                "contig": contig,
                "start": item["start"],
                "end": item["end"],
                "pairs": [item["pair"]],
            }

            continue


        if item["start"] <= current["end"]:

            current["end"] = max(
                current["end"],
                item["end"]
            )

            current["pairs"].append(
                item["pair"]
            )

        else:

            contexts_raw.append(current)

            current = {
                "coassembly": group,
                "contig": contig,
                "start": item["start"],
                "end": item["end"],
                "pairs": [item["pair"]],
            }


    if current is not None:
        contexts_raw.append(current)


contexts_raw.sort(
    key=lambda x: (
        GROUPS.index(x["coassembly"]),
        x["contig"],
        x["start"],
        x["end"]
    )
)


# =====================================================================
# 5. Cargar todas las anotaciones
# =====================================================================

amr_rows, amr_fields = read_tsv(AMR_FILE)

if len(amr_rows) != 565:
    raise RuntimeError(
        f"AMR esperado=565 observado={len(amr_rows)}"
    )


for r in amr_rows:

    r["_start"] = min(
        to_int(r["start"]),
        to_int(r["stop"])
    )

    r["_end"] = max(
        to_int(r["start"]),
        to_int(r["stop"])
    )


amr_by_contig = defaultdict(list)

for r in amr_rows:
    amr_by_contig[
        (
            r["coassembly"],
            r["contig"]
        )
    ].append(r)


is_rows, is_fields = read_tsv(IS_FILE)

if len(is_rows) != 6391:
    raise RuntimeError(
        f"IS esperado=6391 observado={len(is_rows)}"
    )


for r in is_rows:

    r["_start"] = min(
        to_int(r["is_start"]),
        to_int(r["is_end"])
    )

    r["_end"] = max(
        to_int(r["is_start"]),
        to_int(r["is_end"])
    )


is_by_contig = defaultdict(list)

for r in is_rows:
    is_by_contig[
        (
            r["coassembly"],
            r["contig"]
        )
    ].append(r)


int_rows, int_fields = read_tsv(INT_FILE)

if len(int_rows) != 375:
    raise RuntimeError(
        f"Integrones esperado=375 observado={len(int_rows)}"
    )


for r in int_rows:

    r["_start"] = min(
        to_int(r["structure_start"]),
        to_int(r["structure_end"])
    )

    r["_end"] = max(
        to_int(r["structure_start"]),
        to_int(r["structure_end"])
    )


int_by_contig = defaultdict(list)

for r in int_rows:
    int_by_contig[
        (
            r["coassembly"],
            r["contig"]
        )
    ].append(r)


direct_rows, direct_fields = read_tsv(
    DIRECT_FILE
)

direct_ids = {
    r["community_locus_id"]
    for r in direct_rows
}


# =====================================================================
# 6. Construir contextos anotados
# =====================================================================

contexts = []
context_sequences = {}
feature_rows = []
pair_map_rows = []


for idx, ctx in enumerate(contexts_raw, 1):

    context_id = f"MOBCTX{idx:04d}"

    group = ctx["coassembly"]
    contig = ctx["contig"]

    start = ctx["start"]
    end = ctx["end"]

    seq_full = sequences[(group, contig)]

    seq = seq_full[
        start - 1:end
    ]


    all_amr = [
        r
        for r in amr_by_contig.get(
            (group, contig),
            []
        )
        if overlaps(
            start,
            end,
            r["_start"],
            r["_end"]
        )
    ]


    all_is = [
        r
        for r in is_by_contig.get(
            (group, contig),
            []
        )
        if overlaps(
            start,
            end,
            r["_start"],
            r["_end"]
        )
    ]


    all_int = [
        r
        for r in int_by_contig.get(
            (group, contig),
            []
        )
        if overlaps(
            start,
            end,
            r["_start"],
            r["_end"]
        )
    ]


    seed_arg_ids = sorted({
        r["community_locus_id"]
        for r in ctx["pairs"]
    })


    seed_symbols = sorted({
        r["AMR_symbol"]
        for r in ctx["pairs"]
    })


    seed_is_ids = sorted({
        r["IS_id"]
        for r in ctx["pairs"]
    })


    seed_is_families = sorted({
        r["IS_family"]
        for r in ctx["pairs"]
    })


    plasmid_values = {
        str(r["plasmid_candidate"]).strip()
        for r in all_is
    }

    plasmid_candidate = (
        "YES"
        if "YES" in plasmid_values
        else "NO"
    )


    mags = sorted({
        mag
        for r in all_is
        for mag in str(
            r.get(
                "final18_MAGs",
                ""
            )
        ).split(";")
        if mag
    })


    taxa = sorted({
        tax
        for r in all_is
        for tax in str(
            r.get(
                "final18_taxa",
                ""
            )
        ).split(";")
        if tax
    })


    direct_arg_ids = sorted(
        set(seed_arg_ids)
        & direct_ids
    )


    component_signature = (
        "ARG="
        + ",".join(
            sorted({
                r["element_symbol"]
                for r in all_amr
            })
        )
        + "|IS="
        + ",".join(
            sorted({
                f"{r['family']}:{r['type']}"
                for r in all_is
            })
        )
        + "|INT="
        + ",".join(
            sorted({
                r["integron_type"]
                for r in all_int
            })
        )
        + "|PLASMID="
        + plasmid_candidate
    )


    canonical_seq = min(
        seq,
        revcomp(seq)
    )

    canonical_sha = hashlib.sha256(
        canonical_seq.encode()
    ).hexdigest()


    row = {
        "context_id": context_id,
        "coassembly": group,
        "producer": group[0],
        "biological_unit": group[1:],
        "contig": contig,
        "context_start": start,
        "context_end": end,
        "context_length_bp": len(seq),
        "contig_length_bp": len(seq_full),
        "left_edge_clipped":
            "YES" if start == 1 else "NO",
        "right_edge_clipped":
            "YES" if end == len(seq_full) else "NO",
        "seed_pair_count": len(ctx["pairs"]),
        "seed_ARG_count": len(seed_arg_ids),
        "seed_ARG_loci": ";".join(seed_arg_ids),
        "seed_ARG_symbols": ";".join(seed_symbols),
        "seed_IS_count": len(seed_is_ids),
        "seed_IS_ids": ";".join(seed_is_ids),
        "seed_IS_families": ";".join(seed_is_families),
        "all_ARG_count": len(all_amr),
        "all_ARG_symbols": ";".join(
            sorted({
                r["element_symbol"]
                for r in all_amr
            })
        ),
        "all_IS_count": len(all_is),
        "complete_IS_count": sum(
            r["type"] == "c"
            for r in all_is
        ),
        "partial_IS_count": sum(
            r["type"] == "p"
            for r in all_is
        ),
        "all_IS_families": ";".join(
            sorted({
                r["family"]
                for r in all_is
            })
        ),
        "integron_structure_count": len(
            all_int
        ),
        "integron_types": ";".join(
            sorted({
                r["integron_type"]
                for r in all_int
            })
        ),
        "direct_integron_ARG_count": len(
            direct_arg_ids
        ),
        "direct_integron_ARG_loci":
            ";".join(direct_arg_ids),
        "plasmid_candidate":
            plasmid_candidate,
        "n_final18_MAGs": len(mags),
        "final18_MAGs": ";".join(mags),
        "final18_taxa": ";".join(taxa),
        "component_signature":
            component_signature,
        "canonical_sequence_sha256":
            canonical_sha,
    }

    contexts.append(row)

    context_sequences[
        context_id
    ] = seq


    # -------------------------------------------------------------
    # Features
    # -------------------------------------------------------------

    for r in all_amr:

        feature_rows.append({
            "context_id": context_id,
            "feature_type": "ARG",
            "feature_id":
                r["community_locus_id"],
            "feature_label":
                r["element_symbol"],
            "feature_start":
                r["_start"],
            "feature_end":
                r["_end"],
            "feature_strand":
                r.get("strand", ""),
            "feature_subtype":
                r.get("class", ""),
        })


    for r in all_is:

        feature_rows.append({
            "context_id": context_id,
            "feature_type": "IS",
            "feature_id":
                r["IS_id"],
            "feature_label":
                r["family"],
            "feature_start":
                r["_start"],
            "feature_end":
                r["_end"],
            "feature_strand":
                r.get("strand", ""),
            "feature_subtype":
                (
                    "complete"
                    if r["type"] == "c"
                    else "partial"
                ),
        })


    for r in all_int:

        feature_rows.append({
            "context_id": context_id,
            "feature_type":
                "INTEGRON_STRUCTURE",
            "feature_id":
                r["structure_id"],
            "feature_label":
                r["integron_type"],
            "feature_start":
                r["_start"],
            "feature_end":
                r["_end"],
            "feature_strand":
                "",
            "feature_subtype":
                r["integron_type"],
        })


    for r in ctx["pairs"]:

        pair_map_rows.append({
            "pair_id":
                r["_pair_id"],
            "context_id":
                context_id,
            "community_locus_id":
                r["community_locus_id"],
            "AMR_symbol":
                r["AMR_symbol"],
            "IS_id":
                r["IS_id"],
            "IS_family":
                r["IS_family"],
            "IS_type":
                r["IS_type"],
            "distance_bp":
                r["distance_bp"],
            "AMR_IS_relation":
                r["AMR_IS_relation"],
        })


# =====================================================================
# 7. Exact sequence clusters, orientation-insensitive
# =====================================================================

hash_groups = defaultdict(list)

for row in contexts:
    hash_groups[
        row["canonical_sequence_sha256"]
    ].append(
        row["context_id"]
    )


ordered_hash_groups = sorted(
    hash_groups.values(),
    key=lambda members: (
        min(members),
        len(members)
    )
)


exact_cluster_by_context = {}

exact_cluster_rows = []


for idx, members in enumerate(
    ordered_hash_groups,
    1
):

    cluster_id = (
        f"EXACTCTX{idx:04d}"
    )

    members = sorted(members)

    for context_id in members:
        exact_cluster_by_context[
            context_id
        ] = cluster_id

    exact_cluster_rows.append({
        "exact_cluster_id":
            cluster_id,
        "n_contexts":
            len(members),
        "members":
            ";".join(members),
    })


for row in contexts:
    row["exact_cluster_id"] = (
        exact_cluster_by_context[
            row["context_id"]
        ]
    )


# =====================================================================
# 8. Escritura
# =====================================================================

CONTEXT_FIELDS = [
    "context_id",
    "coassembly",
    "producer",
    "biological_unit",
    "contig",
    "context_start",
    "context_end",
    "context_length_bp",
    "contig_length_bp",
    "left_edge_clipped",
    "right_edge_clipped",
    "seed_pair_count",
    "seed_ARG_count",
    "seed_ARG_loci",
    "seed_ARG_symbols",
    "seed_IS_count",
    "seed_IS_ids",
    "seed_IS_families",
    "all_ARG_count",
    "all_ARG_symbols",
    "all_IS_count",
    "complete_IS_count",
    "partial_IS_count",
    "all_IS_families",
    "integron_structure_count",
    "integron_types",
    "direct_integron_ARG_count",
    "direct_integron_ARG_loci",
    "plasmid_candidate",
    "n_final18_MAGs",
    "final18_MAGs",
    "final18_taxa",
    "component_signature",
    "canonical_sequence_sha256",
    "exact_cluster_id",
]

write_tsv(
    OUT / "93D5_local_contexts.tsv",
    contexts,
    CONTEXT_FIELDS,
)


FEATURE_FIELDS = [
    "context_id",
    "feature_type",
    "feature_id",
    "feature_label",
    "feature_start",
    "feature_end",
    "feature_strand",
    "feature_subtype",
]

write_tsv(
    OUT / "93D5_feature_inventory.tsv",
    feature_rows,
    FEATURE_FIELDS,
)


PAIR_MAP_FIELDS = [
    "pair_id",
    "context_id",
    "community_locus_id",
    "AMR_symbol",
    "IS_id",
    "IS_family",
    "IS_type",
    "distance_bp",
    "AMR_IS_relation",
]

write_tsv(
    OUT / "93D5_seed_pair_to_context.tsv",
    pair_map_rows,
    PAIR_MAP_FIELDS,
)


write_tsv(
    OUT / "93D5_exact_context_clusters.tsv",
    exact_cluster_rows,
    [
        "exact_cluster_id",
        "n_contexts",
        "members",
    ],
)


with (
    OUT / "93D5_local_contexts.fna"
).open("w") as fh:

    for row in contexts:

        cid = row["context_id"]

        header = (
            f">{cid}"
            f"|{row['coassembly']}"
            f"|{row['contig']}"
            f":{row['context_start']}"
            f"-{row['context_end']}"
        )

        fh.write(header + "\n")

        seq = context_sequences[cid]

        for i in range(
            0,
            len(seq),
            80
        ):
            fh.write(
                seq[i:i+80]
                + "\n"
            )


with (
    OUT / "93D5_precluster_summary.tsv"
).open("w") as fh:

    fh.write("metric\tvalue\n")
    fh.write(
        f"seed_ARG_IS_pairs\t{len(pairs)}\n"
    )
    fh.write(
        f"unique_seed_ARG_loci\t"
        f"{len(unique_seed_arg)}\n"
    )
    fh.write(
        f"merged_local_contexts\t"
        f"{len(contexts)}\n"
    )
    fh.write(
        f"contexts_on_plasmid_candidate_contigs\t"
        f"{sum(r['plasmid_candidate']=='YES' for r in contexts)}\n"
    )
    fh.write(
        f"contexts_on_final18_MAG_contigs\t"
        f"{sum(r['n_final18_MAGs']>0 for r in contexts)}\n"
    )
    fh.write(
        f"contexts_with_integron_structure\t"
        f"{sum(r['integron_structure_count']>0 for r in contexts)}\n"
    )
    fh.write(
        f"contexts_with_direct_integron_ARG\t"
        f"{sum(r['direct_integron_ARG_count']>0 for r in contexts)}\n"
    )
    fh.write(
        f"orientation_insensitive_exact_clusters\t"
        f"{len(exact_cluster_rows)}\n"
    )
    fh.write(
        f"orientation_insensitive_exact_multimember_clusters\t"
        f"{sum(r['n_contexts']>1 for r in exact_cluster_rows)}\n"
    )


print(
    f"seed_pairs={len(pairs)}"
)

print(
    f"unique_seed_ARG={len(unique_seed_arg)}"
)

print(
    f"merged_contexts={len(contexts)}"
)

print(
    "PRECLUSTER_PASS"
)

PY


FASTA="${OUT}/93D5_local_contexts.fna"
PAF="${OUT}/93D5_context_all_vs_all.paf"


[[ -s "$FASTA" ]] || {
    echo "ERROR: FASTA de contextos vacío" >&2
    exit 5
}


echo
echo "======================================================================"
echo "MINIMAP2 ALL-VS-ALL"
echo "======================================================================"

srun --exact \
    -n 1 \
    -c "$SLURM_CPUS_PER_TASK" \
    "$APPTAINER" exec \
    --bind "${ROOT}:${ROOT}" \
    "$MINIMAP2_SIF" \
    minimap2 \
    -x asm5 \
    -t "$SLURM_CPUS_PER_TASK" \
    -N 1000 \
    --secondary=yes \
    "$FASTA" \
    "$FASTA" \
    > "$PAF"


[[ -s "$PAF" ]] || {
    echo "ERROR: PAF vacío" >&2
    exit 6
}


# ======================================================================
# FASE B
# Clustering operacional de contextos
# ======================================================================

python3 - \
    "$OUT" \
    "$PAF" <<'PY'

import csv
import sys
from collections import defaultdict
from pathlib import Path


OUT = Path(sys.argv[1])
PAF = Path(sys.argv[2])

IDENTITY_THRESHOLD = 0.95
SHORTER_FRACTION_THRESHOLD = 0.80


def read_tsv(path):

    with path.open() as fh:

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
            extrasaction="ignore"
        )

        writer.writeheader()

        writer.writerows(rows)


contexts, context_fields = read_tsv(
    OUT / "93D5_local_contexts.tsv"
)

context_by_id = {
    r["context_id"]: r
    for r in contexts
}

ids = sorted(context_by_id)


# =====================================================================
# Union-find
# =====================================================================

parent = {
    x: x
    for x in ids
}


def find(x):

    while parent[x] != x:

        parent[x] = parent[
            parent[x]
        ]

        x = parent[x]

    return x


def union(a, b):

    ra = find(a)
    rb = find(b)

    if ra == rb:
        return

    if ra < rb:
        parent[rb] = ra
    else:
        parent[ra] = rb


# =====================================================================
# Exact clusters always collapse
# =====================================================================

exact_groups = defaultdict(list)

for r in contexts:
    exact_groups[
        r["exact_cluster_id"]
    ].append(
        r["context_id"]
    )


for members in exact_groups.values():

    if len(members) < 2:
        continue

    first = members[0]

    for x in members[1:]:
        union(first, x)


# =====================================================================
# Parse PAF
# =====================================================================

best_pair = {}


with PAF.open() as fh:

    for line in fh:

        if not line.strip():
            continue

        f = line.rstrip(
            "\n"
        ).split("\t")

        if len(f) < 12:
            continue

        qname = f[0].split(
            "|",
            1
        )[0]

        qlen = int(f[1])
        qstart = int(f[2])
        qend = int(f[3])

        strand = f[4]

        tname = f[5].split(
            "|",
            1
        )[0]

        tlen = int(f[6])
        tstart = int(f[7])
        tend = int(f[8])

        nmatch = int(f[9])
        alnlen = int(f[10])
        mapq = int(f[11])


        if qname == tname:
            continue


        if (
            qname not in context_by_id
            or tname not in context_by_id
        ):
            continue


        identity = (
            nmatch / alnlen
            if alnlen > 0
            else 0.0
        )


        qspan = qend - qstart
        tspan = tend - tstart

        shorter_len = min(
            qlen,
            tlen
        )

        aligned_fraction_shorter = (
            min(
                qspan,
                tspan
            )
            / shorter_len
            if shorter_len > 0
            else 0.0
        )


        key = tuple(
            sorted(
                (
                    qname,
                    tname
                )
            )
        )


        candidate = {
            "context_A":
                key[0],
            "context_B":
                key[1],
            "alignment_identity":
                identity,
            "aligned_fraction_shorter":
                aligned_fraction_shorter,
            "alignment_length":
                alnlen,
            "strand":
                strand,
            "mapq":
                mapq,
        }


        previous = best_pair.get(
            key
        )


        if previous is None:

            best_pair[key] = candidate

        else:

            old_key = (
                previous[
                    "aligned_fraction_shorter"
                ],
                previous[
                    "alignment_identity"
                ],
                previous[
                    "alignment_length"
                ],
            )

            new_key = (
                candidate[
                    "aligned_fraction_shorter"
                ],
                candidate[
                    "alignment_identity"
                ],
                candidate[
                    "alignment_length"
                ],
            )

            if new_key > old_key:
                best_pair[key] = candidate


pair_rows = []

for key in sorted(best_pair):

    r = best_pair[key]

    passed = (
        r["alignment_identity"]
        >= IDENTITY_THRESHOLD
        and
        r[
            "aligned_fraction_shorter"
        ]
        >= SHORTER_FRACTION_THRESHOLD
    )

    r["passes_context_homology_threshold"] = (
        "YES"
        if passed
        else "NO"
    )

    pair_rows.append(r)

    if passed:

        union(
            r["context_A"],
            r["context_B"]
        )


# =====================================================================
# Clusters
# =====================================================================

groups = defaultdict(list)

for context_id in ids:

    groups[
        find(context_id)
    ].append(
        context_id
    )


ordered_clusters = sorted(
    groups.values(),
    key=lambda members: (
        min(members),
        len(members)
    )
)


cluster_rows = []

membership_rows = []


for idx, members in enumerate(
    ordered_clusters,
    1
):

    cluster_id = (
        f"HOMCTX{idx:04d}"
    )

    members = sorted(members)

    rows = [
        context_by_id[x]
        for x in members
    ]


    coassemblies = sorted({
        r["coassembly"]
        for r in rows
    })


    producers = sorted({
        r["producer"]
        for r in rows
    })


    arg_symbols = sorted({
        symbol
        for r in rows
        for symbol in r[
            "all_ARG_symbols"
        ].split(";")
        if symbol
    })


    is_families = sorted({
        family
        for r in rows
        for family in r[
            "all_IS_families"
        ].split(";")
        if family
    })


    int_types = sorted({
        itype
        for r in rows
        for itype in r[
            "integron_types"
        ].split(";")
        if itype
    })


    mags = sorted({
        mag
        for r in rows
        for mag in r[
            "final18_MAGs"
        ].split(";")
        if mag
    })


    taxa = sorted({
        tax
        for r in rows
        for tax in r[
            "final18_taxa"
        ].split(";")
        if tax
    })


    exact_clusters = sorted({
        r["exact_cluster_id"]
        for r in rows
    })


    cluster_rows.append({
        "homology_cluster_id":
            cluster_id,
        "n_contexts":
            len(members),
        "members":
            ";".join(members),
        "n_coassemblies":
            len(coassemblies),
        "coassemblies":
            ";".join(coassemblies),
        "n_producers":
            len(producers),
        "producers":
            ";".join(producers),
        "both_producers":
            (
                "YES"
                if set(producers)
                == {"L", "M"}
                else "NO"
            ),
        "AMR_symbols":
            ";".join(arg_symbols),
        "IS_families":
            ";".join(is_families),
        "integron_types":
            ";".join(int_types),
        "plasmid_candidate_contexts":
            sum(
                r[
                    "plasmid_candidate"
                ] == "YES"
                for r in rows
            ),
        "final18_MAGs":
            ";".join(mags),
        "final18_taxa":
            ";".join(taxa),
        "n_exact_sequence_clusters":
            len(exact_clusters),
        "exact_sequence_clusters":
            ";".join(exact_clusters),
    })


    for context_id in members:

        membership_rows.append({
            "context_id":
                context_id,
            "homology_cluster_id":
                cluster_id,
        })


# =====================================================================
# Recurrence subsets
# =====================================================================

recurrent = [
    r
    for r in cluster_rows
    if int(
        r["n_contexts"]
    ) > 1
]


cross_coassembly = [
    r
    for r in recurrent
    if int(
        r["n_coassemblies"]
    ) > 1
]


both_producers = [
    r
    for r in recurrent
    if r["both_producers"] == "YES"
]


# =====================================================================
# Write
# =====================================================================

PAIR_FIELDS = [
    "context_A",
    "context_B",
    "alignment_identity",
    "aligned_fraction_shorter",
    "alignment_length",
    "strand",
    "mapq",
    "passes_context_homology_threshold",
]

write_tsv(
    OUT / "93D5_context_pairwise_similarity.tsv",
    pair_rows,
    PAIR_FIELDS,
)


CLUSTER_FIELDS = [
    "homology_cluster_id",
    "n_contexts",
    "members",
    "n_coassemblies",
    "coassemblies",
    "n_producers",
    "producers",
    "both_producers",
    "AMR_symbols",
    "IS_families",
    "integron_types",
    "plasmid_candidate_contexts",
    "final18_MAGs",
    "final18_taxa",
    "n_exact_sequence_clusters",
    "exact_sequence_clusters",
]

write_tsv(
    OUT / "93D5_context_homology_clusters.tsv",
    cluster_rows,
    CLUSTER_FIELDS,
)

write_tsv(
    OUT / "93D5_recurrent_context_clusters.tsv",
    recurrent,
    CLUSTER_FIELDS,
)

write_tsv(
    OUT / "93D5_cross_coassembly_context_clusters.tsv",
    cross_coassembly,
    CLUSTER_FIELDS,
)

write_tsv(
    OUT / "93D5_both_producer_context_clusters.tsv",
    both_producers,
    CLUSTER_FIELDS,
)

write_tsv(
    OUT / "93D5_context_cluster_membership.tsv",
    membership_rows,
    [
        "context_id",
        "homology_cluster_id",
    ],
)


exact_rows, _ = read_tsv(
    OUT
    / "93D5_exact_context_clusters.tsv"
)


with (
    OUT
    / "93D5_global_summary.tsv"
).open("w") as fh:

    fh.write("metric\tvalue\n")

    fh.write(
        f"seed_ARG_IS_pairs\t37\n"
    )

    fh.write(
        f"unique_seed_ARG_loci\t33\n"
    )

    fh.write(
        f"merged_local_contexts\t"
        f"{len(contexts)}\n"
    )

    fh.write(
        f"orientation_insensitive_exact_clusters\t"
        f"{len(exact_rows)}\n"
    )

    fh.write(
        f"exact_multimember_clusters\t"
        f"{sum(int(r['n_contexts']) > 1 for r in exact_rows)}\n"
    )

    fh.write(
        f"operational_context_homology_clusters\t"
        f"{len(cluster_rows)}\n"
    )

    fh.write(
        f"homology_multimember_clusters\t"
        f"{len(recurrent)}\n"
    )

    fh.write(
        f"cross_coassembly_homology_clusters\t"
        f"{len(cross_coassembly)}\n"
    )

    fh.write(
        f"both_producer_homology_clusters\t"
        f"{len(both_producers)}\n"
    )

    fh.write(
        f"identity_threshold\t"
        f"{IDENTITY_THRESHOLD}\n"
    )

    fh.write(
        f"aligned_fraction_shorter_threshold\t"
        f"{SHORTER_FRACTION_THRESHOLD}\n"
    )


with (
    OUT
    / "93D5_metadata.tsv"
).open("w") as fh:

    fh.write("field\tvalue\n")

    fh.write(
        "context_definition\t"
        "union_of_ARG_and_IS_pair_plus_2000bp_each_side\n"
    )

    fh.write(
        "overlapping_seed_intervals_same_contig\tMERGED\n"
    )

    fh.write(
        "sequence_comparison_tool\tminimap2_2.28_asm5\n"
    )

    fh.write(
        "context_homology_identity_threshold\t0.95\n"
    )

    fh.write(
        "context_homology_shorter_sequence_fraction_threshold\t0.80\n"
    )

    fh.write(
        "context_homology_cluster_equals_same_mobile_element\tNO\n"
    )

    fh.write(
        "cross_coassembly_recurrence_equals_prevalence\tNO\n"
    )

    fh.write(
        "ARG_IS_colocalization_equals_HGT\tNO\n"
    )

    fh.write(
        "plasmid_candidate_equals_validated_plasmid\tNO\n"
    )

    fh.write(
        "partial_IS_equals_complete_mobile_element\tNO\n"
    )


print(
    f"contexts={len(contexts)}"
)

print(
    f"homology_clusters={len(cluster_rows)}"
)

print(
    f"recurrent_clusters={len(recurrent)}"
)

print(
    f"cross_coassembly={len(cross_coassembly)}"
)

print(
    f"both_producers={len(both_producers)}"
)

print(
    "93D5_CONTEXT_CLUSTERING=PASS"
)

PY


# ======================================================================
# Final QC
# ======================================================================

REQUIRED=(
    "93D5_local_contexts.tsv"
    "93D5_local_contexts.fna"
    "93D5_feature_inventory.tsv"
    "93D5_seed_pair_to_context.tsv"
    "93D5_exact_context_clusters.tsv"
    "93D5_context_pairwise_similarity.tsv"
    "93D5_context_homology_clusters.tsv"
    "93D5_recurrent_context_clusters.tsv"
    "93D5_cross_coassembly_context_clusters.tsv"
    "93D5_both_producer_context_clusters.tsv"
    "93D5_context_cluster_membership.tsv"
    "93D5_global_summary.tsv"
    "93D5_metadata.tsv"
)

for NAME in "${REQUIRED[@]}"
do
    [[ -s "${OUT}/${NAME}" ]] || {
        echo "ERROR: salida faltante/vacía ${OUT}/${NAME}" >&2
        exit 7
    }
done


sha256sum \
    "${OUT}/93D5_local_contexts.tsv" \
    "${OUT}/93D5_local_contexts.fna" \
    "${OUT}/93D5_context_homology_clusters.tsv" \
    > "${OUT}/93D5_output_sha256.tsv"


printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/93D5_COMPLETE.ok"


echo
echo "======================================================================"
echo "93D5 GLOBAL"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93D5_global_summary.tsv"

echo
echo "======================================================================"
echo "RECURRENT CONTEXT CLUSTERS"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93D5_recurrent_context_clusters.tsv"

echo
echo "======================================================================"
echo "93D5 FINALIZO CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "======================================================================"

