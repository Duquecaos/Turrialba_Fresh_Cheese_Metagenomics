#!/bin/bash

#SBATCH --job-name=92f_votu_repair
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=01:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=92f_votu_repair_%j.out
#SBATCH --error=92f_votu_repair_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

INPUT="${ROOT}/89_viral_catalog_hostclean/92E2_hostclean_exact_representatives.fna"

HOSTCLEAN_CLUSTERS="${ROOT}/89_viral_catalog_hostclean/92E2_hostclean_exact_clusters.tsv"
HOSTCLEAN_MEMBERS="${ROOT}/89_viral_catalog_hostclean/92E2_hostclean_exact_members.tsv"

PRECLUSTERS="${ROOT}/86_viral_catalog_preQC/92D_exact_sequence_clusters.tsv"

POSTCHECKV="${ROOT}/90_viral_hostclean_validation/checkv_hostclean/quality_summary.tsv"
PRECHECKV="${ROOT}/88_checkv_viral_catalog/92E1_checkv_integrated.tsv"

OUT="${ROOT}/91_votu_catalog"
WORK="${OUT}/work"

BLASTOUT="${WORK}/92F_all_vs_all_blast.tsv"
ANIOUT="${WORK}/92F_pairwise_ani.tsv"
RAWCLUST="${WORK}/92F_aniclust_raw.tsv"

DONE="${OUT}/92F_COMPLETE.ok"

mkdir -p "$OUT"

echo "======================================================================"
echo "PASO 92F-REPAIR"
echo "Reparación de integración del clustering vOTU"
echo "NO se recalcula BLAST, ANI ni aniclust"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"

# ======================================================================
# 1. VALIDAR QUE EL CLUSTERING PRODUCTIVO YA EXISTE
# ======================================================================

for F in \
    "$INPUT" \
    "$HOSTCLEAN_CLUSTERS" \
    "$HOSTCLEAN_MEMBERS" \
    "$PRECLUSTERS" \
    "$POSTCHECKV" \
    "$PRECHECKV" \
    "$BLASTOUT" \
    "$ANIOUT" \
    "$RAWCLUST"
do
    [[ -s "$F" ]] || {
        echo "ERROR: falta archivo requerido:"
        echo "$F"
        exit 2
    }
done

NSEQ="$(grep -c '^>' "$INPUT")"

echo "VCLEAN sequences = ${NSEQ}"

[[ "$NSEQ" -eq 6748 ]] || {
    echo "ERROR: esperaba 6748 VCLEAN." >&2
    exit 3
}

echo
echo "Archivos productivos existentes:"

du -h \
    "$BLASTOUT" \
    "$ANIOUT" \
    "$RAWCLUST"


# ======================================================================
# 2. PREVIEW DEL FORMATO REAL ANICLUST
# ======================================================================

echo
echo "======================================================================"
echo "FORMATO REAL DE ANICLUST"
echo "======================================================================"

echo "Primeras 10 líneas:"
head -n 10 "$RAWCLUST"

echo
echo "Primeros clusters multimember:"

awk -F '\t' '$2 ~ /,/ {print; n++; if(n==10) exit}' \
    "$RAWCLUST" \
    || true


# ======================================================================
# 3. ELIMINAR SOLO SALIDAS PARCIALES DE LA INTEGRACIÓN
# ======================================================================

rm -f \
    "${OUT}/92F_votu_membership.tsv" \
    "${OUT}/92F_votu_summary.tsv" \
    "${OUT}/92F_votu_representatives.fna" \
    "${OUT}/92F_votu_representative_map.tsv" \
    "${OUT}/92F_global_summary.tsv" \
    "${OUT}/92F_cluster_size_distribution.tsv" \
    "${OUT}/92F_run_metadata.tsv" \
    "${OUT}/92F_sha256.tsv" \
    "$DONE"


# ======================================================================
# 4. INTEGRACIÓN CORREGIDA
# ======================================================================

python3 - <<'PY'

import csv
import os
import statistics

from collections import Counter, defaultdict
from pathlib import Path


root = Path(
    f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"
)

out = root / "91_votu_catalog"
work = out / "work"


fasta_file = (
    root
    / "89_viral_catalog_hostclean"
    / "92E2_hostclean_exact_representatives.fna"
)

hostclean_clusters_file = (
    root
    / "89_viral_catalog_hostclean"
    / "92E2_hostclean_exact_clusters.tsv"
)

hostclean_members_file = (
    root
    / "89_viral_catalog_hostclean"
    / "92E2_hostclean_exact_members.tsv"
)

preclusters_file = (
    root
    / "86_viral_catalog_preQC"
    / "92D_exact_sequence_clusters.tsv"
)

postcheckv_file = (
    root
    / "90_viral_hostclean_validation"
    / "checkv_hostclean"
    / "quality_summary.tsv"
)

precheckv_file = (
    root
    / "88_checkv_viral_catalog"
    / "92E1_checkv_integrated.tsv"
)

rawcluster_file = (
    work
    / "92F_aniclust_raw.tsv"
)


# ======================================================================
# HELPERS
# ======================================================================

def clean(x):

    if x is None:
        return ""

    x = str(x).strip()

    if x.upper() in {
        "",
        "NA",
        "NAN",
        "NONE",
    }:
        return ""

    return x


def num(x):

    x = clean(x)

    if not x:
        return None

    try:
        return float(x)

    except ValueError:
        return None


def splitsemi(x):

    return {
        y.strip()
        for y in clean(x).split(";")
        if y.strip()
    }


def read_fasta(path):

    seqs = {}

    name = None
    pieces = []

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

                    seqs[name] = "".join(pieces)

                name = line[1:].split()[0]
                pieces = []

            else:

                pieces.append(line)

        if name is not None:

            if name in seqs:
                raise RuntimeError(
                    f"FASTA ID duplicado: {name}"
                )

            seqs[name] = "".join(pieces)

    return seqs


def write_fasta(fh, name, seq, width=80):

    fh.write(f">{name}\n")

    for i in range(
        0,
        len(seq),
        width
    ):
        fh.write(
            seq[i:i+width] + "\n"
        )


def read_dict(path, key):

    with path.open() as fh:

        rows = list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )

    result = {}

    for row in rows:

        k = row[key]

        if k in result:
            raise RuntimeError(
                f"Clave duplicada {k} en {path}"
            )

        result[k] = row

    return result


# ======================================================================
# INPUTS
# ======================================================================

seqs = read_fasta(
    fasta_file
)

if len(seqs) != 6748:

    raise RuntimeError(
        f"FASTA VCLEAN: esperaba 6748, "
        f"encontré {len(seqs)}"
    )


hc = read_dict(
    hostclean_clusters_file,
    "hostclean_exact_cluster_id"
)

preclusters = read_dict(
    preclusters_file,
    "exact_cluster_id"
)

postcv = read_dict(
    postcheckv_file,
    "contig_id"
)

precv = read_dict(
    precheckv_file,
    "exact_cluster_id"
)


if set(seqs) != set(hc):

    raise RuntimeError(
        "VCLEAN FASTA != hostclean cluster IDs"
    )


if set(seqs) != set(postcv):

    raise RuntimeError(
        "VCLEAN FASTA != postclean CheckV IDs"
    )


# ======================================================================
# MAPEO VCLEAN -> VEXACT
# ======================================================================

vclean_to_pre = defaultdict(list)

with hostclean_members_file.open() as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in reader:

        vclean_to_pre[
            row["hostclean_exact_cluster_id"]
        ].append(
            row["preclean_exact_cluster_id"]
        )


if set(vclean_to_pre) != set(seqs):

    missing = (
        set(seqs)
        - set(vclean_to_pre)
    )

    raise RuntimeError(
        f"Mapeo VCLEAN -> VEXACT incompleto. "
        f"Faltan {len(missing)}"
    )


# ======================================================================
# PARSEAR ANICLUST
#
# FORMATO REAL:
#
# SEED<TAB>MEMBER1,MEMBER2,MEMBER3,...
#
# Para singletons:
#
# SEED<TAB>SEED
# ======================================================================

seed_to_members = {}

all_members_seen = set()

raw_lines = 0
multimember_raw = 0


with rawcluster_file.open() as fh:

    for lineno, line in enumerate(
        fh,
        start=1
    ):

        line = line.rstrip("\n")

        if not line:
            continue

        raw_lines += 1

        fields = line.split("\t")

        if len(fields) != 2:

            raise RuntimeError(
                f"aniclust línea {lineno}: "
                f"esperaba 2 columnas, "
                f"obtuve {len(fields)}"
            )


        seed = fields[0].strip()

        members = [
            x.strip()
            for x in fields[1].split(",")
            if x.strip()
        ]


        if not members:

            raise RuntimeError(
                f"Cluster vacío en línea {lineno}"
            )


        if seed not in seqs:

            raise RuntimeError(
                f"Seed desconocido: {seed}"
            )


        # El formato observado en el smoke test
        # y en la salida productiva incluye al seed
        # en la lista de miembros.
        if seed not in members:

            raise RuntimeError(
                f"Seed {seed} no aparece en su "
                f"lista de miembros."
            )


        if len(members) != len(set(members)):

            raise RuntimeError(
                f"Miembro repetido dentro del "
                f"cluster {seed}"
            )


        if seed in seed_to_members:

            raise RuntimeError(
                f"Seed repetido: {seed}"
            )


        if len(members) > 1:
            multimember_raw += 1


        for member in members:

            if member not in seqs:

                raise RuntimeError(
                    f"Miembro desconocido: {member}"
                )

            if member in all_members_seen:

                raise RuntimeError(
                    f"VCLEAN asignado a más de "
                    f"un vOTU: {member}"
                )

            all_members_seen.add(
                member
            )


        seed_to_members[
            seed
        ] = members


# ======================================================================
# VALIDACIÓN FUERTE DEL CLUSTERING
# ======================================================================

expected = set(
    seqs
)

if all_members_seen != expected:

    missing = (
        expected
        - all_members_seen
    )

    extra = (
        all_members_seen
        - expected
    )

    raise RuntimeError(
        "aniclust no cubre exactamente los "
        "6748 VCLEAN. "
        f"missing={len(missing)}, "
        f"extra={len(extra)}"
    )


membership_sum = sum(
    len(x)
    for x in seed_to_members.values()
)


if membership_sum != 6748:

    raise RuntimeError(
        f"Membership total={membership_sum}; "
        "esperaba 6748"
    )


print()
print("=" * 78)
print("ANICLUST PARSE QC")
print("=" * 78)

print(
    "raw_cluster_lines",
    raw_lines,
    sep="\t"
)

print(
    "parsed_clusters",
    len(seed_to_members),
    sep="\t"
)

print(
    "multimember_clusters_raw",
    multimember_raw,
    sep="\t"
)

print(
    "total_memberships",
    membership_sum,
    sep="\t"
)

print(
    "unique_members",
    len(all_members_seen),
    sep="\t"
)


# ======================================================================
# IDs VOTU DETERMINISTAS
# ======================================================================

seeds = sorted(
    seed_to_members
)

seed_to_votu = {
    seed: f"VOTU{i:05d}"
    for i, seed in enumerate(
        seeds,
        start=1
    )
}


# ======================================================================
# CRITERIO DE BEST QUALITY MEMBER
# ======================================================================

quality_rank = {
    "Complete": 5,
    "High-quality": 4,
    "Medium-quality": 3,
    "Low-quality": 2,
    "Not-determined": 1,
}


def best_quality_key(member):

    q = postcv[member]

    completeness = num(
        q["completeness"]
    )

    gscore = num(
        hc[member][
            "max_genomad_score"
        ]
    )

    return (
        quality_rank.get(
            q["checkv_quality"],
            0
        ),

        completeness
        if completeness is not None
        else -1,

        len(
            seqs[member]
        ),

        gscore
        if gscore is not None
        else -1,
    )


# ======================================================================
# CONSTRUIR TABLAS
# ======================================================================

membership_rows = []
votu_rows = []


for seed in seeds:

    votu = seed_to_votu[
        seed
    ]

    members = sorted(
        seed_to_members[
            seed
        ]
    )


    # En empates exactos conservamos
    # orden lexicográfico reproducible.
    bestq = sorted(
        members,
        key=lambda m: (
            best_quality_key(m),
            m
        ),
        reverse=True
    )[0]


    all_preclusters = []

    all_coassemblies = set()
    all_producers = set()
    all_taxa = set()
    all_topologies = set()

    original_predictions = 0

    scores = []
    fdrs = []
    hallmarks = []


    for member in members:

        for pre in vclean_to_pre[
            member
        ]:

            all_preclusters.append(
                pre
            )

            p = preclusters[
                pre
            ]

            pcv = precv[
                pre
            ]


            all_coassemblies.update(
                splitsemi(
                    p["coassemblies"]
                )
            )

            all_producers.update(
                splitsemi(
                    p["producers"]
                )
            )

            all_taxa.update(
                splitsemi(
                    p["terminal_taxa"]
                )
            )

            all_topologies.update(
                splitsemi(
                    p["topologies"]
                )
            )


            original_predictions += int(
                p["n_members"]
            )


            score = num(
                pcv[
                    "genomad_max_virus_score"
                ]
            )

            fdr = num(
                pcv[
                    "genomad_min_fdr"
                ]
            )

            hallmark = num(
                pcv[
                    "genomad_max_n_hallmarks"
                ]
            )


            if score is not None:
                scores.append(
                    score
                )

            if fdr is not None:
                fdrs.append(
                    fdr
                )

            if hallmark is not None:
                hallmarks.append(
                    hallmark
                )


    lengths = [
        len(
            seqs[m]
        )
        for m in members
    ]

    seed_length = len(
        seqs[seed]
    )

    max_length = max(
        lengths
    )


    qseed = postcv[
        seed
    ]

    qbest = postcv[
        bestq
    ]


    votu_rows.append({

        "votu_id":
            votu,

        "aniclust_seed":
            seed,

        "n_vclean_members":
            len(members),

        "members":
            ";".join(
                members
            ),

        "seed_length":
            seed_length,

        "max_member_length":
            max_length,

        "seed_is_longest":
            (
                "YES"
                if seed_length
                == max_length
                else "NO"
            ),

        "best_quality_member":
            bestq,

        "seed_checkv_quality":
            qseed[
                "checkv_quality"
            ],

        "seed_completeness":
            clean(
                qseed[
                    "completeness"
                ]
            ),

        "best_checkv_quality":
            qbest[
                "checkv_quality"
            ],

        "best_completeness":
            clean(
                qbest[
                    "completeness"
                ]
            ),

        "n_preclean_exact_clusters":
            len(
                set(
                    all_preclusters
                )
            ),

        "n_original_genomad_predictions":
            original_predictions,

        "n_coassemblies":
            len(
                all_coassemblies
            ),

        "coassemblies":
            ";".join(
                sorted(
                    all_coassemblies
                )
            ),

        "n_producers":
            len(
                all_producers
            ),

        "producers":
            ";".join(
                sorted(
                    all_producers
                )
            ),

        "genomad_topologies":
            ";".join(
                sorted(
                    all_topologies
                )
            ),

        "genomad_terminal_taxa":
            ";".join(
                sorted(
                    all_taxa
                )
            ),

        "max_genomad_virus_score":
            (
                max(scores)
                if scores
                else ""
            ),

        "min_genomad_fdr":
            (
                min(fdrs)
                if fdrs
                else ""
            ),

        "max_genomad_hallmarks":
            (
                int(
                    max(
                        hallmarks
                    )
                )
                if hallmarks
                else ""
            ),
    })


    for member in members:

        q = postcv[
            member
        ]

        membership_rows.append({

            "votu_id":
                votu,

            "aniclust_seed":
                seed,

            "vclean_member":
                member,

            "is_aniclust_seed":
                (
                    "YES"
                    if member == seed
                    else "NO"
                ),

            "is_best_quality_member":
                (
                    "YES"
                    if member == bestq
                    else "NO"
                ),

            "length":
                len(
                    seqs[
                        member
                    ]
                ),

            "checkv_quality":
                q[
                    "checkv_quality"
                ],

            "miuvig_quality":
                q[
                    "miuvig_quality"
                ],

            "completeness":
                clean(
                    q[
                        "completeness"
                    ]
                ),

            "completeness_method":
                clean(
                    q[
                        "completeness_method"
                    ]
                ),

            "viral_genes":
                q[
                    "viral_genes"
                ],

            "host_genes":
                q[
                    "host_genes"
                ],

            "warnings":
                clean(
                    q[
                        "warnings"
                    ]
                ),

            "preclean_exact_clusters":
                ";".join(
                    sorted(
                        vclean_to_pre[
                            member
                        ]
                    )
                ),

            "n_preclean_exact_clusters":
                len(
                    vclean_to_pre[
                        member
                    ]
                ),
        })


# ======================================================================
# VALIDACIONES DE INTEGRIDAD
# ======================================================================

if len(
    membership_rows
) != 6748:

    raise RuntimeError(
        "membership_rows != 6748"
    )


member_ids = {
    x["vclean_member"]
    for x in membership_rows
}


if member_ids != set(
    seqs
):

    raise RuntimeError(
        "Membership final no contiene "
        "exactamente todos los VCLEAN."
    )


if sum(
    int(
        x[
            "n_vclean_members"
        ]
    )
    for x in votu_rows
) != 6748:

    raise RuntimeError(
        "Suma de miembros por vOTU != 6748"
    )


# ======================================================================
# ESCRIBIR MEMBERSHIP
# ======================================================================

membership_file = (
    out
    / "92F_votu_membership.tsv"
)

with membership_file.open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=list(
            membership_rows[
                0
            ].keys()
        ),
        delimiter="\t"
    )

    writer.writeheader()

    writer.writerows(
        membership_rows
    )


# ======================================================================
# ESCRIBIR SUMMARY
# ======================================================================

summary_file = (
    out
    / "92F_votu_summary.tsv"
)

with summary_file.open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=list(
            votu_rows[
                0
            ].keys()
        ),
        delimiter="\t"
    )

    writer.writeheader()

    writer.writerows(
        votu_rows
    )


# ======================================================================
# FASTA CANÓNICO DE vOTUs
#
# Representante = seed original de aniclust.
# ======================================================================

rep_fasta = (
    out
    / "92F_votu_representatives.fna"
)

with rep_fasta.open(
    "w"
) as fh:

    for row in votu_rows:

        votu = row[
            "votu_id"
        ]

        seed = row[
            "aniclust_seed"
        ]

        write_fasta(
            fh,
            votu,
            seqs[
                seed
            ]
        )


# ======================================================================
# MAPA DE REPRESENTANTES
# ======================================================================

repmap = (
    out
    / "92F_votu_representative_map.tsv"
)

with repmap.open(
    "w"
) as fh:

    fh.write(
        "votu_id\t"
        "aniclust_seed\t"
        "best_quality_member\n"
    )

    for row in votu_rows:

        fh.write(
            f"{row['votu_id']}\t"
            f"{row['aniclust_seed']}\t"
            f"{row['best_quality_member']}\n"
        )


# ======================================================================
# RESUMEN GLOBAL
# ======================================================================

n_votu = len(
    votu_rows
)

cluster_sizes = [
    int(
        x[
            "n_vclean_members"
        ]
    )
    for x in votu_rows
]


n_singletons = sum(
    x == 1
    for x in cluster_sizes
)

n_multimember = sum(
    x > 1
    for x in cluster_sizes
)


seed_quality = Counter(
    x[
        "seed_checkv_quality"
    ]
    for x in votu_rows
)

best_quality = Counter(
    x[
        "best_checkv_quality"
    ]
    for x in votu_rows
)


global_summary = [

    (
        "input_hostclean_sequences",
        6748
    ),

    (
        "aniclust_raw_cluster_lines",
        raw_lines
    ),

    (
        "votu_total",
        n_votu
    ),

    (
        "votu_singletons",
        n_singletons
    ),

    (
        "votu_multimember",
        n_multimember
    ),

    (
        "sequences_collapsed_by_votu_clustering",
        6748 - n_votu
    ),

    (
        "largest_votu_members",
        max(
            cluster_sizes
        )
    ),

    (
        "median_votu_members",
        statistics.median(
            cluster_sizes
        )
    ),

    (
        "votu_multiple_coassemblies",
        sum(
            int(
                x[
                    "n_coassemblies"
                ]
            ) > 1
            for x in votu_rows
        )
    ),

    (
        "votu_both_producers",
        sum(
            int(
                x[
                    "n_producers"
                ]
            ) > 1
            for x in votu_rows
        )
    ),

    (
        "aniclust_seed_is_longest",
        sum(
            x[
                "seed_is_longest"
            ]
            == "YES"
            for x in votu_rows
        )
    ),

    (
        "aniclust_seed_not_longest",
        sum(
            x[
                "seed_is_longest"
            ]
            == "NO"
            for x in votu_rows
        )
    ),

    (
        "ANI_threshold_percent",
        95
    ),

    (
        "AF_shorter_threshold_percent",
        85
    ),
]


def slug(x):

    return (
        x.lower()
        .replace("-", "_")
        .replace(" ", "_")
    )


for quality, count in sorted(
    seed_quality.items()
):

    global_summary.append(
        (
            f"seed_quality_{slug(quality)}",
            count
        )
    )


for quality, count in sorted(
    best_quality.items()
):

    global_summary.append(
        (
            f"best_member_quality_{slug(quality)}",
            count
        )
    )


global_file = (
    out
    / "92F_global_summary.tsv"
)

with global_file.open(
    "w"
) as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for key, value in global_summary:

        fh.write(
            f"{key}\t{value}\n"
        )


# ======================================================================
# DISTRIBUCIÓN DE TAMAÑOS
# ======================================================================

size_counts = Counter(
    cluster_sizes
)

size_file = (
    out
    / "92F_cluster_size_distribution.tsv"
)

with size_file.open(
    "w"
) as fh:

    fh.write(
        "vclean_members\tvotu_count\n"
    )

    for size, count in sorted(
        size_counts.items()
    ):

        fh.write(
            f"{size}\t{count}\n"
        )


# ======================================================================
# PRINT
# ======================================================================

print()
print("=" * 78)
print("92F REPAIRED GLOBAL")
print("=" * 78)

for k, v in global_summary:

    print(
        k,
        v,
        sep="\t"
    )


print()
print("=" * 78)
print("20 LARGEST vOTUs")
print("=" * 78)

largest = sorted(
    votu_rows,
    key=lambda x: (
        -int(
            x[
                "n_vclean_members"
            ]
        ),
        x[
            "votu_id"
        ]
    )
)[:20]


for row in largest:

    print(
        row[
            "votu_id"
        ],
        row[
            "n_vclean_members"
        ],
        row[
            "aniclust_seed"
        ],
        row[
            "seed_length"
        ],
        row[
            "best_quality_member"
        ],
        row[
            "best_checkv_quality"
        ],
        row[
            "best_completeness"
        ],
        row[
            "n_coassemblies"
        ],
        row[
            "coassemblies"
        ],
        row[
            "n_producers"
        ],
        row[
            "producers"
        ],
        sep="\t"
    )

PY


# ======================================================================
# 5. CHECKSUMS Y METADATA
# ======================================================================

sha256sum \
    "$INPUT" \
    "$BLASTOUT" \
    "$ANIOUT" \
    "$RAWCLUST" \
    "${OUT}/92F_votu_membership.tsv" \
    "${OUT}/92F_votu_summary.tsv" \
    "${OUT}/92F_votu_representatives.fna" \
    > "${OUT}/92F_sha256.tsv"


{
    echo -e "field\tvalue"

    echo -e \
"analysis\tspecies_rank_vOTU_clustering_repaired"

    echo -e \
"original_clustering_job\t131619"

    echo -e \
"repair_scope\tintegration_only"

    echo -e \
"blast_recomputed\tNO"

    echo -e \
"anicalc_recomputed\tNO"

    echo -e \
"aniclust_recomputed\tNO"

    echo -e \
"aniclust_output_format\tseed_TAB_comma_separated_members"

    echo -e \
"input_sequences\t6748"

    echo -e \
"input_state\tCheckV_host_cleaned_exact_nonredundant"

    echo -e \
"blast\tBLAST_2.16.0"

    echo -e \
"blast_task\tmegablast"

    echo -e \
"max_target_seqs\t10000"

    echo -e \
"min_ani_percent\t95"

    echo -e \
"min_tcov_percent\t85"

    echo -e \
"min_qcov_percent\t0"

    echo -e \
"intended_AF_definition\t85_percent_shorter_sequence"

    echo -e \
"canonical_representative\taniclust_seed"

    echo -e \
"best_quality_member\tstored_separately"

    echo -e \
"quality_filter_before_clustering\tNO"

    echo -e \
"length_filter_before_clustering\tNO"

    echo -e \
"date\t$(date -Iseconds)"

} > "${OUT}/92F_run_metadata.tsv"


touch "$DONE"


# ======================================================================
# 6. VALIDACIÓN FINAL
# ======================================================================

echo
echo "======================================================================"
echo "VALIDACIÓN FINAL"
echo "======================================================================"

for F in \
    "${OUT}/92F_votu_membership.tsv" \
    "${OUT}/92F_votu_summary.tsv" \
    "${OUT}/92F_votu_representatives.fna" \
    "${OUT}/92F_votu_representative_map.tsv" \
    "${OUT}/92F_global_summary.tsv" \
    "${OUT}/92F_cluster_size_distribution.tsv" \
    "${OUT}/92F_run_metadata.tsv" \
    "${OUT}/92F_sha256.tsv" \
    "$DONE"
do

    [[ -s "$F" ]] || {
        echo "ERROR: falta $F" >&2
        exit 4
    }

    echo "OK: $(basename "$F")"

done


NVOTU="$(
    awk 'END{print NR-1}' \
        "${OUT}/92F_votu_summary.tsv"
)"

NMEMBERS="$(
    awk 'END{print NR-1}' \
        "${OUT}/92F_votu_membership.tsv"
)"

NFASTA="$(
    grep -c '^>' \
        "${OUT}/92F_votu_representatives.fna"
)"


echo
echo "vOTUs en summary      = $NVOTU"
echo "VCLEAN memberships    = $NMEMBERS"
echo "vOTUs en FASTA        = $NFASTA"


[[ "$NMEMBERS" -eq 6748 ]] || {
    echo "ERROR: membership != 6748" >&2
    exit 5
}


[[ "$NVOTU" -eq "$NFASTA" ]] || {
    echo "ERROR: summary/FASTA no coinciden." >&2
    exit 6
}


echo
echo "======================================================================"
echo "RESUMEN GLOBAL"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92F_global_summary.tsv"


echo
echo "======================================================================"
echo "DISTRIBUCIÓN DE TAMAÑOS"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92F_cluster_size_distribution.tsv"


echo
echo "======================================================================"
echo "92F-REPAIR FINALIZÓ CORRECTAMENTE."
echo "El clustering de 131619 fue reutilizado."
echo "NO se recalculó BLAST/ANI/aniclust."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "======================================================================"

