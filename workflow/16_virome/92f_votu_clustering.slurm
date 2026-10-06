#!/bin/bash

#SBATCH --job-name=92f_votu
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=92f_votu_%j.out
#SBATCH --error=92f_votu_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPT="${APPT:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"

CHECKV_SIF="${ROOT}/software/checkv/checkv-1.1.1.sif"
BLAST_SIF="${BLAST_SIF:-/opt/ohpc/pub/containers/BIO/blast-2.16.0.sif}"

ANICALC="/usr/local/lib/python3.12/site-packages/checkv/commands/anicalc.py"
ANICLUST="/usr/local/lib/python3.12/site-packages/checkv/commands/aniclust.py"

INPUT="${ROOT}/89_viral_catalog_hostclean/92E2_hostclean_exact_representatives.fna"

HOSTCLEAN_CLUSTERS="${ROOT}/89_viral_catalog_hostclean/92E2_hostclean_exact_clusters.tsv"
HOSTCLEAN_MEMBERS="${ROOT}/89_viral_catalog_hostclean/92E2_hostclean_exact_members.tsv"

PRECLUSTERS="${ROOT}/86_viral_catalog_preQC/92D_exact_sequence_clusters.tsv"
CHECKV="${ROOT}/90_viral_hostclean_validation/checkv_hostclean/quality_summary.tsv"
PRECHECKV="${ROOT}/88_checkv_viral_catalog/92E1_checkv_integrated.tsv"

OUT="${ROOT}/91_votu_catalog"
WORK="${OUT}/work"

mkdir -p "$OUT" "$WORK"

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

DBPREFIX="${WORK}/vclean_db"
BLASTOUT="${WORK}/92F_all_vs_all_blast.tsv"
ANIOUT="${WORK}/92F_pairwise_ani.tsv"
RAWCLUST="${WORK}/92F_aniclust_raw.tsv"

DONE="${OUT}/92F_COMPLETE.ok"

echo "======================================================================"
echo "PASO 92F - CLUSTERING DEFINITIVO EN vOTUs"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "======================================================================"

# ======================================================================
# 1. VALIDACIONES
# ======================================================================

for F in \
    "$CHECKV_SIF" \
    "$BLAST_SIF" \
    "$INPUT" \
    "$HOSTCLEAN_CLUSTERS" \
    "$HOSTCLEAN_MEMBERS" \
    "$PRECLUSTERS" \
    "$CHECKV" \
    "$PRECHECKV"
do
    [[ -s "$F" ]] || {
        echo "ERROR: falta $F" >&2
        exit 2
    }
done

[[ -x "$APPT" ]] || {
    echo "ERROR: falta Apptainer." >&2
    exit 3
}

NSEQ=$(grep -c '^>' "$INPUT")

echo "VCLEAN input sequences: $NSEQ"

[[ "$NSEQ" -eq 6748 ]] || {
    echo "ERROR: esperaba 6748 VCLEAN." >&2
    exit 4
}


# ======================================================================
# 2. Si existe una ejecución FINAL completa, no recalcular
# ======================================================================

if [[ -s "$DONE" ]] \
   && [[ -s "${OUT}/92F_votu_summary.tsv" ]] \
   && [[ -s "${OUT}/92F_votu_membership.tsv" ]] \
   && [[ -s "${OUT}/92F_votu_representatives.fna" ]]
then

    echo "92F ya está completo."
    echo "No se recalculará."

else

    rm -f "$DONE"

    # ==================================================================
    # 3. BLAST DATABASE
    # ==================================================================

    echo
    echo "======================================================================"
    echo "1. MAKEBLASTDB"
    echo "======================================================================"

    rm -f "${DBPREFIX}".*

    srun --exact \
        -n 1 \
        -c 1 \
        "$APPT" exec \
        --bind "${ROOT}:/work" \
        "$BLAST_SIF" \
        makeblastdb \
            -in /work/89_viral_catalog_hostclean/92E2_hostclean_exact_representatives.fna \
            -dbtype nucl \
            -out /work/91_votu_catalog/work/vclean_db


    # ==================================================================
    # 4. ALL-vs-ALL BLASTN
    # ==================================================================

    echo
    echo "======================================================================"
    echo "2. BLASTN ALL-vs-ALL"
    echo "======================================================================"

    rm -f "$BLASTOUT"

    srun --exact \
        -n 1 \
        -c "$SLURM_CPUS_PER_TASK" \
        "$APPT" exec \
        --bind "${ROOT}:/work" \
        "$BLAST_SIF" \
        blastn \
            -task megablast \
            -query /work/89_viral_catalog_hostclean/92E2_hostclean_exact_representatives.fna \
            -db /work/91_votu_catalog/work/vclean_db \
            -outfmt "6 std qlen slen" \
            -max_target_seqs 10000 \
            -num_threads "$SLURM_CPUS_PER_TASK" \
            -out /work/91_votu_catalog/work/92F_all_vs_all_blast.tsv

    [[ -s "$BLASTOUT" ]] || {
        echo "ERROR: BLAST all-vs-all vacío." >&2
        exit 5
    }

    echo "BLAST rows: $(wc -l < "$BLASTOUT")"
    du -h "$BLASTOUT"


    # ==================================================================
    # 5. ANI
    # ==================================================================

    echo
    echo "======================================================================"
    echo "3. ANICALC"
    echo "======================================================================"

    rm -f "$ANIOUT"

    srun --exact \
        -n 1 \
        -c 1 \
        "$APPT" exec \
        --bind "${ROOT}:/work" \
        "$CHECKV_SIF" \
        python3 "$ANICALC" \
            -i /work/91_votu_catalog/work/92F_all_vs_all_blast.tsv \
            -o /work/91_votu_catalog/work/92F_pairwise_ani.tsv

    [[ -s "$ANIOUT" ]] || {
        echo "ERROR: anicalc no produjo salida." >&2
        exit 6
    }

    echo "ANI rows: $(wc -l < "$ANIOUT")"
    du -h "$ANIOUT"


    # ==================================================================
    # 6. vOTU CLUSTERING
    # ==================================================================

    echo
    echo "======================================================================"
    echo "4. ANICLUST 95 ANI / 85 AF"
    echo "======================================================================"

    rm -f "$RAWCLUST"

    srun --exact \
        -n 1 \
        -c 1 \
        "$APPT" exec \
        --bind "${ROOT}:/work" \
        "$CHECKV_SIF" \
        python3 "$ANICLUST" \
            --fna /work/89_viral_catalog_hostclean/92E2_hostclean_exact_representatives.fna \
            --ani /work/91_votu_catalog/work/92F_pairwise_ani.tsv \
            --out /work/91_votu_catalog/work/92F_aniclust_raw.tsv \
            --min_ani 95 \
            --min_tcov 85 \
            --min_qcov 0

    [[ -s "$RAWCLUST" ]] || {
        echo "ERROR: aniclust no produjo salida." >&2
        exit 7
    }


    # ==================================================================
    # 7. INTEGRACIÓN Y QC
    # ==================================================================

    echo
    echo "======================================================================"
    echo "5. INTEGRACION vOTU"
    echo "======================================================================"

    srun --exact -n 1 -c 1 python3 - <<'PY'

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


# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------

def clean(x):
    if x is None:
        return ""
    x = str(x).strip()
    if x.upper() in {"", "NA", "NAN", "NONE"}:
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
                            f"Duplicate FASTA ID: {name}"
                        )
                    seqs[name] = "".join(chunks)

                name = line[1:].split()[0]
                chunks = []

            else:
                chunks.append(line)

        if name is not None:
            if name in seqs:
                raise RuntimeError(
                    f"Duplicate FASTA ID: {name}"
                )
            seqs[name] = "".join(chunks)

    return seqs


def write_fasta_record(fh, name, seq, width=80):

    fh.write(f">{name}\n")

    for i in range(0, len(seq), width):
        fh.write(seq[i:i+width] + "\n")


def read_tsv_dict(path, key):

    with path.open() as fh:

        rows = list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )

    d = {}

    for row in rows:

        k = row[key]

        if k in d:
            raise RuntimeError(
                f"Duplicate key {k} in {path}"
            )

        d[k] = row

    return d


def splitsemi(x):

    return {
        z.strip()
        for z in clean(x).split(";")
        if z.strip()
    }


# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------

seqs = read_fasta(
    fasta_file
)

if len(seqs) != 6748:
    raise RuntimeError(
        f"Expected 6748 sequences, found {len(seqs)}"
    )


hc = read_tsv_dict(
    hostclean_clusters_file,
    "hostclean_exact_cluster_id"
)

preclusters = read_tsv_dict(
    preclusters_file,
    "exact_cluster_id"
)

postcv = read_tsv_dict(
    postcheckv_file,
    "contig_id"
)

precv = read_tsv_dict(
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


# VCLEAN -> VEXACT
vclean_to_pre = defaultdict(list)

with hostclean_members_file.open() as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in r:

        vclean_to_pre[
            row["hostclean_exact_cluster_id"]
        ].append(
            row["preclean_exact_cluster_id"]
        )


if set(vclean_to_pre) != set(seqs):
    raise RuntimeError(
        "Incomplete VCLEAN -> VEXACT mapping"
    )


# ------------------------------------------------------------
# Parse raw aniclust
#
# Formato observado en smoke test:
# <cluster_seed> TAB <member>
# sin encabezado.
# ------------------------------------------------------------

seed_to_members = defaultdict(list)
member_seen = set()

with rawcluster_file.open() as fh:

    for lineno, line in enumerate(
        fh,
        start=1
    ):

        line = line.rstrip("\n")

        if not line:
            continue

        fields = line.split("\t")

        if len(fields) < 2:
            raise RuntimeError(
                f"Invalid aniclust line {lineno}: {line}"
            )

        seed = fields[0]
        member = fields[1]

        if seed not in seqs:
            raise RuntimeError(
                f"Unknown cluster seed: {seed}"
            )

        if member not in seqs:
            raise RuntimeError(
                f"Unknown cluster member: {member}"
            )

        if member in member_seen:
            raise RuntimeError(
                f"Member assigned more than once: {member}"
            )

        member_seen.add(member)
        seed_to_members[seed].append(member)


if member_seen != set(seqs):

    missing = set(seqs) - member_seen

    raise RuntimeError(
        f"Clustering does not contain all sequences; "
        f"missing={len(missing)}"
    )


if sum(
    len(x)
    for x in seed_to_members.values()
) != 6748:

    raise RuntimeError(
        "Unexpected total cluster memberships."
    )


# ------------------------------------------------------------
# Deterministic vOTU IDs
# ------------------------------------------------------------

# Ordenamos por seed para que el ID sea reproducible.
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


# ------------------------------------------------------------
# Quality ranking
# ------------------------------------------------------------

quality_rank = {
    "Complete": 5,
    "High-quality": 4,
    "Medium-quality": 3,
    "Low-quality": 2,
    "Not-determined": 1,
}


def best_quality_key(member):

    q = postcv[member]

    quality = q["checkv_quality"]

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
            quality,
            0
        ),
        (
            completeness
            if completeness is not None
            else -1
        ),
        len(seqs[member]),
        (
            gscore
            if gscore is not None
            else -1
        ),
        member,
    )


# ------------------------------------------------------------
# Build membership and vOTU summaries
# ------------------------------------------------------------

membership = []
votu_summary = []

for seed in seeds:

    votu = seed_to_votu[seed]

    members = sorted(
        seed_to_members[seed]
    )

    bestq = max(
        members,
        key=best_quality_key
    )

    all_coassemblies = set()
    all_producers = set()
    all_taxa = set()
    all_topologies = set()

    original_predictions = 0
    all_preclusters = []

    scores = []
    fdrs = []
    hallmarks = []

    for member in members:

        for pre in vclean_to_pre[member]:

            all_preclusters.append(pre)

            p = preclusters[pre]
            pcv = precv[pre]

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

            s = num(
                pcv[
                    "genomad_max_virus_score"
                ]
            )

            f = num(
                pcv[
                    "genomad_min_fdr"
                ]
            )

            h = num(
                pcv[
                    "genomad_max_n_hallmarks"
                ]
            )

            if s is not None:
                scores.append(s)

            if f is not None:
                fdrs.append(f)

            if h is not None:
                hallmarks.append(h)


    lengths = [
        len(seqs[m])
        for m in members
    ]

    seed_len = len(
        seqs[seed]
    )

    max_len = max(lengths)

    seed_is_longest = (
        seed_len == max_len
    )

    q_seed = postcv[seed]
    q_best = postcv[bestq]

    votu_summary.append({
        "votu_id":
            votu,

        "aniclust_seed":
            seed,

        "n_vclean_members":
            len(members),

        "members":
            ";".join(members),

        "seed_length":
            seed_len,

        "max_member_length":
            max_len,

        "seed_is_longest":
            "YES"
            if seed_is_longest
            else "NO",

        "best_quality_member":
            bestq,

        "seed_checkv_quality":
            q_seed[
                "checkv_quality"
            ],

        "seed_completeness":
            clean(
                q_seed[
                    "completeness"
                ]
            ),

        "best_checkv_quality":
            q_best[
                "checkv_quality"
            ],

        "best_completeness":
            clean(
                q_best[
                    "completeness"
                ]
            ),

        "n_preclean_exact_clusters":
            len(
                set(all_preclusters)
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
            max(scores)
            if scores
            else "",

        "min_genomad_fdr":
            min(fdrs)
            if fdrs
            else "",

        "max_genomad_hallmarks":
            int(max(hallmarks))
            if hallmarks
            else "",
    })


    for member in members:

        q = postcv[member]

        membership.append({
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
                    seqs[member]
                ),

            "checkv_quality":
                q["checkv_quality"],

            "miuvig_quality":
                q["miuvig_quality"],

            "completeness":
                clean(
                    q["completeness"]
                ),

            "completeness_method":
                clean(
                    q[
                        "completeness_method"
                    ]
                ),

            "viral_genes":
                q["viral_genes"],

            "host_genes":
                q["host_genes"],

            "warnings":
                clean(
                    q["warnings"]
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


# ------------------------------------------------------------
# Write membership
# ------------------------------------------------------------

membership_file = (
    out
    / "92F_votu_membership.tsv"
)

with membership_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            membership[0].keys()
        ),
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(membership)


# ------------------------------------------------------------
# Write summary
# ------------------------------------------------------------

summary_file = (
    out
    / "92F_votu_summary.tsv"
)

with summary_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            votu_summary[0].keys()
        ),
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(votu_summary)


# ------------------------------------------------------------
# Canonical representative FASTA
# = aniclust seed
# ------------------------------------------------------------

rep_fasta = (
    out
    / "92F_votu_representatives.fna"
)

with rep_fasta.open("w") as fh:

    for row in votu_summary:

        votu = row["votu_id"]
        seed = row["aniclust_seed"]

        write_fasta_record(
            fh,
            votu,
            seqs[seed]
        )


# ------------------------------------------------------------
# Map canonical VOTU -> VCLEAN seed
# ------------------------------------------------------------

repmap_file = (
    out
    / "92F_votu_representative_map.tsv"
)

with repmap_file.open("w") as fh:

    fh.write(
        "votu_id\taniclust_seed\t"
        "best_quality_member\n"
    )

    for row in votu_summary:

        fh.write(
            f"{row['votu_id']}\t"
            f"{row['aniclust_seed']}\t"
            f"{row['best_quality_member']}\n"
        )


# ------------------------------------------------------------
# Global QC
# ------------------------------------------------------------

n_votu = len(votu_summary)

n_single = sum(
    int(x["n_vclean_members"]) == 1
    for x in votu_summary
)

n_multi = sum(
    int(x["n_vclean_members"]) > 1
    for x in votu_summary
)

n_cross_coassembly = sum(
    int(x["n_coassemblies"]) > 1
    for x in votu_summary
)

n_cross_producer = sum(
    int(x["n_producers"]) > 1
    for x in votu_summary
)

seed_longest = sum(
    x["seed_is_longest"] == "YES"
    for x in votu_summary
)

seed_quality_counts = Counter(
    x["seed_checkv_quality"]
    for x in votu_summary
)

best_quality_counts = Counter(
    x["best_checkv_quality"]
    for x in votu_summary
)

cluster_sizes = [
    int(
        x["n_vclean_members"]
    )
    for x in votu_summary
]


global_summary = [
    (
        "input_hostclean_sequences",
        6748
    ),

    (
        "votu_total",
        n_votu
    ),

    (
        "votu_singletons",
        n_single
    ),

    (
        "votu_multimember",
        n_multi
    ),

    (
        "sequences_collapsed_by_votu_clustering",
        6748 - n_votu
    ),

    (
        "largest_votu_members",
        max(cluster_sizes)
    ),

    (
        "median_votu_members",
        statistics.median(
            cluster_sizes
        )
    ),

    (
        "votu_multiple_coassemblies",
        n_cross_coassembly
    ),

    (
        "votu_both_producers",
        n_cross_producer
    ),

    (
        "aniclust_seed_is_longest",
        seed_longest
    ),

    (
        "aniclust_seed_not_longest",
        n_votu - seed_longest
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


for key, value in sorted(
    seed_quality_counts.items()
):
    global_summary.append(
        (
            f"seed_quality_{slug(key)}",
            value
        )
    )


for key, value in sorted(
    best_quality_counts.items()
):
    global_summary.append(
        (
            f"best_member_quality_{slug(key)}",
            value
        )
    )


global_file = (
    out
    / "92F_global_summary.tsv"
)

with global_file.open("w") as fh:

    fh.write("metric\tvalue\n")

    for key, value in global_summary:
        fh.write(
            f"{key}\t{value}\n"
        )


# ------------------------------------------------------------
# Cluster-size distribution
# ------------------------------------------------------------

size_counts = Counter(
    cluster_sizes
)

size_file = (
    out
    / "92F_cluster_size_distribution.tsv"
)

with size_file.open("w") as fh:

    fh.write(
        "vclean_members\tvotu_count\n"
    )

    for size, count in sorted(
        size_counts.items()
    ):

        fh.write(
            f"{size}\t{count}\n"
        )


# ------------------------------------------------------------
# Strong validation
# ------------------------------------------------------------

if len(membership) != 6748:
    raise RuntimeError(
        f"Membership count={len(membership)}"
    )

if len({
    x["vclean_member"]
    for x in membership
}) != 6748:

    raise RuntimeError(
        "VCLEAN membership not unique."
    )

if len({
    x["votu_id"]
    for x in votu_summary
}) != n_votu:

    raise RuntimeError(
        "Duplicated VOTU IDs."
    )

if sum(
    int(x["n_vclean_members"])
    for x in votu_summary
) != 6748:

    raise RuntimeError(
        "VOTU membership sum != 6748."
    )


print()
print("=" * 78)
print("92F GLOBAL")
print("=" * 78)

for k, v in global_summary:
    print(k, v, sep="\t")

print()
print("=" * 78)
print("LARGEST vOTUs")
print("=" * 78)

for row in sorted(
    votu_summary,
    key=lambda x: (
        -int(
            x["n_vclean_members"]
        ),
        x["votu_id"]
    )
)[:20]:

    print(
        row["votu_id"],
        row["n_vclean_members"],
        row["aniclust_seed"],
        row["seed_length"],
        row["best_checkv_quality"],
        row["best_completeness"],
        row["n_coassemblies"],
        row["coassemblies"],
        row["n_producers"],
        row["producers"],
        sep="\t"
    )

PY


    # ==================================================================
    # 8. CHECKSUMS / METADATA
    # ==================================================================

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
        echo -e "analysis\tspecies_rank_vOTU_clustering"
        echo -e "input_sequences\t6748"
        echo -e "input_state\tCheckV_host_cleaned_exact_nonredundant"
        echo -e "blast\tBLAST_2.16.0"
        echo -e "blast_task\tmegablast"
        echo -e "max_target_seqs\t10000"
        echo -e "anicalc\tCheckV_1.1.1_supporting_code"
        echo -e "aniclust\tCheckV_1.1.1_supporting_code"
        echo -e "min_ani_percent\t95"
        echo -e "min_tcov_percent\t85"
        echo -e "min_qcov_percent\t0"
        echo -e "intended_AF_definition\t85_percent_shorter_sequence"
        echo -e "canonical_representative\taniclust_seed"
        echo -e "best_quality_member\tstored_separately_in_summary"
        echo -e "quality_filter_before_clustering\tNO"
        echo -e "length_filter_before_clustering\tNO"
        echo -e "date\t$(date -Iseconds)"
    } > "${OUT}/92F_run_metadata.tsv"


    touch "$DONE"

fi


# ======================================================================
# 9. VALIDACIÓN FINAL
# ======================================================================

echo
echo "======================================================================"
echo "VALIDACION FINAL"
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
        echo "ERROR: falta salida: $F" >&2
        exit 8
    }

    echo "OK: $(basename "$F")"

done


echo
echo "======================================================================"
echo "RESUMEN GLOBAL"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92F_global_summary.tsv"


echo
echo "======================================================================"
echo "DISTRIBUCION DE TAMAÑOS"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92F_cluster_size_distribution.tsv" \
    | head -n 30


echo
echo "======================================================================"
echo "PASO 92F FINALIZÓ CORRECTAMENTE."
echo "Los vOTUs species-rank quedaron definidos."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "======================================================================"

