#!/bin/bash

#SBATCH --job-name=92e2_viralclean
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=01:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=92e2_viralclean_%j.out
#SBATCH --error=92e2_viralclean_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

PRE="${ROOT}/86_viral_catalog_preQC"
CHECKV="${ROOT}/88_checkv_viral_catalog/checkv_exact_representatives"
INT="${ROOT}/88_checkv_viral_catalog/92E1_checkv_integrated.tsv"

OUT="${ROOT}/89_viral_catalog_hostclean"

mkdir -p "$OUT"

INPUT="${PRE}/92D_exact_representatives.fna"
CONTAM="${CHECKV}/contamination.tsv"

for F in \
    "$INPUT" \
    "$CONTAM" \
    "$INT"
do
    [[ -s "$F" ]] || {
        echo "ERROR: falta $F" >&2
        exit 2
    }
done

echo "======================================================================"
echo "PASO 92E2 - LIMPIEZA DE REGIONES HOSPEDADORAS"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"

srun --exact -n 1 -c 1 python3 - <<'PY'

import csv
import hashlib
import os
import statistics
from collections import Counter, defaultdict
from pathlib import Path

root = Path(
    f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"
)

pre = root / "86_viral_catalog_preQC"
checkv = root / "88_checkv_viral_catalog/checkv_exact_representatives"
integrated_file = root / "88_checkv_viral_catalog/92E1_checkv_integrated.tsv"

out = root / "89_viral_catalog_hostclean"

input_fasta = pre / "92D_exact_representatives.fna"
contam_file = checkv / "contamination.tsv"


# ============================================================
# Helpers
# ============================================================

def clean(x):
    if x is None:
        return ""

    x = str(x).strip()

    if x.upper() in {"", "NA", "NAN", "NONE"}:
        return ""

    return x


def fnum(x):
    x = clean(x)
    return float(x) if x else None


def inum(x):
    x = clean(x)
    return int(float(x)) if x else None


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


def write_record(fh, name, seq, width=80):

    fh.write(f">{name}\n")

    for i in range(0, len(seq), width):
        fh.write(seq[i:i+width] + "\n")


def sha(seq):

    return hashlib.sha256(
        seq.encode("ascii")
    ).hexdigest()


# ============================================================
# 1. Inputs
# ============================================================

seqs = read_fasta(input_fasta)

if len(seqs) != 6751:
    raise RuntimeError(
        f"Esperaba 6751 FASTA; encontré {len(seqs)}"
    )


with integrated_file.open() as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    integrated = {
        x["exact_cluster_id"]: x
        for x in r
    }


if set(seqs) != set(integrated):

    raise RuntimeError(
        "Los IDs del FASTA y 92E1 no coinciden."
    )


with contam_file.open() as fh:

    r = csv.DictReader(
        fh,
        delimiter="\t"
    )

    contamination = {
        x["contig_id"]: x
        for x in r
    }


if set(seqs) != set(contamination):

    raise RuntimeError(
        "Los IDs FASTA/contamination.tsv no coinciden."
    )


# ============================================================
# 2. Reconstrucción de región viral
# ============================================================

cleaned = {}
lineage = []

n_trimmed = 0
total_host_bp_removed = 0


for cid in sorted(seqs):

    original = seqs[cid]

    q = integrated[cid]
    c = contamination[cid]

    checkv_provirus = (
        clean(c.get("provirus")).lower()
        == "yes"
    )

    contamination_pct = fnum(
        q.get("contamination")
    )

    if contamination_pct is None:
        contamination_pct = 0.0

    if not checkv_provirus:

        viral_seq = original

        coords_used = (
            f"1-{len(original)}"
        )

        action = "UNCHANGED"

    else:

        region_types = [
            x.strip().lower()
            for x in clean(
                c.get("region_types")
            ).split(",")
            if x.strip()
        ]

        region_coords = [
            x.strip()
            for x in clean(
                c.get("region_coords_bp")
            ).split(",")
            if x.strip()
        ]

        if not region_types or not region_coords:

            raise RuntimeError(
                f"{cid}: CheckV provirus Yes "
                "pero faltan regiones."
            )

        if len(region_types) != len(region_coords):

            raise RuntimeError(
                f"{cid}: region_types y coordinates "
                "tienen diferente longitud."
            )

        viral_regions = []

        for typ, coords in zip(
            region_types,
            region_coords
        ):

            if typ != "viral":
                continue

            if "-" not in coords:

                raise RuntimeError(
                    f"{cid}: coordenadas inválidas: {coords}"
                )

            start, end = map(
                int,
                coords.split("-", 1)
            )

            if (
                start < 1
                or end > len(original)
                or end < start
            ):

                raise RuntimeError(
                    f"{cid}: coordenadas fuera de rango "
                    f"{start}-{end}/{len(original)}"
                )

            # CheckV usa coordenadas 1-based inclusivas.
            viral_regions.append(
                (
                    start,
                    end,
                    original[
                        start-1:end
                    ]
                )
            )

        if len(viral_regions) != 1:

            raise RuntimeError(
                f"{cid}: esperaba una región viral "
                f"contigua; encontró {len(viral_regions)}"
            )

        start, end, viral_seq = viral_regions[0]

        expected = inum(
            c.get("proviral_length")
        )

        if expected is None:

            raise RuntimeError(
                f"{cid}: proviral_length ausente"
            )

        if len(viral_seq) != expected:

            raise RuntimeError(
                f"{cid}: longitud extraída "
                f"{len(viral_seq)} != CheckV "
                f"proviral_length {expected}"
            )

        coords_used = f"{start}-{end}"
        action = "HOST_TRIMMED"

        n_trimmed += 1

        total_host_bp_removed += (
            len(original)
            - len(viral_seq)
        )


    if not viral_seq:

        raise RuntimeError(
            f"{cid}: secuencia limpia vacía"
        )

    cleaned[cid] = viral_seq

    lineage.append({
        "preclean_exact_cluster_id":
            cid,

        "original_length":
            len(original),

        "cleaned_length":
            len(viral_seq),

        "bp_removed":
            len(original) - len(viral_seq),

        "checkv_provirus":
            q["checkv_provirus"],

        "checkv_contamination":
            q["contamination"],

        "checkv_proviral_length":
            q["checkv_proviral_length"],

        "checkv_region_coords_used":
            coords_used,

        "cleaning_action":
            action,

        "checkv_quality":
            q["checkv_quality"],

        "miuvig_quality":
            q["miuvig_quality"],

        "completeness":
            q["completeness"],

        "completeness_method":
            q["completeness_method"],

        "genomad_topologies":
            q["genomad_topologies"],

        "genomad_max_virus_score":
            q["genomad_max_virus_score"],

        "genomad_min_fdr":
            q["genomad_min_fdr"],

        "genomad_max_n_hallmarks":
            q["genomad_max_n_hallmarks"],

        "sequence_sha256_after_cleaning":
            sha(viral_seq),
    })


if n_trimmed != 186:

    raise RuntimeError(
        f"Esperaba 186 secuencias CheckV provirus; "
        f"se limpiaron {n_trimmed}"
    )


# ============================================================
# 3. Escribir FASTA limpio previo a redereplicación
# ============================================================

all_clean = (
    out
    / "92E2_hostclean_all_6751.fna"
)

with all_clean.open("w") as fh:

    for cid in sorted(cleaned):

        write_record(
            fh,
            cid,
            cleaned[cid]
        )


# ============================================================
# 4. Trazabilidad preclean -> host-clean
# ============================================================

lineage_file = (
    out
    / "92E2_hostclean_lineage.tsv"
)

with lineage_file.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=list(
            lineage[0].keys()
        ),
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(lineage)


# ============================================================
# 5. Redereplicación EXACTA post-limpieza
# ============================================================

by_hash = defaultdict(list)

for row in lineage:

    by_hash[
        row[
            "sequence_sha256_after_cleaning"
        ]
    ].append(row)


cluster_rows = []
member_rows = []

rep_fasta = (
    out
    / "92E2_hostclean_exact_representatives.fna"
)


with rep_fasta.open("w") as fh:

    for i, h in enumerate(
        sorted(by_hash),
        start=1
    ):

        clean_id = (
            f"VCLEAN{i:06d}"
        )

        members = sorted(
            by_hash[h],
            key=lambda x:
                x[
                    "preclean_exact_cluster_id"
                ]
        )

        rep = members[0]

        rep_old = (
            rep[
                "preclean_exact_cluster_id"
            ]
        )

        seq = cleaned[rep_old]

        # SHA collision check
        for m in members:

            mid = (
                m[
                    "preclean_exact_cluster_id"
                ]
            )

            if cleaned[mid] != seq:

                raise RuntimeError(
                    "SHA256 collision"
                )

        write_record(
            fh,
            clean_id,
            seq
        )


        qualities = sorted({
            m["checkv_quality"]
            for m in members
            if m["checkv_quality"]
        })

        topologies = sorted({
            m["genomad_topologies"]
            for m in members
            if m["genomad_topologies"]
        })

        completions = [
            fnum(m["completeness"])
            for m in members
        ]

        completions = [
            x
            for x in completions
            if x is not None
        ]


        cluster_rows.append({
            "hostclean_exact_cluster_id":
                clean_id,

            "sequence_sha256":
                h,

            "length":
                len(seq),

            "n_preclean_clusters":
                len(members),

            "preclean_clusters":
                ";".join(
                    m[
                        "preclean_exact_cluster_id"
                    ]
                    for m in members
                ),

            "representative_preclean_cluster":
                rep_old,

            "n_trimmed_members":
                sum(
                    m["cleaning_action"]
                    == "HOST_TRIMMED"
                    for m in members
                ),

            "checkv_qualities":
                ";".join(qualities),

            "genomad_topologies":
                ";".join(topologies),

            "max_completeness":
                (
                    max(completions)
                    if completions
                    else ""
                ),

            "max_genomad_score":
                max(
                    [
                        fnum(
                            m[
                                "genomad_max_virus_score"
                            ]
                        )
                        for m in members
                        if fnum(
                            m[
                                "genomad_max_virus_score"
                            ]
                        )
                        is not None
                    ]
                    or [0]
                ),

            "min_genomad_fdr":
                min(
                    [
                        fnum(
                            m[
                                "genomad_min_fdr"
                            ]
                        )
                        for m in members
                        if fnum(
                            m[
                                "genomad_min_fdr"
                            ]
                        )
                        is not None
                    ]
                    or [1]
                ),
        })


        for m in members:

            member_rows.append({
                "hostclean_exact_cluster_id":
                    clean_id,

                "preclean_exact_cluster_id":
                    m[
                        "preclean_exact_cluster_id"
                    ],

                "cleaning_action":
                    m["cleaning_action"],

                "original_length":
                    m["original_length"],

                "cleaned_length":
                    m["cleaned_length"],

                "bp_removed":
                    m["bp_removed"],

                "checkv_quality":
                    m["checkv_quality"],

                "completeness":
                    m["completeness"],

                "genomad_topologies":
                    m["genomad_topologies"],
            })


cluster_file = (
    out
    / "92E2_hostclean_exact_clusters.tsv"
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


member_file = (
    out
    / "92E2_hostclean_exact_members.tsv"
)

with member_file.open(
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
# 6. Distribución de longitudes después de limpieza
# ============================================================

length_counts = Counter()

for x in cluster_rows:

    L = int(x["length"])

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

    length_counts[key] += 1


# ============================================================
# 7. Resumen
# ============================================================

trimmed_rows = [
    x
    for x in lineage
    if x["cleaning_action"]
    == "HOST_TRIMMED"
]

removed_values = [
    int(x["bp_removed"])
    for x in trimmed_rows
]


summary = [
    (
        "preclean_exact_sequences",
        6751
    ),

    (
        "checkv_host_trimmed_sequences",
        n_trimmed
    ),

    (
        "unchanged_sequences",
        6751 - n_trimmed
    ),

    (
        "total_host_bp_removed",
        total_host_bp_removed
    ),

    (
        "median_bp_removed_among_trimmed",
        statistics.median(
            removed_values
        )
        if removed_values
        else 0
    ),

    (
        "max_bp_removed",
        max(removed_values)
        if removed_values
        else 0
    ),

    (
        "postclean_exact_clusters",
        len(cluster_rows)
    ),

    (
        "postclean_exact_singletons",
        sum(
            x["n_preclean_clusters"] == 1
            for x in cluster_rows
        )
    ),

    (
        "postclean_exact_multimember_clusters",
        sum(
            x["n_preclean_clusters"] > 1
            for x in cluster_rows
        )
    ),

    (
        "new_exact_collapses_after_trimming",
        6751 - len(cluster_rows)
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
            f"postclean_length_{key}",
            length_counts[key]
        )
    )


summary_file = (
    out
    / "92E2_hostclean_summary.tsv"
)

with summary_file.open("w") as fh:

    fh.write("metric\tvalue\n")

    for k, v in summary:
        fh.write(f"{k}\t{v}\n")


# ============================================================
# 8. QC
# ============================================================

if len(lineage) != 6751:

    raise RuntimeError(
        "Trazabilidad incompleta."
    )


if sum(
    x["n_preclean_clusters"]
    for x in cluster_rows
) != 6751:

    raise RuntimeError(
        "Se perdieron preclean clusters."
    )


if len({
    x["hostclean_exact_cluster_id"]
    for x in cluster_rows
}) != len(cluster_rows):

    raise RuntimeError(
        "IDs VCLEAN duplicados."
    )


print()
print("=" * 78)
print("92E2 SUMMARY")
print("=" * 78)

for k, v in summary:
    print(k, v, sep="\t")

PY


# ======================================================================
# Metadata
# ======================================================================

{
    echo -e "field\tvalue"
    echo -e "analysis\tCheckV_host_region_cleaning"
    echo -e "input_sequences\t6751"
    echo -e "host_boundary_source\tCheckV_1.1.1_contamination.tsv"
    echo -e "coordinate_system\t1_based_inclusive"
    echo -e "sequence_filtering_by_quality\tNO"
    echo -e "sequence_filtering_by_length\tNO"
    echo -e "postclean_dereplication\tSHA256_exact_sequence"
    echo -e "vOTU_clustering\tNOT_PERFORMED"
    echo -e "date\t$(date -Iseconds)"
} > "${OUT}/92E2_run_metadata.tsv"


# ======================================================================
# Validación
# ======================================================================

echo
echo "======================================================================"
echo "VALIDACION"
echo "======================================================================"

for F in \
    "${OUT}/92E2_hostclean_all_6751.fna" \
    "${OUT}/92E2_hostclean_lineage.tsv" \
    "${OUT}/92E2_hostclean_exact_representatives.fna" \
    "${OUT}/92E2_hostclean_exact_clusters.tsv" \
    "${OUT}/92E2_hostclean_exact_members.tsv" \
    "${OUT}/92E2_hostclean_summary.tsv" \
    "${OUT}/92E2_run_metadata.tsv"
do

    [[ -s "$F" ]] || {
        echo "ERROR: falta salida: $F" >&2
        exit 3
    }

    echo "OK: $(basename "$F")"

done


echo
echo "======================================================================"
echo "RESUMEN"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92E2_hostclean_summary.tsv"


echo
echo "======================================================================"
echo "PASO 92E2 FINALIZÓ CORRECTAMENTE."
echo "No se aplicó filtro de longitud ni calidad."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "======================================================================"

