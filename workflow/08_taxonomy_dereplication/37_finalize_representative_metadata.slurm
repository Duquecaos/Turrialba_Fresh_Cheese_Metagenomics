#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

#SBATCH --job-name=finalize_rep18
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

REP="${ROOT}/22_final_representative_mags/fastas"

CDB="${ROOT}/21_drep_final24/drep/data_tables/Cdb.csv"

GTDB="${ROOT}/20_gtdbtk_final24/gtdbtk_final24_taxonomy.tsv"

PASS22="${ROOT}/14_gunc_refined/refined_candidates_90_5_GUNCpass.tsv"

POSTMDM="${ROOT}/17_checkm2_post_mdm/checkm2_original_vs_post_mdm.tsv"

MANIFEST24="${ROOT}/19_final_qc_mags/final_qc_manifest.tsv"

OUT="${ROOT}/22_final_representative_mags"

MASTER="${OUT}/representative_mags_master.tsv"

COMBINED="${OUT}/representative18_traceable.fasta"

HEADERMAP="${OUT}/representative18_contig_header_map.tsv"

CHECKSUMS="${OUT}/representative18_sha256.tsv"


echo "============================================================"
echo "FINALIZACIÓN METADATOS 18 MAGs REPRESENTATIVOS"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "============================================================"
echo


# ============================================================
# VALIDACIONES
# ============================================================

for F in \
    "${CDB}" \
    "${GTDB}" \
    "${PASS22}" \
    "${POSTMDM}" \
    "${MANIFEST24}"
do
    if [[ ! -s "${F}" ]]; then
        echo "ERROR: falta archivo requerido: ${F}" >&2
        exit 1
    fi
done


NREP=$(find "${REP}" \
    -maxdepth 1 \
    -type f \
    -name '*.fa' \
    | wc -l)

echo "MAGs representativos detectados: ${NREP}"

if [[ "${NREP}" -ne 18 ]]; then
    echo "ERROR: se esperaban 18 representantes." >&2
    exit 1
fi


# ============================================================
# GENERAR TABLA MAESTRA + FASTA TRAZABLE
# ============================================================

python3 - \
"${REP}" \
"${CDB}" \
"${GTDB}" \
"${PASS22}" \
"${POSTMDM}" \
"${MANIFEST24}" \
"${MASTER}" \
"${COMBINED}" \
"${HEADERMAP}" <<'PY'

import csv
import sys
from pathlib import Path
from collections import defaultdict

(
    rep_dir,
    cdb_file,
    gtdb_file,
    pass22_file,
    postmdm_file,
    manifest24_file,
    master_file,
    combined_file,
    headermap_file
) = sys.argv[1:10]

rep_dir = Path(rep_dir)

def strip_fasta(name):
    for ext in [".fasta.gz", ".fna.gz", ".fa.gz",
                ".fasta", ".fna", ".fa"]:
        if name.endswith(ext):
            return name[:-len(ext)]
    return name


# ============================================================
# 1. REPRESENTANTES
# ============================================================

rep_files = sorted(rep_dir.glob("*.fa"))

if len(rep_files) != 18:
    raise SystemExit(
        f"ERROR: representantes encontrados = {len(rep_files)}"
    )

representatives = {
    strip_fasta(f.name): f
    for f in rep_files
}


# ============================================================
# 2. CLUSTERS dRep
# ============================================================

clusters = defaultdict(list)

with open(cdb_file, newline="") as fh:
    reader = csv.DictReader(fh)

    for r in reader:
        genome = strip_fasta(r["genome"])
        cluster = r["secondary_cluster"]

        clusters[cluster].append(genome)

if len(clusters) != 18:
    raise SystemExit(
        f"ERROR: clusters secundarios = {len(clusters)}"
    )


rep_to_cluster = {}

for cluster, members in clusters.items():

    reps_here = [
        x for x in members
        if x in representatives
    ]

    if len(reps_here) != 1:
        raise SystemExit(
            f"ERROR cluster {cluster}: "
            f"representantes encontrados = {reps_here}"
        )

    rep_to_cluster[reps_here[0]] = cluster


# ============================================================
# 3. GTDB-Tk
# ============================================================

taxonomy = {}

with open(gtdb_file, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    for r in reader:
        taxonomy[r["user_genome"]] = r


# ============================================================
# 4. CHECKM2
# ============================================================

quality = {}

# 22 originales
with open(pass22_file, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    for r in reader:
        mag = f"{r['coassembly']}__{r['bin']}"

        quality[mag] = {
            "completeness": float(r["completeness"]),
            "contamination": float(r["contamination"])
        }


# Sustituir métricas de los 2 rescatados
rescued_names = {
    "L2_maxbin2.004_sub": "L2__L2_maxbin2.004_sub",
    "concoct_72": "L2__concoct_72",
}

with open(postmdm_file, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    for r in reader:
        if r["MAG"] in rescued_names:

            mag = rescued_names[r["MAG"]]

            quality[mag] = {
                "completeness":
                    float(r["post_completeness"]),
                "contamination":
                    float(r["post_contamination"])
            }


# ============================================================
# 5. ORIGEN
# ============================================================

source_type = {}

with open(manifest24_file, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")

    for r in reader:
        source_type[r["final_id"]] = r["source_type"]


# ============================================================
# 6. TAXONOMÍA POR RANGOS
# ============================================================

def parse_taxonomy(classification):

    ranks = {
        "d": "",
        "p": "",
        "c": "",
        "o": "",
        "f": "",
        "g": "",
        "s": "",
    }

    for item in classification.split(";"):
        if "__" not in item:
            continue

        rank, value = item.split("__", 1)

        if rank in ranks:
            ranks[rank] = value

    return ranks


# ============================================================
# 7. TABLA MAESTRA
# ============================================================

rows = []

for rep in sorted(representatives):

    if rep not in rep_to_cluster:
        raise SystemExit(
            f"ERROR: representante sin cluster: {rep}"
        )

    if rep not in taxonomy:
        raise SystemExit(
            f"ERROR: representante sin GTDB-Tk: {rep}"
        )

    if rep not in quality:
        raise SystemExit(
            f"ERROR: representante sin CheckM2: {rep}"
        )

    cluster = rep_to_cluster[rep]

    members = sorted(clusters[cluster])

    groups = sorted({
        x.split("__", 1)[0]
        for x in members
    })

    tx = taxonomy[rep]
    ranks = parse_taxonomy(tx["classification"])

    rows.append({
        "drep_cluster": cluster,
        "representative_MAG": rep,
        "n_cluster_members": len(members),
        "cluster_members": ";".join(members),
        "coassemblies_represented": ";".join(groups),
        "representative_source":
            source_type.get(rep, "UNKNOWN"),
        "completeness":
            quality[rep]["completeness"],
        "contamination":
            quality[rep]["contamination"],
        "domain": ranks["d"],
        "phylum": ranks["p"],
        "class": ranks["c"],
        "order": ranks["o"],
        "family": ranks["f"],
        "genus": ranks["g"],
        "species": ranks["s"],
        "classification":
            tx["classification"],
        "closest_genome_reference":
            tx.get("closest_genome_reference", ""),
        "closest_genome_ani":
            tx.get("closest_genome_ani", ""),
        "closest_genome_af":
            tx.get("closest_genome_af", ""),
        "classification_method":
            tx.get("classification_method", ""),
        "warnings":
            tx.get("warnings", ""),
    })


fields = [
    "drep_cluster",
    "representative_MAG",
    "n_cluster_members",
    "cluster_members",
    "coassemblies_represented",
    "representative_source",
    "completeness",
    "contamination",
    "domain",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "classification",
    "closest_genome_reference",
    "closest_genome_ani",
    "closest_genome_af",
    "classification_method",
    "warnings",
]


with open(master_file, "w", newline="") as out:
    writer = csv.DictWriter(
        out,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)


# ============================================================
# 8. FASTA CONCATENADO TRAZABLE
# ============================================================

with open(combined_file, "w") as fasta_out, \
     open(headermap_file, "w") as map_out:

    map_out.write(
        "new_contig_id\tMAG\toriginal_header\tlength_bp\n"
    )

    total_contigs = 0
    total_bp = 0

    for rep in sorted(representatives):

        fasta = representatives[rep]

        seq_id = None
        original_header = None
        sequence = []

        counter = 0

        def flush():
            global counter, total_contigs, total_bp

            if original_header is None:
                return

            counter += 1

            seq = "".join(sequence)

            new_id = (
                f"{rep}___contig{counter:06d}"
            )

            fasta_out.write(
                f">{new_id}\n"
            )

            for i in range(0, len(seq), 80):
                fasta_out.write(
                    seq[i:i+80] + "\n"
                )

            map_out.write(
                f"{new_id}\t"
                f"{rep}\t"
                f"{original_header}\t"
                f"{len(seq)}\n"
            )

            total_contigs += 1
            total_bp += len(seq)


        with open(fasta) as fh:

            for line in fh:

                line = line.rstrip("\n")

                if line.startswith(">"):

                    flush()

                    original_header = line[1:]
                    sequence = []

                else:
                    sequence.append(line.strip())

            flush()


print(f"Representantes: {len(rows)}")
print(f"Clusters dRep: {len(clusters)}")
print(f"Contigs concatenados: {total_contigs}")
print(f"bp concatenados: {total_bp}")
PY


# ============================================================
# VALIDAR TABLA
# ============================================================

N_MASTER=$(awk '
NR > 1 && NF > 0 {
    n++
}
END {
    print n+0
}' "${MASTER}")

echo
echo "Filas tabla maestra: ${N_MASTER}"

if [[ "${N_MASTER}" -ne 18 ]]; then
    echo "ERROR: tabla maestra no tiene 18 MAGs." >&2
    exit 1
fi


# ============================================================
# VALIDAR FASTA
# ============================================================

N_CONTIGS=$(grep -c '^>' "${COMBINED}" || true)

echo "Contigs FASTA combinado: ${N_CONTIGS}"

if [[ "${N_CONTIGS}" -lt 18 ]]; then
    echo "ERROR: FASTA combinado inválido." >&2
    exit 1
fi


# ============================================================
# CHECKSUMS
# ============================================================

printf "sha256\tfile\n" > "${CHECKSUMS}"

for F in "${REP}"/*.fa "${COMBINED}"; do
    SUM=$(sha256sum "${F}" | awk '{print $1}')
    printf "%s\t%s\n" "${SUM}" "${F}" \
        >> "${CHECKSUMS}"
done


# ============================================================
# RESUMEN
# ============================================================

echo
echo "===== TABLA MAESTRA ====="

column -t -s $'\t' "${MASTER}"

echo
echo "===== RESUMEN CLUSTERS ====="

awk -F $'\t' '
NR > 1 {
    n++
    members += $3

    if ($15 != "")
        species++
    else if ($14 != "")
        genus_only++
}
END {
    print "Representantes:", n
    print "MAGs originales representados:", members
    print "Con especie GTDB:", species
    print "Sin especie pero con género:", genus_only
}' "${MASTER}"


echo
echo "============================================================"
echo "REPRESENTATIVE MAG METADATA COMPLETADO"
echo "MAGs representativos: 18"
echo "Tabla maestra: ${MASTER}"
echo "FASTA trazable: ${COMBINED}"
echo "Mapa de contigs: ${HEADERMAP}"
echo "Fin: $(date)"
echo "============================================================"
