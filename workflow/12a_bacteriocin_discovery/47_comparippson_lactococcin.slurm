#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

#SBATCH --job-name=comparippson_lcn
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

RESULTS="${ROOT}/26_antismash_rep18/results"
TRUEHITS="${ROOT}/28_antismash_true_hits/antismash_true_bacteriocin_hits.tsv"

DB="${ROOT}/databases/antismash_7.1.0"

OUT="${ROOT}/29_comparippson_lactococcin"

ANTISMASH="${ANTISMASH:-/opt/ohpc/pub/containers/BIO/antismash-7.1.0.sif}"

mkdir -p "${OUT}"

echo "============================================================"
echo "CompaRiPPson manual - candidatos Lactococcin"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "============================================================"

if [[ ! -s "${TRUEHITS}" ]]; then
    echo "ERROR: falta ${TRUEHITS}" >&2
    exit 1
fi

for F in \
    "${DB}/comparippson/mibig/3.1/cores.fa" \
    "${DB}/comparippson/mibig/3.1/metadata.json" \
    "${DB}/comparippson/asdb/4.0/cores.fa" \
    "${DB}/comparippson/asdb/4.0/metadata.json"
do
    if [[ ! -s "${F}" ]]; then
        echo "ERROR: falta ${F}" >&2
        exit 1
    fi
done


# ============================================================
# EXTRAER LOS CDS Lactococcin Y PROPONER CORE TRAS GG
# ============================================================

echo
echo "============================================================"
echo "EXTRAYENDO CANDIDATOS"
echo "============================================================"

srun --exact -n 1 -c 1 \
apptainer exec \
    --cleanenv \
    --bind "${ROOT}:/project:ro" \
    --bind "${OUT}:/out" \
    "${ANTISMASH}" \
    python3 - <<'PY'

import csv
from pathlib import Path
from Bio import SeqIO

ROOT = Path("/project")
TRUEHITS = (
    ROOT /
    "28_antismash_true_hits" /
    "antismash_true_bacteriocin_hits.tsv"
)
RESULTS = ROOT / "26_antismash_rep18" / "results"
OUT = Path("/out")

candidates = []

with TRUEHITS.open(newline="") as fh:

    reader = csv.DictReader(
        fh,
        delimiter="\t"
    )

    for row in reader:

        if row["detected_profile"] != "Lactococcin":
            continue

        candidates.append(row)

if len(candidates) != 5:
    raise SystemExit(
        f"ERROR: se esperaban 5 hits Lactococcin "
        f"y se encontraron {len(candidates)}"
    )

full_records = []
core_records = []
table = []

for row in candidates:

    mag = row["MAG"]
    locus = row["locus_tag"]
    region_file = row["region_file"]

    gbk = RESULTS / mag / region_file

    record = SeqIO.read(
        str(gbk),
        "genbank"
    )

    feature = None

    for f in record.features:

        if f.type != "CDS":
            continue

        tags = f.qualifiers.get(
            "locus_tag",
            []
        )

        if locus in tags:
            feature = f
            break

    if feature is None:
        raise RuntimeError(
            f"No se encontró {locus} en {gbk}"
        )

    translation = "".join(
        feature.qualifiers.get(
            "translation",
            []
        )
    )

    if not translation:
        raise RuntimeError(
            f"{mag} {locus} sin traducción"
        )

    # Buscar el primer GG compatible con un leader:
    # no demasiado cerca del inicio y dejando un core >= 8 aa.
    gg_pos = None

    for i in range(
        8,
        len(translation) - 9
    ):

        if translation[i:i+2] == "GG":
            gg_pos = i
            break

    if gg_pos is not None:

        leader = translation[:gg_pos + 2]
        core = translation[gg_pos + 2:]
        method = "after_first_GG"

    else:

        leader = ""
        core = translation
        method = "full_ORF_no_suitable_GG"

    query_id = (
        f"{mag}__{locus}"
    )

    full_records.append(
        f">{query_id}\n"
        f"{translation}\n"
    )

    core_records.append(
        f">{query_id}\n"
        f"{core}\n"
    )

    table.append({
        "query_id": query_id,
        "MAG": mag,
        "genus": row["genus"],
        "species": row["species"],
        "contig": row["contig"],
        "region": row["region"],
        "locus_tag": locus,
        "contig_edge": row["contig_edge"],
        "full_translation": translation,
        "full_length_aa": len(translation),
        "core_method": method,
        "leader_candidate": leader,
        "leader_length_aa": len(leader),
        "core_candidate": core,
        "core_length_aa": len(core),
    })


with (
    OUT /
    "lactococcin_full_orfs.faa"
).open("w") as fh:
    fh.writelines(full_records)


with (
    OUT /
    "lactococcin_core_candidates.faa"
).open("w") as fh:
    fh.writelines(core_records)


fields = [
    "query_id",
    "MAG",
    "genus",
    "species",
    "contig",
    "region",
    "locus_tag",
    "contig_edge",
    "full_translation",
    "full_length_aa",
    "core_method",
    "leader_candidate",
    "leader_length_aa",
    "core_candidate",
    "core_length_aa",
]

with (
    OUT /
    "lactococcin_candidate_sequences.tsv"
).open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(table)

print("Candidatos Lactococcin:", len(table))

for r in table:
    print(
        r["query_id"],
        "full=", r["full_translation"],
        "leader=", r["leader_candidate"],
        "core=", r["core_candidate"],
        "method=", r["core_method"],
        sep="\t"
    )

PY


# ============================================================
# BLASTP CON EL MISMO FORMATO DE CompaRiPPson 7.1
# ============================================================

BLASTFMT="6 qacc sacc nident qseq qstart qend qlen sseq sstart send slen"

echo
echo "============================================================"
echo "BLASTP - CORES"
echo "============================================================"

for DATABASE in mibig asdb; do

    if [[ "${DATABASE}" == "mibig" ]]; then
        DBPATH="/db/comparippson/mibig/3.1/cores.fa"
    else
        DBPATH="/db/comparippson/asdb/4.0/cores.fa"
    fi

    echo "Base: ${DATABASE}"

    srun --exact -n 1 -c 1 \
    apptainer exec \
        --cleanenv \
        --bind "${DB}:/db:ro" \
        --bind "${OUT}:/out" \
        "${ANTISMASH}" \
        blastp \
            -query /out/lactococcin_core_candidates.faa \
            -db "${DBPATH}" \
            -outfmt "${BLASTFMT}" \
        > "${OUT}/${DATABASE}_cores.raw.tsv"

done


echo
echo "============================================================"
echo "BLASTP - ORF COMPLETO, CONTROL EXPLORATORIO"
echo "============================================================"

for DATABASE in mibig asdb; do

    if [[ "${DATABASE}" == "mibig" ]]; then
        DBPATH="/db/comparippson/mibig/3.1/cores.fa"
    else
        DBPATH="/db/comparippson/asdb/4.0/cores.fa"
    fi

    srun --exact -n 1 -c 1 \
    apptainer exec \
        --cleanenv \
        --bind "${DB}:/db:ro" \
        --bind "${OUT}:/out" \
        "${ANTISMASH}" \
        blastp \
            -query /out/lactococcin_full_orfs.faa \
            -db "${DBPATH}" \
            -outfmt "${BLASTFMT}" \
        > "${OUT}/${DATABASE}_full.raw.tsv"

done


# ============================================================
# REPRODUCIR SCORING CompaRiPPson + METADATOS
# ============================================================

echo
echo "============================================================"
echo "INTERPRETANDO RESULTADOS"
echo "============================================================"

srun --exact -n 1 -c 1 \
apptainer exec \
    --cleanenv \
    --bind "${DB}:/db:ro" \
    --bind "${OUT}:/out" \
    "${ANTISMASH}" \
    python3 - <<'PY'

import csv
import json
from pathlib import Path
from collections import defaultdict

OUT = Path("/out")
DBROOT = Path("/db")

DATABASES = {
    "mibig": {
        "version": "3.1",
        "metadata":
            DBROOT /
            "comparippson/mibig/3.1/metadata.json",
    },
    "asdb": {
        "version": "4.0",
        "metadata":
            DBROOT /
            "comparippson/asdb/4.0/metadata.json",
    },
}


def build_text(template, fields, field_names):

    result = template

    for field in field_names:
        result = result.replace(
            f"@{field}@",
            str(fields.get(field, ""))
        )

    return result


all_rows = []

for dbname, info in DATABASES.items():

    with info["metadata"].open() as fh:
        meta = json.load(fh)

    entries = meta["entries"]
    fields = meta["fields"]

    for mode in ["cores", "full"]:

        raw = (
            OUT /
            f"{dbname}_{mode}.raw.tsv"
        )

        if not raw.exists():
            continue

        with raw.open() as fh:

            for line in fh:

                line = line.rstrip("\n")

                if not line:
                    continue

                parts = line.split("\t")

                if len(parts) != 11:
                    raise RuntimeError(
                        f"Línea BLAST inesperada: {line}"
                    )

                (
                    qacc,
                    sacc,
                    nident,
                    qseq,
                    qstart,
                    qend,
                    qlen,
                    sseq,
                    sstart,
                    send,
                    slen,
                ) = parts

                nident = int(nident)
                qlen = int(qlen)
                slen = int(slen)

                qstart = int(qstart)
                qend = int(qend)
                sstart = int(sstart)
                send = int(send)

                unique_id = (
                    sacc.split("|", 1)[0]
                )

                ref = entries.get(
                    unique_id,
                    {}
                )

                # Métrica EXACTA usada por CompaRiPPson
                similarity = (
                    nident /
                    max(qlen, slen)
                )

                # antiSMASH default:
                # hit.similarity > 0.05
                if similarity <= 0.05:
                    continue

                aln_len = len(
                    qseq.replace("-", "")
                )

                identity_aln = (
                    100 * nident /
                    max(1, len(qseq))
                )

                qcov = (
                    100 *
                    (qend - qstart + 1) /
                    qlen
                )

                scov = (
                    100 *
                    (send - sstart + 1) /
                    slen
                )

                identifier = build_text(
                    meta.get(
                        "id_format",
                        unique_id
                    ),
                    ref,
                    fields
                )

                description = build_text(
                    meta.get(
                        "description_format",
                        ""
                    ),
                    ref,
                    fields
                )

                all_rows.append({
                    "query": qacc,
                    "query_mode": mode,
                    "database":
                        meta.get(
                            "name",
                            dbname
                        ),
                    "database_version":
                        meta.get(
                            "version",
                            info["version"]
                        ),
                    "subject_unique_id":
                        unique_id,
                    "reference_identifier":
                        identifier,
                    "reference_description":
                        description,
                    "comparippson_similarity_pct":
                        round(
                            similarity * 100,
                            2
                        ),
                    "n_identical":
                        nident,
                    "query_length":
                        qlen,
                    "reference_length":
                        slen,
                    "alignment_identity_pct":
                        round(
                            identity_aln,
                            2
                        ),
                    "query_coverage_pct":
                        round(
                            qcov,
                            2
                        ),
                    "reference_coverage_pct":
                        round(
                            scov,
                            2
                        ),
                    "query_aligned":
                        qseq,
                    "reference_aligned":
                        sseq,
                    "reference_metadata":
                        json.dumps(
                            ref,
                            sort_keys=True
                        ),
                })


all_rows.sort(
    key=lambda x: (
        x["query"],
        x["query_mode"],
        x["database"],
        -x["comparippson_similarity_pct"],
        x["reference_identifier"],
    )
)


fields = [
    "query",
    "query_mode",
    "database",
    "database_version",
    "subject_unique_id",
    "reference_identifier",
    "reference_description",
    "comparippson_similarity_pct",
    "n_identical",
    "query_length",
    "reference_length",
    "alignment_identity_pct",
    "query_coverage_pct",
    "reference_coverage_pct",
    "query_aligned",
    "reference_aligned",
    "reference_metadata",
]


with (
    OUT /
    "comparippson_lactococcin_all_hits.tsv"
).open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(all_rows)


# Top 10 por query / modo / base
grouped = defaultdict(list)

for row in all_rows:
    grouped[
        (
            row["query"],
            row["query_mode"],
            row["database"],
        )
    ].append(row)


top_rows = []

for key in sorted(grouped):

    rows = grouped[key]

    rows.sort(
        key=lambda x:
        -x["comparippson_similarity_pct"]
    )

    top_rows.extend(
        rows[:10]
    )


with (
    OUT /
    "comparippson_lactococcin_top10.tsv"
).open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(top_rows)


# Resumen máximo
summary_fields = [
    "query",
    "query_mode",
    "database",
    "best_reference",
    "best_description",
    "best_similarity_pct",
    "alignment_identity_pct",
    "query_coverage_pct",
    "reference_coverage_pct",
]

summary = []

for key in sorted(grouped):

    best = max(
        grouped[key],
        key=lambda x:
        x["comparippson_similarity_pct"]
    )

    summary.append({
        "query": best["query"],
        "query_mode":
            best["query_mode"],
        "database":
            best["database"],
        "best_reference":
            best["reference_identifier"],
        "best_description":
            best["reference_description"],
        "best_similarity_pct":
            best["comparippson_similarity_pct"],
        "alignment_identity_pct":
            best["alignment_identity_pct"],
        "query_coverage_pct":
            best["query_coverage_pct"],
        "reference_coverage_pct":
            best["reference_coverage_pct"],
    })


with (
    OUT /
    "comparippson_lactococcin_best_hits.tsv"
).open("w", newline="") as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=summary_fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(summary)


print(
    "Hits CompaRiPPson >5%:",
    len(all_rows)
)

print(
    "Grupos query/modo/base:",
    len(grouped)
)

PY


echo
echo "============================================================"
echo "SECUENCIAS CANDIDATAS"
echo "============================================================"

column -t -s $'\t' \
"${OUT}/lactococcin_candidate_sequences.tsv"

echo
echo "============================================================"
echo "MEJORES HITS CompaRiPPson"
echo "============================================================"

if [[ $(wc -l < \
"${OUT}/comparippson_lactococcin_best_hits.tsv") -gt 1 ]]; then

    column -t -s $'\t' \
    "${OUT}/comparippson_lactococcin_best_hits.tsv"

else
    echo "No se encontraron hits CompaRiPPson."
fi

echo
echo "============================================================"
echo "COMPARIPPSON LACTOCOCCIN COMPLETADO"
echo "Resultados: ${OUT}"
echo "Fin: $(date)"
echo "============================================================"
