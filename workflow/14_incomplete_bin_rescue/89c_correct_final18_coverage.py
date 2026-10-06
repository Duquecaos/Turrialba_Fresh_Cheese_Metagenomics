#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict, Counter
import csv
import os

USER = os.environ["USER"]
ROOT = Path(f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba")

OUT89 = ROOT / "75_incomplete_functional_sequence_redundancy"
OUT88 = ROOT / "74_incomplete_functional_evidence_context"

SELECTED = OUT89 / "89A_selected_rescue_functional_genes.tsv"
GENE2MAG = OUT89 / "final18_gene_to_MAG.tsv"

RAW = OUT89 / "89B_vs_final18.coverage_corrected.raw.tsv"
CLUSTERS = OUT89 / "89E_rescue_near_identity_cluster_members.tsv"

CEP = OUT88 / "88D_priority_CEP_like_candidates.tsv"
BA = OUT88 / "88E_context_supported_biogenic_amine_candidates.tsv"

BEST_OUT = OUT89 / "89K_best_final18_protein_match_corrected.tsv"
CEP_OUT = OUT89 / "89L_priority_CEP_sequence_redundancy_corrected.tsv"
BA_OUT = OUT89 / "89M_context_supported_BA_sequence_redundancy_corrected.tsv"
BA_STRONG_OUT = OUT89 / "89N_strong_BA_complete_lowcontam_corrected.tsv"
SUMMARY_OUT = OUT89 / "89O_final18_match_summary_corrected.tsv"


def read_tsv(path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


selected_rows = read_tsv(SELECTED)
selected = {r["gene_id"]: r for r in selected_rows}

gene2mag = {}
for r in read_tsv(GENE2MAG):
    gene2mag[r["gene_id"]] = r["MAG"]

cluster_meta = {}
for r in read_tsv(CLUSTERS):
    cluster_meta[r["gene_id"]] = {
        "rescue_cluster_id": r["cluster_id"],
        "rescue_cluster_n_members": r["cluster_n_members"],
    }


rank = {
    "A_exact_full_length": 6,
    "B_near_identical_full_length": 5,
    "C_near_identical_query_fragment": 4,
    "D_close_homolog": 3,
    "E_homologous_or_domain_level": 2,
    "F_weak_similarity": 1,
    "G_no_detected_similarity": 0,
}


def classify(pident, qcov, scov):

    if (
        pident >= 99.999
        and qcov >= 99.5
        and scov >= 99.5
    ):
        return "A_exact_full_length"

    if (
        pident >= 95
        and qcov >= 90
        and scov >= 90
    ):
        return "B_near_identical_full_length"

    if (
        pident >= 95
        and qcov >= 90
        and scov < 90
    ):
        return "C_near_identical_query_fragment"

    if (
        pident >= 70
        and qcov >= 80
        and scov >= 70
    ):
        return "D_close_homolog"

    if (
        pident >= 40
        and qcov >= 70
    ):
        return "E_homologous_or_domain_level"

    return "F_weak_similarity"


hits = defaultdict(list)

with RAW.open() as fh:

    r = csv.reader(fh, delimiter="\t")

    for row in r:

        if not row:
            continue

        (
            qseqid,
            sseqid,
            pident,
            alnlen,
            qlen,
            slen,
            qstart,
            qend,
            sstart,
            send,
            evalue,
            bitscore,
        ) = row

        pident = float(pident)
        qlen = int(qlen)
        slen = int(slen)

        qstart = int(qstart)
        qend = int(qend)
        sstart = int(sstart)
        send = int(send)

        qspan = abs(qend - qstart) + 1
        sspan = abs(send - sstart) + 1

        qcov = 100.0 * qspan / qlen
        scov = 100.0 * sspan / slen

        if qcov > 100.000001 or scov > 100.000001:
            raise RuntimeError(
                f"Cobertura imposible: "
                f"{qseqid} -> {sseqid}; "
                f"qcov={qcov}, scov={scov}"
            )

        cls = classify(
            pident,
            qcov,
            scov,
        )

        hits[qseqid].append({
            "subject_gene_id": sseqid,
            "subject_MAG": gene2mag.get(sseqid, ""),
            "pident": pident,
            "qcov_pct": qcov,
            "scov_pct": scov,
            "evalue": evalue,
            "bitscore": float(bitscore),
            "match_class": cls,
        })


best = {}

for gid, meta in selected.items():

    exact_n = int(
        meta.get(
            "exact_final18_protein_match_n",
            "0"
        ) or 0
    )

    if exact_n > 0:

        ids = [
            x
            for x in meta.get(
                "exact_final18_gene_ids",
                ""
            ).split(";")
            if x
        ]

        sid = ids[0] if ids else ""

        best[gid] = {
            "subject_gene_id": sid,
            "subject_MAG": gene2mag.get(sid, ""),
            "pident": 100.0,
            "qcov_pct": 100.0,
            "scov_pct": 100.0,
            "evalue": "0",
            "bitscore": "",
            "match_class": "A_exact_full_length",
        }

        continue

    hs = hits.get(gid, [])

    if not hs:

        best[gid] = {
            "subject_gene_id": "",
            "subject_MAG": "",
            "pident": "",
            "qcov_pct": "",
            "scov_pct": "",
            "evalue": "",
            "bitscore": "",
            "match_class": "G_no_detected_similarity",
        }

        continue

    hs.sort(
        key=lambda x: (
            rank[x["match_class"]],
            x["bitscore"],
            x["pident"],
            x["qcov_pct"],
            x["scov_pct"],
        ),
        reverse=True,
    )

    best[gid] = hs[0]


best_fields = [
    "gene_id",
    "contig_id",
    "aa_length",
    "partial",
    "gene_evidence_scope",
    "strongest_incomplete_context",
    "coassemblies",
    "producers",
    "evidence_tables",
    "functional_labels",
    "subject_gene_id",
    "subject_MAG",
    "pident",
    "qcov_pct",
    "scov_pct",
    "evalue",
    "bitscore",
    "match_class",
]


with BEST_OUT.open("w", newline="") as out:

    w = csv.DictWriter(
        out,
        delimiter="\t",
        fieldnames=best_fields,
        lineterminator="\n",
    )

    w.writeheader()

    for gid in sorted(selected):

        rec = selected[gid].copy()
        rec.update(best[gid])

        w.writerow({
            k: rec.get(k, "")
            for k in best_fields
        })


extra = [
    "rescue_cluster_id",
    "rescue_cluster_n_members",
    "final18_match_class",
    "final18_subject_gene",
    "final18_subject_MAG",
    "final18_pident",
    "final18_qcov_pct",
    "final18_scov_pct",
]


def add_sequence_info(row):

    gid = row["gene_id"]

    return {
        **row,
        **cluster_meta[gid],
        "final18_match_class":
            best[gid]["match_class"],
        "final18_subject_gene":
            best[gid]["subject_gene_id"],
        "final18_subject_MAG":
            best[gid]["subject_MAG"],
        "final18_pident":
            best[gid]["pident"],
        "final18_qcov_pct":
            best[gid]["qcov_pct"],
        "final18_scov_pct":
            best[gid]["scov_pct"],
    }


def write_joined(rows, outfile):

    joined = [
        add_sequence_info(r)
        for r in rows
    ]

    if not joined:
        raise RuntimeError(
            f"Sin registros para {outfile}"
        )

    fields = list(rows[0].keys()) + extra

    with outfile.open("w", newline="") as out:

        w = csv.DictWriter(
            out,
            delimiter="\t",
            fieldnames=fields,
            lineterminator="\n",
        )

        w.writeheader()

        for row in joined:
            w.writerow(row)

    return joined


cep_rows = read_tsv(CEP)
ba_rows = read_tsv(BA)

cep_joined = write_joined(
    cep_rows,
    CEP_OUT,
)

ba_joined = write_joined(
    ba_rows,
    BA_OUT,
)


strong_ba = [
    r
    for r in ba_joined
    if (
        r.get("prediction_integrity")
        == "complete_predicted_CDS"
        and
        r.get("context_support_class")
        in {
            "B50_89_contam_le5",
            "B50_89_contam_gt5_le10",
        }
    )
]

if not strong_ba:
    raise RuntimeError(
        "El subconjunto BA fuerte quedó vacío."
    )


with BA_STRONG_OUT.open("w", newline="") as out:

    fields = list(strong_ba[0].keys())

    w = csv.DictWriter(
        out,
        delimiter="\t",
        fieldnames=fields,
        lineterminator="\n",
    )

    w.writeheader()

    for row in strong_ba:
        w.writerow(row)


counts = Counter(
    x["match_class"]
    for x in best.values()
)


with SUMMARY_OUT.open("w", newline="") as out:

    w = csv.writer(
        out,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writerow(["metric", "value"])

    w.writerow([
        "selected_rescue_functional_genes",
        len(selected),
    ])

    for cls in sorted(
        rank,
        key=rank.get,
        reverse=True,
    ):
        w.writerow([
            f"final18_match__{cls}",
            counts.get(cls, 0),
        ])

    w.writerow([
        "priority_CEP_candidates",
        len(cep_joined),
    ])

    w.writerow([
        "context_supported_BA_candidates",
        len(ba_joined),
    ])

    w.writerow([
        "strong_BA_complete_B50_89_contam_le10",
        len(strong_ba),
    ])

    w.writerow([
        "strong_BA_unique_sequence_clusters",
        len({
            r["rescue_cluster_id"]
            for r in strong_ba
        }),
    ])

    w.writerow([
        "priority_CEP_unique_sequence_clusters",
        len({
            r["rescue_cluster_id"]
            for r in cep_joined
        }),
    ])


print("=" * 72)
print("PASO 89c - CORRECCION DE COBERTURA CONTRA FINAL18")
print("=" * 72)

for cls in sorted(
    rank,
    key=rank.get,
    reverse=True,
):
    print(
        f"{cls:40s} "
        f"{counts.get(cls, 0)}"
    )

print()
print(
    "CEP prioritarios: "
    f"{len(cep_joined)} genes / "
    f"{len(set(r['rescue_cluster_id'] for r in cep_joined))} clusters"
)

print(
    "BA fuertes: "
    f"{len(strong_ba)} genes / "
    f"{len(set(r['rescue_cluster_id'] for r in strong_ba))} clusters"
)

print()
print("PASO 89c FINALIZÓ CORRECTAMENTE.")
print("ES SEGURO SALIR.")
