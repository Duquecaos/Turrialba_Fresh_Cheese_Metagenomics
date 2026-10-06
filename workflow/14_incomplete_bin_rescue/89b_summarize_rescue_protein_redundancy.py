#!/usr/bin/env python3

from pathlib import Path
from collections import defaultdict, Counter
import csv
import os


USER = os.environ["USER"]

ROOT = Path(f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba")
OUT = ROOT / "75_incomplete_functional_sequence_redundancy"

SELECTED = OUT / "89A_selected_rescue_functional_genes.tsv"
GENE2MAG = OUT / "final18_gene_to_MAG.tsv"

VS_FINAL = OUT / "89B_vs_final18.raw.tsv"
SELF = OUT / "89C_rescue_self_95id_90cov.raw.tsv"

OUT_BEST = OUT / "89D_best_final18_protein_match.tsv"
OUT_CLMAP = OUT / "89E_rescue_near_identity_cluster_members.tsv"
OUT_CLSUM = OUT / "89F_rescue_near_identity_cluster_summary.tsv"
OUT_CTX = OUT / "89G_match_class_by_context.tsv"

OUT_CEP = OUT / "89H_priority_CEP_sequence_redundancy.tsv"
OUT_BA = OUT / "89I_context_supported_BA_sequence_redundancy.tsv"
OUT_BA_STRONG = OUT / "89J_strong_BA_complete_lowcontam.tsv"

OUT_SUM = OUT / "protein_redundancy_summary.tsv"

BA88 = (
    ROOT
    / "74_incomplete_functional_evidence_context"
    / "88E_context_supported_biogenic_amine_candidates.tsv"
)

CEP88 = (
    ROOT
    / "74_incomplete_functional_evidence_context"
    / "88D_priority_CEP_like_candidates.tsv"
)


# ------------------------------------------------------------------
# Utilidades
# ------------------------------------------------------------------

def read_dict(path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def safe_float(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default


def classify_hit(pident, qcov, scov):
    """
    Clases operacionales de similitud contra final18.
    No equivalen a identidad taxonómica ni a identidad funcional probada.
    """

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


rank = {
    "A_exact_full_length": 6,
    "B_near_identical_full_length": 5,
    "C_near_identical_query_fragment": 4,
    "D_close_homolog": 3,
    "E_homologous_or_domain_level": 2,
    "F_weak_similarity": 1,
    "G_no_detected_similarity": 0,
}


# ------------------------------------------------------------------
# Metadatos seleccionados
# ------------------------------------------------------------------

selected_rows = read_dict(SELECTED)
selected = {
    r["gene_id"]: r
    for r in selected_rows
}

genes = set(selected)


# ------------------------------------------------------------------
# Mapa final18
# ------------------------------------------------------------------

gene2mag = {}

with GENE2MAG.open() as fh:
    r = csv.DictReader(fh, delimiter="\t")

    for row in r:
        gene2mag[row["gene_id"]] = row["MAG"]


# ------------------------------------------------------------------
# Hits contra final18
# ------------------------------------------------------------------

hits_by_query = defaultdict(list)

with VS_FINAL.open() as fh:

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
            evalue,
            bitscore,
        ) = row

        pident = float(pident)
        alnlen = float(alnlen)
        qlen = float(qlen)
        slen = float(slen)
        bitscore = float(bitscore)

        qcov = 100.0 * alnlen / qlen if qlen else 0.0
        scov = 100.0 * alnlen / slen if slen else 0.0

        cls = classify_hit(
            pident,
            qcov,
            scov,
        )

        hits_by_query[qseqid].append({
            "subject_gene_id": sseqid,
            "subject_MAG": gene2mag.get(sseqid, ""),
            "pident": pident,
            "alignment_length": int(alnlen),
            "qlen": int(qlen),
            "slen": int(slen),
            "qcov_pct": qcov,
            "scov_pct": scov,
            "evalue": evalue,
            "bitscore": bitscore,
            "match_class": cls,
        })


best_matches = {}

for gid in sorted(genes):

    exact_hash_n = int(
        selected[gid].get(
            "exact_final18_protein_match_n",
            "0"
        ) or 0
    )

    hs = hits_by_query.get(gid, [])

    if exact_hash_n > 0:

        exact_ids = [
            x for x in
            selected[gid]
            .get("exact_final18_gene_ids", "")
            .split(";")
            if x
        ]

        subject = exact_ids[0] if exact_ids else ""

        best_matches[gid] = {
            "subject_gene_id": subject,
            "subject_MAG": gene2mag.get(subject, ""),
            "pident": 100.0,
            "alignment_length": "",
            "qlen": "",
            "slen": "",
            "qcov_pct": 100.0,
            "scov_pct": 100.0,
            "evalue": "0",
            "bitscore": "",
            "match_class": "A_exact_full_length",
            "exact_sequence_hash_match": 1,
        }

        continue

    if not hs:

        best_matches[gid] = {
            "subject_gene_id": "",
            "subject_MAG": "",
            "pident": "",
            "alignment_length": "",
            "qlen": "",
            "slen": "",
            "qcov_pct": "",
            "scov_pct": "",
            "evalue": "",
            "bitscore": "",
            "match_class": "G_no_detected_similarity",
            "exact_sequence_hash_match": 0,
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

    best = hs[0].copy()
    best["exact_sequence_hash_match"] = 0
    best_matches[gid] = best


# ------------------------------------------------------------------
# Escribir mejores matches
# ------------------------------------------------------------------

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
    "exact_sequence_hash_match",
]

with OUT_BEST.open("w", newline="") as out:

    w = csv.DictWriter(
        out,
        delimiter="\t",
        fieldnames=best_fields,
        lineterminator="\n",
    )

    w.writeheader()

    for gid in sorted(genes):

        rec = selected[gid].copy()
        rec.update(best_matches[gid])

        w.writerow({
            k: rec.get(k, "")
            for k in best_fields
        })


# ------------------------------------------------------------------
# Unión de proteínas rescatadas casi idénticas
# 95% identidad + 90% cobertura de ambas proteínas
# ------------------------------------------------------------------

parent = {
    g: g
    for g in genes
}


def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union(a, b):
    ra = find(a)
    rb = find(b)

    if ra != rb:
        if ra < rb:
            parent[rb] = ra
        else:
            parent[ra] = rb


with SELF.open() as fh:

    r = csv.reader(fh, delimiter="\t")

    for row in r:

        if not row:
            continue

        (
            q,
            s,
            pident,
            alnlen,
            qlen,
            slen,
            evalue,
            bitscore,
        ) = row

        if q == s:
            continue

        if q not in genes or s not in genes:
            continue

        pident = float(pident)
        alnlen = float(alnlen)
        qlen = float(qlen)
        slen = float(slen)

        qcov = 100.0 * alnlen / qlen if qlen else 0
        scov = 100.0 * alnlen / slen if slen else 0

        if (
            pident >= 95
            and qcov >= 90
            and scov >= 90
        ):
            union(q, s)


clusters = defaultdict(list)

for gid in sorted(genes):
    clusters[find(gid)].append(gid)


cluster_items = sorted(
    clusters.values(),
    key=lambda xs: xs[0],
)

gene_cluster = {}

for i, members in enumerate(
    cluster_items,
    start=1,
):
    cid = f"RCLUST{i:05d}"

    for g in members:
        gene_cluster[g] = cid


# ------------------------------------------------------------------
# Tabla miembros
# ------------------------------------------------------------------

with OUT_CLMAP.open("w", newline="") as out:

    fields = [
        "cluster_id",
        "gene_id",
        "cluster_n_members",
        "contig_id",
        "producers",
        "coassemblies",
        "partial",
        "gene_evidence_scope",
        "evidence_tables",
        "functional_labels",
        "final18_match_class",
        "final18_subject_gene",
        "final18_subject_MAG",
    ]

    w = csv.DictWriter(
        out,
        delimiter="\t",
        fieldnames=fields,
        lineterminator="\n",
    )

    w.writeheader()

    cluster_size = {
        gene_cluster[g]: len(clusters[find(g)])
        for g in genes
    }

    for gid in sorted(genes):

        meta = selected[gid]
        best = best_matches[gid]

        w.writerow({
            "cluster_id": gene_cluster[gid],
            "gene_id": gid,
            "cluster_n_members":
                cluster_size[gene_cluster[gid]],
            "contig_id":
                meta.get("contig_id", ""),
            "producers":
                meta.get("producers", ""),
            "coassemblies":
                meta.get("coassemblies", ""),
            "partial":
                meta.get("partial", ""),
            "gene_evidence_scope":
                meta.get("gene_evidence_scope", ""),
            "evidence_tables":
                meta.get("evidence_tables", ""),
            "functional_labels":
                meta.get("functional_labels", ""),
            "final18_match_class":
                best["match_class"],
            "final18_subject_gene":
                best["subject_gene_id"],
            "final18_subject_MAG":
                best["subject_MAG"],
        })


# ------------------------------------------------------------------
# Resumen de clusters
# ------------------------------------------------------------------

cluster_summary_rows = []

for members in cluster_items:

    cid = gene_cluster[members[0]]

    contigs = sorted({
        selected[g].get("contig_id", "")
        for g in members
    })

    producers = sorted({
        x
        for g in members
        for x in selected[g]
        .get("producers", "")
        .split(";")
        if x
    })

    coassemblies = sorted({
        x
        for g in members
        for x in selected[g]
        .get("coassemblies", "")
        .split(";")
        if x
    })

    tables = sorted({
        x
        for g in members
        for x in selected[g]
        .get("evidence_tables", "")
        .split(";")
        if x
    })

    labels = sorted({
        x
        for g in members
        for x in selected[g]
        .get("functional_labels", "")
        .split(";")
        if x
    })

    classes = [
        best_matches[g]["match_class"]
        for g in members
    ]

    best_class = max(
        classes,
        key=lambda x: rank[x],
    )

    cluster_summary_rows.append({
        "cluster_id": cid,
        "n_members": len(members),
        "gene_ids": ";".join(sorted(members)),
        "n_contigs": len(contigs),
        "contigs": ";".join(contigs),
        "producers": ";".join(producers),
        "coassemblies": ";".join(coassemblies),
        "evidence_tables": ";".join(tables),
        "functional_labels": ";".join(labels),
        "best_final18_match_class": best_class,
        "n_members_exact_final18":
            sum(
                c == "A_exact_full_length"
                for c in classes
            ),
        "n_members_near_identical_final18":
            sum(
                c in {
                    "A_exact_full_length",
                    "B_near_identical_full_length",
                }
                for c in classes
            ),
    })


with OUT_CLSUM.open("w", newline="") as out:

    fields = list(
        cluster_summary_rows[0].keys()
    ) if cluster_summary_rows else [
        "cluster_id"
    ]

    w = csv.DictWriter(
        out,
        delimiter="\t",
        fieldnames=fields,
        lineterminator="\n",
    )

    w.writeheader()

    for row in cluster_summary_rows:
        w.writerow(row)


# ------------------------------------------------------------------
# Match class por contexto
# ------------------------------------------------------------------

ctx_counts = Counter()

for gid in genes:

    meta = selected[gid]
    best = best_matches[gid]

    ctx_counts[(
        meta.get("gene_evidence_scope", ""),
        meta.get("strongest_incomplete_context", ""),
        best["match_class"],
    )] += 1


with OUT_CTX.open("w", newline="") as out:

    w = csv.writer(
        out,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writerow([
        "gene_evidence_scope",
        "strongest_incomplete_context",
        "final18_match_class",
        "n_genes",
    ])

    for key, n in sorted(ctx_counts.items()):
        w.writerow([*key, n])


# ------------------------------------------------------------------
# CEP prioritarios
# ------------------------------------------------------------------

cep_rows = read_dict(CEP88)

cep_extra = [
    "rescue_cluster_id",
    "rescue_cluster_n_members",
    "final18_match_class",
    "final18_subject_gene",
    "final18_subject_MAG",
    "final18_pident",
    "final18_qcov_pct",
    "final18_scov_pct",
]

if cep_rows:

    fields = list(cep_rows[0].keys()) + cep_extra

    with OUT_CEP.open("w", newline="") as out:

        w = csv.DictWriter(
            out,
            delimiter="\t",
            fieldnames=fields,
            lineterminator="\n",
        )

        w.writeheader()

        for row in cep_rows:

            gid = row["gene_id"]
            best = best_matches[gid]

            cid = gene_cluster[gid]

            row.update({
                "rescue_cluster_id": cid,
                "rescue_cluster_n_members":
                    len(clusters[find(gid)]),
                "final18_match_class":
                    best["match_class"],
                "final18_subject_gene":
                    best["subject_gene_id"],
                "final18_subject_MAG":
                    best["subject_MAG"],
                "final18_pident":
                    best["pident"],
                "final18_qcov_pct":
                    best["qcov_pct"],
                "final18_scov_pct":
                    best["scov_pct"],
            })

            w.writerow(row)


# ------------------------------------------------------------------
# Aminas context-supported
# ------------------------------------------------------------------

ba_rows = read_dict(BA88)

if ba_rows:

    fields = list(ba_rows[0].keys()) + cep_extra

    with OUT_BA.open("w", newline="") as out:

        w = csv.DictWriter(
            out,
            delimiter="\t",
            fieldnames=fields,
            lineterminator="\n",
        )

        w.writeheader()

        for row in ba_rows:

            gid = row["gene_id"]
            best = best_matches[gid]

            cid = gene_cluster[gid]

            row.update({
                "rescue_cluster_id": cid,
                "rescue_cluster_n_members":
                    len(clusters[find(gid)]),
                "final18_match_class":
                    best["match_class"],
                "final18_subject_gene":
                    best["subject_gene_id"],
                "final18_subject_MAG":
                    best["subject_MAG"],
                "final18_pident":
                    best["pident"],
                "final18_qcov_pct":
                    best["qcov_pct"],
                "final18_scov_pct":
                    best["scov_pct"],
            })

            w.writerow(row)


# ------------------------------------------------------------------
# Subconjunto BA más defendible:
# CDS completa + B50-89 + contaminación <=10
# ------------------------------------------------------------------

strong_ba = []

for row in ba_rows:

    if row.get("prediction_integrity") != "complete_predicted_CDS":
        continue

    cls = row.get("context_support_class", "")

    if cls not in {
        "B50_89_contam_le5",
        "B50_89_contam_gt5_le10",
    }:
        continue

    gid = row["gene_id"]
    best = best_matches[gid]
    cid = gene_cluster[gid]

    rec = row.copy()

    rec.update({
        "rescue_cluster_id": cid,
        "rescue_cluster_n_members":
            len(clusters[find(gid)]),
        "final18_match_class":
            best["match_class"],
        "final18_subject_gene":
            best["subject_gene_id"],
        "final18_subject_MAG":
            best["subject_MAG"],
        "final18_pident":
            best["pident"],
        "final18_qcov_pct":
            best["qcov_pct"],
        "final18_scov_pct":
            best["scov_pct"],
    })

    strong_ba.append(rec)


if strong_ba:

    fields = list(strong_ba[0].keys())

    with OUT_BA_STRONG.open("w", newline="") as out:

        w = csv.DictWriter(
            out,
            delimiter="\t",
            fieldnames=fields,
            lineterminator="\n",
        )

        w.writeheader()

        for row in strong_ba:
            w.writerow(row)


# ------------------------------------------------------------------
# Resumen
# ------------------------------------------------------------------

class_counts = Counter(
    best_matches[g]["match_class"]
    for g in genes
)

multi_clusters = sum(
    len(x) > 1
    for x in cluster_items
)

singleton_clusters = sum(
    len(x) == 1
    for x in cluster_items
)

clusters_near_final = 0

for row in cluster_summary_rows:
    if row["n_members_near_identical_final18"] > 0:
        clusters_near_final += 1


with OUT_SUM.open("w", newline="") as out:

    w = csv.writer(
        out,
        delimiter="\t",
        lineterminator="\n",
    )

    w.writerow(["metric", "value"])

    w.writerow([
        "selected_rescue_functional_genes",
        len(genes),
    ])

    for cls in sorted(
        rank,
        key=lambda x: rank[x],
        reverse=True,
    ):
        w.writerow([
            f"final18_match__{cls}",
            class_counts.get(cls, 0),
        ])

    w.writerow([
        "rescue_near_identity_clusters_95id_90cov",
        len(cluster_items),
    ])

    w.writerow([
        "rescue_singleton_clusters",
        singleton_clusters,
    ])

    w.writerow([
        "rescue_multi_member_clusters",
        multi_clusters,
    ])

    w.writerow([
        "clusters_with_exact_or_near_identical_final18_member",
        clusters_near_final,
    ])

    w.writerow([
        "priority_CEP_candidates",
        len(cep_rows),
    ])

    w.writerow([
        "context_supported_BA_candidates",
        len(ba_rows),
    ])

    w.writerow([
        "strong_BA_complete_B50_89_contam_le10",
        len(strong_ba),
    ])


print("=" * 72)
print("PASO 89b - RESUMEN DE REDUNDANCIA PROTEICA")
print("=" * 72)

print(f"Genes funcionales rescatados : {len(genes)}")

for cls in sorted(
    rank,
    key=lambda x: rank[x],
    reverse=True,
):
    print(
        f"{cls:40s} "
        f"{class_counts.get(cls, 0)}"
    )

print()
print(
    "Clusters rescate 95% id / 90% cobertura : "
    f"{len(cluster_items)}"
)
print(
    f"Clusters singleton                       : {singleton_clusters}"
)
print(
    f"Clusters con >1 miembro                  : {multi_clusters}"
)
print(
    "Clusters con miembro exacto/casi idéntico final18: "
    f"{clusters_near_final}"
)
print()
print(f"CEP prioritarios                         : {len(cep_rows)}")
print(f"BA con soporte contextual                : {len(ba_rows)}")
print(
    "BA fuertes: CDS completa + B50-89 + contam <=10: "
    f"{len(strong_ba)}"
)
print()
print("PASO 89b FINALIZÓ CORRECTAMENTE.")
print("ES SEGURO SALIR.")
