#!/usr/bin/env python3

import csv
import os
import re
from collections import defaultdict

USER = os.environ["USER"]

ROOT = f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"

GENES = os.path.join(
    ROOT,
    "53_metabolic_annotation_consolidated",
    "gene_function_inventory_all_predicted.tsv"
)

FAA = os.path.join(
    ROOT,
    "23_eggnog_rep18",
    "results",
    "representative18.emapper.genepred.fasta"
)

MASTER = os.path.join(
    ROOT,
    "22_final_representative_mags",
    "representative_mags_master.tsv"
)

OUTDIR = os.path.join(
    ROOT,
    "57_proteolytic_system_reconstruction"
)

os.makedirs(OUTDIR, exist_ok=True)


# ============================================================
# IO
# ============================================================

def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(rows)


genes = read_tsv(GENES)
master = read_tsv(MASTER)

MAG_ORDER = [
    x["representative_MAG"]
    for x in master
]

MAG_SET = set(MAG_ORDER)

meta = {
    x["representative_MAG"]: x
    for x in master
}


# ============================================================
# PROTEIN LENGTHS + SEQUENCES
# ============================================================

protein_seq = {}
protein_len = {}

current = None
seq = []

with open(FAA) as f:
    for line in f:

        line = line.rstrip()

        if line.startswith(">"):

            if current is not None:
                s = "".join(seq)
                protein_seq[current] = s
                protein_len[current] = len(s)

            current = line[1:].split()[0]
            seq = []

        else:
            seq.append(line.strip())

    if current is not None:
        s = "".join(seq)
        protein_seq[current] = s
        protein_len[current] = len(s)


# ============================================================
# HELPERS
# ============================================================

def norm(x):
    return (x or "").strip().lower()


def split_values(x):
    if not x:
        return set()

    return {
        z.strip()
        for z in re.split(r"[,;]", x)
        if z.strip()
    }


def get_contig_index(gene_id):

    if "___" not in gene_id:
        return "", None

    right = gene_id.split("___", 1)[1]

    m = re.match(r"(.+)_([0-9]+)$", right)

    if not m:
        return right, None

    return m.group(1), int(m.group(2))


def normalize_named_component(name, prefix):

    n = norm(name)

    m = re.fullmatch(
        rf"{prefix}([abcdf])(?:[_-]?[0-9]+)?",
        n
    )

    if not m:
        return None

    return m.group(1).upper()


# ============================================================
# INDEX GENES
# ============================================================

rows_by_mag = defaultdict(list)

for row in genes:

    mag = row.get("MAG", "")

    if mag not in MAG_SET:
        continue

    contig, idx = get_contig_index(
        row.get("gene_id", "")
    )

    row["_contig"] = contig
    row["_index"] = idx
    row["_name"] = norm(
        row.get("Preferred_name", "")
    )
    row["_desc"] = norm(
        row.get("Description", "")
    )
    row["_pfam"] = norm(
        row.get("PFAMs", "")
    )

    rows_by_mag[mag].append(row)


# ============================================================
# SURFACE / EXTRACELLULAR PROTEINASE AUDIT
# ============================================================

surface_candidates = []

for mag in MAG_ORDER:

    for row in rows_by_mag[mag]:

        gid = row["gene_id"]
        n = row["_name"]
        d = row["_desc"]
        p = row["_pfam"]

        length = protein_len.get(gid, 0)
        seq = protein_seq.get(gid, "")

        explicit = False
        reason = []

        if re.fullmatch(
            r"prt[pbh](?:[_-]?[0-9]+)?",
            n
        ):
            explicit = True
            reason.append("explicit_Prt_name")

        if (
            "cell envelope proteinase" in d
            or "cell-envelope proteinase" in d
            or "lactocepin" in d
            or "caseinase" in d
            or "caseinolytic proteinase" in d
            or "extracellular proteinase" in d
            or "extracellular protease" in d
        ):
            explicit = True
            reason.append(
                "explicit_extracellular_description"
            )

        s8 = (
            "peptidase_s8" in p
            or "peptidase_s53" in p
            or "subtilisin" in d
            or "subtilase" in d
        )

        if s8:
            reason.append("S8_S53_or_subtilase")

        large_s8 = s8 and length >= 1200

        if large_s8:
            reason.append("large_>=1200aa")

        # Búsqueda informativa de motivo LPXTG-like
        tail = seq[-120:] if seq else ""

        lpxtg_like = bool(
            re.search(
                r"LP.TG",
                tail,
                re.IGNORECASE
            )
        )

        if lpxtg_like:
            reason.append("C_terminal_LPXTG_like")

        if explicit or s8:

            if explicit:
                tier = "A_explicit_surface_proteinase"

            elif large_s8:
                tier = "B_large_S8_S53_candidate"

            else:
                tier = "C_other_S8_S53_candidate"

            surface_candidates.append({
                "MAG":
                    mag,

                "genus":
                    meta[mag].get(
                        "genus",
                        ""
                    ),

                "species":
                    meta[mag].get(
                        "species",
                        ""
                    ),

                "gene_id":
                    gid,

                "contig":
                    row["_contig"],

                "gene_index":
                    row["_index"],

                "aa_length":
                    length,

                "Preferred_name":
                    row.get(
                        "Preferred_name",
                        ""
                    ),

                "Description":
                    row.get(
                        "Description",
                        ""
                    ),

                "KEGG_ko":
                    row.get(
                        "KEGG_ko",
                        ""
                    ),

                "PFAMs":
                    row.get(
                        "PFAMs",
                        ""
                    ),

                "LPXTG_like_Cterminal":
                    int(lpxtg_like),

                "candidate_tier":
                    tier,

                "candidate_reason":
                    ";".join(reason)
            })


write_tsv(
    os.path.join(
        OUTDIR,
        "surface_proteinase_candidate_audit.tsv"
    ),
    surface_candidates,
    [
        "MAG",
        "genus",
        "species",
        "gene_id",
        "contig",
        "gene_index",
        "aa_length",
        "Preferred_name",
        "Description",
        "KEGG_ko",
        "PFAMs",
        "LPXTG_like_Cterminal",
        "candidate_tier",
        "candidate_reason"
    ]
)


# ============================================================
# NAMED OPP / DPP COMPONENTS
# ============================================================

transport_components = []

for mag in MAG_ORDER:

    for row in rows_by_mag[mag]:

        for system in ["opp", "dpp"]:

            comp = normalize_named_component(
                row.get(
                    "Preferred_name",
                    ""
                ),
                system
            )

            if comp is None:
                continue

            transport_components.append({
                "MAG":
                    mag,

                "system":
                    system.upper(),

                "component":
                    comp,

                "gene_id":
                    row["gene_id"],

                "contig":
                    row["_contig"],

                "gene_index":
                    row["_index"],

                "Preferred_name":
                    row.get(
                        "Preferred_name",
                        ""
                    ),

                "Description":
                    row.get(
                        "Description",
                        ""
                    ),

                "KEGG_ko":
                    row.get(
                        "KEGG_ko",
                        ""
                    )
            })


write_tsv(
    os.path.join(
        OUTDIR,
        "named_Opp_Dpp_components.tsv"
    ),
    transport_components,
    [
        "MAG",
        "system",
        "component",
        "gene_id",
        "contig",
        "gene_index",
        "Preferred_name",
        "Description",
        "KEGG_ko"
    ]
)


# ============================================================
# IDENTIFICAR CLUSTERS CONTIGUOS
# ============================================================

cluster_rows = []

required = {
    "A",
    "B",
    "C",
    "D",
    "F"
}

for mag in MAG_ORDER:

    for system in ["OPP", "DPP"]:

        by_contig = defaultdict(list)

        for x in transport_components:

            if (
                x["MAG"] == mag
                and x["system"] == system
                and x["gene_index"] is not None
            ):
                by_contig[
                    x["contig"]
                ].append(x)

        for contig, rr in by_contig.items():

            rr.sort(
                key=lambda x: x["gene_index"]
            )

            components = {
                x["component"]
                for x in rr
            }

            indices = [
                x["gene_index"]
                for x in rr
            ]

            span = (
                max(indices) - min(indices)
                if indices
                else None
            )

            complete = (
                required.issubset(
                    components
                )
                and span is not None
                and span <= 20
            )

            if complete:
                status = (
                    "complete_named_cluster"
                )

            elif (
                len(
                    required.intersection(
                        components
                    )
                ) >= 3
                and span is not None
                and span <= 20
            ):
                status = (
                    "partial_named_cluster"
                )

            else:
                status = (
                    "dispersed_or_insufficient"
                )

            cluster_rows.append({
                "MAG":
                    mag,

                "system":
                    system,

                "contig":
                    contig,

                "components":
                    ",".join(
                        sorted(
                            components
                        )
                    ),

                "n_required_components":
                    len(
                        required.intersection(
                            components
                        )
                    ),

                "min_gene_index":
                    min(indices)
                    if indices
                    else "",

                "max_gene_index":
                    max(indices)
                    if indices
                    else "",

                "ORF_span":
                    span
                    if span is not None
                    else "",

                "cluster_status":
                    status,

                "genes":
                    ";".join(
                        x["gene_id"]
                        for x in rr
                    )
            })


write_tsv(
    os.path.join(
        OUTDIR,
        "named_peptide_transport_clusters.tsv"
    ),
    cluster_rows,
    [
        "MAG",
        "system",
        "contig",
        "components",
        "n_required_components",
        "min_gene_index",
        "max_gene_index",
        "ORF_span",
        "cluster_status",
        "genes"
    ]
)


# ============================================================
# POT / DTP TRANSPORTERS
# ============================================================

pot_rows = []

for mag in MAG_ORDER:

    for row in rows_by_mag[mag]:

        n = row["_name"]

        kos = split_values(
            row.get(
                "KEGG_ko",
                ""
            )
        )

        name_match = bool(
            re.fullmatch(
                r"dtp[abdt](?:[_-]?[0-9]+)?",
                n
            )
        )

        yjdl_match = (
            n == "yjdl"
            and "ko:K03305" in kos
        )

        if name_match or yjdl_match:

            pot_rows.append({
                "MAG":
                    mag,

                "gene_id":
                    row["gene_id"],

                "contig":
                    row["_contig"],

                "gene_index":
                    row["_index"],

                "Preferred_name":
                    row.get(
                        "Preferred_name",
                        ""
                    ),

                "Description":
                    row.get(
                        "Description",
                        ""
                    ),

                "KEGG_ko":
                    row.get(
                        "KEGG_ko",
                        ""
                    )
            })


write_tsv(
    os.path.join(
        OUTDIR,
        "POT_di_tripeptide_transporters.tsv"
    ),
    pot_rows,
    [
        "MAG",
        "gene_id",
        "contig",
        "gene_index",
        "Preferred_name",
        "Description",
        "KEGG_ko"
    ]
)


# ============================================================
# PEPTIDASAS CLAVE
# ============================================================

PEPTIDASE_TYPES = {
    "pepn": "PepN",
    "pepc": "PepC",
    "pepx": "PepX",
    "pepo": "PepO",
    "pepf": "PepF",
    "pepq": "PepQ",
    "pepv": "PepV",
    "pept": "PepT",
    "pepd": "PepD",
    "pepp": "PepP",
    "pepa": "PepA",
    "pepda": "PepDA",
    "pepb": "PepB",
    "pepe": "PepE",
    "peps": "PepS"
}

peptidase_rows = []

for mag in MAG_ORDER:

    seen_gene = set()

    for row in rows_by_mag[mag]:

        n = row["_name"]

        base = re.sub(
            r"[_-]?[0-9]+$",
            "",
            n
        )

        if base not in PEPTIDASE_TYPES:
            continue

        gid = row["gene_id"]

        if gid in seen_gene:
            continue

        seen_gene.add(gid)

        peptidase_rows.append({
            "MAG":
                mag,

            "peptidase_type":
                PEPTIDASE_TYPES[
                    base
                ],

            "gene_id":
                gid,

            "contig":
                row["_contig"],

            "gene_index":
                row["_index"],

            "Preferred_name":
                row.get(
                    "Preferred_name",
                    ""
                ),

            "Description":
                row.get(
                    "Description",
                    ""
                ),

            "KEGG_ko":
                row.get(
                    "KEGG_ko",
                    ""
                ),

            "EC":
                row.get(
                    "EC",
                    ""
                ),

            "PFAMs":
                row.get(
                    "PFAMs",
                    ""
                )
        })


write_tsv(
    os.path.join(
        OUTDIR,
        "key_intracellular_peptidases.tsv"
    ),
    peptidase_rows,
    [
        "MAG",
        "peptidase_type",
        "gene_id",
        "contig",
        "gene_index",
        "Preferred_name",
        "Description",
        "KEGG_ko",
        "EC",
        "PFAMs"
    ]
)


# ============================================================
# MATRIZ FUNCIONAL REFINADA
# ============================================================

result = []

for mag in MAG_ORDER:

    explicit_surface = [
        x
        for x in surface_candidates
        if (
            x["MAG"] == mag
            and x["candidate_tier"]
            == "A_explicit_surface_proteinase"
        )
    ]

    large_s8 = [
        x
        for x in surface_candidates
        if (
            x["MAG"] == mag
            and x["candidate_tier"]
            == "B_large_S8_S53_candidate"
        )
    ]

    opp_complete = [
        x
        for x in cluster_rows
        if (
            x["MAG"] == mag
            and x["system"] == "OPP"
            and x["cluster_status"]
            == "complete_named_cluster"
        )
    ]

    opp_partial = [
        x
        for x in cluster_rows
        if (
            x["MAG"] == mag
            and x["system"] == "OPP"
            and x["cluster_status"]
            == "partial_named_cluster"
        )
    ]

    dpp_complete = [
        x
        for x in cluster_rows
        if (
            x["MAG"] == mag
            and x["system"] == "DPP"
            and x["cluster_status"]
            == "complete_named_cluster"
        )
    ]

    dpp_partial = [
        x
        for x in cluster_rows
        if (
            x["MAG"] == mag
            and x["system"] == "DPP"
            and x["cluster_status"]
            == "partial_named_cluster"
        )
    ]

    pots = [
        x
        for x in pot_rows
        if x["MAG"] == mag
    ]

    peps = [
        x
        for x in peptidase_rows
        if x["MAG"] == mag
    ]

    pep_types = sorted({
        x["peptidase_type"]
        for x in peps
    })

    # --------------------------------------------------------
    # CASEIN HYDROLYSIS:
    # no inferir ausencia biológica; sólo evidencia recuperada
    # --------------------------------------------------------

    if explicit_surface:
        casein_status = (
            "explicit_surface_proteinase_recovered"
        )

    elif large_s8:
        casein_status = (
            "large_S8_S53_candidate_requires_validation"
        )

    else:
        casein_status = (
            "no_surface_proteinase_recovered_in_MAG"
        )

    # --------------------------------------------------------
    # PEPTIDE UTILIZATION:
    # clasificación operacional transparente
    # --------------------------------------------------------

    has_complete_transport = bool(
        opp_complete
        or dpp_complete
    )

    has_transport = bool(
        has_complete_transport
        or opp_partial
        or dpp_partial
        or pots
    )

    if (
        has_complete_transport
        and len(pep_types) >= 4
    ):
        peptide_status = (
            "strong_complete_transport_plus_peptidases"
        )

    elif (
        has_transport
        and len(pep_types) >= 3
    ):
        peptide_status = (
            "supported_transport_plus_peptidases"
        )

    elif len(pep_types) >= 3:
        peptide_status = (
            "peptidase_repertoire_without_resolved_transport"
        )

    else:
        peptide_status = (
            "limited_recovered_evidence"
        )

    m = meta[mag]

    result.append({
        "MAG":
            mag,

        "genus":
            m.get(
                "genus",
                ""
            ),

        "species":
            m.get(
                "species",
                ""
            ),

        "completeness":
            m.get(
                "completeness",
                ""
            ),

        "surface_proteinase_explicit_n":
            len(explicit_surface),

        "large_S8_S53_candidate_n":
            len(large_s8),

        "casein_hydrolysis_evidence":
            casein_status,

        "Opp_complete_clusters":
            len(opp_complete),

        "Opp_partial_clusters":
            len(opp_partial),

        "Dpp_complete_clusters":
            len(dpp_complete),

        "Dpp_partial_clusters":
            len(dpp_partial),

        "POT_Dtp_genes":
            len(pots),

        "n_key_peptidase_types":
            len(pep_types),

        "key_peptidase_types":
            ";".join(
                pep_types
            ),

        "peptide_utilization_evidence":
            peptide_status
    })


write_tsv(
    os.path.join(
        OUTDIR,
        "proteolytic_capability_matrix_18MAGs.tsv"
    ),
    result,
    [
        "MAG",
        "genus",
        "species",
        "completeness",
        "surface_proteinase_explicit_n",
        "large_S8_S53_candidate_n",
        "casein_hydrolysis_evidence",
        "Opp_complete_clusters",
        "Opp_partial_clusters",
        "Dpp_complete_clusters",
        "Dpp_partial_clusters",
        "POT_Dtp_genes",
        "n_key_peptidase_types",
        "key_peptidase_types",
        "peptide_utilization_evidence"
    ]
)


# ============================================================
# MAGs DE INTERÉS PARA BACTERIOCINAS
# ============================================================

FOCUS = [
    "L2__L2_maxbin2.004_sub",
    "L3__concoct_29",
    "M2__M2_maxbin2.004",
    "L2__L2_maxbin2.011_sub"
]

focus = [
    x
    for x in result
    if x["MAG"] in FOCUS
]


write_tsv(
    os.path.join(
        OUTDIR,
        "proteolytic_capability_bacteriocin_focus_MAGs.tsv"
    ),
    focus,
    [
        "MAG",
        "genus",
        "species",
        "completeness",
        "surface_proteinase_explicit_n",
        "large_S8_S53_candidate_n",
        "casein_hydrolysis_evidence",
        "Opp_complete_clusters",
        "Opp_partial_clusters",
        "Dpp_complete_clusters",
        "Dpp_partial_clusters",
        "POT_Dtp_genes",
        "n_key_peptidase_types",
        "key_peptidase_types",
        "peptide_utilization_evidence"
    ]
)


# ============================================================
# CONSOLA
# ============================================================

print("=" * 60)
print("PASO 78b COMPLETADO")
print("=" * 60)

print(
    f"MAGs analizados: {len(MAG_ORDER)}"
)

print(
    "Proteinasa superficial explícita: "
    f"{sum(x['surface_proteinase_explicit_n'] > 0 for x in result)}/18"
)

print(
    "Candidatos grandes S8/S53: "
    f"{sum(x['large_S8_S53_candidate_n'] > 0 for x in result)}/18"
)

print(
    "Opp completo: "
    f"{sum(x['Opp_complete_clusters'] > 0 for x in result)}/18"
)

print(
    "Dpp completo: "
    f"{sum(x['Dpp_complete_clusters'] > 0 for x in result)}/18"
)

print(
    "POT/Dtp: "
    f"{sum(x['POT_Dtp_genes'] > 0 for x in result)}/18"
)

print()

for x in result:
    print(
        f"{x['MAG']:28s} "
        f"OppC={x['Opp_complete_clusters']} "
        f"DppC={x['Dpp_complete_clusters']} "
        f"POT={x['POT_Dtp_genes']} "
        f"PepTypes={x['n_key_peptidase_types']:2d} "
        f"{x['peptide_utilization_evidence']}"
    )

print()
print(
    f"Salida: {OUTDIR}"
)

print(
    "PASO 78b FINALIZÓ CORRECTAMENTE; "
    "ES SEGURO SALIR."
)
