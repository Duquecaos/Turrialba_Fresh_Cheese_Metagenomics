#!/bin/bash

#SBATCH --job-name=93d4_IS_context
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=01:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=93d4_IS_context_%j.out
#SBATCH --error=93d4_IS_context_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

ISROOT="${ROOT}/95_mobilome/93D3_ISEScan"

AMR="${ROOT}/80_resistome_community_consolidated/91E_community_AMR_loci.tsv"

INTEGRONS="${ROOT}/95_mobilome/93D2_integron_ARG_plasmid/93D2_integron_structures_integrated.tsv"

DIRECT_INT_AMR="${ROOT}/95_mobilome/93D2_integron_ARG_plasmid/93D2_AMR_integron_direct_links.tsv"

GENOMAD="${ROOT}/85_genomad_coassemblies"

MAGDIR="${ROOT}/22_final_representative_mags/fastas"

MAGTAX="${ROOT}/43_MAG_abundance_taxonomy/representative_MAG_taxonomy.tsv"

OUT="${ROOT}/95_mobilome/93D4_IS_ARG_integron_plasmid"

mkdir -p "$OUT"


echo "======================================================================"
echo "93D4 - IS + ARG + INTEGRON + PLASMID + MAG"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"


for F in \
    "$AMR" \
    "$INTEGRONS" \
    "$DIRECT_INT_AMR" \
    "$MAGTAX"
do
    [[ -s "$F" ]] || {
        echo "ERROR: falta entrada $F" >&2
        exit 2
    }
done


[[ -d "$ISROOT" ]] || {
    echo "ERROR: falta $ISROOT" >&2
    exit 3
}

[[ -d "$GENOMAD" ]] || {
    echo "ERROR: falta $GENOMAD" >&2
    exit 4
}

[[ -d "$MAGDIR" ]] || {
    echo "ERROR: falta $MAGDIR" >&2
    exit 5
}


python3 - \
    "$ISROOT" \
    "$AMR" \
    "$INTEGRONS" \
    "$DIRECT_INT_AMR" \
    "$GENOMAD" \
    "$MAGDIR" \
    "$MAGTAX" \
    "$OUT" <<'PY'

import csv
import hashlib
import sys
from collections import Counter, defaultdict
from pathlib import Path


ISROOT = Path(sys.argv[1])
AMR_FILE = Path(sys.argv[2])
INTEGRON_FILE = Path(sys.argv[3])
DIRECT_INT_AMR_FILE = Path(sys.argv[4])
GENOMAD_ROOT = Path(sys.argv[5])
MAGDIR = Path(sys.argv[6])
MAGTAX_FILE = Path(sys.argv[7])
OUT = Path(sys.argv[8])

GROUPS = ["L1", "L2", "L3", "M1", "M2", "M3"]


# =====================================================================
# Utilidades
# =====================================================================

def read_tsv(path):
    with path.open(errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
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


def to_int(value):
    value = str(value).strip()
    if value in {"", "NA", "None", "nan"}:
        return None
    return int(float(value))


def to_float(value):
    value = str(value).strip()
    if value in {"", "NA", "None", "nan"}:
        return None
    return float(value)


def pick_field(fields, candidates, required=True):
    exact = {x: x for x in fields}
    lower = {x.lower(): x for x in fields}

    for candidate in candidates:
        if candidate in exact:
            return exact[candidate]

    for candidate in candidates:
        if candidate.lower() in lower:
            return lower[candidate.lower()]

    if required:
        raise RuntimeError(
            "No se encontró columna requerida. "
            f"Candidatas={candidates}; columnas={fields}"
        )

    return None


def interval_relation(a_start, a_end, b_start, b_end):
    """
    A = ARG
    B = IS
    """

    a_start, a_end = sorted((a_start, a_end))
    b_start, b_end = sorted((b_start, b_end))

    overlap = (
        max(a_start, b_start)
        <= min(a_end, b_end)
    )

    if overlap:

        if (
            a_start >= b_start
            and a_end <= b_end
        ):
            relation = "ARG_CONTAINED_IN_IS"

        elif (
            b_start >= a_start
            and b_end <= a_end
        ):
            relation = "IS_CONTAINED_IN_ARG"

        else:
            relation = "PARTIAL_COORDINATE_OVERLAP"

        return 0, relation

    if a_end < b_start:
        distance = b_start - a_end

    else:
        distance = a_start - b_end

    if distance <= 1000:
        relation = "WITHIN_1KB"

    elif distance <= 5000:
        relation = "WITHIN_5KB"

    else:
        relation = "SAME_CONTIG_GT5KB"

    return distance, relation


def simple_distance(a_start, a_end, b_start, b_end):

    a_start, a_end = sorted((a_start, a_end))
    b_start, b_end = sorted((b_start, b_end))

    if max(a_start, b_start) <= min(a_end, b_end):
        return 0, True

    if a_end < b_start:
        return b_start - a_end, False

    return a_start - b_end, False


def sha256_file(path):
    h = hashlib.sha256()

    with path.open("rb") as fh:
        for block in iter(
            lambda: fh.read(1024 * 1024),
            b""
        ):
            h.update(block)

    return h.hexdigest()


# =====================================================================
# 1. ISEScan
# =====================================================================

is_rows = []

for group in GROUPS:

    result_pointer = (
        ISROOT
        / group
        / "93D3_result_file.tsv"
    )

    if not result_pointer.is_file():
        raise RuntimeError(
            f"Falta {result_pointer}"
        )

    line = result_pointer.read_text().strip()

    parts = line.split("\t")

    if (
        len(parts) != 2
        or parts[0] != "result_tsv"
    ):
        raise RuntimeError(
            f"Formato inesperado en {result_pointer}: {line}"
        )

    primary = Path(parts[1])

    if not primary.is_file():
        raise RuntimeError(
            f"Falta TSV ISEScan: {primary}"
        )

    rows, fields = read_tsv(primary)

    required = {
        "seqID",
        "family",
        "cluster",
        "isBegin",
        "isEnd",
        "isLen",
        "orfBegin",
        "orfEnd",
        "strand",
        "type",
    }

    missing = required - set(fields)

    if missing:
        raise RuntimeError(
            f"{group}: faltan columnas ISEScan "
            f"{sorted(missing)}"
        )

    for r in rows:

        start = to_int(r["isBegin"])
        end = to_int(r["isEnd"])

        if start is None or end is None:
            raise RuntimeError(
                f"{group}: coordenada IS vacía"
            )

        start, end = sorted((start, end))

        typ = str(r["type"]).strip()

        if typ not in {"c", "p"}:
            raise RuntimeError(
                f"{group}: type inesperado={typ}"
            )

        is_rows.append({
            "coassembly": group,
            "producer": group[0],
            "biological_unit": group[1:],
            "contig": str(r["seqID"]).strip(),
            "family": str(r["family"]).strip(),
            "cluster": str(r["cluster"]).strip(),
            "is_start": start,
            "is_end": end,
            "is_length_bp": to_int(r["isLen"]),
            "orf_start": to_int(r["orfBegin"]),
            "orf_end": to_int(r["orfEnd"]),
            "strand": str(r["strand"]).strip(),
            "type": typ,
            "IS_completeness": (
                "complete"
                if typ == "c"
                else "partial"
            ),
            "score": str(r.get("score", "")).strip(),
            "ir_identity": str(r.get("irId", "")).strip(),
            "ir_length": str(r.get("irLen", "")).strip(),
            "ncopy4is": str(r.get("ncopy4is", "")).strip(),
            "tir": str(r.get("tir", "")).strip(),
        })


is_rows.sort(
    key=lambda x: (
        GROUPS.index(x["coassembly"]),
        x["contig"],
        x["is_start"],
        x["is_end"],
        x["family"],
        x["type"],
    )
)


for i, row in enumerate(is_rows, 1):
    row["IS_id"] = f"ISLOC{i:06d}"


# QC conocido de 93D3

if len(is_rows) != 6391:
    raise RuntimeError(
        f"ISEScan esperado=6391 observado={len(is_rows)}"
    )

n_complete = sum(
    r["type"] == "c"
    for r in is_rows
)

n_partial = sum(
    r["type"] == "p"
    for r in is_rows
)

if n_complete != 1796:
    raise RuntimeError(
        f"IS completos esperado=1796 observado={n_complete}"
    )

if n_partial != 4595:
    raise RuntimeError(
        f"IS parciales esperado=4595 observado={n_partial}"
    )


is_by_contig = defaultdict(list)

for r in is_rows:
    is_by_contig[
        (r["coassembly"], r["contig"])
    ].append(r)


# =====================================================================
# 2. ARG comunitarios
# =====================================================================

amr_rows, amr_fields = read_tsv(AMR_FILE)

required_amr = {
    "community_locus_id",
    "coassembly",
    "contig",
    "start",
    "stop",
    "protein_id",
    "element_symbol",
    "element_name",
    "class",
    "method",
}

missing = required_amr - set(amr_fields)

if missing:
    raise RuntimeError(
        f"AMR: faltan columnas {sorted(missing)}"
    )

if len(amr_rows) != 565:
    raise RuntimeError(
        f"ARG esperado=565 observado={len(amr_rows)}"
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


# =====================================================================
# 3. IntegronFinder integrado
# =====================================================================

integrons, int_fields = read_tsv(INTEGRON_FILE)

required_int = {
    "structure_id",
    "coassembly",
    "contig",
    "integron_type",
    "structure_start",
    "structure_end",
}

missing = required_int - set(int_fields)

if missing:
    raise RuntimeError(
        f"Integrones: faltan columnas {sorted(missing)}"
    )

if len(integrons) != 375:
    raise RuntimeError(
        f"Integrones esperado=375 observado={len(integrons)}"
    )


integrons_by_contig = defaultdict(list)

for r in integrons:

    r["_start"] = min(
        to_int(r["structure_start"]),
        to_int(r["structure_end"])
    )

    r["_end"] = max(
        to_int(r["structure_start"]),
        to_int(r["structure_end"])
    )

    integrons_by_contig[
        (r["coassembly"], r["contig"])
    ].append(r)


# =====================================================================
# 4. ARG directamente asociados a estructuras IntegronFinder
# =====================================================================

direct_int_amr, direct_fields = read_tsv(
    DIRECT_INT_AMR_FILE
)

required_direct = {
    "community_locus_id",
    "coassembly",
    "contig",
    "AMR_start",
    "AMR_end",
    "AMR_symbol",
    "structure_id",
    "integron_type",
    "link_evidence",
}

missing = required_direct - set(direct_fields)

if missing:
    raise RuntimeError(
        f"Direct integron-AMR: faltan {sorted(missing)}"
    )


direct_unique_amr = {
    r["community_locus_id"]
    for r in direct_int_amr
}

if len(direct_unique_amr) != 6:
    raise RuntimeError(
        "Esperaba 6 ARG únicos directamente "
        f"asociados a integrones; observé "
        f"{len(direct_unique_amr)}"
    )


# =====================================================================
# 5. geNomad plasmid candidates
# =====================================================================

plasmid_rows = []
plasmids_by_contig = defaultdict(list)


for group in GROUPS:

    gdir = GENOMAD_ROOT / group

    files = sorted(
        gdir.rglob("*plasmid_summary.tsv")
    )

    if len(files) != 1:
        raise RuntimeError(
            f"{group}: esperaba exactamente 1 "
            f"plasmid_summary.tsv y encontré {len(files)}: "
            f"{files}"
        )

    rows, fields = read_tsv(files[0])

    required = {
        "seq_name",
        "plasmid_score",
        "topology",
        "conjugation_genes",
        "amr_genes",
    }

    missing = required - set(fields)

    if missing:
        raise RuntimeError(
            f"{group}: geNomad faltan "
            f"{sorted(missing)}"
        )

    for r in rows:

        seq_name = str(
            r["seq_name"]
        ).strip()

        prefix = f"{group}__"

        contig = (
            seq_name[len(prefix):]
            if seq_name.startswith(prefix)
            else seq_name
        )

        x = {
            "coassembly": group,
            "contig": contig,
            "plasmid_seq_name": seq_name,
            "plasmid_score": str(
                r["plasmid_score"]
            ).strip(),
            "plasmid_topology": str(
                r["topology"]
            ).strip(),
            "plasmid_conjugation_genes": str(
                r["conjugation_genes"]
            ).strip(),
            "plasmid_amr_genes": str(
                r["amr_genes"]
            ).strip(),
        }

        plasmid_rows.append(x)

        plasmids_by_contig[
            (group, contig)
        ].append(x)


if len(plasmid_rows) != 13599:
    raise RuntimeError(
        f"Plásmidos esperado=13599 observado="
        f"{len(plasmid_rows)}"
    )


# =====================================================================
# 6. MAG final18: contig original -> MAG
# =====================================================================

mag_files = sorted(
    p
    for p in MAGDIR.iterdir()
    if (
        p.is_file()
        and p.suffix.lower()
        in {".fa", ".fna", ".fasta", ".fas"}
    )
)

if len(mag_files) != 18:
    raise RuntimeError(
        f"Esperaba 18 FASTA MAG; encontré {len(mag_files)}"
    )


contig_to_mags = defaultdict(set)
mag_ids = []
mag_contig_count = 0


for fasta in mag_files:

    mag = fasta.stem

    if "__" not in mag:
        raise RuntimeError(
            f"No puedo inferir coassembly de MAG: {mag}"
        )

    group = mag.split("__", 1)[0]

    if group not in GROUPS:
        raise RuntimeError(
            f"Grupo inesperado para MAG {mag}: {group}"
        )

    mag_ids.append(mag)

    with fasta.open() as fh:

        for line in fh:

            if not line.startswith(">"):
                continue

            contig = (
                line[1:]
                .strip()
                .split()[0]
            )

            contig_to_mags[
                (group, contig)
            ].add(mag)

            mag_contig_count += 1


if mag_contig_count != 2463:
    raise RuntimeError(
        f"Contigs final18 esperado=2463 "
        f"observado={mag_contig_count}"
    )


ambiguous_mag_contigs = {
    k: v
    for k, v in contig_to_mags.items()
    if len(v) > 1
}

if ambiguous_mag_contigs:
    raise RuntimeError(
        "Hay contigs asignados a >1 MAG final18: "
        f"{list(ambiguous_mag_contigs.items())[:5]}"
    )


# =====================================================================
# 7. Taxonomía MAG
# =====================================================================

tax_rows, tax_fields = read_tsv(MAGTAX_FILE)

mag_col = pick_field(
    tax_fields,
    [
        "MAG",
        "mag",
        "mag_id",
        "MAG_ID",
        "representative_MAG",
        "representative_mag",
        "bin",
        "bin_id",
    ]
)

tax_col = pick_field(
    tax_fields,
    [
        "taxon_label",
        "taxonomy",
        "gtdb_taxonomy",
        "classification",
        "taxon",
    ]
)


tax_by_mag = {}

for r in tax_rows:

    mag = str(
        r[mag_col]
    ).strip()

    tax = str(
        r[tax_col]
    ).strip()

    tax_by_mag[mag] = tax


matched_taxa = sum(
    mag in tax_by_mag
    for mag in mag_ids
)

if matched_taxa != 18:
    raise RuntimeError(
        f"Taxonomía final18: esperaba 18 MAG "
        f"con taxonomía, observé {matched_taxa}"
    )


# =====================================================================
# 8. Integrar IS con plasmid, integron, MAG
# =====================================================================

integrated_is = []

for isr in is_rows:

    key = (
        isr["coassembly"],
        isr["contig"]
    )

    plasmids = plasmids_by_contig.get(
        key,
        []
    )

    ints = integrons_by_contig.get(
        key,
        []
    )

    mags = sorted(
        contig_to_mags.get(
            key,
            set()
        )
    )

    taxa = sorted({
        tax_by_mag.get(m, "")
        for m in mags
        if tax_by_mag.get(m, "")
    })


    int_overlaps = []
    nearest_int_distance = None
    nearest_int_ids = []

    for intr in ints:

        dist, overlap = simple_distance(
            isr["is_start"],
            isr["is_end"],
            intr["_start"],
            intr["_end"],
        )

        if overlap:
            int_overlaps.append(intr)

        if (
            nearest_int_distance is None
            or dist < nearest_int_distance
        ):
            nearest_int_distance = dist
            nearest_int_ids = [
                intr["structure_id"]
            ]

        elif dist == nearest_int_distance:
            nearest_int_ids.append(
                intr["structure_id"]
            )


    x = dict(isr)

    x.update({
        "plasmid_candidate": (
            "YES" if plasmids else "NO"
        ),
        "plasmid_seq_names": ";".join(
            sorted({
                p["plasmid_seq_name"]
                for p in plasmids
            })
        ),
        "plasmid_max_score": (
            max(
                (
                    to_float(
                        p["plasmid_score"]
                    )
                    for p in plasmids
                    if to_float(
                        p["plasmid_score"]
                    ) is not None
                ),
                default=None
            )
        ),
        "plasmid_conjugation_genes": ";".join(
            sorted({
                p["plasmid_conjugation_genes"]
                for p in plasmids
                if p["plasmid_conjugation_genes"]
                not in {"", "NA", "None"}
            })
        ),
        "n_integron_structures_same_contig": len(
            ints
        ),
        "n_integron_structures_overlapping_IS": len(
            int_overlaps
        ),
        "integron_overlap_ids": ";".join(
            sorted({
                x["structure_id"]
                for x in int_overlaps
            })
        ),
        "integron_overlap_types": ";".join(
            sorted({
                x["integron_type"]
                for x in int_overlaps
            })
        ),
        "nearest_integron_distance_bp": (
            ""
            if nearest_int_distance is None
            else nearest_int_distance
        ),
        "nearest_integron_ids": ";".join(
            sorted(set(nearest_int_ids))
        ),
        "n_final18_MAGs": len(mags),
        "final18_MAGs": ";".join(mags),
        "final18_taxa": ";".join(taxa),
    })

    integrated_is.append(x)


# =====================================================================
# 9. ARG ↔ IS
# =====================================================================

amr_nearest = []
amr_is_pairs_5kb = []


for amr in amr_rows:

    key = (
        amr["coassembly"],
        amr["contig"]
    )

    is_candidates = is_by_contig.get(
        key,
        []
    )

    nearest_distance = None
    nearest_rows = []

    for isr in is_candidates:

        distance, relation = interval_relation(
            amr["_start"],
            amr["_end"],
            isr["is_start"],
            isr["is_end"],
        )

        if (
            nearest_distance is None
            or distance < nearest_distance
        ):
            nearest_distance = distance
            nearest_rows = [
                (isr, relation)
            ]

        elif distance == nearest_distance:
            nearest_rows.append(
                (isr, relation)
            )


        if distance <= 5000:

            plasmid = bool(
                plasmids_by_contig.get(
                    key,
                    []
                )
            )

            mags = sorted(
                contig_to_mags.get(
                    key,
                    set()
                )
            )

            taxa = sorted({
                tax_by_mag.get(m, "")
                for m in mags
                if tax_by_mag.get(m, "")
            })


            amr_is_pairs_5kb.append({
                "community_locus_id":
                    amr["community_locus_id"],
                "coassembly":
                    amr["coassembly"],
                "producer":
                    amr["producer"],
                "biological_unit":
                    amr["biological_unit"],
                "contig":
                    amr["contig"],
                "AMR_start":
                    amr["_start"],
                "AMR_end":
                    amr["_end"],
                "protein_id":
                    amr["protein_id"],
                "AMR_symbol":
                    amr["element_symbol"],
                "AMR_name":
                    amr["element_name"],
                "AMR_class":
                    amr["class"],
                "AMR_method":
                    amr["method"],
                "IS_id":
                    isr["IS_id"],
                "IS_family":
                    isr["family"],
                "IS_type":
                    isr["type"],
                "IS_completeness":
                    isr["IS_completeness"],
                "IS_start":
                    isr["is_start"],
                "IS_end":
                    isr["is_end"],
                "distance_bp":
                    distance,
                "AMR_IS_relation":
                    relation,
                "plasmid_candidate":
                    "YES" if plasmid else "NO",
                "n_final18_MAGs":
                    len(mags),
                "final18_MAGs":
                    ";".join(mags),
                "final18_taxa":
                    ";".join(taxa),
            })


    plasmid = bool(
        plasmids_by_contig.get(
            key,
            []
        )
    )

    mags = sorted(
        contig_to_mags.get(
            key,
            set()
        )
    )

    taxa = sorted({
        tax_by_mag.get(m, "")
        for m in mags
        if tax_by_mag.get(m, "")
    })


    if nearest_rows:

        nearest_ids = ";".join(
            sorted({
                x[0]["IS_id"]
                for x in nearest_rows
            })
        )

        nearest_families = ";".join(
            sorted({
                x[0]["family"]
                for x in nearest_rows
            })
        )

        nearest_types = ";".join(
            sorted({
                x[0]["type"]
                for x in nearest_rows
            })
        )

        nearest_relations = ";".join(
            sorted({
                x[1]
                for x in nearest_rows
            })
        )

    else:

        nearest_ids = ""
        nearest_families = ""
        nearest_types = ""
        nearest_relations = ""


    amr_nearest.append({
        "community_locus_id":
            amr["community_locus_id"],
        "coassembly":
            amr["coassembly"],
        "producer":
            amr["producer"],
        "biological_unit":
            amr["biological_unit"],
        "contig":
            amr["contig"],
        "AMR_start":
            amr["_start"],
        "AMR_end":
            amr["_end"],
        "protein_id":
            amr["protein_id"],
        "AMR_symbol":
            amr["element_symbol"],
        "AMR_name":
            amr["element_name"],
        "AMR_class":
            amr["class"],
        "AMR_method":
            amr["method"],
        "n_IS_same_contig":
            len(is_candidates),
        "nearest_IS_distance_bp": (
            ""
            if nearest_distance is None
            else nearest_distance
        ),
        "nearest_IS_ids":
            nearest_ids,
        "nearest_IS_families":
            nearest_families,
        "nearest_IS_types":
            nearest_types,
        "nearest_IS_relations":
            nearest_relations,
        "plasmid_candidate":
            "YES" if plasmid else "NO",
        "n_final18_MAGs":
            len(mags),
        "final18_MAGs":
            ";".join(mags),
        "final18_taxa":
            ";".join(taxa),
    })


# =====================================================================
# 10. IS ↔ IntegronFinder
# =====================================================================

is_integron_pairs_5kb = []


for isr in is_rows:

    key = (
        isr["coassembly"],
        isr["contig"]
    )

    for intr in integrons_by_contig.get(
        key,
        []
    ):

        distance, overlap = simple_distance(
            isr["is_start"],
            isr["is_end"],
            intr["_start"],
            intr["_end"],
        )

        if distance > 5000:
            continue

        if overlap:
            relation = "IS_OVERLAPS_INTEGRON_STRUCTURE"

        elif distance <= 1000:
            relation = "IS_WITHIN_1KB_OF_INTEGRON"

        else:
            relation = "IS_WITHIN_5KB_OF_INTEGRON"


        is_integron_pairs_5kb.append({
            "coassembly":
                isr["coassembly"],
            "producer":
                isr["producer"],
            "biological_unit":
                isr["biological_unit"],
            "contig":
                isr["contig"],
            "IS_id":
                isr["IS_id"],
            "IS_family":
                isr["family"],
            "IS_type":
                isr["type"],
            "IS_start":
                isr["is_start"],
            "IS_end":
                isr["is_end"],
            "structure_id":
                intr["structure_id"],
            "integron_type":
                intr["integron_type"],
            "structure_start":
                intr["_start"],
            "structure_end":
                intr["_end"],
            "distance_bp":
                distance,
            "IS_integron_relation":
                relation,
            "plasmid_candidate":
                "YES"
                if plasmids_by_contig.get(
                    key,
                    []
                )
                else "NO",
        })


# =====================================================================
# 11. ARG directamente asociado a IntegronFinder ↔ IS
# =====================================================================

direct_int_amr_is = []


for r in direct_int_amr:

    key = (
        r["coassembly"],
        r["contig"]
    )

    astart = min(
        to_int(r["AMR_start"]),
        to_int(r["AMR_end"])
    )

    aend = max(
        to_int(r["AMR_start"]),
        to_int(r["AMR_end"])
    )

    candidates = is_by_contig.get(
        key,
        []
    )

    nearest_distance = None
    nearest = []


    for isr in candidates:

        distance, relation = interval_relation(
            astart,
            aend,
            isr["is_start"],
            isr["is_end"],
        )

        if (
            nearest_distance is None
            or distance < nearest_distance
        ):
            nearest_distance = distance
            nearest = [
                (isr, relation)
            ]

        elif distance == nearest_distance:
            nearest.append(
                (isr, relation)
            )


    if nearest:

        ids = ";".join(
            sorted({
                x[0]["IS_id"]
                for x in nearest
            })
        )

        families = ";".join(
            sorted({
                x[0]["family"]
                for x in nearest
            })
        )

        types = ";".join(
            sorted({
                x[0]["type"]
                for x in nearest
            })
        )

        relations = ";".join(
            sorted({
                x[1]
                for x in nearest
            })
        )

    else:

        ids = ""
        families = ""
        types = ""
        relations = ""


    x = dict(r)

    x.update({
        "nearest_IS_distance_bp": (
            ""
            if nearest_distance is None
            else nearest_distance
        ),
        "nearest_IS_ids":
            ids,
        "nearest_IS_families":
            families,
        "nearest_IS_types":
            types,
        "nearest_IS_relations":
            relations,
        "IS_within_1kb": (
            "YES"
            if (
                nearest_distance is not None
                and nearest_distance <= 1000
            )
            else "NO"
        ),
        "IS_within_5kb": (
            "YES"
            if (
                nearest_distance is not None
                and nearest_distance <= 5000
            )
            else "NO"
        ),
    })

    direct_int_amr_is.append(x)


# =====================================================================
# 12. Resumen familias IS con contexto
# =====================================================================

is_ids_amr_overlap = {
    r["IS_id"]
    for r in amr_is_pairs_5kb
    if r["AMR_IS_relation"]
    in {
        "ARG_CONTAINED_IN_IS",
        "IS_CONTAINED_IN_ARG",
        "PARTIAL_COORDINATE_OVERLAP",
    }
}

is_ids_amr_1kb = {
    r["IS_id"]
    for r in amr_is_pairs_5kb
    if r["distance_bp"] <= 1000
}

is_ids_int_overlap = {
    r["IS_id"]
    for r in is_integron_pairs_5kb
    if r["distance_bp"] == 0
}


family_summary = []

families = sorted({
    r["family"]
    for r in is_rows
})


for family in families:

    fam_rows = [
        r
        for r in integrated_is
        if r["family"] == family
    ]

    family_summary.append({
        "family": family,
        "IS_predictions": len(fam_rows),
        "complete_IS": sum(
            r["type"] == "c"
            for r in fam_rows
        ),
        "partial_IS": sum(
            r["type"] == "p"
            for r in fam_rows
        ),
        "unique_contigs": len({
            (
                r["coassembly"],
                r["contig"]
            )
            for r in fam_rows
        }),
        "IS_on_plasmid_candidate_contigs": sum(
            r["plasmid_candidate"] == "YES"
            for r in fam_rows
        ),
        "IS_on_final18_MAG_contigs": sum(
            r["n_final18_MAGs"] > 0
            for r in fam_rows
        ),
        "IS_overlapping_AMR": sum(
            r["IS_id"] in is_ids_amr_overlap
            for r in fam_rows
        ),
        "IS_within_1kb_AMR": sum(
            r["IS_id"] in is_ids_amr_1kb
            for r in fam_rows
        ),
        "IS_overlapping_integron_structure": sum(
            r["IS_id"] in is_ids_int_overlap
            for r in fam_rows
        ),
    })


family_summary.sort(
    key=lambda x: (
        -x["IS_predictions"],
        x["family"]
    )
)


# =====================================================================
# 13. Resumen global
# =====================================================================

amr_ids_same_contig = {
    r["community_locus_id"]
    for r in amr_nearest
    if r["n_IS_same_contig"] > 0
}

amr_ids_overlap = {
    r["community_locus_id"]
    for r in amr_is_pairs_5kb
    if r["AMR_IS_relation"]
    in {
        "ARG_CONTAINED_IN_IS",
        "IS_CONTAINED_IN_ARG",
        "PARTIAL_COORDINATE_OVERLAP",
    }
}

amr_ids_contained = {
    r["community_locus_id"]
    for r in amr_is_pairs_5kb
    if r["AMR_IS_relation"]
    == "ARG_CONTAINED_IN_IS"
}

amr_ids_1kb = {
    r["community_locus_id"]
    for r in amr_is_pairs_5kb
    if r["distance_bp"] <= 1000
}

amr_ids_5kb = {
    r["community_locus_id"]
    for r in amr_is_pairs_5kb
}

amr_ids_overlap_plasmid = {
    r["community_locus_id"]
    for r in amr_is_pairs_5kb
    if (
        r["community_locus_id"]
        in amr_ids_overlap
        and r["plasmid_candidate"] == "YES"
    )
}

amr_ids_1kb_plasmid = {
    r["community_locus_id"]
    for r in amr_is_pairs_5kb
    if (
        r["distance_bp"] <= 1000
        and r["plasmid_candidate"] == "YES"
    )
}


direct_with_is_overlap = {
    r["community_locus_id"]
    for r in direct_int_amr_is
    if (
        str(
            r["nearest_IS_distance_bp"]
        ).strip() != ""
        and int(
            r["nearest_IS_distance_bp"]
        ) == 0
    )
}

direct_with_is_1kb = {
    r["community_locus_id"]
    for r in direct_int_amr_is
    if r["IS_within_1kb"] == "YES"
}

direct_with_is_5kb = {
    r["community_locus_id"]
    for r in direct_int_amr_is
    if r["IS_within_5kb"] == "YES"
}


global_metrics = [
    ("ISEScan_predictions_total", len(is_rows)),
    ("complete_IS_c", n_complete),
    ("partial_IS_p", n_partial),
    (
        "unique_contigs_with_IS",
        len(is_by_contig)
    ),
    (
        "IS_predictions_on_plasmid_candidate_contigs",
        sum(
            r["plasmid_candidate"] == "YES"
            for r in integrated_is
        )
    ),
    (
        "unique_IS_contigs_plasmid_candidate",
        len({
            (
                r["coassembly"],
                r["contig"]
            )
            for r in integrated_is
            if r["plasmid_candidate"] == "YES"
        })
    ),
    (
        "IS_predictions_on_final18_MAG_contigs",
        sum(
            r["n_final18_MAGs"] > 0
            for r in integrated_is
        )
    ),
    (
        "unique_final18_MAGs_with_IS",
        len({
            mag
            for r in integrated_is
            for mag in r["final18_MAGs"].split(";")
            if mag
        })
    ),
    ("community_AMR_loci_total", len(amr_rows)),
    (
        "AMR_loci_on_contigs_with_IS",
        len(amr_ids_same_contig)
    ),
    (
        "AMR_loci_with_coordinate_overlap_IS",
        len(amr_ids_overlap)
    ),
    (
        "AMR_loci_fully_contained_in_IS",
        len(amr_ids_contained)
    ),
    (
        "AMR_loci_with_IS_within_1kb_including_overlap",
        len(amr_ids_1kb)
    ),
    (
        "AMR_loci_with_IS_within_5kb_including_overlap",
        len(amr_ids_5kb)
    ),
    (
        "AMR_IS_pairs_within_5kb",
        len(amr_is_pairs_5kb)
    ),
    (
        "AMR_overlap_IS_on_plasmid_candidate_contig",
        len(amr_ids_overlap_plasmid)
    ),
    (
        "AMR_with_IS_within_1kb_on_plasmid_candidate_contig",
        len(amr_ids_1kb_plasmid)
    ),
    ("integron_structures_total", len(integrons)),
    (
        "IS_predictions_overlapping_integron_structure",
        len(is_ids_int_overlap)
    ),
    (
        "IS_integron_pairs_within_5kb",
        len(is_integron_pairs_5kb)
    ),
    (
        "direct_integron_AMR_unique_loci",
        len(direct_unique_amr)
    ),
    (
        "direct_integron_AMR_with_IS_overlap",
        len(direct_with_is_overlap)
    ),
    (
        "direct_integron_AMR_with_IS_within_1kb",
        len(direct_with_is_1kb)
    ),
    (
        "direct_integron_AMR_with_IS_within_5kb",
        len(direct_with_is_5kb)
    ),
    ("genomad_plasmid_candidates", len(plasmid_rows)),
    ("final18_MAGs", len(mag_files)),
    ("final18_contigs", mag_contig_count),
]


# =====================================================================
# 14. Escribir salidas
# =====================================================================

IS_FIELDS = [
    "IS_id",
    "coassembly",
    "producer",
    "biological_unit",
    "contig",
    "family",
    "cluster",
    "type",
    "IS_completeness",
    "is_start",
    "is_end",
    "is_length_bp",
    "orf_start",
    "orf_end",
    "strand",
    "score",
    "ir_identity",
    "ir_length",
    "ncopy4is",
    "tir",
    "plasmid_candidate",
    "plasmid_seq_names",
    "plasmid_max_score",
    "plasmid_conjugation_genes",
    "n_integron_structures_same_contig",
    "n_integron_structures_overlapping_IS",
    "integron_overlap_ids",
    "integron_overlap_types",
    "nearest_integron_distance_bp",
    "nearest_integron_ids",
    "n_final18_MAGs",
    "final18_MAGs",
    "final18_taxa",
]

write_tsv(
    OUT / "93D4_IS_integrated.tsv",
    integrated_is,
    IS_FIELDS,
)


AMR_NEAREST_FIELDS = [
    "community_locus_id",
    "coassembly",
    "producer",
    "biological_unit",
    "contig",
    "AMR_start",
    "AMR_end",
    "protein_id",
    "AMR_symbol",
    "AMR_name",
    "AMR_class",
    "AMR_method",
    "n_IS_same_contig",
    "nearest_IS_distance_bp",
    "nearest_IS_ids",
    "nearest_IS_families",
    "nearest_IS_types",
    "nearest_IS_relations",
    "plasmid_candidate",
    "n_final18_MAGs",
    "final18_MAGs",
    "final18_taxa",
]

write_tsv(
    OUT / "93D4_AMR_nearest_IS.tsv",
    amr_nearest,
    AMR_NEAREST_FIELDS,
)


AMR_IS_FIELDS = [
    "community_locus_id",
    "coassembly",
    "producer",
    "biological_unit",
    "contig",
    "AMR_start",
    "AMR_end",
    "protein_id",
    "AMR_symbol",
    "AMR_name",
    "AMR_class",
    "AMR_method",
    "IS_id",
    "IS_family",
    "IS_type",
    "IS_completeness",
    "IS_start",
    "IS_end",
    "distance_bp",
    "AMR_IS_relation",
    "plasmid_candidate",
    "n_final18_MAGs",
    "final18_MAGs",
    "final18_taxa",
]

write_tsv(
    OUT / "93D4_AMR_IS_pairs_within5kb.tsv",
    amr_is_pairs_5kb,
    AMR_IS_FIELDS,
)


INT_IS_FIELDS = [
    "coassembly",
    "producer",
    "biological_unit",
    "contig",
    "IS_id",
    "IS_family",
    "IS_type",
    "IS_start",
    "IS_end",
    "structure_id",
    "integron_type",
    "structure_start",
    "structure_end",
    "distance_bp",
    "IS_integron_relation",
    "plasmid_candidate",
]

write_tsv(
    OUT / "93D4_IS_integron_pairs_within5kb.tsv",
    is_integron_pairs_5kb,
    INT_IS_FIELDS,
)


DIRECT_FIELDS = (
    direct_fields
    + [
        "nearest_IS_distance_bp",
        "nearest_IS_ids",
        "nearest_IS_families",
        "nearest_IS_types",
        "nearest_IS_relations",
        "IS_within_1kb",
        "IS_within_5kb",
    ]
)

write_tsv(
    OUT / "93D4_direct_integron_AMR_with_IS.tsv",
    direct_int_amr_is,
    DIRECT_FIELDS,
)


FAMILY_FIELDS = [
    "family",
    "IS_predictions",
    "complete_IS",
    "partial_IS",
    "unique_contigs",
    "IS_on_plasmid_candidate_contigs",
    "IS_on_final18_MAG_contigs",
    "IS_overlapping_AMR",
    "IS_within_1kb_AMR",
    "IS_overlapping_integron_structure",
]

write_tsv(
    OUT / "93D4_IS_family_context_summary.tsv",
    family_summary,
    FAMILY_FIELDS,
)


with (
    OUT / "93D4_global_summary.tsv"
).open("w") as fh:

    fh.write("metric\tvalue\n")

    for key, value in global_metrics:
        fh.write(
            f"{key}\t{value}\n"
        )


# =====================================================================
# 15. High-context subsets
# =====================================================================

high_context = [
    r
    for r in amr_is_pairs_5kb
    if (
        r["distance_bp"] <= 1000
        and r["plasmid_candidate"] == "YES"
    )
]

write_tsv(
    OUT / "93D4_AMR_IS_plasmid_high_context.tsv",
    high_context,
    AMR_IS_FIELDS,
)


direct_multicontext = [
    r
    for r in direct_int_amr_is
    if (
        r["IS_within_5kb"] == "YES"
        and r.get(
            "plasmid_candidate",
            "NO"
        ) == "YES"
    )
]

write_tsv(
    OUT / "93D4_integron_AMR_IS_plasmid_contexts.tsv",
    direct_multicontext,
    DIRECT_FIELDS,
)


# =====================================================================
# 16. Provenance/QC
# =====================================================================

with (
    OUT / "93D4_metadata.tsv"
).open("w") as fh:

    fh.write("field\tvalue\n")
    fh.write(
        "IS_tool\tISEScan_1.7.3\n"
    )
    fh.write(
        "IS_complete_code\tc\n"
    )
    fh.write(
        "IS_partial_code\tp\n"
    )
    fh.write(
        "AMR_source\tAMRFinderPlus_community_loci\n"
    )
    fh.write(
        "plasmid_source\tgeNomad_candidate_predictions\n"
    )
    fh.write(
        "integron_source\tIntegronFinder_2.0.6\n"
    )
    fh.write(
        "coordinate_overlap_definition\tinclusive_reported_intervals\n"
    )
    fh.write(
        "IS_proximity_threshold_1\t1000_bp\n"
    )
    fh.write(
        "IS_proximity_threshold_2\t5000_bp\n"
    )
    fh.write(
        "IS_prediction_equals_active_transposition\tNO\n"
    )
    fh.write(
        "same_contig_colocalization_equals_HGT\tNO\n"
    )
    fh.write(
        "plasmid_candidate_equals_validated_plasmid\tNO\n"
    )
    fh.write(
        "cross_coassembly_detection_equals_prevalence\tNO\n"
    )
    fh.write(
        "MAG_association_basis\texact_original_contig_membership\n"
    )


input_files = [
    AMR_FILE,
    INTEGRON_FILE,
    DIRECT_INT_AMR_FILE,
    MAGTAX_FILE,
]

with (
    OUT / "93D4_input_sha256.tsv"
).open("w") as fh:

    fh.write("file\tsha256\n")

    for path in input_files:
        fh.write(
            f"{path}\t{sha256_file(path)}\n"
        )


required_outputs = [
    "93D4_IS_integrated.tsv",
    "93D4_AMR_nearest_IS.tsv",
    "93D4_AMR_IS_pairs_within5kb.tsv",
    "93D4_IS_integron_pairs_within5kb.tsv",
    "93D4_direct_integron_AMR_with_IS.tsv",
    "93D4_IS_family_context_summary.tsv",
    "93D4_global_summary.tsv",
    "93D4_AMR_IS_plasmid_high_context.tsv",
    "93D4_integron_AMR_IS_plasmid_contexts.tsv",
    "93D4_metadata.tsv",
    "93D4_input_sha256.tsv",
]


for name in required_outputs:

    path = OUT / name

    if not path.is_file():
        raise RuntimeError(
            f"Salida faltante: {path}"
        )


print("============================================================")
print("93D4 GLOBAL")
print("============================================================")

for key, value in global_metrics:
    print(f"{key}\t{value}")


print()
print(
    "93D4_INTEGRATION=PASS"
)

PY


printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/93D4_COMPLETE.ok"


echo
echo "======================================================================"
echo "93D4 GLOBAL SUMMARY"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93D4_global_summary.tsv"

echo
echo "======================================================================"
echo "93D4 FINALIZO CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "======================================================================"

