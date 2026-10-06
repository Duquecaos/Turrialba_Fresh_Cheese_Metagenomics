#!/bin/bash

#SBATCH --job-name=92h0b_directhost
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=92h0b_directhost_%j.out
#SBATCH --error=92h0b_directhost_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

OUT="${ROOT}/94_votu_host_taxonomy/audit"
MAGDIR="${ROOT}/22_final_representative_mags/fastas"

VIRAL="${ROOT}/86_viral_catalog_preQC/92D_viral_catalog_all.tsv"
EXACTMEM="${ROOT}/86_viral_catalog_preQC/92D_exact_sequence_members.tsv"
CLEANMEM="${ROOT}/89_viral_catalog_hostclean/92E2_hostclean_exact_members.tsv"
VOTUMEM="${ROOT}/91_votu_catalog/92F_votu_membership.tsv"
VOTUSUM="${ROOT}/91_votu_catalog/92F_votu_summary.tsv"

TRACE="${ROOT}/22_final_representative_mags/representative18_traceable.fasta"

mkdir -p "$OUT"

echo "======================================================================"
echo "92H0b - AUDITORIA HOST/TAXONOMIA + ENLACE DIRECTO PROVIRUS-MAG"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"

for F in \
    "$VIRAL" \
    "$EXACTMEM" \
    "$CLEANMEM" \
    "$VOTUMEM" \
    "$VOTUSUM"
do
    [[ -s "$F" ]] || {
        echo "ERROR: falta $F" >&2
        exit 2
    }
done

[[ -d "$MAGDIR" ]] || {
    echo "ERROR: falta $MAGDIR" >&2
    exit 3
}


# ======================================================================
# 1. TRACEABLE PREVIEW SIN PIPE grep|head
# ======================================================================

if [[ -s "$TRACE" ]]
then

    awk '
        /^>/ {
            print
            n++
            if (n == 30)
                exit
        }
    ' "$TRACE" \
    > "${OUT}/92H0b_traceable_headers_preview.txt"

    awk '
        /^>/ {n++}
        END {print n+0}
    ' "$TRACE" \
    > "${OUT}/92H0b_traceable_contig_count.txt"

else

    echo "TRACEABLE_FASTA_MISSING" \
    > "${OUT}/92H0b_traceable_headers_preview.txt"

    echo "0" \
    > "${OUT}/92H0b_traceable_contig_count.txt"

fi


# ======================================================================
# 2. DIRECT HOST LINK:
#
# geNomad provirus
# -> exact cluster
# -> hostclean cluster
# -> vOTU
# -> exact original contig match inside final18 MAG
# ======================================================================

python3 - \
    "$VIRAL" \
    "$EXACTMEM" \
    "$CLEANMEM" \
    "$VOTUMEM" \
    "$MAGDIR" \
    "$OUT" <<'PY'

import csv
import sys

from collections import Counter, defaultdict
from pathlib import Path


(
    viral_file,
    exactmem_file,
    cleanmem_file,
    votumem_file,
    magdir,
    outdir,
) = map(Path, sys.argv[1:])


# ============================================================
# Helpers
# ============================================================

def read_tsv(path):

    with path.open() as fh:

        return list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


def clean(x):

    if x is None:
        return ""

    return str(x).strip()


# ============================================================
# 18 MAGs:
# (coassembly, original k141 contig) -> MAG(s)
# ============================================================

fasta_files = sorted(
    list(
        magdir.glob("*.fa")
    )
    + list(
        magdir.glob("*.fna")
    )
    + list(
        magdir.glob("*.fasta")
    )
)

if len(fasta_files) != 18:

    raise RuntimeError(
        f"Esperaba 18 MAG FASTA; "
        f"encontré {len(fasta_files)}"
    )


contig_to_mags = defaultdict(set)

mag_contig_counts = Counter()


for fasta in fasta_files:

    mag = fasta.stem

    if "__" not in mag:

        raise RuntimeError(
            f"No puedo inferir coassembly de {mag}"
        )

    coassembly = mag.split(
        "__",
        1
    )[0]


    with fasta.open() as fh:

        for line in fh:

            if not line.startswith(">"):
                continue

            contig = (
                line[1:]
                .strip()
                .split()[0]
            )

            key = (
                coassembly,
                contig,
            )

            contig_to_mags[
                key
            ].add(
                mag
            )

            mag_contig_counts[
                mag
            ] += 1


# ============================================================
# seq_name -> VEXACT
# ============================================================

seq_to_vexact = {}

for row in read_tsv(
    exactmem_file
):

    seq = row[
        "seq_name"
    ]

    vex = row[
        "exact_cluster_id"
    ]

    if seq in seq_to_vexact:

        raise RuntimeError(
            f"seq_name duplicado en exact members: {seq}"
        )

    seq_to_vexact[
        seq
    ] = vex


# ============================================================
# VEXACT -> VCLEAN
# ============================================================

vexact_to_vclean = {}

for row in read_tsv(
    cleanmem_file
):

    vex = row[
        "preclean_exact_cluster_id"
    ]

    vclean = row[
        "hostclean_exact_cluster_id"
    ]

    if vex in vexact_to_vclean:

        raise RuntimeError(
            f"VEXACT duplicado en hostclean members: {vex}"
        )

    vexact_to_vclean[
        vex
    ] = vclean


# ============================================================
# VCLEAN -> VOTU
# ============================================================

vclean_to_votu = {}

for row in read_tsv(
    votumem_file
):

    vclean = row[
        "vclean_member"
    ]

    votu = row[
        "votu_id"
    ]

    if vclean in vclean_to_votu:

        raise RuntimeError(
            f"VCLEAN duplicado en vOTU membership: {vclean}"
        )

    vclean_to_votu[
        vclean
    ] = votu


# ============================================================
# Viral catalog
# ============================================================

viral_rows = read_tsv(
    viral_file
)

if len(viral_rows) != 6888:

    raise RuntimeError(
        f"viral rows={len(viral_rows)}; esperaba 6888"
    )


context_counts = Counter()
topology_counts = Counter()

provirus_rows = []

direct_rows = []

unmapped_chain = 0
ambiguous_mag_hits = 0


for row in viral_rows:

    seq = row[
        "seq_name"
    ]

    coassembly = row[
        "coassembly"
    ]

    contig = row[
        "parent_contig_original"
    ]

    topology = clean(
        row[
            "topology"
        ]
    )

    context = clean(
        row[
            "viral_context"
        ]
    )


    topology_counts[
        topology
    ] += 1

    context_counts[
        context
    ] += 1


    is_provirus = (
        "PROVIRUS"
        in context.upper()
        or topology.upper()
        == "PROVIRUS"
    )


    if not is_provirus:
        continue


    vex = seq_to_vexact.get(
        seq
    )

    vclean = (
        vexact_to_vclean.get(
            vex
        )
        if vex
        else None
    )

    votu = (
        vclean_to_votu.get(
            vclean
        )
        if vclean
        else None
    )


    if not (
        vex
        and vclean
        and votu
    ):

        unmapped_chain += 1


    mags = sorted(
        contig_to_mags.get(
            (
                coassembly,
                contig,
            ),
            set(),
        )
    )


    if len(mags) == 0:

        status = (
            "NO_FINAL18_MAG_MATCH"
        )

    elif len(mags) == 1:

        status = (
            "FINAL18_DIRECT"
        )

    else:

        status = (
            "AMBIGUOUS_FINAL18"
        )

        ambiguous_mag_hits += 1


    base = {

        "seq_name":
            seq,

        "coassembly":
            coassembly,

        "producer":
            row[
                "producer"
            ],

        "biological_unit":
            row[
                "biological_unit"
            ],

        "parent_contig_original":
            contig,

        "coordinates":
            row[
                "coordinates"
            ],

        "viral_length":
            row[
                "length"
            ],

        "viral_context":
            context,

        "topology":
            topology,

        "taxonomy":
            row[
                "taxonomy"
            ],

        "terminal_taxon":
            row[
                "terminal_taxon"
            ],

        "virus_score":
            row[
                "virus_score"
            ],

        "fdr":
            row[
                "fdr"
            ],

        "n_hallmarks":
            row[
                "n_hallmarks"
            ],

        "exact_cluster_id":
            vex or "",

        "hostclean_cluster_id":
            vclean or "",

        "votu_id":
            votu or "",

        "host_link_status":
            status,

        "final18_MAGs":
            ";".join(
                mags
            ),
    }


    provirus_rows.append(
        base
    )


    for mag in mags:

        d = dict(
            base
        )

        d[
            "final18_MAG"
        ] = mag

        direct_rows.append(
            d
        )


# ============================================================
# Write all provirus audit
# ============================================================

all_file = (
    outdir
    / "92H0b_all_provirus_to_final18_audit.tsv"
)

with all_file.open(
    "w",
    newline=""
) as fh:

    fields = list(
        provirus_rows[
            0
        ].keys()
    )

    w = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    w.writeheader()
    w.writerows(
        provirus_rows
    )


# ============================================================
# Direct exact links
# ============================================================

direct_file = (
    outdir
    / "92H0b_direct_provirus_MAG_links.tsv"
)

if direct_rows:

    with direct_file.open(
        "w",
        newline=""
    ) as fh:

        fields = list(
            direct_rows[
                0
            ].keys()
        )

        w = csv.DictWriter(
            fh,
            fieldnames=fields,
            delimiter="\t"
        )

        w.writeheader()
        w.writerows(
            direct_rows
        )

else:

    direct_file.write_text(
        "seq_name\tvotu_id\tfinal18_MAG\n"
    )


# ============================================================
# Collapse to VOTU-MAG association
# ============================================================

votu_mag = defaultdict(
    lambda: {
        "n_provirus_predictions": 0,
        "seq_names": set(),
        "parent_contigs": set(),
        "coassemblies": set(),
        "terminal_taxa": set(),
    }
)


for row in direct_rows:

    votu = row[
        "votu_id"
    ]

    mag = row[
        "final18_MAG"
    ]

    if not votu:
        continue

    key = (
        votu,
        mag,
    )

    x = votu_mag[
        key
    ]

    x[
        "n_provirus_predictions"
    ] += 1

    x[
        "seq_names"
    ].add(
        row[
            "seq_name"
        ]
    )

    x[
        "parent_contigs"
    ].add(
        row[
            "parent_contig_original"
        ]
    )

    x[
        "coassemblies"
    ].add(
        row[
            "coassembly"
        ]
    )

    x[
        "terminal_taxa"
    ].add(
        row[
            "terminal_taxon"
        ]
    )


collapsed_file = (
    outdir
    / "92H0b_direct_votu_MAG_links.tsv"
)

with collapsed_file.open(
    "w",
    newline=""
) as fh:

    fields = [
        "votu_id",
        "final18_MAG",
        "host_evidence",
        "n_provirus_predictions",
        "n_parent_contigs",
        "parent_contigs",
        "seq_names",
        "coassemblies",
        "viral_terminal_taxa",
    ]

    w = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    w.writeheader()


    for (
        votu,
        mag
    ), x in sorted(
        votu_mag.items()
    ):

        w.writerow({

            "votu_id":
                votu,

            "final18_MAG":
                mag,

            "host_evidence":
                "DIRECT_PROVIRUS_CONTIG_IN_MAG",

            "n_provirus_predictions":
                x[
                    "n_provirus_predictions"
                ],

            "n_parent_contigs":
                len(
                    x[
                        "parent_contigs"
                    ]
                ),

            "parent_contigs":
                ";".join(
                    sorted(
                        x[
                            "parent_contigs"
                        ]
                    )
                ),

            "seq_names":
                ";".join(
                    sorted(
                        x[
                            "seq_names"
                        ]
                    )
                ),

            "coassemblies":
                ";".join(
                    sorted(
                        x[
                            "coassemblies"
                        ]
                    )
                ),

            "viral_terminal_taxa":
                ";".join(
                    sorted(
                        x[
                            "terminal_taxa"
                        ]
                    )
                ),
        })


# ============================================================
# Context/topology counts
# ============================================================

context_file = (
    outdir
    / "92H0b_viral_context_counts.tsv"
)

with context_file.open(
    "w"
) as fh:

    fh.write(
        "viral_context\tcount\n"
    )

    for key, value in sorted(
        context_counts.items()
    ):

        fh.write(
            f"{key}\t{value}\n"
        )


topology_file = (
    outdir
    / "92H0b_topology_counts.tsv"
)

with topology_file.open(
    "w"
) as fh:

    fh.write(
        "topology\tcount\n"
    )

    for key, value in sorted(
        topology_counts.items()
    ):

        fh.write(
            f"{key}\t{value}\n"
        )


# ============================================================
# Summary
# ============================================================

unique_direct_votus = {
    row[
        "votu_id"
    ]
    for row in direct_rows
    if row[
        "votu_id"
    ]
}

unique_direct_mags = {
    row[
        "final18_MAG"
    ]
    for row in direct_rows
}

unique_direct_contigs = {
    (
        row[
            "coassembly"
        ],
        row[
            "parent_contig_original"
        ],
    )
    for row in direct_rows
}


summary = [
    (
        "viral_predictions_total",
        len(
            viral_rows
        ),
    ),

    (
        "provirus_predictions_detected",
        len(
            provirus_rows
        ),
    ),

    (
        "provirus_chain_unmapped_to_vOTU",
        unmapped_chain,
    ),

    (
        "provirus_predictions_matching_final18_MAG_contig",
        sum(
            x[
                "host_link_status"
            ]
            == "FINAL18_DIRECT"
            for x in provirus_rows
        ),
    ),

    (
        "provirus_predictions_no_final18_match",
        sum(
            x[
                "host_link_status"
            ]
            == "NO_FINAL18_MAG_MATCH"
            for x in provirus_rows
        ),
    ),

    (
        "ambiguous_final18_contig_matches",
        ambiguous_mag_hits,
    ),

    (
        "unique_direct_vOTUs",
        len(
            unique_direct_votus
        ),
    ),

    (
        "unique_direct_MAGs",
        len(
            unique_direct_mags
        ),
    ),

    (
        "unique_direct_host_contigs",
        len(
            unique_direct_contigs
        ),
    ),

    (
        "unique_vOTU_MAG_pairs",
        len(
            votu_mag
        ),
    ),

    (
        "final18_MAGs_total",
        len(
            fasta_files
        ),
    ),

    (
        "final18_contigs_total",
        sum(
            mag_contig_counts.values()
        ),
    ),
]


summary_file = (
    outdir
    / "92H0b_direct_host_summary.tsv"
)

with summary_file.open(
    "w"
) as fh:

    fh.write(
        "metric\tvalue\n"
    )

    for key, value in summary:

        fh.write(
            f"{key}\t{value}\n"
        )


print()
print("=" * 78)
print("DIRECT PROVIRUS -> FINAL18 MAG")
print("=" * 78)

for key, value in summary:

    print(
        key,
        value,
        sep="\t"
    )

PY


# ======================================================================
# 3. TAXONOMY CANDIDATES
#
# safe_find evita que permisos parciales aborten el job.
# ======================================================================

safe_find () {
    find "$@" 2>/dev/null || true
}


safe_find "$ROOT" \
    -type f \
    \( \
        -iname '*gtdb*.tsv' \
        -o -iname '*gtdb*.txt' \
        -o -iname '*bac120*summary*' \
        -o -iname '*ar53*summary*' \
        -o -iname '*taxonomy*.tsv' \
        -o -iname '*classification*.tsv' \
    \) \
    -printf '%s\t%p\n' \
| sort -nr \
> "${OUT}/92H0b_taxonomy_candidates.tsv"


python3 - \
    "${OUT}/92H0b_taxonomy_candidates.tsv" \
    "${OUT}/92H0b_taxonomy_candidate_preview.tsv" <<'PY'

import sys
from pathlib import Path

source = Path(sys.argv[1])
dest = Path(sys.argv[2])

rows = []

if source.exists():

    for line in source.read_text().splitlines():

        if not line.strip():
            continue

        try:
            size, path = line.split(
                "\t",
                1
            )
        except ValueError:
            continue

        p = Path(path)

        preview = ""

        try:

            with p.open(
                errors="replace"
            ) as fh:

                preview = (
                    fh.readline()
                    .strip()
                )

        except Exception as e:

            preview = (
                f"UNREADABLE:{type(e).__name__}"
            )


        if len(preview) > 250:
            preview = (
                preview[:247]
                + "..."
            )


        rows.append(
            (
                size,
                path,
                preview,
            )
        )


with dest.open("w") as fh:

    fh.write(
        "bytes\tpath\tfirst_line\n"
    )

    for row in rows:

        fh.write(
            "\t".join(row)
            + "\n"
        )

PY


# ======================================================================
# 4. HOST TOOLS
# ======================================================================

{
    echo -e "tool\tPATH_status\tpath"

    for TOOL in \
        iphop \
        minced \
        blastn \
        makeblastdb \
        crisprcasfinder \
        CRISPRCasFinder.pl \
        crispridentify \
        spacerblast \
        prodigal
    do

        P="$(
            command -v "$TOOL" \
            2>/dev/null \
            || true
        )"

        if [[ -n "$P" ]]
        then

            echo -e \
"${TOOL}\tFOUND\t${P}"

        else

            echo -e \
"${TOOL}\tNOT_FOUND\t"

        fi

    done

} > "${OUT}/92H0b_host_tools_PATH.tsv"


safe_find "${CONTAINER_ROOT:-/opt/ohpc/pub/containers/BIO}" \
    -maxdepth 1 \
    -type f \
    \( \
        -iname '*iphop*' \
        -o -iname '*crispr*' \
        -o -iname '*minced*' \
        -o -iname '*blast*' \
        -o -iname '*prodigal*' \
        -o -iname '*genomad*' \
    \) \
    -printf '%f\t%p\n' \
| sort \
> "${OUT}/92H0b_host_container_candidates.tsv"


# ======================================================================
# 5. HOST DATABASE CANDIDATES
# ======================================================================

if [[ -d "${ROOT}/databases" ]]
then

    safe_find "${ROOT}/databases" \
        -maxdepth 4 \
        -type d \
        \( \
            -iname '*iphop*' \
            -o -iname '*crispr*' \
            -o -iname '*virus*' \
            -o -iname '*phage*' \
            -o -iname '*gtdb*' \
        \) \
        -print \
    | sort \
    > "${OUT}/92H0b_host_database_candidates.txt"

else

    : > "${OUT}/92H0b_host_database_candidates.txt"

fi


# ======================================================================
# 6. geNomad FILE INVENTORY
# ======================================================================

if [[ -d "${ROOT}/85_genomad_coassemblies" ]]
then

    safe_find "${ROOT}/85_genomad_coassemblies" \
        -type f \
        \( \
            -name '*virus_summary.tsv' \
            -o -name '*virus_taxonomy.tsv' \
            -o -name '*virus.fna' \
            -o -name '*provirus*.tsv' \
        \) \
        -printf '%s\t%p\n' \
    | sort -nr \
    > "${OUT}/92H0b_genomad_files.tsv"

else

    : > "${OUT}/92H0b_genomad_files.tsv"

fi


# ======================================================================
# 7. GLOBAL SUMMARY
# ======================================================================

{
    echo -e "metric\tvalue"

    awk -F '\t' \
        'NR>1 {print $1"\t"$2}' \
        "${OUT}/92H0b_direct_host_summary.tsv"

    echo -e \
"taxonomy_candidate_files\t$(
    awk 'END{print NR+0}' \
    "${OUT}/92H0b_taxonomy_candidates.tsv"
)"

    echo -e \
"host_container_candidates\t$(
    awk 'END{print NR+0}' \
    "${OUT}/92H0b_host_container_candidates.tsv"
)"

    echo -e \
"host_database_candidates\t$(
    awk 'END{print NR+0}' \
    "${OUT}/92H0b_host_database_candidates.txt"
)"

} > "${OUT}/92H0b_summary.tsv"


# ======================================================================
# 8. FINAL VALIDATION
# ======================================================================

for F in \
    "${OUT}/92H0b_direct_host_summary.tsv" \
    "${OUT}/92H0b_all_provirus_to_final18_audit.tsv" \
    "${OUT}/92H0b_direct_provirus_MAG_links.tsv" \
    "${OUT}/92H0b_direct_votu_MAG_links.tsv" \
    "${OUT}/92H0b_viral_context_counts.tsv" \
    "${OUT}/92H0b_topology_counts.tsv" \
    "${OUT}/92H0b_taxonomy_candidates.tsv" \
    "${OUT}/92H0b_taxonomy_candidate_preview.tsv" \
    "${OUT}/92H0b_host_tools_PATH.tsv" \
    "${OUT}/92H0b_host_container_candidates.tsv" \
    "${OUT}/92H0b_host_database_candidates.txt" \
    "${OUT}/92H0b_genomad_files.tsv" \
    "${OUT}/92H0b_summary.tsv"
do

    [[ -e "$F" ]] || {
        echo "ERROR: falta $F" >&2
        exit 4
    }

done


printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"$SLURM_JOB_ID" \
> "${OUT}/92H0b_COMPLETE.ok"


echo
echo "======================================================================"
echo "GLOBAL SUMMARY"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92H0b_summary.tsv"


echo
echo "======================================================================"
echo "VIRAL CONTEXTS"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92H0b_viral_context_counts.tsv"


echo
echo "======================================================================"
echo "TOPOLOGIES"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92H0b_topology_counts.tsv"


echo
echo "======================================================================"
echo "DIRECT vOTU -> MAG LINKS"
echo "======================================================================"

head -n 30 \
    "${OUT}/92H0b_direct_votu_MAG_links.tsv" \
| column -t -s $'\t'


echo
echo "======================================================================"
echo "TAXONOMY CANDIDATES"
echo "======================================================================"

cat \
    "${OUT}/92H0b_taxonomy_candidate_preview.tsv"


echo
echo "======================================================================"
echo "HOST TOOLS"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/92H0b_host_tools_PATH.tsv"

cat \
    "${OUT}/92H0b_host_container_candidates.tsv"


echo
echo "======================================================================"
echo "92H0b FINALIZO CORRECTAMENTE."
echo "No se modificaron FASTA ni resultados previos."
echo "ES SEGURO SALIR."
echo "======================================================================"

