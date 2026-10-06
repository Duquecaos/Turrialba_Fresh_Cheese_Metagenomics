#!/bin/bash

# Resolve helper scripts relative to this staged module.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


set -euo pipefail

export LC_ALL=C
export LANG=C
export LANGUAGE=C

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

OUT="${ROOT}/98_read_taxonomy/97C2B_reference_panels"
CACHE="${OUT}/cache"
GENOMES="${OUT}/genomes"
PANELS="${OUT}/panels"
INDEX="${OUT}/bowtie2_indexes"

mkdir -p \
    "$OUT" \
    "$CACHE" \
    "$GENOMES" \
    "$PANELS" \
    "$INDEX"


APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"
BOWTIE2_SIF="${BOWTIE2_SIF:-/opt/ohpc/pub/containers/BIO/bowtie2-2.5.4.sif}"

THREADS="${SLURM_CPUS_PER_TASK:-4}"


echo "============================================================"
echo "97C2B - COMPETITIVE REFSEQ REFERENCE PANELS"
echo "START: $(date -Iseconds)"
echo "============================================================"


# ============================================================
# 1. Preconditions
# ============================================================

command -v curl >/dev/null 2>&1 || {
    echo "ERROR: curl not found" >&2
    exit 10
}

[[ -x "$APPTAINER" ]] || {
    echo "ERROR: Apptainer missing: $APPTAINER" >&2
    exit 11
}

[[ -s "$BOWTIE2_SIF" ]] || {
    echo "ERROR: Bowtie2 container missing: $BOWTIE2_SIF" >&2
    exit 12
}


# ============================================================
# 2. Requested competitive panel definition
#
# role:
# target   = species being validated
# neighbor = competing close/co-occurring taxon
# ============================================================

cat > "${OUT}/97C2B_requested_species.tsv" <<'EOF'
panel	role	species
salmonella	target	Salmonella enterica
salmonella	neighbor	Salmonella bongori
salmonella	neighbor	Escherichia coli
salmonella	neighbor	Citrobacter freundii
salmonella	neighbor	Enterobacter hormaechei
salmonella	neighbor	Klebsiella pneumoniae
salmonella	neighbor	Pantoea agglomerans
salmonella	neighbor	Serratia marcescens
listeria	target	Listeria monocytogenes
listeria	neighbor	Listeria innocua
listeria	neighbor	Listeria ivanovii
listeria	neighbor	Listeria seeligeri
listeria	neighbor	Listeria welshimeri
staph_aureus	target	Staphylococcus aureus
staph_aureus	neighbor	Staphylococcus argenteus
staph_aureus	neighbor	Staphylococcus schweitzeri
staph_aureus	neighbor	Staphylococcus epidermidis
staph_aureus	neighbor	Staphylococcus haemolyticus
staph_aureus	neighbor	Staphylococcus saprophyticus
staph_aureus	neighbor	Mammaliicoccus sciuri
ecoli	target	Escherichia coli
ecoli	neighbor	Escherichia albertii
ecoli	neighbor	Escherichia fergusonii
ecoli	neighbor	Shigella flexneri
ecoli	neighbor	Shigella sonnei
ecoli	neighbor	Shigella dysenteriae
ecoli	neighbor	Shigella boydii
bcereus	target	Bacillus cereus
bcereus	neighbor	Bacillus thuringiensis
bcereus	neighbor	Bacillus anthracis
bcereus	neighbor	Bacillus mycoides
bcereus	neighbor	Bacillus toyonensis
bcereus	neighbor	Bacillus wiedmannii
campylobacter	target	Campylobacter jejuni
campylobacter	neighbor	Campylobacter coli
campylobacter	neighbor	Campylobacter lari
campylobacter	neighbor	Campylobacter upsaliensis
campylobacter	neighbor	Campylobacter fetus
cperfringens	target	Clostridium perfringens
cperfringens	neighbor	Clostridium baratii
cperfringens	neighbor	Clostridium septicum
cperfringens	neighbor	Clostridium novyi
cperfringens	neighbor	Clostridium botulinum
cronobacter	target	Cronobacter sakazakii
cronobacter	neighbor	Cronobacter malonaticus
cronobacter	neighbor	Cronobacter turicensis
cronobacter	neighbor	Cronobacter dublinensis
cronobacter	neighbor	Cronobacter universalis
cronobacter	neighbor	Cronobacter muytjensii
yersinia	target	Yersinia enterocolitica
yersinia	neighbor	Yersinia pseudotuberculosis
yersinia	neighbor	Yersinia pestis
yersinia	neighbor	Yersinia frederiksenii
yersinia	neighbor	Yersinia intermedia
yersinia	neighbor	Yersinia kristensenii
EOF


# ============================================================
# 3. Download current RefSeq assembly summary
# ============================================================

ASSEMBLY_URL="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
ASSEMBLY="${CACHE}/assembly_summary_refseq.txt"

curl \
    -L \
    --fail \
    --retry 5 \
    --retry-delay 5 \
    --connect-timeout 30 \
    --max-time 1800 \
    -o "${ASSEMBLY}.tmp" \
    "$ASSEMBLY_URL"

mv \
    "${ASSEMBLY}.tmp" \
    "$ASSEMBLY"

[[ -s "$ASSEMBLY" ]] || {
    echo "ERROR: empty assembly summary" >&2
    exit 20
}

sha256sum \
    "$ASSEMBLY" \
    > "${OUT}/97C2B_assembly_summary.sha256"


# ============================================================
# 4. Select deterministic RefSeq assemblies
# ============================================================

REF_SELECTOR="${SCRIPT_DIR}/97c2b_select_refseq_references.py"

if [[ ! -x "$REF_SELECTOR" ]]
then
    echo "ERROR: missing executable RefSeq selector: $REF_SELECTOR" >&2
    exit 22
fi


python3 \
    "$REF_SELECTOR" \
    "$ASSEMBLY" \
    "${OUT}/97C2B_requested_species.tsv" \
    "${OUT}/97C2B_selected_references.tsv" \
    "${OUT}/97C2B_missing_requested_species.tsv" \
    "${OUT}/97C2B_assembly_parser_audit.tsv"


if [[ ! -s "${OUT}/97C2B_selected_references.tsv" ]]
then
    echo "ERROR: no selected references produced" >&2
    exit 23
fi


if [[ ! -s "${OUT}/97C2B_assembly_parser_audit.tsv" ]]
then
    echo "ERROR: parser audit missing" >&2
    exit 24
fi
# 5. Download selected genomes
# ============================================================

printf \
"panel\trole\tspecies\tassembly_accession\tgz_path\tstatus\n" \
> "${OUT}/97C2B_download_status.tsv"


tail -n +2 \
"${OUT}/97C2B_selected_references.tsv" \
| while IFS=$'\t' read -r \
    PANEL \
    ROLE \
    SPECIES \
    ACCESSION \
    REFSEQ_CATEGORY \
    ORGANISM \
    TAXID \
    SPECIES_TAXID \
    LEVEL \
    RELEASE_DATE \
    FTP_PATH \
    URL

do

    SLUG=$(
        printf '%s' "$SPECIES" \
        | tr ' /()[]:' '_' \
        | tr -cd 'A-Za-z0-9_.-'
    )

    GZ="${GENOMES}/${SLUG}__${ACCESSION}.genomic.fna.gz"


    if [[ ! -s "$GZ" ]]
    then

        echo "Downloading $SPECIES $ACCESSION"

        curl \
            -L \
            --fail \
            --retry 5 \
            --retry-delay 5 \
            --connect-timeout 30 \
            --max-time 1800 \
            -o "${GZ}.tmp" \
            "$URL"

        mv \
            "${GZ}.tmp" \
            "$GZ"

    fi


    if gzip -t "$GZ"
    then
        STATUS="PASS"
    else
        STATUS="FAIL"
    fi


    printf \
"%s\t%s\t%s\t%s\t%s\t%s\n" \
"$PANEL" \
"$ROLE" \
"$SPECIES" \
"$ACCESSION" \
"$GZ" \
"$STATUS" \
>> "${OUT}/97C2B_download_status.tsv"

done


DOWNLOAD_FAILS=$(
    awk -F'\t' '
        NR > 1 && $6 != "PASS" {
            n++
        }
        END {
            print n+0
        }
    ' "${OUT}/97C2B_download_status.tsv"
)

if [[ "$DOWNLOAD_FAILS" -ne 0 ]]
then
    echo "ERROR: genome download validation failures=$DOWNLOAD_FAILS" >&2
    exit 30
fi


# ============================================================
# 6. Build chromosome-focused per-species FASTA
#
# Keep replicons >= 1,000,000 bp.
# This deliberately removes most plasmids/small replicons.
# ============================================================

python3 - \
"${OUT}/97C2B_selected_references.tsv" \
"$GENOMES" \
"${OUT}/97C2B_chromosome_reference_summary.tsv" <<'PY'
import csv
import gzip
import re
import sys
from pathlib import Path


manifest = Path(sys.argv[1])
genome_dir = Path(sys.argv[2])
summary_out = Path(sys.argv[3])

MIN_LEN = 1_000_000


def slug(s):

    x = re.sub(
        r"[^A-Za-z0-9_.-]+",
        "_",
        s
    )

    return x.strip("_")


def fasta_iter(handle):

    header = None
    seq = []

    for line in handle:

        line = line.rstrip(
            "\n"
        )

        if line.startswith(">"):

            if header is not None:
                yield header, "".join(seq)

            header = line[1:]
            seq = []

        else:

            seq.append(
                line.strip()
            )

    if header is not None:
        yield header, "".join(seq)


with manifest.open() as fh:

    rows = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )


# Unique species/accession combinations
unique = {}

for r in rows:

    key = (
        r["requested_species"],
        r["assembly_accession"],
    )

    unique[key] = r


summary = []


for (
    species,
    accession
), r in sorted(
    unique.items()
):

    species_slug = slug(
        species
    )

    gz = (
        genome_dir
        / f"{species_slug}__{accession}.genomic.fna.gz"
    )

    out = (
        genome_dir
        / f"{species_slug}__{accession}.chromosome.fna"
    )


    kept = 0
    kept_bp = 0
    total_seq = 0
    total_bp = 0


    with gzip.open(
        gz,
        "rt",
        errors="replace"
    ) as ih, out.open(
        "w"
    ) as oh:

        for header, seq in fasta_iter(
            ih
        ):

            total_seq += 1
            total_bp += len(seq)

            if len(seq) < MIN_LEN:
                continue

            kept += 1
            kept_bp += len(seq)

            orig_id = header.split()[0]

            new_header = (
                f"{species_slug}"
                f"|{accession}"
                f"|{orig_id}"
            )

            oh.write(
                f">{new_header}\n"
            )

            for i in range(
                0,
                len(seq),
                80
            ):

                oh.write(
                    seq[
                        i:i+80
                    ]
                    + "\n"
                )


    if kept == 0:

        raise RuntimeError(
            f"No >=1 Mb replicon retained for {species} {accession}"
        )


    summary.append({
        "species":
            species,

        "species_slug":
            species_slug,

        "assembly_accession":
            accession,

        "all_replicons":
            total_seq,

        "all_bp":
            total_bp,

        "retained_replicons_ge1Mb":
            kept,

        "retained_bp":
            kept_bp,

        "chromosome_fasta":
            str(out),
    })


fields = list(
    summary[0].keys()
)


with summary_out.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    w.writeheader()
    w.writerows(summary)


print(
    f"UNIQUE_REFERENCES={len(summary)}"
)
PY


# ============================================================
# 7. Build each competitive panel FASTA
# ============================================================

python3 - \
"${OUT}/97C2B_selected_references.tsv" \
"${OUT}/97C2B_chromosome_reference_summary.tsv" \
"$PANELS" \
"${OUT}/97C2B_panel_summary.tsv" <<'PY'
import csv
import shutil
import sys
from collections import defaultdict
from pathlib import Path


manifest = Path(sys.argv[1])
chrom_summary = Path(sys.argv[2])
panel_dir = Path(sys.argv[3])
summary_out = Path(sys.argv[4])


with manifest.open() as fh:

    rows = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )


with chrom_summary.open() as fh:

    chrom_rows = list(
        csv.DictReader(
            fh,
            delimiter="\t"
        )
    )


chrom_lookup = {
    (
        r["species"],
        r["assembly_accession"]
    ):
        r

    for r in chrom_rows
}


by_panel = defaultdict(list)

for r in rows:
    by_panel[
        r["panel"]
    ].append(r)


summary = []


for panel, rr in sorted(
    by_panel.items()
):

    targets = [
        r
        for r in rr
        if r["role"] == "target"
    ]

    if len(targets) != 1:

        raise RuntimeError(
            f"{panel}: expected exactly one target, found {len(targets)}"
        )


    if len(rr) < 3:

        raise RuntimeError(
            f"{panel}: fewer than 3 reference species selected"
        )


    outfile = (
        panel_dir
        / f"{panel}.competitive.fna"
    )


    total_bp = 0
    total_replicons = 0


    with outfile.open(
        "wb"
    ) as oh:

        for r in sorted(
            rr,
            key=lambda x: (
                x["role"] != "target",
                x["requested_species"]
            )
        ):

            key = (
                r["requested_species"],
                r["assembly_accession"]
            )

            c = chrom_lookup[
                key
            ]

            fp = Path(
                c[
                    "chromosome_fasta"
                ]
            )

            with fp.open(
                "rb"
            ) as ih:

                shutil.copyfileobj(
                    ih,
                    oh
                )

            total_bp += int(
                c[
                    "retained_bp"
                ]
            )

            total_replicons += int(
                c[
                    "retained_replicons_ge1Mb"
                ]
            )


    summary.append({
        "panel":
            panel,

        "target_species":
            targets[0][
                "requested_species"
            ],

        "species_in_panel":
            len(rr),

        "neighbor_species":
            len(rr) - 1,

        "retained_replicons":
            total_replicons,

        "total_reference_bp":
            total_bp,

        "panel_fasta":
            str(outfile),
    })


fields = list(
    summary[0].keys()
)


with summary_out.open(
    "w",
    newline=""
) as fh:

    w = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    w.writeheader()
    w.writerows(summary)


if len(summary) != 9:

    raise RuntimeError(
        f"Expected 9 panels; found {len(summary)}"
    )


print(
    f"PANELS={len(summary)}"
)
PY


# ============================================================
# 8. Build Bowtie2 indexes
# ============================================================

tail -n +2 \
"${OUT}/97C2B_panel_summary.tsv" \
| while IFS=$'\t' read -r \
    PANEL \
    TARGET \
    SPECIES_N \
    NEIGHBOR_N \
    REPLICONS \
    TOTAL_BP \
    PANEL_FASTA

do

    PREFIX="${INDEX}/${PANEL}"


    if [[ ! -s "${PREFIX}.1.bt2" && ! -s "${PREFIX}.1.bt2l" ]]
    then

        echo "Building Bowtie2 index: $PANEL"

        "$APPTAINER" exec \
            --bind "${ROOT}:${ROOT}" \
            "$BOWTIE2_SIF" \
            bowtie2-build \
            --threads "$THREADS" \
            "$PANEL_FASTA" \
            "$PREFIX"

    fi

done


# ============================================================
# 9. Validate indexes
# ============================================================

printf \
"panel\tindex_prefix\tstatus\n" \
> "${OUT}/97C2B_index_status.tsv"


for PANEL in \
    salmonella \
    listeria \
    staph_aureus \
    ecoli \
    bcereus \
    campylobacter \
    cperfringens \
    cronobacter \
    yersinia

do

    PREFIX="${INDEX}/${PANEL}"

    SHORT_OK=1
    LONG_OK=1


    for EXT in \
        1.bt2 \
        2.bt2 \
        3.bt2 \
        4.bt2 \
        rev.1.bt2 \
        rev.2.bt2
    do

        if [[ ! -s "${PREFIX}.${EXT}" ]]
        then
            SHORT_OK=0
        fi

    done


    for EXT in \
        1.bt2l \
        2.bt2l \
        3.bt2l \
        4.bt2l \
        rev.1.bt2l \
        rev.2.bt2l
    do

        if [[ ! -s "${PREFIX}.${EXT}" ]]
        then
            LONG_OK=0
        fi

    done


    if [[ "$SHORT_OK" -eq 1 || "$LONG_OK" -eq 1 ]]
    then
        STATUS="PASS"
    else
        STATUS="FAIL"
    fi


    printf \
"%s\t%s\t%s\n" \
"$PANEL" \
"$PREFIX" \
"$STATUS" \
>> "${OUT}/97C2B_index_status.tsv"

done


INDEX_FAILS=$(
    awk -F'\t' '
        NR > 1 && $3 != "PASS" {
            n++
        }
        END {
            print n+0
        }
    ' "${OUT}/97C2B_index_status.tsv"
)


if [[ "$INDEX_FAILS" -ne 0 ]]
then
    echo "ERROR: Bowtie2 index failures=$INDEX_FAILS" >&2
    exit 40
fi


# ============================================================
# 10. Checksums
# ============================================================

find \
    "$GENOMES" \
    "$PANELS" \
    -type f \
    \( \
        -name '*.fna' \
        -o -name '*.fna.gz' \
    \) \
    -print0 \
| sort -z \
| xargs -0 sha256sum \
> "${OUT}/97C2B_reference_files.sha256"


# ============================================================
# 11. Global summary
# ============================================================

REQUESTED=$(
    tail -n +2 \
    "${OUT}/97C2B_requested_species.tsv" \
    | wc -l
)


SELECTED=$(
    tail -n +2 \
    "${OUT}/97C2B_selected_references.tsv" \
    | wc -l
)


MISSING=$(
    tail -n +2 \
    "${OUT}/97C2B_missing_requested_species.tsv" \
    | wc -l
)


TOTAL_BP=$(
    awk -F'\t' '
        NR > 1 {
            x += $6
        }
        END {
            print x+0
        }
    ' "${OUT}/97C2B_panel_summary.tsv"
)


{
    printf "metric\tvalue\n"

    printf "requested_panel_species_rows\t%s\n" "$REQUESTED"
    printf "selected_reference_rows\t%s\n" "$SELECTED"
    printf "missing_optional_reference_rows\t%s\n" "$MISSING"

    printf "competitive_panels\t9\n"
    printf "target_species\t9\n"

    printf "reference_selection\tone_latest_complete_RefSeq_genome_per_species\n"
    printf "selection_priority\treference_then_representative_then_other_complete\n"

    printf "minimum_replicon_length_bp\t1000000\n"
    printf "small_replicons_plasmids_primary_panel\tEXCLUDED_BY_SIZE\n"

    printf "sum_panel_reference_bp\t%s\n" "$TOTAL_BP"

    printf "bowtie2_indexes\t9\n"
    printf "index_failures\t%s\n" "$INDEX_FAILS"

    printf "reads_mapped\tNO\n"
    printf "pathogen_confirmation\tNO\n"

    printf "next_step\t97C2C_competitive_mapping\n"

} > "${OUT}/97C2B_global_summary.tsv"


# ============================================================
# 12. Methodological scope
# ============================================================

cat > "${OUT}/97C2B_methodological_scope.tsv" <<'EOF'
field	value
analysis	competitive_reference_panel_preparation
source	NCBI_RefSeq
reference_selection	one_complete_genome_per_species
reference_priority	reference_genome_then_representative_genome_then_other_complete_genome
competitor_design	target_plus_close_and_or_common_neighbor_species
small_replicons	removed_by_ge1Mb_filter
plasmid_based_validation	NO
mapping_performed	NO
target_detection_claim	NO
limitation	single_reference_per_species_does_not_capture_all_intraspecies_diversity
Ecoli_Shigella	known_short_read_resolution_limitation
Bacillus_cereus_group	known_species_resolution_limitation
competitive_mapping_goal	high_MAPQ_target_specific_chromosomal_support
amplicon_data_used	NO
EOF


# ============================================================
# 13. Final validation
# ============================================================

for F in \
    "${OUT}/97C2B_requested_species.tsv" \
    "${OUT}/97C2B_selected_references.tsv" \
    "${OUT}/97C2B_missing_requested_species.tsv" \
    "${OUT}/97C2B_download_status.tsv" \
    "${OUT}/97C2B_chromosome_reference_summary.tsv" \
    "${OUT}/97C2B_panel_summary.tsv" \
    "${OUT}/97C2B_index_status.tsv" \
    "${OUT}/97C2B_reference_files.sha256" \
    "${OUT}/97C2B_global_summary.tsv" \
    "${OUT}/97C2B_methodological_scope.tsv"
do

    if [[ ! -s "$F" ]]
    then
        echo "ERROR: missing/empty $F" >&2
        exit 50
    fi

done


printf \
"COMPLETED\t%s\n" \
"$(date -Iseconds)" \
> "${OUT}/97C2B_COMPLETE.ok"


echo
echo "============================================================"
echo "97C2B GLOBAL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OUT}/97C2B_global_summary.tsv"


echo
echo "============================================================"
echo "97C2B FINALIZO CORRECTAMENTE"
echo "REFERENCE PANELS + BOWTIE2 INDEXES READY"
echo "ES SEGURO SALIR"
echo "============================================================"
