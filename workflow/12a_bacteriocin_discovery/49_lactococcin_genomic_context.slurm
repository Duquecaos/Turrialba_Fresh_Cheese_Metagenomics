#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

#SBATCH --job-name=lactococcin_context
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

OUT="${ROOT}/31_lactococcin_context"

ANTISMASH="${ANTISMASH:-/opt/ohpc/pub/containers/BIO/antismash-7.1.0.sif}"

mkdir -p "${OUT}"

echo "============================================================"
echo "CONTEXTO GENÓMICO DE CANDIDATOS LACTOCOCCIN"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "============================================================"


# ============================================================
# EXTRAER CONTIGS Y ANALIZAR STOP DE lcnB-LIKE
# ============================================================

srun --exact -n 1 -c 1 \
apptainer exec \
    --cleanenv \
    --bind "${ROOT}:/project:ro" \
    --bind "${OUT}:/out" \
    "${ANTISMASH}" \
python3 - <<'PY'

from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path

ROOT = Path("/project")
OUT = Path("/out")

LACTIS_GBK = (
    ROOT /
    "26_antismash_rep18/results/"
    "M2__M2_maxbin2.004/"
    "k141_64662.region001.gbk"
)

ATL_GBK = (
    ROOT /
    "26_antismash_rep18/results/"
    "L2__L2_maxbin2.011_sub/"
    "k141_100305.region001.gbk"
)

lactis = SeqIO.read(
    str(LACTIS_GBK),
    "genbank"
)

atl = SeqIO.read(
    str(ATL_GBK),
    "genbank"
)

print("===== CONTIGS =====")
print(
    "Lactococcus lactis:",
    lactis.id,
    len(lactis.seq),
    "bp"
)

print(
    "Atlantibacter:",
    atl.id,
    len(atl.seq),
    "bp"
)


# ============================================================
# Guardar contigs
# ============================================================

with (
    OUT / "lactis_k141_64662.fasta"
).open("w") as fh:

    fh.write(
        f">{lactis.id}\n"
        f"{str(lactis.seq)}\n"
    )

with (
    OUT / "atlantibacter_k141_100305.fasta"
).open("w") as fh:

    fh.write(
        f">{atl.id}\n"
        f"{str(atl.seq)}\n"
    )


# ============================================================
# Encontrar ctg152_6
# ============================================================

target = None

for feature in lactis.features:

    if feature.type != "CDS":
        continue

    if (
        "ctg152_6"
        in feature.qualifiers.get(
            "locus_tag",
            []
        )
    ):
        target = feature
        break

if target is None:
    raise RuntimeError(
        "No se encontró ctg152_6"
    )

if target.location.strand != 1:
    raise RuntimeError(
        "ctg152_6 no está en hebra + "
        "como se esperaba"
    )

start = int(
    target.location.start
)

end = int(
    target.location.end
)


# ============================================================
# Mismo marco: 240 nt upstream + CDS
# ============================================================

window_start = max(
    0,
    start - 240
)

# Alinear el inicio de la ventana
# con el frame de ctg152_6
delta = start - window_start
window_start += delta % 3

window = lactis.seq[
    window_start:end
]

aa = str(
    window.translate(
        to_stop=False
    )
)

cds_aa_start = (
    start - window_start
) // 3

upstream_aa = aa[
    :cds_aa_start
]

print()
print("============================================================")
print("ANÁLISIS DEL STOP EN lcnB-LIKE")
print("============================================================")

print(
    "Inicio ctg152_6:",
    start + 1
)

print(
    "Secuencia AA upstream:"
)

print(upstream_aa)

print()
print(
    "Últimos 45 aa upstream:"
)

print(
    upstream_aa[-45:]
)


# ============================================================
# Buscar patrón MKN*LNF
# ============================================================

pattern = "MKN*LNF"

idx = upstream_aa.find(
    pattern
)

if idx == -1:

    print()
    print(
        "ADVERTENCIA: no se encontró "
        "MKN*LNF automáticamente."
    )

else:

    stop_aa_index = (
        idx + 3
    )

    stop_nt_start = (
        window_start +
        stop_aa_index * 3
    )

    stop_codon = str(
        lactis.seq[
            stop_nt_start:
            stop_nt_start + 3
        ]
    )

    stop_aa = str(
        Seq(stop_codon).translate()
    )

    print()
    print(
        "Patrón encontrado:",
        pattern
    )

    print(
        "Codón STOP:",
        stop_codon
    )

    print(
        "Traducción:",
        stop_aa
    )

    print(
        "Coordenadas 1-based:",
        f"{stop_nt_start + 1}-"
        f"{stop_nt_start + 3}"
    )

    print()
    print(
        "Contexto de codones:"
    )

    left = max(
        window_start,
        stop_nt_start - 30
    )

    right = min(
        len(lactis.seq),
        stop_nt_start + 33
    )

    # mantener frame respecto al STOP
    left -= (
        stop_nt_start - left
    ) % 3

    for pos in range(
        left,
        right,
        3
    ):

        codon = str(
            lactis.seq[
                pos:pos + 3
            ]
        )

        if len(codon) < 3:
            continue

        amino = str(
            Seq(codon).translate()
        )

        marker = (
            "  <<< STOP CANDIDATO"
            if pos == stop_nt_start
            else ""
        )

        print(
            f"{pos+1:6d}-"
            f"{pos+3:<6d} "
            f"{codon} "
            f"{amino}"
            f"{marker}"
        )

    with (
        OUT /
        "lcnB_like_stop_position.tsv"
    ).open("w") as fh:

        fh.write(
            "contig\t"
            "position_start\t"
            "position_end\t"
            "codon\t"
            "translation\n"
        )

        fh.write(
            f"{lactis.id}\t"
            f"{stop_nt_start+1}\t"
            f"{stop_nt_start+3}\t"
            f"{stop_codon}\t"
            f"{stop_aa}\n"
        )


# ============================================================
# Localizar ctg152_9 y ctg195_2
# ============================================================

def find_feature(record, locus):

    for feature in record.features:

        if feature.type != "CDS":
            continue

        if locus in feature.qualifiers.get(
            "locus_tag",
            []
        ):
            return feature

    raise RuntimeError(
        f"No se encontró {locus}"
    )


f_lactis = find_feature(
    lactis,
    "ctg152_9"
)

f_atl = find_feature(
    atl,
    "ctg195_2"
)

print()
print("============================================================")
print("POSICIONES DEL GEN IDÉNTICO")
print("============================================================")

print(
    "Lactococcus ctg152_9:",
    int(f_lactis.location.start) + 1,
    "-",
    int(f_lactis.location.end),
    "strand",
    f_lactis.location.strand
)

print(
    "Atlantibacter ctg195_2:",
    int(f_atl.location.start) + 1,
    "-",
    int(f_atl.location.end),
    "strand",
    f_atl.location.strand
)


# ============================================================
# Extraer ventanas ±1000 bp alrededor del gen
# ============================================================

def write_window(
    record,
    feature,
    outfile,
    flank=1000
):

    start = max(
        0,
        int(feature.location.start)
        - flank
    )

    end = min(
        len(record.seq),
        int(feature.location.end)
        + flank
    )

    seq = record.seq[
        start:end
    ]

    with outfile.open("w") as fh:

        fh.write(
            f">{record.id}_"
            f"{start+1}_{end}\n"
            f"{str(seq)}\n"
        )

    return start + 1, end


l1, l2 = write_window(
    lactis,
    f_lactis,
    OUT /
    "lactis_ctg152_9_flank1000.fasta"
)

a1, a2 = write_window(
    atl,
    f_atl,
    OUT /
    "atlantibacter_ctg195_2_flank1000.fasta"
)

print()
print(
    "Ventana Lactococcus:",
    f"{l1}-{l2}"
)

print(
    "Ventana Atlantibacter:",
    f"{a1}-{a2}"
)

PY


# ============================================================
# BLASTn: contig completo
# ============================================================

echo
echo "============================================================"
echo "BLASTn ENTRE CONTIGS COMPLETOS"
echo "============================================================"

srun --exact -n 1 -c 1 \
apptainer exec \
    --cleanenv \
    --bind "${OUT}:/out" \
    "${ANTISMASH}" \
    blastn \
        -query /out/atlantibacter_k141_100305.fasta \
        -subject /out/lactis_k141_64662.fasta \
        -dust no \
        -word_size 7 \
        -evalue 1e-5 \
        -outfmt \
'6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore' \
        -out /out/whole_contig_blastn.tsv


# ============================================================
# BLASTn: ventanas ±1 kb
# ============================================================

echo
echo "============================================================"
echo "BLASTn ENTRE VENTANAS ±1 kb"
echo "============================================================"

srun --exact -n 1 -c 1 \
apptainer exec \
    --cleanenv \
    --bind "${OUT}:/out" \
    "${ANTISMASH}" \
    blastn \
        -query /out/atlantibacter_ctg195_2_flank1000.fasta \
        -subject /out/lactis_ctg152_9_flank1000.fasta \
        -dust no \
        -word_size 7 \
        -evalue 1e-5 \
        -outfmt \
'6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore' \
        -out /out/flank1000_blastn.tsv


echo
echo "============================================================"
echo "MEJORES ALINEAMIENTOS - CONTIG COMPLETO"
echo "============================================================"

if [[ -s "${OUT}/whole_contig_blastn.tsv" ]]; then

    head -n 20 \
    "${OUT}/whole_contig_blastn.tsv" \
    | column -t -s $'\t'

else
    echo "Sin alineamientos BLASTn."
fi


echo
echo "============================================================"
echo "MEJORES ALINEAMIENTOS - ±1 kb"
echo "============================================================"

if [[ -s "${OUT}/flank1000_blastn.tsv" ]]; then

    head -n 20 \
    "${OUT}/flank1000_blastn.tsv" \
    | column -t -s $'\t'

else
    echo "Sin alineamientos BLASTn."
fi


echo
echo "============================================================"
echo "CONTEXTO LACTOCOCCIN COMPLETADO"
echo "Resultados: ${OUT}"
echo "Fin: $(date)"
echo "============================================================"
