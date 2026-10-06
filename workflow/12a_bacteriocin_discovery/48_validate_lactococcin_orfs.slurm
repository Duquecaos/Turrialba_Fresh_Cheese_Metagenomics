#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00

#SBATCH --job-name=lactococcin_validate
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"
OUT="${ROOT}/30_lactococcin_validation"

ANTISMASH="${ANTISMASH:-/opt/ohpc/pub/containers/BIO/antismash-7.1.0.sif}"

mkdir -p "${OUT}"

echo "============================================================"
echo "VALIDACIÓN ORFs LACTOCOCCIN-LIKE"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "============================================================"

srun --exact -n 1 -c 1 \
apptainer exec \
    --cleanenv \
    --bind "${ROOT}:/project:ro" \
    --bind "${OUT}:/out" \
    "${ANTISMASH}" \
python3 - <<'PY'

from Bio import SeqIO
from Bio import pairwise2
from pathlib import Path

ROOT = Path("/project")
OUT = Path("/out")

# ============================================================
# Referencia lactococcina B
# UniProt P35518
# ============================================================

LCNB_PRECURSOR = (
    "MKNQLNFNIVSDEELAEVNGG"
    "SLQYVMSAGPYTWYKDTRTGKTICKQTIDTASYTFGVMAEGWGKTFH"
)

LCNB_MATURE = (
    "SLQYVMSAGPYTWYKDTRTGKTICKQTIDTASYTFGVMAEGWGKTFH"
)

LCNB_TRUNC_7_47 = (
    "SAGPYTWYKDTRTGKTICKQTIDTASYTFGVMAEGWGKTFH"
)


# ============================================================
# Archivos candidatos
# ============================================================

candidates = {
    "petauri_ctg14_2": (
        ROOT /
        "26_antismash_rep18/results/"
        "L2__L2_maxbin2.004_sub/"
        "k141_179388.region001.gbk",
        "ctg14_2"
    ),

    "atlantibacter_ctg195_2": (
        ROOT /
        "26_antismash_rep18/results/"
        "L2__L2_maxbin2.011_sub/"
        "k141_100305.region001.gbk",
        "ctg195_2"
    ),

    "lactis_ctg152_6": (
        ROOT /
        "26_antismash_rep18/results/"
        "M2__M2_maxbin2.004/"
        "k141_64662.region001.gbk",
        "ctg152_6"
    ),

    "lactis_ctg152_9": (
        ROOT /
        "26_antismash_rep18/results/"
        "M2__M2_maxbin2.004/"
        "k141_64662.region001.gbk",
        "ctg152_9"
    ),

    "lactis_ctg25_2": (
        ROOT /
        "26_antismash_rep18/results/"
        "M2__M2_maxbin2.004/"
        "k141_78158.region001.gbk",
        "ctg25_2"
    ),
}


def get_cds(gbk, locus):

    rec = SeqIO.read(str(gbk), "genbank")

    for f in rec.features:

        if f.type != "CDS":
            continue

        if locus in f.qualifiers.get("locus_tag", []):
            return rec, f

    raise RuntimeError(
        f"No se encontró {locus} en {gbk}"
    )


seqs = {}
ntseqs = {}


# ============================================================
# Candidatos
# ============================================================

print("===== CINCO CANDIDATOS =====")

for name, (gbk, locus) in candidates.items():

    rec, f = get_cds(gbk, locus)

    aa = "".join(
        f.qualifiers.get("translation", [])
    )

    nt = str(
        f.extract(rec.seq)
    )

    seqs[name] = aa
    ntseqs[name] = nt

    print()
    print(name)
    print("AA:", aa)
    print("Longitud AA:", len(aa))
    print("NT:", nt)
    print("Longitud NT:", len(nt))


# ============================================================
# ctg152_6 versus lactococcina B
# ============================================================

print()
print("=" * 80)
print("ctg152_6 vs LACTOCOCCINA B")
print("=" * 80)

query = seqs["lactis_ctg152_6"]

print("Referencia precursor:")
print(LCNB_PRECURSOR)

print()
print("Referencia madura:")
print(LCNB_MATURE)

print()
print("ctg152_6:")
print(query)

print()
print(
    "ctg152_6[1:] == LcnB madura residuos 7-47:",
    query[1:] == LCNB_TRUNC_7_47
)

print(
    "Secuencia después de Met inicial:",
    query[1:]
)


# ============================================================
# Examinar 180 nt upstream del inicio de ctg152_6
# manteniendo el mismo marco de lectura
# ============================================================

gbk, locus = candidates["lactis_ctg152_6"]
rec, feature = get_cds(gbk, locus)

start = int(feature.location.start)
end = int(feature.location.end)

UPSTREAM = 180

window_start = max(
    0,
    start - UPSTREAM
)

# Ajustar para que la distancia hasta el CDS
# sea múltiplo de 3
delta = start - window_start
window_start += delta % 3

window = rec.seq[
    window_start:end
]

translated = str(
    window.translate(to_stop=False)
)

cds_aa_start = (
    start - window_start
) // 3

upstream_aa = translated[:cds_aa_start]
cds_plus_upstream = translated

print()
print("=" * 80)
print("CONTEXTO EN MARCO DE ctg152_6")
print("=" * 80)

print(
    "Coordenada inicio CDS:",
    start + 1
)

print(
    "Ventana upstream inicia:",
    window_start + 1
)

print(
    "AA upstream mismo marco:"
)

print(upstream_aa)

print()
print(
    "Últimos 30 aa antes del inicio:"
)

print(
    upstream_aa[-30:]
)

print()
print(
    "Traducción upstream + CDS:"
)

print(cds_plus_upstream)

last_stop = upstream_aa.rfind("*")

if last_stop >= 0:
    since_stop = upstream_aa[
        last_stop + 1:
    ]
else:
    since_stop = upstream_aa

print()
print(
    "AA desde último STOP antes del CDS:"
)

print(since_stop)

print()
print(
    "¿Aparece precursor completo LcnB en ventana?:",
    LCNB_PRECURSOR in cds_plus_upstream
)

print(
    "¿Aparece motivo líder+mature "
    "'AEVNGGSLQYV'?:",
    "AEVNGGSLQYV" in cds_plus_upstream
)

print(
    "¿Aparece 'GGSLQYV'?:",
    "GGSLQYV" in cds_plus_upstream
)


# ============================================================
# ctg152_9 vs ctg195_2
# ============================================================

print()
print("=" * 80)
print("ctg152_9 vs ctg195_2")
print("=" * 80)

a = seqs["lactis_ctg152_9"]
b = seqs["atlantibacter_ctg195_2"]

na = ntseqs["lactis_ctg152_9"]
nb = ntseqs["atlantibacter_ctg195_2"]

print(
    "Proteína exactamente idéntica:",
    a == b
)

print(
    "ADN exactamente idéntico:",
    na == nb
)

print()
print("Lactococcus AA:")
print(a)

print("Atlantibacter AA:")
print(b)

print()
print("Lactococcus NT:")
print(na)

print("Atlantibacter NT:")
print(nb)


# ============================================================
# Comparación por pares de proteínas
# ============================================================

print()
print("=" * 80)
print("IDENTIDAD GLOBAL ENTRE LOS CINCO ORFs")
print("=" * 80)

names = list(seqs)

for i in range(len(names)):

    for j in range(i + 1, len(names)):

        n1 = names[i]
        n2 = names[j]

        s1 = seqs[n1]
        s2 = seqs[n2]

        aln = pairwise2.align.globalxx(
            s1,
            s2,
            one_alignment_only=True
        )[0]

        aligned1 = aln.seqA
        aligned2 = aln.seqB

        matches = sum(
            x == y and x != "-"
            for x, y
            in zip(aligned1, aligned2)
        )

        denom = max(
            len(s1),
            len(s2)
        )

        identity = (
            100 * matches / denom
        )

        print(
            f"{n1}\t{n2}\t"
            f"{matches}/{denom}\t"
            f"{identity:.2f}%"
        )


# ============================================================
# Guardar resumen
# ============================================================

with (
    OUT /
    "lactococcin_sequences.tsv"
).open("w") as fh:

    fh.write(
        "candidate\tAA_sequence\t"
        "length_aa\tNT_sequence\n"
    )

    for name in names:

        fh.write(
            f"{name}\t"
            f"{seqs[name]}\t"
            f"{len(seqs[name])}\t"
            f"{ntseqs[name]}\n"
        )

PY

echo
echo "============================================================"
echo "VALIDACIÓN LACTOCOCCIN COMPLETADA"
echo "Resultados: ${OUT}"
echo "Fin: $(date)"
echo "============================================================"
