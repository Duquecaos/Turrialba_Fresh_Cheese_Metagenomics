#!/bin/bash

#SBATCH --job-name=92c_genomad
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=23:30:00
#SBATCH --array=0-5%2

#SBATCH --mail-type=END,FAIL

#SBATCH --output=92c_genomad_%A_%a.out
#SBATCH --error=92c_genomad_%A_%a.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"

SIF="${ROOT}/software/genomad/genomad-1.11.0.sif"
DB="${ROOT}/databases/genomad/genomad_db"

OUT="${ROOT}/85_genomad_coassemblies"
INPUTDIR="${OUT}/inputs"
MAPDIR="${OUT}/contig_maps"

mkdir -p \
    "$OUT" \
    "$INPUTDIR" \
    "$MAPDIR"

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

# IMPORTANTE:
# no usar GROUPS: es variable especial de Bash
COASSEMBLIES=(
    L1
    L2
    L3
    M1
    M2
    M3
)

G="${COASSEMBLIES[$SLURM_ARRAY_TASK_ID]}"

SRC="${ROOT}/05_assembly/${G}/final.contigs.fa"

STAGED="${INPUTDIR}/${G}.fna"
CONTIG_MAP="${MAPDIR}/${G}.contig_map.tsv"

GOUT="${OUT}/${G}"

THREADS="$SLURM_CPUS_PER_TASK"
SPLITS=16

echo "======================================================================"
echo "PASO 92C - geNomad"
echo "Coensamblaje: ${G}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Job: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Inicio: $(date)"
echo "Threads: ${THREADS}"
echo "Splits: ${SPLITS}"
echo "======================================================================"

# ======================================================================
# 1. VALIDACIONES
# ======================================================================

[[ -x "$APPTAINER" ]] || {
    echo "ERROR: falta Apptainer." >&2
    exit 2
}

[[ -s "$SIF" ]] || {
    echo "ERROR: falta geNomad SIF: $SIF" >&2
    exit 3
}

[[ -d "$DB" ]] || {
    echo "ERROR: falta geNomad DB: $DB" >&2
    exit 4
}

[[ -s "$SRC" ]] || {
    echo "ERROR: falta assembly: $SRC" >&2
    exit 5
}


# ======================================================================
# 2. CREAR FASTA CON IDs GLOBALMENTE UNICOS
# ======================================================================

echo
echo "======================================================================"
echo "1. PREPARANDO FASTA CON IDs UNICOS"
echo "======================================================================"

TMP_FASTA="${STAGED}.tmp.${SLURM_JOB_ID}"

python3 - \
    "$SRC" \
    "$TMP_FASTA" \
    "$CONTIG_MAP" \
    "$G" <<'PY'

import sys
from pathlib import Path

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
mapfile = Path(sys.argv[3])
group = sys.argv[4]

seen = set()
n_contigs = 0
total_bp = 0

with (
    src.open() as inp,
    dst.open("w") as out,
    mapfile.open("w") as mp
):

    mp.write(
        "coassembly\tprefixed_contig\toriginal_contig\n"
    )

    for line in inp:

        if line.startswith(">"):

            original = line[1:].strip().split()[0]

            if original in seen:
                raise RuntimeError(
                    f"Contig ID duplicado en {group}: {original}"
                )

            seen.add(original)

            prefixed = f"{group}__{original}"

            out.write(f">{prefixed}\n")

            mp.write(
                f"{group}\t{prefixed}\t{original}\n"
            )

            n_contigs += 1

        else:

            seq = line.strip()

            if seq:
                total_bp += len(seq)
                out.write(seq + "\n")

print("coassembly:", group)
print("contigs:", n_contigs)
print("total_bp:", total_bp)

if n_contigs == 0:
    raise RuntimeError("FASTA sin contigs")

PY

mv "$TMP_FASTA" "$STAGED"

[[ -s "$STAGED" ]] || {
    echo "ERROR: FASTA staged vacío." >&2
    exit 6
}

[[ -s "$CONTIG_MAP" ]] || {
    echo "ERROR: mapa de IDs vacío." >&2
    exit 7
}

echo
echo "SHA256 original:"
sha256sum "$SRC"

echo
echo "SHA256 staged:"
sha256sum "$STAGED"

echo
echo "Primeros IDs:"
grep '^>' "$STAGED" \
    | head -n 5 \
    || true


# ======================================================================
# 3. EVITAR SOBREESCRIBIR UNA EJECUCION COMPLETA
# ======================================================================

VSUM="${GOUT}/${G}_summary/${G}_virus_summary.tsv"
PSUM="${GOUT}/${G}_summary/${G}_plasmid_summary.tsv"

if [[ -s "$VSUM" && -s "$PSUM" ]]
then

    echo
    echo "======================================================================"
    echo "RESULTADO COMPLETO PREEXISTENTE"
    echo "No se recalculará ${G}."
    echo "======================================================================"

else

    if [[ -d "$GOUT" ]]
    then
        echo "Eliminando salida incompleta previa:"
        echo "$GOUT"
        rm -rf "$GOUT"
    fi

    # ==================================================================
    # 4. geNomad
    # ==================================================================

    echo
    echo "======================================================================"
    echo "2. EJECUTANDO geNomad END-TO-END"
    echo "======================================================================"

    # No --cleanup:
    # conservamos resultados intermedios para análisis posteriores.
    #
    # --enable-score-calibration:
    # permite estimación de FDR.
    #
    # --splits:
    # reduce uso de memoria en búsqueda MMseqs2.

    srun --exact \
        -n 1 \
        -c "$SLURM_CPUS_PER_TASK" \
        "$APPTAINER" run \
        --bind "${ROOT}:/work" \
        "$SIF" \
        end-to-end \
        --threads "$THREADS" \
        --splits "$SPLITS" \
        --enable-score-calibration \
        "/work/85_genomad_coassemblies/inputs/${G}.fna" \
        "/work/85_genomad_coassemblies/${G}" \
        "/work/databases/genomad/genomad_db"

fi


# ======================================================================
# 5. VALIDAR SALIDAS
# ======================================================================

echo
echo "======================================================================"
echo "3. VALIDACION DE SALIDAS"
echo "======================================================================"

VSUM="${GOUT}/${G}_summary/${G}_virus_summary.tsv"
VFNA="${GOUT}/${G}_summary/${G}_virus.fna"
VGEN="${GOUT}/${G}_summary/${G}_virus_genes.tsv"
VFAA="${GOUT}/${G}_summary/${G}_virus_proteins.faa"

PSUM="${GOUT}/${G}_summary/${G}_plasmid_summary.tsv"
PFNA="${GOUT}/${G}_summary/${G}_plasmid.fna"

for F in \
    "$VSUM" \
    "$VFNA" \
    "$VGEN" \
    "$VFAA" \
    "$PSUM" \
    "$PFNA"
do
    [[ -e "$F" ]] || {
        echo "ERROR: falta salida: $F" >&2
        exit 8
    }

    echo "OK: $F"
done


# ======================================================================
# 6. RESUMEN DEL COENSAMBLAJE
# ======================================================================

echo
echo "======================================================================"
echo "4. RESUMEN ${G}"
echo "======================================================================"

python3 - \
    "$VSUM" \
    "$PSUM" \
    "$G" \
    "$OUT" \
    "$SRC" \
    "$STAGED" \
    "$SIF" <<'PY'

import csv
import hashlib
import sys
from pathlib import Path
from collections import Counter

vfile = Path(sys.argv[1])
pfile = Path(sys.argv[2])
group = sys.argv[3]
out = Path(sys.argv[4])
src = Path(sys.argv[5])
staged = Path(sys.argv[6])
sif = Path(sys.argv[7])


def sha256(path):
    h = hashlib.sha256()

    with path.open("rb") as fh:
        for block in iter(
            lambda: fh.read(1024 * 1024),
            b""
        ):
            h.update(block)

    return h.hexdigest()


def read_tsv(path):
    with path.open() as fh:
        return list(
            csv.DictReader(
                fh,
                delimiter="\t"
            )
        )


viruses = read_tsv(vfile)
plasmids = read_tsv(pfile)

topologies = Counter(
    x.get("topology", "NA")
    for x in viruses
)

n_provirus = sum(
    x.get("topology") == "Provirus"
    for x in viruses
)

n_non_provirus = (
    len(viruses) - n_provirus
)

virus_lengths = [
    int(x["length"])
    for x in viruses
    if x.get("length")
]

plasmid_lengths = [
    int(x["length"])
    for x in plasmids
    if x.get("length")
]


summary = out / f"{group}.92C_summary.tsv"

with summary.open("w") as fh:

    fh.write("metric\tvalue\n")

    values = [
        ("coassembly", group),
        ("viral_predictions", len(viruses)),
        ("provirus_predictions", n_provirus),
        (
            "non_provirus_viral_predictions",
            n_non_provirus
        ),
        ("plasmid_predictions", len(plasmids)),
        (
            "viral_total_bp",
            sum(virus_lengths)
        ),
        (
            "viral_max_bp",
            max(virus_lengths)
            if virus_lengths
            else 0
        ),
        (
            "plasmid_total_bp",
            sum(plasmid_lengths)
        ),
        (
            "input_sha256",
            sha256(src)
        ),
        (
            "staged_input_sha256",
            sha256(staged)
        ),
        (
            "container_sha256",
            sha256(sif)
        ),
        ("threads", 16),
        ("splits", 16),
        (
            "score_calibration",
            "ENABLED"
        ),
        (
            "classification_filtering",
            "genomad_default"
        ),
        (
            "cleanup",
            "NO"
        ),
    ]

    for k, v in values:
        fh.write(f"{k}\t{v}\n")

    for topology, n in sorted(
        topologies.items()
    ):
        fh.write(
            f"topology_{topology}\t{n}\n"
        )


print("coassembly:", group)
print("viral_predictions:", len(viruses))
print("proviruses:", n_provirus)
print(
    "non_provirus_viruses:",
    n_non_provirus
)
print("plasmid_predictions:", len(plasmids))

print()
print("Topologies:")

for k, v in sorted(
    topologies.items()
):
    print(k, v, sep="\t")

PY


# ======================================================================
# 7. INVENTARIO
# ======================================================================

find "$GOUT" \
    -type f \
    -printf '%P\t%s\n' \
    | sort \
    > "${OUT}/${G}.92C_file_inventory.tsv"


echo
echo "===== RESUMEN ====="

column -t -s $'\t' \
    "${OUT}/${G}.92C_summary.tsv"


echo
echo "======================================================================"
echo "PASO 92C ${G} FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "======================================================================"
