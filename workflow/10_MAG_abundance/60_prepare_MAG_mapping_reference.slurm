#!/bin/bash

#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00

#SBATCH --job-name=prepare_MAG_ref
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

#SBATCH --mail-type=END,FAIL

set -euo pipefail


# ============================================================
# CONFIGURACIÓN
# ============================================================

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

MAGDIR="${ROOT}/22_final_representative_mags/fastas"
READDIR="${ROOT}/04_host_removed"

BASE="${ROOT}/39_mag_abundance"
REFDIR="${BASE}/reference"
INDEXDIR="${REFDIR}/bowtie2_index"

REF="${REFDIR}/representative_18.prefixed.fna"
STB="${REFDIR}/representative_18.stb"
ORIGINAL_MAP="${REFDIR}/contig_name_map.tsv"
MAG_MANIFEST="${REFDIR}/MAG_manifest.tsv"
SAMPLE_MANIFEST="${BASE}/samples_hostfree.tsv"
STATS="${REFDIR}/reference_stats.tsv"

INDEX_PREFIX="${INDEXDIR}/representative_18"

mkdir -p "${REFDIR}" "${INDEXDIR}"


echo "============================================================"
echo "PREPARACIÓN DE REFERENCIA PARA ABUNDANCIA DE MAGs"
echo "Inicio: $(date)"
echo "Nodo: $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "============================================================"


# ============================================================
# PREFLIGHT: 18 MAGs
# ============================================================

echo
echo "============================================================"
echo "PREFLIGHT: MAGs REPRESENTATIVOS"
echo "============================================================"

mapfile -t MAGFILES < <(
    find "${MAGDIR}" \
        -maxdepth 1 \
        -type f \
        \( \
            -name '*.fa' \
            -o -name '*.fna' \
            -o -name '*.fasta' \
        \) \
        | sort
)

N_MAGS="${#MAGFILES[@]}"

echo "MAGs encontrados: ${N_MAGS}"

if [[ "${N_MAGS}" -ne 18 ]]; then
    echo "ERROR: se esperaban exactamente 18 MAGs representativos." >&2
    exit 1
fi

for F in "${MAGFILES[@]}"
do
    test -s "${F}" || {
        echo "ERROR: MAG vacío: ${F}" >&2
        exit 1
    }
done

echo "MAGs representativos: OK"


# ============================================================
# CONSTRUIR REFERENCIA CON HEADERS ÚNICOS
# ============================================================

echo
echo "============================================================"
echo "CONSTRUYENDO REFERENCIA PREFIJADA"
echo "============================================================"

: > "${REF}"

printf \
'prefixed_contig\tMAG\n' \
> "${STB}"

printf \
'prefixed_contig\toriginal_contig\tMAG\n' \
> "${ORIGINAL_MAP}"

printf \
'MAG\tsource_file\tcontigs\tbp\n' \
> "${MAG_MANIFEST}"


for F in "${MAGFILES[@]}"
do
    NAME=$(basename "${F}")

    case "${NAME}" in
        *.fasta)
            MAG="${NAME%.fasta}"
            ;;
        *.fna)
            MAG="${NAME%.fna}"
            ;;
        *.fa)
            MAG="${NAME%.fa}"
            ;;
        *)
            echo "ERROR: extensión inesperada: ${NAME}" >&2
            exit 1
            ;;
    esac

    echo "Procesando: ${MAG}"

    N_CONTIGS=$(
        awk '
            /^>/ {n++}
            END {print n+0}
        ' "${F}"
    )

    N_BP=$(
        awk '
            !/^>/ {
                gsub(/[[:space:]]/, "")
                bp += length($0)
            }
            END {print bp+0}
        ' "${F}"
    )

    printf \
        '%s\t%s\t%s\t%s\n' \
        "${MAG}" \
        "${F}" \
        "${N_CONTIGS}" \
        "${N_BP}" \
        >> "${MAG_MANIFEST}"

    awk \
        -v mag="${MAG}" \
        -v stb="${STB}" \
        -v cmap="${ORIGINAL_MAP}" \
        '
        /^>/ {
            original=$1
            sub(/^>/, "", original)

            prefixed=mag "__" original

            print ">" prefixed

            print prefixed "\t" mag >> stb
            print prefixed "\t" original "\t" mag >> cmap

            next
        }

        {
            print
        }
        ' \
        "${F}" \
        >> "${REF}"
done


# ============================================================
# VALIDACIONES DE REFERENCIA
# ============================================================

echo
echo "============================================================"
echo "VALIDANDO REFERENCIA"
echo "============================================================"

test -s "${REF}" || {
    echo "ERROR: referencia vacía." >&2
    exit 1
}

N_REF_CONTIGS=$(
    grep -c '^>' "${REF}"
)

N_STB=$(
    awk 'END {print NR-1}' "${STB}"
)

N_MAP=$(
    awk 'END {print NR-1}' "${ORIGINAL_MAP}"
)

N_UNIQUE=$(
    grep '^>' "${REF}" \
        | sed 's/^>//' \
        | sort -u \
        | wc -l
)

TOTAL_BP=$(
    awk '
        !/^>/ {
            gsub(/[[:space:]]/, "")
            bp += length($0)
        }
        END {print bp+0}
    ' "${REF}"
)

echo "MAGs: ${N_MAGS}"
echo "Contigs referencia: ${N_REF_CONTIGS}"
echo "Entradas STB: ${N_STB}"
echo "Entradas mapa: ${N_MAP}"
echo "Headers únicos: ${N_UNIQUE}"
echo "Total bp: ${TOTAL_BP}"


if [[ "${N_REF_CONTIGS}" -ne "${N_STB}" ]]; then
    echo "ERROR: referencia y STB tienen diferente número de contigs." >&2
    exit 1
fi

if [[ "${N_REF_CONTIGS}" -ne "${N_MAP}" ]]; then
    echo "ERROR: referencia y mapa tienen diferente número de contigs." >&2
    exit 1
fi

if [[ "${N_REF_CONTIGS}" -ne "${N_UNIQUE}" ]]; then
    echo "ERROR: todavía existen headers duplicados." >&2
    exit 1
fi


# ============================================================
# MANIFIESTO DE LAS 18 MUESTRAS HOST-FREE
# ============================================================

echo
echo "============================================================"
echo "CONSTRUYENDO MANIFIESTO HOST-FREE"
echo "============================================================"

printf \
'sample\tproducer\tbiological_unit\tweek\tR1\tR2\n' \
> "${SAMPLE_MANIFEST}"

N_SAMPLES=0

for PRODUCER in L M
do
    for UNIT in 1 2 3
    do
        for TIMEPOINT in 1 2 3
        do
            SAMPLE="${PRODUCER}${UNIT}_${TIMEPOINT}"
            WEEK=$((TIMEPOINT - 1))

            R1="${READDIR}/${SAMPLE}/${SAMPLE}_1.hostfree.fastq.gz"
            R2="${READDIR}/${SAMPLE}/${SAMPLE}_2.hostfree.fastq.gz"

            test -s "${R1}" || {
                echo "ERROR: falta R1: ${R1}" >&2
                exit 1
            }

            test -s "${R2}" || {
                echo "ERROR: falta R2: ${R2}" >&2
                exit 1
            }

            printf \
                '%s\t%s\t%s\t%s\t%s\t%s\n' \
                "${SAMPLE}" \
                "${PRODUCER}" \
                "${UNIT}" \
                "${WEEK}" \
                "${R1}" \
                "${R2}" \
                >> "${SAMPLE_MANIFEST}"

            N_SAMPLES=$((N_SAMPLES + 1))
        done
    done
done

if [[ "${N_SAMPLES}" -ne 18 ]]; then
    echo "ERROR: se esperaban 18 muestras." >&2
    exit 1
fi

echo "Muestras host-free: ${N_SAMPLES}"


# ============================================================
# STATS REPRODUCIBLES
# ============================================================

printf \
'metric\tvalue\n' \
> "${STATS}"

printf 'representative_MAGs\t%s\n' "${N_MAGS}" >> "${STATS}"
printf 'reference_contigs\t%s\n' "${N_REF_CONTIGS}" >> "${STATS}"
printf 'unique_contig_headers\t%s\n' "${N_UNIQUE}" >> "${STATS}"
printf 'reference_bp\t%s\n' "${TOTAL_BP}" >> "${STATS}"
printf 'hostfree_samples\t%s\n' "${N_SAMPLES}" >> "${STATS}"

sha256sum "${REF}" \
    > "${REF}.sha256"


# ============================================================
# LOCALIZAR BOWTIE2
# ============================================================

echo
echo "============================================================"
echo "LOCALIZANDO BOWTIE2"
echo "============================================================"

CONTAINER_ROOT="${CONTAINER_ROOT:-/opt/ohpc/pub/containers/BIO}"

BOWTIE2_SIF=$(
    find "${CONTAINER_ROOT}" \
        -maxdepth 1 \
        -type f \
        -iname '*bowtie2*.sif' \
        | sort -V \
        | tail -n 1
)

if [[ -z "${BOWTIE2_SIF}" ]]; then
    echo "ERROR: no se encontró contenedor Bowtie2." >&2
    exit 1
fi

if command -v apptainer >/dev/null 2>&1; then
    CONTAINER_CMD="apptainer"
elif command -v singularity >/dev/null 2>&1; then
    CONTAINER_CMD="singularity"
else
    echo "ERROR: no se encontró apptainer ni singularity." >&2
    exit 1
fi

echo "Contenedor: ${BOWTIE2_SIF}"
echo "Runtime: ${CONTAINER_CMD}"

"${CONTAINER_CMD}" exec \
    "${BOWTIE2_SIF}" \
    bowtie2-build --version


# ============================================================
# CONSTRUIR ÍNDICE BOWTIE2
# ============================================================

echo
echo "============================================================"
echo "CONSTRUYENDO ÍNDICE BOWTIE2"
echo "============================================================"

rm -f \
    "${INDEX_PREFIX}".*.bt2 \
    "${INDEX_PREFIX}".*.bt2l

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    "${CONTAINER_CMD}" exec \
        "${BOWTIE2_SIF}" \
        bowtie2-build \
            --threads "${SLURM_CPUS_PER_TASK}" \
            "${REF}" \
            "${INDEX_PREFIX}"


# ============================================================
# VALIDAR ÍNDICE
# ============================================================

N_INDEX=$(
    find "${INDEXDIR}" \
        -maxdepth 1 \
        -type f \
        \( \
            -name 'representative_18*.bt2' \
            -o -name 'representative_18*.bt2l' \
        \) \
        | wc -l
)

echo
echo "Archivos de índice Bowtie2: ${N_INDEX}"

if [[ "${N_INDEX}" -ne 6 ]]; then
    echo "ERROR: se esperaban 6 archivos del índice Bowtie2." >&2

    find "${INDEXDIR}" \
        -maxdepth 1 \
        -type f \
        -printf '%f\t%s bytes\n' \
        | sort

    exit 1
fi


# ============================================================
# RESUMEN
# ============================================================

echo
echo "============================================================"
echo "ESTADÍSTICAS DE REFERENCIA"
echo "============================================================"

column -t -s $'\t' "${STATS}"

echo
echo "============================================================"
echo "MAGs"
echo "============================================================"

column -t -s $'\t' "${MAG_MANIFEST}"

echo
echo "============================================================"
echo "MUESTRAS HOST-FREE"
echo "============================================================"

column -t -s $'\t' "${SAMPLE_MANIFEST}"

echo
echo "============================================================"
echo "ÍNDICE"
echo "============================================================"

find "${INDEXDIR}" \
    -maxdepth 1 \
    -type f \
    -printf '%f\t%s bytes\n' \
    | sort


echo
echo "============================================================"
echo "PASO 60 COMPLETADO: REFERENCIA DE 18 MAGs E ÍNDICE BOWTIE2"
echo "Referencia: ${REF}"
echo "STB: ${STB}"
echo "Manifiesto: ${SAMPLE_MANIFEST}"
echo "Índice: ${INDEX_PREFIX}"
echo "Fin: $(date)"
echo "============================================================"
