#!/bin/bash

#SBATCH --job-name=92g1_votumap
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --array=0-17%2

#SBATCH --mail-type=END,FAIL

#SBATCH --output=92g1_votumap_%A_%a.out
#SBATCH --error=92g1_votumap_%A_%a.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"

BOWTIE2_SIF="${BOWTIE2_SIF:-/opt/ohpc/pub/containers/BIO/bowtie2-2.5.4.sif}"
SAMTOOLS_SIF="${SAMTOOLS_SIF:-/opt/ohpc/pub/containers/BIO/samtools-1.21.sif}"

INDEX="${ROOT}/92_votu_abundance/index/votu_4491"
REFERENCE="${ROOT}/92_votu_abundance/reference/votu_4491.fna"

OUT="${ROOT}/92_votu_abundance/mapping"

mkdir -p "$OUT"

export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

# ----------------------------------------------------------------------
# 12 threads Bowtie2 + 4 threads samtools sort = 16 CPUs allocated.
# ----------------------------------------------------------------------

BT2_THREADS=12
SORT_THREADS=4

SAMPLES=(
    L1_1
    L1_2
    L1_3
    L2_1
    L2_2
    L2_3
    L3_1
    L3_2
    L3_3
    M1_1
    M1_2
    M1_3
    M2_1
    M2_2
    M2_3
    M3_1
    M3_2
    M3_3
)

SAMPLE="${SAMPLES[$SLURM_ARRAY_TASK_ID]}"

R1="${ROOT}/04_host_removed/${SAMPLE}/${SAMPLE}_1.hostfree.fastq.gz"
R2="${ROOT}/04_host_removed/${SAMPLE}/${SAMPLE}_2.hostfree.fastq.gz"

SOUT="${OUT}/${SAMPLE}"

mkdir -p "$SOUT"

BAM="${SOUT}/${SAMPLE}.votu.sorted.bam"
BAM_TMP="${SOUT}/${SAMPLE}.votu.sorted.tmp.bam"

BT2LOG="${SOUT}/${SAMPLE}.bowtie2.log"
FLAGSTAT="${SOUT}/${SAMPLE}.flagstat.txt"
IDXSTATS="${SOUT}/${SAMPLE}.idxstats.tsv"
METRICS="${SOUT}/${SAMPLE}.mapping_metrics.tsv"
META="${SOUT}/${SAMPLE}.run_metadata.tsv"

DONE="${SOUT}/${SAMPLE}.92G1_COMPLETE.ok"


echo "======================================================================"
echo "PASO 92G1 - MAPEO DE READS CONTRA 4491 vOTUs"
echo "Muestra: ${SAMPLE}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Job: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Inicio: $(date)"
echo "======================================================================"


# ======================================================================
# 1. VALIDACIONES
# ======================================================================

for F in \
    "$BOWTIE2_SIF" \
    "$SAMTOOLS_SIF" \
    "$REFERENCE" \
    "$R1" \
    "$R2"
do
    [[ -s "$F" ]] || {
        echo "ERROR: falta archivo:"
        echo "$F"
        exit 2
    }
done


[[ -x "$APPTAINER" ]] || {
    echo "ERROR: falta Apptainer." >&2
    exit 3
}


# Validar índice Bowtie2
for SUFFIX in \
    1 \
    2 \
    3 \
    4 \
    rev.1 \
    rev.2
do

    if [[ ! -s "${INDEX}.${SUFFIX}.bt2" \
       && ! -s "${INDEX}.${SUFFIX}.bt2l" ]]
    then

        echo "ERROR: índice incompleto: ${SUFFIX}" >&2
        exit 4

    fi
done


# ======================================================================
# 2. EVITAR RECALCULAR MUESTRAS YA COMPLETAS
# ======================================================================

if [[ -s "$DONE" \
   && -s "$BAM" \
   && -s "${BAM}.bai" \
   && -s "$METRICS" ]]
then

    echo
    echo "Muestra ${SAMPLE} ya completada."
    echo "No se recalculará."

else

    rm -f \
        "$DONE" \
        "$BAM_TMP" \
        "$BAM" \
        "${BAM}.bai" \
        "$BT2LOG" \
        "$FLAGSTAT" \
        "$IDXSTATS" \
        "$METRICS" \
        "$META"


    # ==================================================================
    # 3. MAPEO
    # ==================================================================

    echo
    echo "======================================================================"
    echo "1. BOWTIE2 --very-sensitive"
    echo "======================================================================"

    export \
        APPTAINER \
        ROOT \
        BOWTIE2_SIF \
        SAMTOOLS_SIF \
        INDEX \
        R1 \
        R2 \
        BAM_TMP \
        BT2LOG \
        SAMPLE \
        BT2_THREADS \
        SORT_THREADS


    srun --exact \
        -n 1 \
        -c "$SLURM_CPUS_PER_TASK" \
        bash -c '
            set -euo pipefail

            "$APPTAINER" exec \
                --bind "${ROOT}:${ROOT}" \
                "$BOWTIE2_SIF" \
                bowtie2 \
                    --very-sensitive \
                    --no-unal \
                    -p "$BT2_THREADS" \
                    --rg-id "$SAMPLE" \
                    --rg "SM:${SAMPLE}" \
                    -x "$INDEX" \
                    -1 "$R1" \
                    -2 "$R2" \
                    2> "$BT2LOG" \
            | \
            "$APPTAINER" exec \
                --bind "${ROOT}:${ROOT}" \
                "$SAMTOOLS_SIF" \
                samtools sort \
                    -@ "$SORT_THREADS" \
                    -m 2G \
                    -o "$BAM_TMP" \
                    -
        '


    [[ -s "$BAM_TMP" ]] || {
        echo "ERROR: BAM temporal vacío." >&2
        exit 5
    }


    # ==================================================================
    # 4. VALIDAR BAM
    # ==================================================================

    echo
    echo "======================================================================"
    echo "2. VALIDANDO BAM"
    echo "======================================================================"

    "$APPTAINER" exec \
        --bind "${ROOT}:${ROOT}" \
        "$SAMTOOLS_SIF" \
        samtools quickcheck \
            -v \
            "$BAM_TMP"


    mv \
        "$BAM_TMP" \
        "$BAM"


    # ==================================================================
    # 5. INDEX
    # ==================================================================

    srun --exact \
        -n 1 \
        -c 2 \
        "$APPTAINER" exec \
        --bind "${ROOT}:${ROOT}" \
        "$SAMTOOLS_SIF" \
        samtools index \
            -@ 2 \
            "$BAM"


    [[ -s "${BAM}.bai" ]] || {
        echo "ERROR: no se generó BAI." >&2
        exit 6
    }


    # ==================================================================
    # 6. FLAGSTAT + IDXSTATS
    # ==================================================================

    "$APPTAINER" exec \
        --bind "${ROOT}:${ROOT}" \
        "$SAMTOOLS_SIF" \
        samtools flagstat \
            "$BAM" \
        > "$FLAGSTAT"


    "$APPTAINER" exec \
        --bind "${ROOT}:${ROOT}" \
        "$SAMTOOLS_SIF" \
        samtools idxstats \
            "$BAM" \
        > "$IDXSTATS"


    # ==================================================================
    # 7. MÉTRICAS
    # ==================================================================

    PRIMARY_MAPPED="$(
        "$APPTAINER" exec \
            --bind "${ROOT}:${ROOT}" \
            "$SAMTOOLS_SIF" \
            samtools view \
                -c \
                -F 3844 \
                "$BAM"
    )"


    MAPQ20_MAPPED="$(
        "$APPTAINER" exec \
            --bind "${ROOT}:${ROOT}" \
            "$SAMTOOLS_SIF" \
            samtools view \
                -c \
                -F 3844 \
                -q 20 \
                "$BAM"
    )"


    VOTU_WITH_READS="$(
        awk \
            '$1!="*" && $3>0 {n++}
             END{print n+0}' \
            "$IDXSTATS"
    )"


    REF_ROWS="$(
        awk \
            '$1!="*" {n++}
             END{print n+0}' \
            "$IDXSTATS"
    )"


    [[ "$REF_ROWS" -eq 4491 ]] || {
        echo "ERROR: BAM tiene ${REF_ROWS} referencias; esperaba 4491." >&2
        exit 7
    }


    OVERALL_RATE="$(
        awk \
            '/overall alignment rate/{
                print $1
            }' \
            "$BT2LOG" \
        | tail -n 1
    )"


    if [[ -z "$OVERALL_RATE" ]]
    then
        OVERALL_RATE="NA"
    fi


    {
        echo -e "metric\tvalue"

        echo -e \
"sample\t${SAMPLE}"

        echo -e \
"primary_mapped_alignment_records\t${PRIMARY_MAPPED}"

        echo -e \
"mapq20_primary_mapped_alignment_records\t${MAPQ20_MAPPED}"

        echo -e \
"vOTUs_with_at_least_one_mapped_record\t${VOTU_WITH_READS}"

        echo -e \
"reference_vOTUs\t${REF_ROWS}"

        echo -e \
"bowtie2_overall_alignment_rate\t${OVERALL_RATE}"

        echo -e \
"bam_bytes\t$(stat -c %s "$BAM")"

    } > "$METRICS"


    # ==================================================================
    # 8. METADATA
    # ==================================================================

    {
        echo -e "field\tvalue"

        echo -e \
"analysis\tvOTU_read_mapping"

        echo -e \
"sample\t${SAMPLE}"

        echo -e \
"reference\t92_votu_abundance/reference/votu_4491.fna"

        echo -e \
"reference_vOTUs\t4491"

        echo -e \
"reference_sha256\t$(sha256sum "$REFERENCE" | awk '{print $1}')"

        echo -e \
"bowtie2_version\t2.5.4"

        echo -e \
"bowtie2_preset\tvery-sensitive"

        echo -e \
"alignment_mode\tend-to-end"

        echo -e \
"unaligned_records_in_BAM\tNO"

        echo -e \
"reporting_mode\tBowtie2_default_best_alignment"

        echo -e \
"mixed_alignment_policy\tBowtie2_default"

        echo -e \
"discordant_alignment_policy\tBowtie2_default"

        echo -e \
"explicit_read_identity_filter\tNO"

        echo -e \
"explicit_MAPQ_filter_during_mapping\tNO"

        echo -e \
"downstream_MAPQ20_sensitivity\tPLANNED"

        echo -e \
"downstream_identity_sensitivity\tPLANNED"

        echo -e \
"biological_presence_threshold\tNOT_DEFINED"

        echo -e \
"R1\t${R1}"

        echo -e \
"R2\t${R2}"

        echo -e \
"threads_bowtie2\t${BT2_THREADS}"

        echo -e \
"threads_samtools_sort\t${SORT_THREADS}"

        echo -e \
"date\t$(date -Iseconds)"

    } > "$META"


    # ==================================================================
    # 9. MARCADOR
    # ==================================================================

    printf \
        "COMPLETED\t%s\tjob=%s\tsample=%s\n" \
        "$(date -Iseconds)" \
        "$SLURM_JOB_ID" \
        "$SAMPLE" \
        > "$DONE"

fi


# ======================================================================
# 10. VALIDACIÓN FINAL
# ======================================================================

echo
echo "======================================================================"
echo "3. VALIDACIÓN FINAL ${SAMPLE}"
echo "======================================================================"

for F in \
    "$BAM" \
    "${BAM}.bai" \
    "$BT2LOG" \
    "$FLAGSTAT" \
    "$IDXSTATS" \
    "$METRICS" \
    "$META" \
    "$DONE"
do

    [[ -s "$F" ]] || {
        echo "ERROR: salida ausente/vacía:"
        echo "$F"
        exit 8
    }

    echo "OK: $(basename "$F")"

done


echo
echo "===== MÉTRICAS ====="

column -t -s $'\t' \
    "$METRICS"


echo
echo "===== BOWTIE2 FINAL ====="

tail -n 15 \
    "$BT2LOG"


echo
echo "======================================================================"
echo "PASO 92G1 ${SAMPLE} FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "======================================================================"

