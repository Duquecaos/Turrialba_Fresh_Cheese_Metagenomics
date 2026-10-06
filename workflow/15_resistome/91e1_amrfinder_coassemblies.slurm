#!/bin/bash

#SBATCH --job-name=91e1_amrcomm
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --array=0-5%2

#SBATCH --mail-type=END,FAIL

#SBATCH --output=91e1_amrcomm_%A_%a.out
#SBATCH --error=91e1_amrcomm_%A_%a.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"
OUTROOT="${ROOT}/79_resistome_coassemblies"

APPTAINER="${APPTAINER:-/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer}"
AMR_SIF="${AMR_SIF:-/opt/ohpc/pub/containers/BIO/ncbi-amrfinderplus-4.0.3.sif}"

DB_CONT="/work/databases/amrfinderplus/latest"

COASSEMBLIES=(
    L1
    L2
    L3
    M1
    M2
    M3
)

GROUP="${COASSEMBLIES[$SLURM_ARRAY_TASK_ID]}"

ASSEMBLY="${ROOT}/05_assembly/${GROUP}/final.contigs.fa"

OUT="${OUTROOT}/${GROUP}"

FAA="${OUT}/${GROUP}.prodigal.faa"
RAWGFF="${OUT}/${GROUP}.prodigal.raw.gff3"
STDGFF="${OUT}/${GROUP}.prodigal.amrfinder.gff3"

AMRF="${OUT}/${GROUP}.amrfinder_plus.tsv"
AMR_ONLY="${OUT}/${GROUP}.AMR.tsv"
STRESS_ONLY="${OUT}/${GROUP}.STRESS.tsv"
NON_AMR="${OUT}/${GROUP}.non_AMR.tsv"

mkdir -p "$OUT"

echo "======================================================================"
echo "PASO 91E1 - RESISTOMA COENSAMBLAJE"
echo "Grupo: ${GROUP}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Job: ${SLURM_JOB_ID}"
echo "Host: $(hostname)"
echo "Inicio: $(date)"
echo "======================================================================"

# ----------------------------------------------------------------------
# Validaciones
# ----------------------------------------------------------------------

[[ -x "$APPTAINER" ]] || {
    echo "ERROR: Apptainer no encontrado." >&2
    exit 2
}

[[ -s "$AMR_SIF" ]] || {
    echo "ERROR: contenedor AMRFinderPlus no encontrado." >&2
    exit 3
}

[[ -s "$ASSEMBLY" ]] || {
    echo "ERROR: falta assembly: $ASSEMBLY" >&2
    exit 4
}

[[ -d "${ROOT}/databases/amrfinderplus/2025-07-16.1" ]] || {
    echo "ERROR: falta database AMRFinderPlus." >&2
    exit 5
}

echo
echo "Assembly:"
ls -lh "$ASSEMBLY"

echo
echo "======================================================================"
echo "1. LOCALIZANDO PRODIGAL"
echo "======================================================================"

PRODIGAL_SIF="$(
    find "${CONTAINER_ROOT:-/opt/ohpc/pub/containers/BIO}" \
        -maxdepth 1 \
        -type f \
        -iname 'prodigal*.sif' \
        2>/dev/null \
    | sort -V \
    | tail -n 1
)"

USE_AMR_CONTAINER_FOR_PRODIGAL=0

if [[ -n "${PRODIGAL_SIF}" && -s "${PRODIGAL_SIF}" ]]; then

    echo "Prodigal container:"
    echo "$PRODIGAL_SIF"

else

    echo "No se encontró contenedor Prodigal independiente."
    echo "Probando dentro de AMRFinderPlus..."

    if srun --exact \
        -n 1 \
        -c 1 \
        "$APPTAINER" exec \
        "$AMR_SIF" \
        sh -lc 'command -v prodigal >/dev/null 2>&1'
    then
        PRODIGAL_SIF="$AMR_SIF"
        USE_AMR_CONTAINER_FOR_PRODIGAL=1
        echo "Prodigal disponible dentro de AMRFinderPlus."
    else
        echo "ERROR: no se encontró Prodigal disponible." >&2
        echo "Busque con:" >&2
        echo "find ${CONTAINER_ROOT:-/opt/ohpc/pub/containers/BIO} -iname '*prodigal*.sif'" >&2
        exit 6
    fi
fi

echo
echo "======================================================================"
echo "2. VERSIONES"
echo "======================================================================"

srun --exact \
    -n 1 \
    -c 1 \
    "$APPTAINER" exec \
    "$PRODIGAL_SIF" \
    prodigal -v \
    2>&1 \
    | head -n 10

"$APPTAINER" exec \
    --bind "${ROOT}:/work" \
    "$AMR_SIF" \
    amrfinder \
        -d "$DB_CONT" \
        -V

echo
echo "======================================================================"
echo "3. PRODIGAL - MODO METAGENOMICO"
echo "======================================================================"

srun --exact \
    -n 1 \
    -c 1 \
    "$APPTAINER" exec \
    --bind "${ROOT}:/work" \
    "$PRODIGAL_SIF" \
    prodigal \
        -i "/work/05_assembly/${GROUP}/final.contigs.fa" \
        -a "/work/79_resistome_coassemblies/${GROUP}/${GROUP}.prodigal.faa" \
        -o "/work/79_resistome_coassemblies/${GROUP}/${GROUP}.prodigal.raw.gff3" \
        -f gff \
        -p meta \
        -m \
        -q

[[ -s "$FAA" ]] || {
    echo "ERROR: Prodigal no produjo FAA." >&2
    exit 7
}

[[ -s "$RAWGFF" ]] || {
    echo "ERROR: Prodigal no produjo GFF." >&2
    exit 8
}

echo -n "Proteínas predichas: "
grep -c '^>' "$FAA"

echo
echo "======================================================================"
echo "4. GFF ESTANDAR PARA AMRFINDERPLUS"
echo "======================================================================"

GROUP="$GROUP" python3 - <<'PY'
import os
import re
from pathlib import Path

group = os.environ["GROUP"]

root = Path(
    f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"
)

assembly = (
    root
    / "05_assembly"
    / group
    / "final.contigs.fa"
)

faa = (
    root
    / "79_resistome_coassemblies"
    / group
    / f"{group}.prodigal.faa"
)

gff = (
    root
    / "79_resistome_coassemblies"
    / group
    / f"{group}.prodigal.amrfinder.gff3"
)

# ----------------------------------------------------------
# Contigs válidos
# ----------------------------------------------------------

contigs = set()

with assembly.open() as fh:
    for line in fh:
        if line.startswith(">"):
            contigs.add(
                line[1:].split()[0]
            )

if not contigs:
    raise RuntimeError(
        "Assembly sin contigs."
    )

# Prodigal FASTA:
# >protein_id # start # end # strand # metadata

pattern = re.compile(
    r"^>(\S+)\s+#\s+"
    r"(\d+)\s+#\s+"
    r"(\d+)\s+#\s+"
    r"(-?1)\s+#\s*(.*)$"
)

n = 0
bad = []
missing = []
seen = set()

with faa.open() as inp, gff.open("w") as out:

    out.write("##gff-version 3\n")

    for line in inp:

        if not line.startswith(">"):
            continue

        n += 1

        line = line.rstrip()

        m = pattern.match(line)

        if not m:
            bad.append(line)
            continue

        pid = m.group(1)
        start = int(m.group(2))
        end = int(m.group(3))
        strand_code = int(m.group(4))

        if pid in seen:
            raise RuntimeError(
                f"Protein ID duplicado: {pid}"
            )

        seen.add(pid)

        # Prodigal agrega _N al contig.
        if "_" not in pid:
            raise RuntimeError(
                f"No se puede derivar contig: {pid}"
            )

        contig = pid.rsplit("_", 1)[0]

        if contig not in contigs:
            missing.append(
                (pid, contig)
            )
            continue

        lo = min(start, end)
        hi = max(start, end)

        strand = (
            "+"
            if strand_code == 1
            else "-"
        )

        attrs = (
            f"ID={pid};"
            f"Name={pid}"
        )

        out.write(
            "\t".join([
                contig,
                "Prodigal",
                "CDS",
                str(lo),
                str(hi),
                ".",
                strand,
                "0",
                attrs,
            ])
            + "\n"
        )

print("group:", group)
print("assembly_contigs:", len(contigs))
print("proteins:", n)
print("bad_headers:", len(bad))
print("missing_contigs:", len(missing))

if bad:
    for x in bad[:10]:
        print("BAD_HEADER:", x)

if missing:
    for x in missing[:10]:
        print(
            "MISSING_CONTIG:",
            *x,
            sep="\t"
        )

if bad:
    raise RuntimeError(
        "Cabeceras Prodigal no interpretables."
    )

if missing:
    raise RuntimeError(
        "Proteínas con contigs ausentes del assembly."
    )

if n != len(seen):
    raise RuntimeError(
        "Proteínas != IDs únicos."
    )

print("GFF OK:", gff)
PY

[[ -s "$STDGFF" ]] || {
    echo "ERROR: GFF estándar vacío." >&2
    exit 9
}

echo -n "CDS GFF: "
awk '$3=="CDS"{n++} END{print n+0}' "$STDGFF"

echo
echo "======================================================================"
echo "5. AMRFINDERPLUS COMBINADO"
echo "======================================================================"

srun --exact \
    -n 1 \
    -c "${SLURM_CPUS_PER_TASK}" \
    "$APPTAINER" exec \
    --bind "${ROOT}:/work" \
    "$AMR_SIF" \
    amrfinder \
        --protein "/work/79_resistome_coassemblies/${GROUP}/${GROUP}.prodigal.faa" \
        --nucleotide "/work/05_assembly/${GROUP}/final.contigs.fa" \
        --gff "/work/79_resistome_coassemblies/${GROUP}/${GROUP}.prodigal.amrfinder.gff3" \
        --annotation_format standard \
        --plus \
        --database "$DB_CONT" \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --output "/work/79_resistome_coassemblies/${GROUP}/${GROUP}.amrfinder_plus.tsv"

[[ -s "$AMRF" ]] || {
    echo "ERROR: AMRFinderPlus no produjo salida." >&2
    exit 10
}

echo
echo "======================================================================"
echo "6. SEPARAR AMR / STRESS / OTROS"
echo "======================================================================"

GROUP="$GROUP" python3 - <<'PY'
import csv
import os
from pathlib import Path
from collections import Counter

group = os.environ["GROUP"]

root = Path(
    f"/scratch/global/{os.environ['USER']}/Shotgun_MAGs_Turrialba"
)

base = (
    root
    / "79_resistome_coassemblies"
    / group
)

src = base / f"{group}.amrfinder_plus.tsv"

with src.open(newline="") as fh:
    r = csv.DictReader(
        fh,
        delimiter="\t"
    )
    rows = list(r)
    fields = r.fieldnames

if not fields:
    raise RuntimeError(
        "AMRFinderPlus sin encabezado."
    )

def write_subset(path, subset):

    with path.open(
        "w",
        newline=""
    ) as out:

        w = csv.DictWriter(
            out,
            fieldnames=fields,
            delimiter="\t"
        )

        w.writeheader()
        w.writerows(subset)

amr = [
    x for x in rows
    if x.get("Type", "").strip() == "AMR"
]

stress = [
    x for x in rows
    if x.get("Type", "").strip() == "STRESS"
]

non_amr = [
    x for x in rows
    if x.get("Type", "").strip() != "AMR"
]

write_subset(
    base / f"{group}.AMR.tsv",
    amr
)

write_subset(
    base / f"{group}.STRESS.tsv",
    stress
)

write_subset(
    base / f"{group}.non_AMR.tsv",
    non_amr
)

print("group:", group)
print("all_records:", len(rows))
print("AMR:", len(amr))
print("STRESS:", len(stress))
print("non_AMR:", len(non_amr))

print()
print("AMR classes:")

for k, n in Counter(
    x.get("Class", "")
    for x in amr
).most_common():

    print(k, n, sep="\t")

print()
print("AMR methods:")

for k, n in Counter(
    x.get("Method", "")
    for x in amr
).most_common():

    print(k, n, sep="\t")
PY

echo
echo "======================================================================"
echo "7. METADATA"
echo "======================================================================"

{
    echo -e "field\tvalue"
    echo -e "coassembly\t${GROUP}"
    echo -e "analysis\tAMRFinderPlus_community_coassembly"
    echo -e "AMRFinderPlus_version\t4.0.3"
    echo -e "AMRFinderPlus_database\t2025-07-16.1"
    echo -e "plus_mode\ttrue"
    echo -e "organism_option\tNOT_USED"
    echo -e "Prodigal_mode\tmeta"
    echo -e "Prodigal_container\t${PRODIGAL_SIF}"
    echo -e "assembly\t05_assembly/${GROUP}/final.contigs.fa"
    echo -e "date\t$(date -Iseconds)"
} > "${OUT}/${GROUP}.91E1_metadata.tsv"

sha256sum \
    "$ASSEMBLY" \
    "$FAA" \
    "$STDGFF" \
    "$AMRF" \
    > "${OUT}/${GROUP}.91E1_sha256.tsv"

echo
echo "======================================================================"
echo "PASO 91E1 ${GROUP} FINALIZÓ CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "Fin: $(date)"
echo "======================================================================"
