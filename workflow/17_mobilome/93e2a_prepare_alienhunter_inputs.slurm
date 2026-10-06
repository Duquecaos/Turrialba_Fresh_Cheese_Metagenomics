#!/bin/bash

#SBATCH --job-name=93e2a_AH_prep
#SBATCH --partition=shared
#SBATCH --qos=sm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00

#SBATCH --mail-type=END,FAIL

#SBATCH --output=93e2a_AH_prep_%j.out
#SBATCH --error=93e2a_AH_prep_%j.err

set -euo pipefail

ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

MAGDIR="${ROOT}/22_final_representative_mags/fastas"
OUT="${ROOT}/95_mobilome/93E2_alienhunter_final18"

INPUTDIR="${OUT}/inputs"

rm -rf "$OUT"
mkdir -p "$INPUTDIR"

echo "======================================================================"
echo "93E2A - PREPARAR CONTIGS >=50 KB PARA ALIENHUNTER"
echo "Inicio: $(date)"
echo "Job: ${SLURM_JOB_ID}"
echo "======================================================================"

[[ -d "$MAGDIR" ]] || {
    echo "ERROR: falta $MAGDIR" >&2
    exit 2
}


python3 - \
    "$MAGDIR" \
    "$INPUTDIR" \
    "$OUT" <<'PY'

import csv
import hashlib
import sys
from pathlib import Path


MAGDIR = Path(sys.argv[1])
INPUTDIR = Path(sys.argv[2])
OUT = Path(sys.argv[3])

MINLEN = 50000


def fasta_records(path):

    name = None
    seq = []

    with path.open() as fh:

        for line in fh:

            if line.startswith(">"):

                if name is not None:
                    yield name, "".join(seq).upper()

                name = (
                    line[1:]
                    .strip()
                    .split()[0]
                )

                seq = []

            else:

                seq.append(
                    line.strip()
                )

        if name is not None:
            yield name, "".join(seq).upper()


fastas = sorted(
    p
    for p in MAGDIR.iterdir()
    if (
        p.is_file()
        and p.suffix.lower()
        in {".fa", ".fna", ".fasta", ".fas"}
    )
)

if len(fastas) != 18:
    raise RuntimeError(
        f"Esperaba 18 MAG FASTA; encontré {len(fastas)}"
    )


rows = []

seen_contigs = set()

for fasta in fastas:

    mag = fasta.stem

    if "__" not in mag:
        raise RuntimeError(
            f"No puedo inferir coassembly del MAG {mag}"
        )

    coassembly = mag.split("__", 1)[0]

    for contig, seq in fasta_records(fasta):

        if len(seq) < MINLEN:
            continue

        key = (
            coassembly,
            contig
        )

        if key in seen_contigs:
            raise RuntimeError(
                f"Contig duplicado dentro del conjunto final18: {key}"
            )

        seen_contigs.add(key)

        idx = len(rows)

        input_id = (
            f"AHCONTIG{idx:04d}"
        )

        outfile = (
            INPUTDIR
            / f"{input_id}.fna"
        )

        with outfile.open("w") as fh:

            fh.write(
                f">{contig}\n"
            )

            for i in range(
                0,
                len(seq),
                80
            ):
                fh.write(
                    seq[i:i+80]
                    + "\n"
                )

        sha = hashlib.sha256(
            seq.encode()
        ).hexdigest()

        rows.append({
            "array_index": idx,
            "input_id": input_id,
            "coassembly": coassembly,
            "MAG": mag,
            "contig": contig,
            "contig_length_bp": len(seq),
            "input_fasta": str(outfile),
            "sequence_sha256": sha,
        })


if len(rows) != 271:
    raise RuntimeError(
        f"Esperaba 271 contigs >=50 kb; observé {len(rows)}"
    )


manifest = (
    OUT
    / "93E2A_manifest.tsv"
)

fields = [
    "array_index",
    "input_id",
    "coassembly",
    "MAG",
    "contig",
    "contig_length_bp",
    "input_fasta",
    "sequence_sha256",
]

with manifest.open(
    "w",
    newline=""
) as fh:

    writer = csv.DictWriter(
        fh,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)


with (
    OUT
    / "93E2A_summary.tsv"
).open("w") as fh:

    fh.write(
        "metric\tvalue\n"
    )

    fh.write(
        "final18_MAGs\t18\n"
    )

    fh.write(
        f"contigs_ge_50kb\t{len(rows)}\n"
    )

    fh.write(
        f"total_bp_screened\t"
        f"{sum(r['contig_length_bp'] for r in rows)}\n"
    )

    fh.write(
        f"minimum_contig_bp\t"
        f"{min(r['contig_length_bp'] for r in rows)}\n"
    )

    fh.write(
        f"maximum_contig_bp\t"
        f"{max(r['contig_length_bp'] for r in rows)}\n"
    )

    fh.write(
        f"MAGs_represented\t"
        f"{len({r['MAG'] for r in rows})}\n"
    )


print(
    f"CONTIGS_GE50KB={len(rows)}"
)

print(
    f"MAGS_REPRESENTED={len({r['MAG'] for r in rows})}"
)

print(
    "93E2A_PREP=PASS"
)

PY


# ======================================================================
# Validaciones
# ======================================================================

MANIFEST="${OUT}/93E2A_manifest.tsv"

[[ -s "$MANIFEST" ]] || {
    echo "ERROR: manifest vacío" >&2
    exit 3
}


N_MANIFEST=$(
    tail -n +2 "$MANIFEST" |
    wc -l
)

N_FASTA=$(
    find "$INPUTDIR" \
        -maxdepth 1 \
        -type f \
        -name 'AHCONTIG*.fna' \
        | wc -l
)

if [[ "$N_MANIFEST" -ne 271 ]]
then
    echo "ERROR: manifest tiene $N_MANIFEST entradas" >&2
    exit 4
fi

if [[ "$N_FASTA" -ne 271 ]]
then
    echo "ERROR: hay $N_FASTA FASTA individuales" >&2
    exit 5
fi


sha256sum \
    "$MANIFEST" \
    > "${OUT}/93E2A_manifest_sha256.tsv"


{
    echo -e "field\tvalue"
    echo -e "scope\tfinal18_MAG_contigs"
    echo -e "minimum_contig_length_bp\t50000"
    echo -e "expected_contigs\t271"
    echo -e "AlienHunter_prediction_term\tcomposition_based_candidate_GI_region"
    echo -e "negative_prediction_equals_GI_absence\tNO"
    echo -e "fragmented_MAG_negative_result_equals_absence\tNO"
} > "${OUT}/93E2A_metadata.tsv"


printf \
"COMPLETED\t%s\tjob=%s\n" \
"$(date -Iseconds)" \
"${SLURM_JOB_ID}" \
> "${OUT}/93E2A_COMPLETE.ok"


echo
echo "======================================================================"
echo "93E2A SUMMARY"
echo "======================================================================"

column -t -s $'\t' \
    "${OUT}/93E2A_summary.tsv"

echo
echo "Manifest rows: $N_MANIFEST"
echo "Input FASTA:   $N_FASTA"

echo
echo "======================================================================"
echo "93E2A FINALIZO CORRECTAMENTE."
echo "ES SEGURO SALIR."
echo "======================================================================"

