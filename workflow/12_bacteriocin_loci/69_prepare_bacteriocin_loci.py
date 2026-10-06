#!/usr/bin/env python3

# ============================================================
# PASO 69
#
# Reconstrucción reproducible de candidatos bacteriocínicos:
#
# 1. Localiza exactamente las proteínas canónicas en sus MAGs.
# 2. Recupera segmento nucleotídico codificante.
# 3. Busca cada proteína exacta en los 18 MAGs representativos.
# 4. Extrae contextos +/- 5 kb.
# 5. Fusiona contextos solapados del mismo MAG/contig.
# 6. Detecta secuencias nucleotídicas idénticas entre MAGs.
#
# IMPORTANTE:
# - El candidato lcnB-like interrumpido se conserva como
#   evidencia de locus, NO como precursor funcional.
# - Los contextos se mantienen en orientación genómica.
# - Los segmentos codificantes se escriben en orientación
#   5'->3' del ORF.
# ============================================================

from pathlib import Path
from collections import defaultdict
import csv
import hashlib
import os
import re
import sys


# ============================================================
# RUTAS
# ============================================================

USER = os.environ["USER"]

ROOT = Path(
    f"/scratch/global/{USER}/Shotgun_MAGs_Turrialba"
)

MASTER = (
    ROOT
    / "38_bacteriocin_master"
)

METADATA = (
    MASTER
    / "structural_candidates_metadata.tsv"
)

PROTEINS = (
    MASTER
    / "structural_candidates_protein.faa"
)

MAGDIR = (
    ROOT
    / "22_final_representative_mags"
    / "fastas"
)

OUT = (
    ROOT
    / "47_bacteriocin_locus_reference"
)

OUT.mkdir(
    parents=True,
    exist_ok=True
)

FLANK = 5000


# ============================================================
# CÓDIGO GENÉTICO
# ============================================================

CODON_TABLE = {
    # T
    "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*",
    "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",

    # C
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",

    # A
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",

    # G
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
}

STOP_CODONS = {
    "TAA",
    "TAG",
    "TGA",
}

# Código bacteriano 11 puede usar estos codones como inicio.
# En una proteína anotada se representan como M.
BACTERIAL_START_CODONS = {
    "ATG",
    "GTG",
    "TTG",
    "CTG",
    "ATT",
    "ATC",
    "ATA",
}


# ============================================================
# HELPERS
# ============================================================

def revcomp(seq):
    table = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan"
    )

    return seq.translate(table)[::-1]


def translate(seq):
    seq = seq.upper()

    aa = []

    for i in range(
        0,
        len(seq) - 2,
        3
    ):
        codon = seq[i:i+3]

        aa.append(
            CODON_TABLE.get(
                codon,
                "X"
            )
        )

    return "".join(aa)


def read_tsv(path):
    with open(
        path,
        newline=""
    ) as handle:

        return list(
            csv.DictReader(
                handle,
                delimiter="\t"
            )
        )


def write_tsv(
    path,
    rows,
    fields
):
    with open(
        path,
        "w",
        newline=""
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path):
    records = {}

    name = None
    parts = []

    with open(path) as handle:

        for raw in handle:

            line = raw.strip()

            if not line:
                continue

            if line.startswith(">"):

                if name is not None:

                    if name in records:
                        raise RuntimeError(
                            f"FASTA ID duplicado: {name}"
                        )

                    records[name] = (
                        "".join(parts)
                        .upper()
                    )

                name = (
                    line[1:]
                    .split()[0]
                )

                parts = []

            else:
                parts.append(line)

    if name is not None:

        if name in records:
            raise RuntimeError(
                f"FASTA ID duplicado: {name}"
            )

        records[name] = (
            "".join(parts)
            .upper()
        )

    return records


def read_candidate_proteins(path):
    raw = read_fasta(path)

    result = {}

    full_headers = {}

    for header, seq in raw.items():

        candidate = header.split("|")[0]

        if candidate in result:
            raise RuntimeError(
                f"Candidato duplicado en FAA: "
                f"{candidate}"
            )

        result[candidate] = seq

        full_headers[candidate] = header

    return result, full_headers


def write_fasta(
    path,
    records,
    width=80
):
    with open(path, "w") as out:

        for header, seq in records:

            out.write(
                f">{header}\n"
            )

            for i in range(
                0,
                len(seq),
                width
            ):
                out.write(
                    seq[i:i+width]
                    + "\n"
                )


def sha256_text(seq):
    return hashlib.sha256(
        seq.encode()
    ).hexdigest()


def find_all(text, pattern):
    start = 0

    while True:

        pos = text.find(
            pattern,
            start
        )

        if pos < 0:
            return

        yield pos

        start = pos + 1


# ============================================================
# BUSCAR UNA PROTEÍNA EN LAS 6 FASES
#
# Acepta:
# 1) coincidencia exacta de traducción;
# 2) primer aa M codificado por start bacteriano alternativo.
# ============================================================

def find_peptide_matches(
    dna,
    peptide
):
    dna = dna.upper()

    n = len(dna)

    found = {}

    for strand, oriented in [
        ("+", dna),
        ("-", revcomp(dna))
    ]:

        for frame in range(3):

            usable = (
                len(oriented) - frame
            )

            usable -= (
                usable % 3
            )

            if usable <= 0:
                continue

            coding = oriented[
                frame:
                frame + usable
            ]

            aa = translate(coding)

            # ----------------------------------------
            # Coincidencia exacta
            # ----------------------------------------

            for aa_pos in find_all(
                aa,
                peptide
            ):

                coding_start = (
                    frame
                    + aa_pos * 3
                )

                coding_end = (
                    coding_start
                    + len(peptide) * 3
                )

                coding_nt = oriented[
                    coding_start:
                    coding_end
                ]

                if strand == "+":

                    orig_start0 = (
                        coding_start
                    )

                    orig_end0 = (
                        coding_end
                    )

                else:

                    orig_start0 = (
                        n - coding_end
                    )

                    orig_end0 = (
                        n - coding_start
                    )

                terminal_stop = ""

                if (
                    coding_end + 3
                    <= len(oriented)
                ):

                    next_codon = oriented[
                        coding_end:
                        coding_end + 3
                    ]

                    if (
                        next_codon
                        in STOP_CODONS
                    ):
                        terminal_stop = (
                            next_codon
                        )

                key = (
                    orig_start0,
                    orig_end0,
                    strand
                )

                found[key] = {
                    "start0":
                        orig_start0,

                    "end0":
                        orig_end0,

                    "strand":
                        strand,

                    "coding_nt":
                        coding_nt,

                    "match_mode":
                        "exact_translation",

                    "start_codon":
                        coding_nt[:3],

                    "terminal_stop_codon":
                        terminal_stop,
                }


            # ----------------------------------------
            # Start bacteriano alternativo
            #
            # Solo cuando la proteína anotada empieza M.
            # ----------------------------------------

            if (
                peptide.startswith("M")
                and len(peptide) > 1
            ):

                suffix = peptide[1:]

                for suffix_pos in find_all(
                    aa,
                    suffix
                ):

                    aa_pos = (
                        suffix_pos - 1
                    )

                    if aa_pos < 0:
                        continue

                    coding_start = (
                        frame
                        + aa_pos * 3
                    )

                    coding_end = (
                        coding_start
                        + len(peptide) * 3
                    )

                    if (
                        coding_end
                        > len(oriented)
                    ):
                        continue

                    start_codon = oriented[
                        coding_start:
                        coding_start + 3
                    ]

                    if (
                        start_codon
                        not in
                        BACTERIAL_START_CODONS
                    ):
                        continue

                    coding_nt = oriented[
                        coding_start:
                        coding_end
                    ]

                    # Confirmar el resto exacto.
                    translated = translate(
                        coding_nt
                    )

                    if len(
                        translated
                    ) != len(
                        peptide
                    ):
                        continue

                    if (
                        translated[1:]
                        != peptide[1:]
                    ):
                        continue

                    if strand == "+":

                        orig_start0 = (
                            coding_start
                        )

                        orig_end0 = (
                            coding_end
                        )

                    else:

                        orig_start0 = (
                            n - coding_end
                        )

                        orig_end0 = (
                            n - coding_start
                        )

                    terminal_stop = ""

                    if (
                        coding_end + 3
                        <= len(oriented)
                    ):

                        next_codon = oriented[
                            coding_end:
                            coding_end + 3
                        ]

                        if (
                            next_codon
                            in STOP_CODONS
                        ):
                            terminal_stop = (
                                next_codon
                            )

                    key = (
                        orig_start0,
                        orig_end0,
                        strand
                    )

                    if key not in found:

                        found[key] = {
                            "start0":
                                orig_start0,

                            "end0":
                                orig_end0,

                            "strand":
                                strand,

                            "coding_nt":
                                coding_nt,

                            "match_mode":
                                "bacterial_start_override",

                            "start_codon":
                                start_codon,

                            "terminal_stop_codon":
                                terminal_stop,
                        }

    return list(
        found.values()
    )


# ============================================================
# MERGE DE CONTEXTOS SOLAPADOS
# ============================================================

def merge_contexts(
    interval_rows,
    genome_sequences,
    prefix
):
    grouped = defaultdict(list)

    for row in interval_rows:

        grouped[
            (
                row["MAG"],
                row["contig"]
            )
        ].append(row)

    merged = []

    for (
        mag,
        contig
    ), rows in sorted(
        grouped.items()
    ):

        rows = sorted(
            rows,
            key=lambda x: (
                x["context_start0"],
                x["context_end0"]
            )
        )

        current = None

        for row in rows:

            if current is None:

                current = {
                    "MAG":
                        mag,

                    "contig":
                        contig,

                    "context_start0":
                        row[
                            "context_start0"
                        ],

                    "context_end0":
                        row[
                            "context_end0"
                        ],

                    "candidate_set":
                        {
                            row["candidate"]
                        },

                    "hit_rows":
                        [row],
                }

                continue

            if (
                row["context_start0"]
                <= current[
                    "context_end0"
                ]
            ):

                current[
                    "context_end0"
                ] = max(
                    current[
                        "context_end0"
                    ],
                    row[
                        "context_end0"
                    ]
                )

                current[
                    "candidate_set"
                ].add(
                    row["candidate"]
                )

                current[
                    "hit_rows"
                ].append(row)

            else:

                merged.append(current)

                current = {
                    "MAG":
                        mag,

                    "contig":
                        contig,

                    "context_start0":
                        row[
                            "context_start0"
                        ],

                    "context_end0":
                        row[
                            "context_end0"
                        ],

                    "candidate_set":
                        {
                            row["candidate"]
                        },

                    "hit_rows":
                        [row],
                }

        if current is not None:
            merged.append(current)

    # ----------------------------------------
    # Asignar IDs y secuencias
    # ----------------------------------------

    result = []

    fasta = []

    for idx, locus in enumerate(
        merged,
        start=1
    ):

        locus_id = (
            f"{prefix}{idx:03d}"
        )

        mag = locus["MAG"]
        contig = locus["contig"]

        seq = (
            genome_sequences[
                mag
            ][
                contig
            ]
        )

        s0 = locus[
            "context_start0"
        ]

        e0 = locus[
            "context_end0"
        ]

        context_seq = seq[s0:e0]

        candidates = sorted(
            locus[
                "candidate_set"
            ]
        )

        hit_rows = locus[
            "hit_rows"
        ]

        has_interrupted = any(
            (
                h.get(
                    "evidence_level",
                    ""
                )
                == "interrupted"
            )
            or (
                "premature_stop"
                in h.get(
                    "status",
                    ""
                )
            )
            for h in hit_rows
        )

        result.append({
            "locus_id":
                locus_id,

            "MAG":
                mag,

            "contig":
                contig,

            "context_start":
                s0 + 1,

            "context_end":
                e0,

            "context_length":
                len(context_seq),

            "contig_length":
                len(seq),

            "left_edge_truncated":
                int(s0 == 0),

            "right_edge_truncated":
                int(
                    e0 == len(seq)
                ),

            "n_candidate_hits":
                len(hit_rows),

            "n_unique_candidates":
                len(candidates),

            "candidates":
                ";".join(candidates),

            "contains_interrupted_candidate":
                int(has_interrupted),

            "context_sha256":
                sha256_text(
                    context_seq
                ),
        })

        fasta.append(
            (
                locus_id,
                context_seq
            )
        )

    return result, fasta


# ============================================================
# CARGAR METADATA Y PROTEÍNAS
# ============================================================

metadata = read_tsv(
    METADATA
)

proteins, protein_headers = (
    read_candidate_proteins(
        PROTEINS
    )
)


if len(metadata) != 8:
    raise RuntimeError(
        f"Esperábamos 8 candidatos en metadata; "
        f"hay {len(metadata)}"
    )

if len(proteins) != 8:
    raise RuntimeError(
        f"Esperábamos 8 proteínas; "
        f"hay {len(proteins)}"
    )


meta_by_candidate = {
    row["candidate"]: row
    for row in metadata
}

if (
    set(meta_by_candidate)
    != set(proteins)
):
    raise RuntimeError(
        "Los candidatos del metadata y FAA "
        "no coinciden exactamente."
    )


candidate_order = [
    row["candidate"]
    for row in metadata
]


# ============================================================
# EXTRAER CONTIG FUENTE DEL ID DEL CANDIDATO
# ============================================================

for candidate in candidate_order:

    match = re.search(
        r"(k141_\d+)",
        candidate
    )

    if not match:
        raise RuntimeError(
            f"No puedo inferir contig de "
            f"{candidate}"
        )

    meta_by_candidate[
        candidate
    ][
        "source_contig"
    ] = match.group(1)


# ============================================================
# CARGAR 18 MAGs
# ============================================================

mag_files = sorted(
    list(
        MAGDIR.glob("*.fa")
    )
    + list(
        MAGDIR.glob("*.fna")
    )
    + list(
        MAGDIR.glob("*.fasta")
    )
)


if len(mag_files) != 18:
    raise RuntimeError(
        f"Esperábamos 18 MAG FASTA; "
        f"hay {len(mag_files)}"
    )


genomes = {}

for path in mag_files:

    mag = path.name

    for suffix in [
        ".fasta",
        ".fna",
        ".fa"
    ]:
        if mag.endswith(suffix):
            mag = mag[
                :-len(suffix)
            ]
            break

    genomes[mag] = read_fasta(
        path
    )


if len(genomes) != 18:
    raise RuntimeError(
        "Los MAG IDs no son 18 únicos."
    )


for candidate in candidate_order:

    mag = meta_by_candidate[
        candidate
    ][
        "MAG"
    ]

    contig = meta_by_candidate[
        candidate
    ][
        "source_contig"
    ]

    if mag not in genomes:
        raise RuntimeError(
            f"MAG fuente no encontrado: {mag}"
        )

    if (
        contig
        not in genomes[mag]
    ):
        raise RuntimeError(
            f"Contig fuente no encontrado: "
            f"{mag}::{contig}"
        )


# ============================================================
# BÚSQUEDA DE LOS 8 CANDIDATOS EN LOS 18 MAGs
# ============================================================

all_hits = []


for mag in sorted(genomes):

    print(
        f"Buscando candidatos en {mag}...",
        flush=True
    )

    for contig, dna in genomes[
        mag
    ].items():

        for candidate in candidate_order:

            peptide = proteins[
                candidate
            ]

            hits = find_peptide_matches(
                dna,
                peptide
            )

            for hit in hits:

                meta = meta_by_candidate[
                    candidate
                ]

                start0 = hit[
                    "start0"
                ]

                end0 = hit[
                    "end0"
                ]

                context_start0 = max(
                    0,
                    start0 - FLANK
                )

                context_end0 = min(
                    len(dna),
                    end0 + FLANK
                )

                all_hits.append({
                    "candidate":
                        candidate,

                    "query_source_MAG":
                        meta["MAG"],

                    "query_source_contig":
                        meta[
                            "source_contig"
                        ],

                    "taxon":
                        meta["taxon"],

                    "evidence_level":
                        meta[
                            "evidence_level"
                        ],

                    "status":
                        meta["status"],

                    "matched_MAG":
                        mag,

                    "MAG":
                        mag,

                    "contig":
                        contig,

                    "start0":
                        start0,

                    "end0":
                        end0,

                    "start":
                        start0 + 1,

                    "end":
                        end0,

                    "strand":
                        hit["strand"],

                    "peptide_length":
                        len(peptide),

                    "coding_nt_length":
                        len(
                            hit[
                                "coding_nt"
                            ]
                        ),

                    "match_mode":
                        hit[
                            "match_mode"
                        ],

                    "start_codon":
                        hit[
                            "start_codon"
                        ],

                    "terminal_stop_codon":
                        hit[
                            "terminal_stop_codon"
                        ],

                    "contig_length":
                        len(dna),

                    "context_start0":
                        context_start0,

                    "context_end0":
                        context_end0,

                    "context_start":
                        context_start0 + 1,

                    "context_end":
                        context_end0,

                    "context_length":
                        (
                            context_end0
                            - context_start0
                        ),

                    "coding_nt":
                        hit[
                            "coding_nt"
                        ],

                    "coding_nt_sha256":
                        sha256_text(
                            hit[
                                "coding_nt"
                            ]
                        ),
                })


# ============================================================
# VALIDAR EXACTAMENTE UNA LOCALIZACIÓN FUENTE POR CANDIDATO
# ============================================================

source_hits = []

for candidate in candidate_order:

    meta = meta_by_candidate[
        candidate
    ]

    intended = [
        h
        for h in all_hits
        if (
            h["candidate"]
            == candidate
            and h["matched_MAG"]
            == meta["MAG"]
            and h["contig"]
            == meta[
                "source_contig"
            ]
        )
    ]

    if len(intended) != 1:

        raise RuntimeError(
            f"{candidate}: esperaba 1 hit "
            f"en {meta['MAG']}::"
            f"{meta['source_contig']}; "
            f"encontré {len(intended)}"
        )

    hit = dict(
        intended[0]
    )

    hit[
        "is_source_location"
    ] = 1

    source_hits.append(hit)


source_keys = {
    (
        h["candidate"],
        h["matched_MAG"],
        h["contig"],
        h["start"],
        h["end"],
        h["strand"]
    )
    for h in source_hits
}


for hit in all_hits:

    key = (
        hit["candidate"],
        hit["matched_MAG"],
        hit["contig"],
        hit["start"],
        hit["end"],
        hit["strand"]
    )

    hit[
        "is_source_location"
    ] = int(
        key in source_keys
    )


# ============================================================
# TABLA: LOCALIZACIONES FUENTE
# ============================================================

source_location_rows = []

source_fasta = []

for hit in source_hits:

    source_location_rows.append({
        "candidate":
            hit["candidate"],

        "MAG":
            hit["matched_MAG"],

        "taxon":
            hit["taxon"],

        "evidence_level":
            hit["evidence_level"],

        "status":
            hit["status"],

        "contig":
            hit["contig"],

        "start":
            hit["start"],

        "end":
            hit["end"],

        "strand":
            hit["strand"],

        "peptide_length":
            hit["peptide_length"],

        "coding_nt_length":
            hit["coding_nt_length"],

        "match_mode":
            hit["match_mode"],

        "start_codon":
            hit["start_codon"],

        "terminal_stop_codon":
            hit[
                "terminal_stop_codon"
            ],

        "contig_length":
            hit["contig_length"],

        "coding_nt_sha256":
            hit[
                "coding_nt_sha256"
            ],
    })

    source_fasta.append(
        (
            hit["candidate"],
            hit["coding_nt"]
        )
    )


write_tsv(
    OUT
    / "source_candidate_locations.tsv",

    source_location_rows,

    [
        "candidate",
        "MAG",
        "taxon",
        "evidence_level",
        "status",
        "contig",
        "start",
        "end",
        "strand",
        "peptide_length",
        "coding_nt_length",
        "match_mode",
        "start_codon",
        "terminal_stop_codon",
        "contig_length",
        "coding_nt_sha256",
    ]
)


write_fasta(
    OUT
    / "source_candidate_coding_segments.fna",

    source_fasta
)


# ============================================================
# TABLA: TODOS LOS MATCHES EXACTOS ENTRE 18 MAGs
# ============================================================

all_hit_rows = []

all_hit_segment_fasta = []


for idx, hit in enumerate(
    sorted(
        all_hits,
        key=lambda h: (
            h["candidate"],
            h["matched_MAG"],
            h["contig"],
            h["start"],
            h["strand"]
        )
    ),
    start=1
):

    match_id = (
        f"MATCH{idx:04d}"
    )

    row = {
        "match_id":
            match_id,

        "candidate":
            hit["candidate"],

        "query_source_MAG":
            hit[
                "query_source_MAG"
            ],

        "query_source_contig":
            hit[
                "query_source_contig"
            ],

        "matched_MAG":
            hit["matched_MAG"],

        "contig":
            hit["contig"],

        "start":
            hit["start"],

        "end":
            hit["end"],

        "strand":
            hit["strand"],

        "peptide_length":
            hit[
                "peptide_length"
            ],

        "coding_nt_length":
            hit[
                "coding_nt_length"
            ],

        "match_mode":
            hit["match_mode"],

        "start_codon":
            hit["start_codon"],

        "terminal_stop_codon":
            hit[
                "terminal_stop_codon"
            ],

        "coding_nt_sha256":
            hit[
                "coding_nt_sha256"
            ],

        "context_start":
            hit[
                "context_start"
            ],

        "context_end":
            hit[
                "context_end"
            ],

        "context_length":
            hit[
                "context_length"
            ],

        "is_source_location":
            hit[
                "is_source_location"
            ],

        "evidence_level":
            hit[
                "evidence_level"
            ],

        "status":
            hit["status"],
    }

    all_hit_rows.append(row)

    all_hit_segment_fasta.append(
        (
            match_id,
            hit["coding_nt"]
        )
    )


write_tsv(
    OUT
    / "all_exact_protein_matches.tsv",

    all_hit_rows,

    [
        "match_id",
        "candidate",
        "query_source_MAG",
        "query_source_contig",
        "matched_MAG",
        "contig",
        "start",
        "end",
        "strand",
        "peptide_length",
        "coding_nt_length",
        "match_mode",
        "start_codon",
        "terminal_stop_codon",
        "coding_nt_sha256",
        "context_start",
        "context_end",
        "context_length",
        "is_source_location",
        "evidence_level",
        "status",
    ]
)


write_fasta(
    OUT
    / "all_exact_protein_match_segments.fna",

    all_hit_segment_fasta
)


# ============================================================
# RESUMEN CROSS-MAG POR CANDIDATO
# ============================================================

cross_summary = []

for candidate in candidate_order:

    meta = meta_by_candidate[
        candidate
    ]

    hits = [
        h
        for h in all_hits
        if h["candidate"]
        == candidate
    ]

    matched_mags = sorted({
        h["matched_MAG"]
        for h in hits
    })

    non_source_mags = sorted({
        h["matched_MAG"]
        for h in hits
        if h["matched_MAG"]
        != meta["MAG"]
    })

    cross_summary.append({
        "candidate":
            candidate,

        "source_MAG":
            meta["MAG"],

        "source_contig":
            meta["source_contig"],

        "evidence_level":
            meta["evidence_level"],

        "status":
            meta["status"],

        "exact_match_count":
            len(hits),

        "unique_MAG_count":
            len(matched_mags),

        "matched_MAGs":
            ";".join(
                matched_mags
            ),

        "non_source_MAG_count":
            len(
                non_source_mags
            ),

        "non_source_MAGs":
            ";".join(
                non_source_mags
            ),

        "has_non_source_exact_match":
            int(
                len(
                    non_source_mags
                ) > 0
            ),
    })


write_tsv(
    OUT
    / "candidate_crossMAG_summary.tsv",

    cross_summary,

    [
        "candidate",
        "source_MAG",
        "source_contig",
        "evidence_level",
        "status",
        "exact_match_count",
        "unique_MAG_count",
        "matched_MAGs",
        "non_source_MAG_count",
        "non_source_MAGs",
        "has_non_source_exact_match",
    ]
)


# ============================================================
# GRUPOS DE SEGMENTOS NUCLEOTÍDICOS IDÉNTICOS
# ============================================================

nt_groups = defaultdict(list)

for row in all_hit_rows:

    nt_groups[
        row[
            "coding_nt_sha256"
        ]
    ].append(row)


duplicate_rows = []

group_number = 0

for sha, rows in sorted(
    nt_groups.items()
):

    unique_locations = {
        (
            r["matched_MAG"],
            r["contig"],
            r["start"],
            r["end"],
            r["strand"]
        )
        for r in rows
    }

    if len(
        unique_locations
    ) < 2:
        continue

    group_number += 1

    duplicate_rows.append({
        "duplicate_group":
            f"NTDUP{group_number:03d}",

        "sha256":
            sha,

        "nt_length":
            rows[0][
                "coding_nt_length"
            ],

        "n_hits":
            len(rows),

        "n_unique_locations":
            len(
                unique_locations
            ),

        "candidates":
            ";".join(
                sorted({
                    r["candidate"]
                    for r in rows
                })
            ),

        "MAGs":
            ";".join(
                sorted({
                    r["matched_MAG"]
                    for r in rows
                })
            ),

        "locations":
            ";".join(
                sorted(
                    f"{mag}:{contig}:"
                    f"{start}-{end}:"
                    f"{strand}"
                    for (
                        mag,
                        contig,
                        start,
                        end,
                        strand
                    )
                    in unique_locations
                )
            ),
    })


write_tsv(
    OUT
    / "exact_nucleotide_duplicate_groups.tsv",

    duplicate_rows,

    [
        "duplicate_group",
        "sha256",
        "nt_length",
        "n_hits",
        "n_unique_locations",
        "candidates",
        "MAGs",
        "locations",
    ]
)


# ============================================================
# CONTEXTOS FUENTE +/- 5 kb
#
# IMPORTANTE:
# Fusionar contextos que se solapan en el mismo MAG/contig.
# ============================================================

source_interval_rows = []

for hit in source_hits:

    row = dict(hit)

    row[
        "context_start0"
    ] = max(
        0,
        hit["start0"] - FLANK
    )

    row[
        "context_end0"
    ] = min(
        hit["contig_length"],
        hit["end0"] + FLANK
    )

    source_interval_rows.append(row)


source_loci, source_loci_fasta = (
    merge_contexts(
        source_interval_rows,
        genomes,
        "SRCLOC"
    )
)


write_tsv(
    OUT
    / "source_locus_metadata.tsv",

    source_loci,

    [
        "locus_id",
        "MAG",
        "contig",
        "context_start",
        "context_end",
        "context_length",
        "contig_length",
        "left_edge_truncated",
        "right_edge_truncated",
        "n_candidate_hits",
        "n_unique_candidates",
        "candidates",
        "contains_interrupted_candidate",
        "context_sha256",
    ]
)


write_fasta(
    OUT
    / "source_locus_contexts_5kb.fna",

    source_loci_fasta
)


# ============================================================
# TODOS LOS LOCI DONDE HAY PROTEÍNAS EXACTAMENTE IGUALES
#
# Esto crea la referencia de atribución competitiva.
# ============================================================

all_match_loci, all_match_loci_fasta = (
    merge_contexts(
        all_hits,
        genomes,
        "ATTRLOC"
    )
)


write_tsv(
    OUT
    / "all_match_locus_metadata.tsv",

    all_match_loci,

    [
        "locus_id",
        "MAG",
        "contig",
        "context_start",
        "context_end",
        "context_length",
        "contig_length",
        "left_edge_truncated",
        "right_edge_truncated",
        "n_candidate_hits",
        "n_unique_candidates",
        "candidates",
        "contains_interrupted_candidate",
        "context_sha256",
    ]
)


write_fasta(
    OUT
    / "all_match_locus_contexts_5kb.fna",

    all_match_loci_fasta
)


# ============================================================
# README
# ============================================================

README = OUT / "README_locus_reference.txt"

with open(
    README,
    "w"
) as out:

    out.write(
        "Bacteriocin/RiPP candidate locus reference\n"
    )

    out.write(
        "=========================================\n\n"
    )

    out.write(
        f"Candidates: {len(candidate_order)}\n"
    )

    out.write(
        f"Representative MAGs searched: {len(genomes)}\n"
    )

    out.write(
        f"Context flank: +/- {FLANK} bp\n\n"
    )

    out.write(
        "source_candidate_coding_segments.fna:\n"
        "  Exact nucleotide segment encoding each canonical "
        "candidate peptide in coding orientation.\n\n"
    )

    out.write(
        "source_locus_contexts_5kb.fna:\n"
        "  Genomic contexts around source candidates. "
        "Overlapping +/-5 kb contexts on the same contig "
        "are merged into one locus to prevent artificial "
        "multimapping between candidates in the same locus.\n\n"
    )

    out.write(
        "all_exact_protein_matches.tsv:\n"
        "  Exact occurrences of each canonical peptide across "
        "all 18 representative MAGs.\n\n"
    )

    out.write(
        "all_match_locus_contexts_5kb.fna:\n"
        "  Contexts surrounding all exact peptide matches "
        "across representative MAGs. This is intended as an "
        "attribution/ambiguity-control reference.\n\n"
    )

    out.write(
        "exact_nucleotide_duplicate_groups.tsv:\n"
        "  Coding segments that are nucleotide-identical at "
        "more than one genomic location.\n\n"
    )

    out.write(
        "Coordinates in TSV files are 1-based inclusive. "
        "Internal start0/end0 coordinates are not exported.\n\n"
    )

    out.write(
        "Coding FASTA sequences exclude terminal stop codons. "
        "terminal_stop_codon is reported separately when the "
        "immediately downstream codon in coding orientation "
        "is TAA/TAG/TGA.\n\n"
    )

    out.write(
        "Alternative bacterial initiation codons are allowed "
        "when reproducing an annotated N-terminal methionine; "
        "such matches are marked bacterial_start_override.\n\n"
    )

    out.write(
        "IMPORTANT: the premature-stop lcnB-like candidate is "
        "retained for locus tracking only and must NOT be "
        "interpreted as an intact functional precursor.\n\n"
    )

    out.write(
        "Exact peptide or nucleotide identity across MAGs does "
        "not establish horizontal transfer or biological "
        "activity.\n"
    )


# ============================================================
# FINAL QC
# ============================================================

if len(source_location_rows) != 8:
    raise RuntimeError(
        "No se obtuvieron las 8 "
        "localizaciones fuente."
    )

if len(cross_summary) != 8:
    raise RuntimeError(
        "No se generaron 8 filas "
        "de resumen cross-MAG."
    )

if len(source_loci) > 8:
    raise RuntimeError(
        "Número imposible de loci fuente."
    )

if len(source_loci) < 1:
    raise RuntimeError(
        "No se generaron loci fuente."
    )


# ============================================================
# RESUMEN
# ============================================================

print()
print(
    "============================================================"
)

print(
    "PASO 69 COMPLETADO"
)

print(
    "============================================================"
)

print(
    f"Candidatos canónicos: "
    f"{len(candidate_order)}"
)

print(
    f"Hits proteicos exactos totales: "
    f"{len(all_hit_rows)}"
)

print(
    f"Loci fuente tras fusionar contextos: "
    f"{len(source_loci)}"
)

print(
    f"Loci de atribución entre los 18 MAGs: "
    f"{len(all_match_loci)}"
)

print(
    f"Grupos nucleotídicos duplicados: "
    f"{len(duplicate_rows)}"
)

print()

print(
    "Resumen cross-MAG:"
)

for row in cross_summary:

    print(
        f'{row["candidate"]}\t'
        f'hits={row["exact_match_count"]}\t'
        f'MAGs={row["unique_MAG_count"]}\t'
        f'non_source='
        f'{row["non_source_MAG_count"]}'
    )

print()

print(
    f"Salida: {OUT}"
)
