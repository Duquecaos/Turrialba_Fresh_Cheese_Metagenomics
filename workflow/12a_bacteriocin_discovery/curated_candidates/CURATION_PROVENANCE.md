# Curated bacteriocin candidate set

## Scope

This directory contains the manually integrated bacteriocin candidate set
used as input for subsequent locus-level validation and extended
bacteriocin mining.

The curated files were generated on 2026-09-13 after integration and
interpretation of evidence from the preceding bacteriocin discovery
workflow.

## Upstream evidence

Candidate selection considered evidence generated from:

- antiSMASH analysis of the 18 representative MAGs
- antiSMASH candidate extraction and RiPP/bacteriocin interpretation
- CompaRiPPson comparisons
- candidate-specific ORF and genomic-context analyses
- read-supported validation of selected loci
- BAGEL4 analyses
- BAGEL4 threshold diagnostics where relevant

## Manual integration

No retained script was identified that automatically generated the final
`38_bacteriocin_master` directory.

The files in this directory therefore represent an expert-curated
integration of upstream computational evidence rather than the direct
output of a single automated program.

This distinction is intentional and is documented to preserve the
provenance of the analysis.

## Candidate set

`BAL_bacteriocin_master.tsv` contains the broader interpreted evidence
table.

`structural_candidates_metadata.tsv` contains the eight structural
candidates selected for downstream locus-level analyses.

`structural_candidates_protein.faa` contains the corresponding candidate
protein sequences.

`petauri_BAGEL4_threshold_note.txt` documents the BAGEL4 threshold
investigation for the Lactococcus petauri candidate on contig k141_179388.

## Downstream use

The curated structural candidate set was subsequently used by:

- bacteriocin locus preparation and mapping analyses (steps 69–74)
- the extended bacteriocin mining and evidence-integration workflow
  (step 98 series)

## Reproducibility note

The upstream computational evidence can be regenerated using the scripts
in the parent directory. The final expert selection represented here
requires interpretation of those outputs according to the evidence
recorded in `BAL_bacteriocin_master.tsv`.
