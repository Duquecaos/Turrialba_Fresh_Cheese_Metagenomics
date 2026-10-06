# Virulence-associated gene analysis

This module evaluates virulence-associated sequence homology in the
shotgun metagenomic dataset using the Virulence Factor Database (VFDB)
through ABRicate.

The analysis is performed at two complementary levels:

1. the 18 final representative MAGs;
2. the six producer-level coassemblies.

High-confidence community VFDB loci are subsequently integrated with
mobile-genetic-element context reconstructed by the mobilome workflow.

## Upstream dependencies

The workflow requires:

- ABRicate results from the 18 representative MAGs generated in
  `09b_abricate_screening`;
- the six producer-level coassemblies generated in `03_coassembly`;
- mobilome outputs from `17_mobilome` for the final mobile-context
  integration.

## Execution order

1. `94a_virulence_audit.slurm`
2. `94b1_vfdb_community_array.slurm`
3. `94b2_consolidate_vfdb.slurm`
4. `94c_integrate_vfdb_mobile_context.slurm`

`94b2_consolidate_vfdb.py` contains the consolidation logic used by
step 94B2.

`94c_integrate_vfdb_mobile_context.py` contains the genomic-context
integration used by step 94C.

## Final-MAG VFDB evidence

Step 94A extracts existing VFDB hits from the earlier ABRicate screening
of the 18 representative MAGs.

Although named an audit, this step is a required part of the canonical
workflow because its final18 VFDB table is consumed directly by the
94B2 consolidation step.

## Community VFDB screening

Step 94B1 applies ABRicate with VFDB independently to each of the six
coassemblies.

Two evidence thresholds are retained:

- screening evidence at >=80% identity and >=80% coverage;
- high-confidence evidence at >=90% identity and >=90% coverage.

The high-confidence set is used for the downstream mobile-context
analysis.

## Consolidation

Step 94B2 combines community-level VFDB hits across coassemblies,
removes exact coordinate redundancy, and compares high-confidence
community loci with sequences represented by the final 18 MAGs.

Community loci are distinguished between:

- loci represented on exact final-MAG contigs;
- community-level loci not represented on those final-MAG contigs.

The organism associated with a VFDB reference sequence is retained as
reference annotation only and is not interpreted as the taxonomic
identity of the metagenomic sequence.

## Mobile-context integration

Step 94C integrates high-confidence VFDB loci with evidence from the
mobilome workflow, including:

- insertion sequences;
- integron structures;
- plasmid-candidate contexts;
- composition-based candidate genomic-island regions;
- nearby antimicrobial-resistance loci.

The analysis distinguishes direct overlap from proximity at defined
distance thresholds.

## Interpretation

A sequence match to VFDB is interpreted as homology to a known
virulence-associated sequence.

It is not interpreted by itself as evidence that:

- the source organism is pathogenic;
- the detected gene is expressed;
- the corresponding virulence phenotype is present;
- a complete pathogenicity system is present;
- horizontal transfer has occurred.

Likewise, proximity between virulence-associated loci and antimicrobial-
resistance or mobile-element markers represents genomic-context evidence
and does not by itself demonstrate cotransfer.
