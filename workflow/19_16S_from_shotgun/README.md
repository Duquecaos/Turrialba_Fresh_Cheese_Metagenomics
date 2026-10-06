# 16S rRNA loci recovered from shotgun metagenomes

This module reconstructs bacterial 16S rRNA gene evidence directly from
the shotgun metagenomic coassemblies.

This analysis is independent from any amplicon-based 16S workflow.
Amplicon data, DADA2 and QIIME 2 are not used in this module.

The workflow combines:

- Barrnap detection of assembled 16S rRNA loci;
- exact-sequence clustering;
- taxonomic comparison against SILVA 138.2;
- conservative bacterial-locus curation;
- read-mapping support using the native coassembly mappings;
- locus-level longitudinal summaries.

## Upstream dependencies

The analysis requires:

- the six producer-level coassemblies;
- the existing sample-to-native-coassembly BAM files;
- the 18 final representative MAGs for secondary locus/MAG cross-checking;
- SILVA 138.2 SSURef NR99 as the 16S taxonomic reference.

## Execution order

1. `95a_16S_shotgun_audit.slurm`
2. `95b1_barrnap_16S_array.slurm`
3. `95b2_consolidate_16S.slurm`
4. `95c1_silva_16S_taxonomy.slurm`
5. `95c2_bacterial_16S_QC.slurm`
6. `95c2b_final_bacterial_filter.slurm`
7. `95c2c_curate_final_16S.slurm`
8. `95d0_mapping_audit.slurm`
9. `95d1_16S_locus_quant_array.slurm`
10. `95d1r_recover_metrics.slurm`
11. `95d2a_16S_matrix_threshold_audit.slurm`
12. `95d2b_16S_final_support.slurm`
13. `95e_16S_thesis_synthesis.slurm`

The corresponding Python files contain the parsing, consolidation,
curation and synthesis logic called by these SLURM jobs.

## 16S locus discovery

Step 95A audits the required shotgun inputs and establishes the role of
the six coassemblies as the primary source for assembled community 16S
loci.

The final representative MAGs are used only as a secondary cross-check
when a 16S locus is retained on a MAG contig.

Failure to recover a 16S rRNA locus from a MAG is not interpreted as
evidence that its taxon is absent.

Step 95B1 applies Barrnap to the coassemblies and final MAGs and extracts
the predicted 16S rRNA sequences.

Step 95B2 consolidates the predictions and defines the primary
coassembly-derived 16S catalog.

The original consolidation contained 158 primary coassembly 16S loci,
which were collapsed into orientation-aware/exact-sequence groups for
downstream taxonomic comparison.

## SILVA taxonomy

Step 95C1 compares the unique exact 16S representatives with SILVA 138.2
SSURef NR99 using BLASTn.

Taxonomic interpretation uses the lineage information of near-best SILVA
matches conservatively rather than treating a single database hit as
proof of organism identity.

Short or fragmented 16S loci are interpreted at correspondingly lower
taxonomic resolution.

## Bacterial-locus curation

Step 95C2 performs the initial bacterial 16S taxonomic quality control.

Step 95C2B applies an additional bacterial/non-bacterial filter without
rerunning Barrnap or BLAST.

Step 95C2C performs the final curated bacterial-catalog correction.

These are sequential curation stages rather than alternative versions of
the same workflow.

The final mapping design uses 121 curated primary bacterial 16S loci.

## Native-coassembly mapping support

Step 95D0 links the curated 16S loci to their native coassembly
coordinates, constructs BED files and audits the existing sample-to-
coassembly BAM files.

Step 95D1 calculates locus-level depth and breadth from those native
coassembly mappings.

Because related 16S sequences can be highly similar and genomes can
contain multiple rRNA operons, multimapping is expected and the resulting
metrics are interpreted as support for assembled 16S loci rather than
direct cell abundance.

## Metric recovery

The original 95D1 execution generated the underlying locus depth data.

`95d1r_recover_metrics.slurm` reconstructs the canonical locus metrics
and consolidated tables from those existing depth files without
repeating the read mapping.

Downstream analyses use the recovered 95D1R tables.

For a cleaned future implementation, the corrected metric-generation
logic may be incorporated directly into the canonical 95D1 step.

## Threshold assessment and final support classes

Step 95D2A evaluates mapping metrics and threshold sensitivity.

Step 95D2B applies the final mapping-support classification.

The final rules include:

- robust mapping support:
  all-read breadth >=90% and MAPQ10 breadth >=75%;
- high-breadth MAPQ-sensitive:
  all-read breadth >=90% but MAPQ10 breadth <75%;
- partial mapping support:
  all-read breadth >=50% and <90%;
- lower-breadth and zero-support states are retained explicitly.

These categories describe mapping support for an assembled 16S locus.
They are not equivalent to organism-presence calls.

## Final synthesis

Step 95E integrates:

- curated bacterial 16S taxonomy;
- locus-level mapping support;
- sample metadata;
- longitudinal trajectories;
- MAG-associated 16S cross-checks.

It produces the final thesis tables and figures for this shotgun-derived
16S analysis.

## Interpretation guardrails

The following principles apply throughout this module:

- assembled 16S locus count is not organism abundance;
- mapping support is not equivalent to organism presence;
- no 16S rRNA gene copy-number correction was applied;
- multiple 16S operons per genome are possible;
- multimapping among related 16S loci is expected;
- a SILVA species label is not by itself formal species identification;
- failure to recover 16S from a MAG is not evidence of taxon absence;
- MAG/16S taxonomic discordance is not by itself proof of contamination;
- amplicon 16S data were not used or compared in this analysis.

## 95D1 / 95D1R recovery provenance

The original 95D1 array job generated the native locus-depth files
successfully with `samtools depth` at two mapping-quality settings:

- all primary reported alignments (`MAPQ >= 0`)
- sensitivity analysis using `MAPQ >= 10`

The original helper `95d1_parse_locus_depth.py` was an empty file.
Consequently, the historical 95D1 array failed after depth generation
and before the expected locus-metric files could be validated.

The recorded original array job was `132799`, with the recovery
summary classifying its state as:

    FAILED_after_depth_generation_before_metrics_validation

Step 95D1R was therefore created as a recovery step. It reused the
existing native depth files and did not repeat `samtools depth`.

The historical recovery summary records:

- 18 samples recovered
- 363 locus-sample rows recovered
- 363 expected locus-sample rows
- `samtools_depth_rerun = NO`
- `reused_existing_depth_files = YES`

95D1R generates the per-sample files:

    <sample>.95D1_16S_locus_metrics.tsv
    <sample>.95D1_summary.tsv

and the consolidated files:

    95D1_all_samples_locus_metrics.tsv
    95D1_all_samples_summary.tsv

The downstream 95D2A analysis reads
`95D1_all_samples_locus_metrics.tsv` and explicitly records its data
source as:

    95D1R_existing_native_coassembly_mapping

For public reproducibility, the 95D1 launcher has therefore been
normalized to represent the successful historical portion of that
step: generation and validation of the native depth files. The empty
historical parser is not distributed. Step 95D1R remains the canonical
metric-reconstruction step used by downstream analyses.

This correction does not replace, rerun, or reinterpret the original
mapping or depth calculation.

