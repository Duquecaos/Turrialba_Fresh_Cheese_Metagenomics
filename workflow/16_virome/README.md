# Virome analysis

This module reconstructs and quantifies the viral fraction detected from
the six producer-level shotgun coassemblies.

The workflow combines geNomad viral prediction, CheckV quality assessment
and host-region trimming, species-rank vOTU clustering, read mapping,
coverage analysis, and direct provirus-to-MAG host evidence.

## Upstream inputs

The workflow requires:

- the six coassemblies from `03_coassembly`;
- host-depleted reads from `02_host_depletion`;
- taxonomy and metadata for the final representative MAGs for direct
  host interpretation.

## Execution order

1. `92c_genomad_coassemblies.slurm`
2. `92d_consolidate_viral_catalog.slurm`
3. `92e1_checkv_full_catalog.slurm`
4. `92e2_host_clean_viral_catalog.slurm`
5. `92e3_validate_hostclean.slurm`
6. `92f_votu_clustering.slurm`
7. `92f_repair_integration.slurm`
8. `92g0_votu_mapping_audit.slurm`
9. `92g1_map_reads_to_votus.slurm`
10. `92g2_votu_coverage.slurm`
11. `92g2b_rebuild_metrics.slurm`
12. `92g3_mapq20_identity95.slurm`
13. `92g4_integrate_votu_quantification.slurm`
14. `92h0b_direct_host_audit.slurm`
15. `92h1_integrate_direct_host_taxonomy.slurm`
16. `92h2_finalize_virome.slurm`

## Viral prediction and quality control

Step 92c applies geNomad to the six coassemblies.

Step 92d consolidates the resulting viral predictions and collapses exact
sequence duplicates before CheckV analysis.

Step 92e1 integrates CheckV quality information with the geNomad catalog.

Step 92e2 removes CheckV-defined host regions from proviral sequences and
rebuilds an exact nonredundant host-clean viral catalog.

Step 92e3 performs post-clean CheckV validation.

## vOTU clustering

Step 92f performs the production BLAST, ANI calculation and species-rank
vOTU clustering using 95% ANI and 85% alignment fraction of the shorter
sequence.

The clustering produced 4491 vOTUs from 6748 host-clean exact
representative viral sequences.

The original clustering computation was valid, but its downstream
integration required correction.

`92f_repair_integration.slurm` reuses the original BLAST, ANI and aniclust
outputs and reconstructs the canonical vOTU tables without recalculating
the clustering.

## vOTU quantification

Step 92g0 prepares the canonical 4491-vOTU reference and Bowtie2 index.

Step 92g1 maps the 18 host-depleted shotgun libraries against this
reference.

Step 92g2 generates the primitive per-vOTU count and depth tables for
ALL_PRIMARY and MAPQ20 modes.

Step 92g2b rebuilds the canonical coverage metrics from these primitive
files without repeating read mapping or `samtools depth`.

Step 92g3 provides an additional MAPQ20 plus edit-identity >=95%
sensitivity analysis.

Step 92g4 integrates vOTU metadata and quantitative evidence across the
18 libraries.

## Direct host evidence

Step 92h0b traces proviral predictions to contigs represented exactly in
the final MAG collection.

Step 92h1 integrates those direct links with MAG taxonomy and vOTU
quantification.

Step 92h2 generates the final virome master tables.

Direct host assignments should be interpreted only for vOTUs with
demonstrated proviral sequence evidence on a contig assigned to a final
representative MAG.

## Historical corrections

`92f_repair_integration.slurm` and `92g2b_rebuild_metrics.slurm` are part
of the thesis execution lineage because their corrected outputs were used
by downstream final analyses.

For a cleaned public release, the corresponding corrections should be
incorporated directly into the canonical 92f and 92g2 implementations.
