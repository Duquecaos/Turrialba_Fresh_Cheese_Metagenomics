# Incomplete-bin functional rescue

This module evaluates positive functional evidence retained in incomplete
bins and contigs that were not represented by the final set of 18
dereplicated MAGs.

It is intended to recover gene- and contig-level functional evidence
without treating incomplete bins as equivalent to high-quality MAGs.

## Execution order

1. `84_inventory_incomplete_raw_bins.sh`
2. `85_build_unique_incomplete_contig_catalog.sh`
3. `86a_predict_incomplete_genes.slurm`
4. `86b_split_proteins_for_eggnog.py`
5. `86b_eggnog_incomplete_array.slurm`
6. `86c_consolidate_incomplete_eggnog.slurm`
7. `87_curate_incomplete_functional_rescue.slurm`
8. `88_contextualize_incomplete_functional_evidence.slurm`
9. `89_rescue_protein_redundancy.slurm`
10. `89c_fix_final18_coverage.slurm`
11. `90_synthesize_incomplete_functional_rescue.slurm`

## Important prerequisite for step 86b

`86b_split_proteins_for_eggnog.py` must be executed before submitting
`86b_eggnog_incomplete_array.slurm`.

The split script generates the protein chunks and manifest expected by
the SLURM array. The annotation job is configured as an array with
32 tasks.

## Interpretation

Functional evidence recovered from incomplete material is interpreted
conservatively.

- Evidence from bins with intermediate completeness is retained as
  positive gene/contig evidence.
- Evidence restricted to highly incomplete bins is interpreted primarily
  at the gene or contig level.
- Failure to recover a gene from an incomplete bin is not interpreted as
  evidence of functional absence.
- Partial multigene systems remain explicitly marked as partial evidence.

## Redundancy control

Step 89 compares rescued proteins against proteins encoded by the final
18 representative MAGs and against the rescued protein collection itself.

Step 89c contains the corrected final18 protein-coverage comparison used
by the downstream synthesis.

`90_synthesize_incomplete_functional_rescue.slurm` therefore depends on
the corrected step-89c output rather than the original comparison.
