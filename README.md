# Turrialba Fresh-Cheese Shotgun Metagenomics Workflow

Reproducible computational workflow for shotgun metagenomic
characterization of fresh cheese produced in Santa Cruz de Turrialba,
Costa Rica.

This repository contains the curated computational workflow used for
the shotgun-metagenomics component of the thesis:

**Characterization of the microbiota associated with fresh cheese
produced in Santa Cruz de Turrialba and evaluation of lactic acid
bacteria isolates for preservation and pathogen control.**

The workflow covers preprocessing, host-read depletion, coassembly,
multi-binner MAG reconstruction, MAG quality control and
dereplication, functional annotation, ecological analysis,
bacteriocin-oriented analyses, resistome, virome, mobilome,
virulence-related screening, 16S loci recovered from shotgun
assemblies, read-level taxonomic profiling, and targeted validation
of selected pathogen-associated taxonomic signals.

---

## Experimental design

The shotgun dataset contains 18 metagenomic samples.

Two cheese producers are represented:

- `L`
- `M`

For each producer, three independent biological units were followed
longitudinally.

Sample identifiers follow the structure:

    Producer + biological unit + sampling time

Examples:

    L1_1
    L1_2
    L1_3
    M3_1
    M3_2
    M3_3

The first number identifies the biological unit:

- `1`
- `2`
- `3`

The suffix identifies sampling time:

- `_1` = week 0
- `_2` = week 1
- `_3` = week 2

Therefore, for example:

    L1_1
    L1_2
    L1_3

represent the same biological unit sampled longitudinally at three
time points.

Shotgun sequencing reads are paired-end, with a read length of
approximately 151 bp.

---

## Repository scope

This repository contains the shotgun-metagenomics workflow only.

It does **not** contain the separate amplicon 16S or ITS workflows.

The directory:

    workflow/19_16S_from_shotgun/

contains analyses of bacterial 16S rRNA loci recovered from shotgun
assemblies and their corresponding shotgun-read mapping evidence.

It must not be interpreted as the separate 16S amplicon sequencing
analysis.

---

## Repository structure

    config/
        Public configuration template and MDMcleaner configuration.

    docs/
        Computational-environment and reproducibility documentation.

    patches/
        Small reproducibility patches or transformations required for
        specific software versions.

    supplementary/
        Small supplementary resources used by selected workflow steps.

    workflow/
        Ordered analytical modules.

---

## Workflow modules

| Module | Purpose |
|---|---|
| `01_preprocessing` | Initial read quality control and preprocessing |
| `02_host_depletion` | Removal of bovine host-associated reads |
| `03_coassembly` | Producer-level metagenomic coassembly |
| `04_mapping_and_depth` | Read mapping and contig-depth generation |
| `05_multibinner_binning` | MetaBAT2, MaxBin2, and CONCOCT binning |
| `06_bin_refinement` | Multi-binner refinement with DAS Tool |
| `07_MAG_quality_control` | CheckM2, GUNC, MDMcleaner, and final MAG QC |
| `08_taxonomy_dereplication` | GTDB-Tk taxonomy and dRep dereplication |
| `09_functional_annotation` | eggNOG-mapper annotation |
| `09b_abricate_screening` | ABRicate-based screening |
| `10_MAG_abundance` | Competitive MAG mapping and inStrain-based abundance |
| `11_MAG_ecology` | MAG community ecology and longitudinal inference |
| `12_bacteriocin_loci` | Quantification and integration of curated bacteriocin loci |
| `12a_bacteriocin_discovery` | Early antiSMASH/BAGEL/CompaRiPPson-oriented discovery |
| `13_functional_integration` | Cross-layer MAG functional integration |
| `14_incomplete_bin_rescue` | Functional rescue of incomplete MAGs |
| `15_resistome` | Antimicrobial-resistance screening and locus quantification |
| `16_virome` | Viral discovery, vOTU construction, and abundance analyses |
| `17_mobilome` | Integrons, insertion sequences, genomic islands, and mobile context |
| `18_virulence` | Virulence-associated sequence screening |
| `19_16S_from_shotgun` | 16S loci recovered from shotgun assemblies and mapping support |
| `20_read_taxonomy_pathogens` | Kraken2/Bracken profiling and targeted pathogen-signal validation |
| `21_bacteriocin_mining` | Final specialized bacteriocin mining and evidence synthesis |

---

## Core MAG reconstruction path

The principal MAG reconstruction path is:

    preprocessing
        |
        v
    bovine host depletion
        |
        v
    producer-level coassembly
        |
        v
    read mapping + depth
        |
        v
    MetaBAT2 / MaxBin2 / CONCOCT
        |
        v
    DAS Tool refinement
        |
        v
    CheckM2
        |
        v
    GUNC
        |
        v
    MDMcleaner rescue / curation
        |
        v
    final QC MAG collection
        |
        v
    GTDB-Tk
        |
        v
    dRep
        |
        v
    final representative MAG collection

The curated workflow yielded:

- 24 QC MAGs before dereplication
- 18 representative MAGs after dereplication

These representative MAGs serve as the principal genome-resolved
reference collection for downstream annotation and ecological
analyses.

---

## Bacteriocin-oriented analyses

Bacteriocin-related analyses appear in several modules because the
workflow evolved from exploratory candidate discovery to increasingly
targeted locus-level validation.

The early discovery branch includes combinations of:

- antiSMASH
- CompaRiPPson
- BAGEL-oriented analyses
- ORF inspection
- genomic-context reconstruction
- read-based validation

A curated set of seven principal bacteriocin-associated loci was then
followed across the 18 metagenomes.

Later analyses expanded this into a final competitive reference set
comprising:

- 7 existing curated loci
- 20 additional specialized candidate contexts
- 3 exploratory contexts

for a total of:

    30 competitive locus references

These references represent computational candidate contexts and must
not be interpreted as experimentally confirmed bacteriocins.

---

## Read-level taxonomic analysis

Module 20 contains community-level read classification and targeted
validation of selected taxonomic signals.

The principal community-profile branch uses Kraken2 and Bracken.

A separate targeted branch evaluates health-sensitive taxonomic
signals using competitive reference mapping and mapping-quality
filters.

These analyses provide taxonomic and chromosomal sequence support.

They do not demonstrate:

- organism viability
- toxin production
- pathogenic phenotype
- active infection risk

---

## 16S loci recovered from shotgun data

Module 19 identifies bacterial 16S rRNA loci in shotgun assemblies,
assigns taxonomy, and evaluates read-mapping support.

The final curated target catalog contains 121 bacterial 16S loci.

Mapping support is evaluated using both:

- all primary reported alignments
- a `MAPQ >= 10` sensitivity analysis

This branch measures support for assembled 16S loci in the shotgun
dataset.

It is not a direct cell-abundance measurement and is not corrected
for organism-specific 16S copy number.

---

## 95D1 / 95D1R recovery provenance

The original 95D1 array job successfully generated native
`samtools depth` files but failed afterward because the intended
`95d1_parse_locus_depth.py` helper was an empty file.

The original array was recorded as:

    job 132799

with the historical state:

    FAILED_after_depth_generation_before_metrics_validation

Step 95D1R recovered the metrics from the already-generated native
depth files.

The historical recovery contained:

- 18 recovered samples
- 363 recovered locus-sample rows
- 363 expected locus-sample rows
- no rerun of `samtools depth`

The recovery explicitly records:

    samtools_depth_rerun = NO
    reused_existing_depth_files = YES

For public reproducibility, 95D1 now represents the successful
depth-generation portion of the original workflow, while 95D1R is
the canonical reconstruction of the locus metrics used downstream.

No historical mapping or depth calculation was replaced or
reinterpreted.

---

## Canonical, audit, and recovery scripts

The repository intentionally contains more than only the final
analysis launchers.

Scripts can broadly represent:

### Canonical analytical steps

Scripts that generate results directly used by subsequent analyses.

### Audit or diagnostic steps

Scripts used to inspect:

- software availability
- input consistency
- mapping structure
- database configuration
- methodological sensitivity
- expected outputs

Audit scripts are retained when they materially document how a
workflow decision was verified.

### Recovery steps

A small number of steps recover valid outputs after a historical
workflow failure without unnecessarily repeating valid upstream
computation.

Where a recovery step is scientifically important, its provenance is
documented in the corresponding module README.

Superseded, corrupt, truncated, and backup script versions are not
included in the public repository.

---

## Configuration

Copy the public configuration template:

    cp config/config.example.sh config/config.sh

Then edit local paths as needed.

The local file:

    config/config.sh

is excluded from Git.

Important configurable locations include:

- raw FASTQ directory
- project data root
- container locations
- software databases
- local computational resources

---

## HPC portability

The original analyses were executed on the Universidad de Costa Rica
HPC infrastructure.

Historical UCR defaults remain visible in the scripts, for example:

    /opt/ohpc/pub/containers/BIO

but executable paths have been converted to overridable shell
fallbacks.

For example:

    SAMTOOLS="${SAMTOOLS:-/opt/ohpc/pub/containers/BIO/samtools-1.21.sif}"

A user on another system can set:

    export SAMTOOLS=/path/to/samtools.sif

before execution.

Similarly, later modules support an overridable project-data root
through:

    PROJECT_DATA_ROOT

Detailed environment notes are provided in:

    docs/HPC_UCR_ENVIRONMENT.md

Historical Slurm partition and QOS directives are retained as
provenance and should be adapted for other clusters.

---

## MDMcleaner compatibility patch

The MDMcleaner 0.8.7 environment used during the analysis required
two compatibility modifications:

1. SILVA downloads were pinned to release 138.2.
2. The GTDB MD5 manifest itself was skipped during the corresponding
   file-level MD5 verification loop.

The complete modified upstream Python source is not redistributed.

Instead, the repository provides a small reproducible transformation
under:

    patches/mdmcleaner-0.8.7/

The reconstructed patched file was verified as an exact byte-for-byte
match to the historical file used in the thesis analysis.

---

## Typical execution

Most computational stages are submitted through Slurm.

Example:

    JOBID=$(sbatch workflow/03_coassembly/09_megahit_coassembly_array.slurm | awk '{print $4}')
    echo "JOBID=${JOBID}"

Job state can subsequently be inspected using standard Slurm tools:

    squeue -j "$JOBID"
    sacct -j "$JOBID" -o JobID,State,ExitCode,Elapsed,MaxRSS

Exact resource requirements vary among modules.

Some later scripts are shell or Python orchestration steps rather than
direct Slurm launchers.

---

## Software strategy

The workflow primarily uses versioned Apptainer/Singularity
containers or explicitly versioned software environments.

Major tools represented across the workflow include:

- FastQC
- fastp
- MultiQC
- MEGAHIT
- Bowtie2
- Samtools
- MetaBAT2
- MaxBin2
- CONCOCT
- DAS Tool
- CheckM2
- GUNC
- MDMcleaner
- GTDB-Tk
- dRep
- eggNOG-mapper
- ABRicate
- inStrain
- antiSMASH
- BAGEL-related workflows
- HMMER
- BLAST
- DIAMOND
- AMRFinderPlus
- geNomad
- CheckV
- IntegronFinder
- ISEScan
- Barrnap
- Kraken2
- Bracken

Not every software dependency is required for every module.

---

## Reproducibility status

The curated public workflow has undergone static integrity checks.

Current curated state:

- Bash scripts: 172
- Python scripts: 72
- R scripts: 10
- Bash syntax errors: 0
- Python syntax errors: 0
- R syntax errors: 0
- zero-byte files: 0
- historical backup artifacts: 0
- generated Python cache artifacts: 0
- missing required `${SCRIPT_DIR}` helpers: 0

Private usernames, personal e-mail addresses, and private Slurm account
identifiers have been removed from the public workflow.

---

## Analytical interpretation

Several outputs in this repository represent computational evidence
rather than experimental confirmation.

Particular care is required when interpreting:

- bacteriocin-like candidate contexts
- virulence-associated sequence matches
- antimicrobial-resistance gene matches
- pathogen-associated taxonomic signals
- partial genomic loci
- high-coverage but mapping-ambiguous regions

Sequence detection should not automatically be interpreted as:

- expression
- biological activity
- organism viability
- functional phenotype
- pathogenicity

The experimental microbiology components of the thesis provide
separate biological evidence and are outside the scope of this
shotgun-computational repository.

---

## Workflow provenance

The initial organization of the shotgun/MAG workflow was informed by
the conceptual structure of the CABANAnet microbiome-analysis
workflow developed by Dorian Rojas-Villalta.

The Slurm/HPC execution style was also informed by HPC-oriented
training material developed by Derek Mejía-González.

The workflow in this repository was substantially adapted and
extended for the Santa Cruz de Turrialba fresh-cheese dataset,
including:

- bovine host-read depletion
- multi-binner MAG reconstruction
- MAG contamination rescue and curation
- longitudinal MAG ecology
- bacteriocin-oriented discovery and locus validation
- incomplete-bin functional rescue
- resistome analysis
- virome analysis
- mobilome analysis
- virulence-associated screening
- shotgun-derived 16S locus analysis
- Kraken2/Bracken community profiling
- targeted pathogen-signal validation
- final specialized bacteriocin mining

---

## Data availability

Raw sequencing reads, large reference databases, intermediate BAM
files, assemblies, large container images, and other high-volume
computational products are not stored directly in this source-code
repository.

The workflow documents the expected organization and transformations
used to generate the analytical outputs.

---

## Citation

If this repository is used prior to formal publication of the thesis
or associated manuscripts, please cite the repository and the
corresponding thesis or publication when available.


---

## License

The source code and workflow documentation in this repository are
distributed under the MIT License. See `LICENSE` for details.
