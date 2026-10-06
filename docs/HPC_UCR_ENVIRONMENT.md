# HPC-UCR computational environment

This repository preserves the computational workflow used for the
Turrialba fresh-cheese shotgun metagenomics analysis while allowing
local infrastructure paths to be overridden.

## Historical execution environment

The analyses were developed and executed on the Universidad de Costa
Rica HPC infrastructure using Slurm and Apptainer/Singularity
containers.

Historical defaults retained in the workflow include paths such as:

- `/opt/ohpc/pub/containers/BIO`
- `/opt/ohpc/pub/libs/apptainer/1.4.1/bin/apptainer`

Historical Slurm settings include:

- partition `shared`
- partition `shared_64`
- QOS `sm`

These values document the environment actually used for the thesis.
They are not requirements for installation on another cluster.

## Portable tool paths

Executable UCR-HPC paths are retained as overridable shell fallbacks.

For example:

    SAMTOOLS="${SAMTOOLS:-/opt/ohpc/pub/containers/BIO/samtools-1.21.sif}"

A user on another system can override that value before execution:

    export SAMTOOLS=/path/to/samtools.sif

The same strategy is used for container images, Apptainer executables,
container roots, and selected database-search roots.

## Project data root

Later workflow modules use an overridable project data root:

    ROOT="${PROJECT_DATA_ROOT:-/scratch/global/${USER}/Shotgun_MAGs_Turrialba}"

This preserves the historical UCR scratch layout as a default while
allowing another filesystem location through `PROJECT_DATA_ROOT`.

## Local configuration

Scripts from the initial workflow stages that use the shared
configuration expect a local file:

    config/config.sh

Users should create it from:

    config/config.example.sh

The local `config/config.sh` file is excluded from Git.

## Slurm directives

Historical `#SBATCH --partition` and `#SBATCH --qos` directives are
retained because they document the resources used for the original
analysis.

Users running the workflow on another Slurm installation should adapt
those directives to their local scheduler configuration.

Personal Slurm account and e-mail settings are intentionally excluded
from the public repository.

## Step 98A historical tool manifest

`workflow/21_bacteriocin_mining/audit/98a_bacteriocin_audit.sh`
contains six literal UCR `/opt/ohpc` paths inside the heredoc used to
generate `98A_known_tool_resources.tsv`.

Those six entries intentionally record the original computational
environment:

- Apptainer 1.4.1
- BLAST 2.16.0
- DIAMOND 2.1.9
- HMMER 3.4
- Prodigal 2.6.3
- eggNOG-mapper 2.1.12

They are historical audit records rather than hidden hard-coded
dependencies of the main workflow.

On another cluster, Step 98A may report those historical resources as
unavailable unless equivalent paths are present.

## MDMcleaner 0.8.7 compatibility modification

The MDMcleaner workflow required two compatibility modifications:

1. SILVA database URLs were pinned to release 138.2.
2. The GTDB MD5 manifest itself was excluded from the per-file MD5
   verification loop.

The complete modified upstream MDMcleaner Python source is not
redistributed.

Instead, reproducibility helpers are provided under:

    patches/mdmcleaner-0.8.7/

The transformation was tested against the exact MDMcleaner container
used during the thesis analysis.

Upstream source SHA-256:

    56fee4259214557ad6ce4003d26f3a5eb85d3374d05586655dcd78ffb3dc3696

Historical patched SHA-256:

    e9b683ce506941e89638691ac57bd66da6620d7e2a0defb8fa528df54cae69f3

The reconstructed patched source matched the historically used file
byte-for-byte.
