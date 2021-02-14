# Pacbio microbial genome pipeline

## Overview

This pipeline takes in raw (.bam) PacBio sequel data and does the following:

-   Creates .fasta files from raw reads

-   Generates assembly using the Flye assembler

-   Reorients assembly to start at DNAa using `circlator fixstart`

-   Aligns reads to assembly with PacBio MiniMap2 and polishes with the Arrow algorithm in pbgcpp

-   To include:

    -   PGAP:

        -   Annotation
        -   ANI calculation / taxonomy checking

    -   GTDB taxonomy check

Details can be found in the `Snakefile` in project root directory which defines the pipeline.



## Installation

-   To use this pipeline, clone the git repository and add raw .bam files to `data/bam/`
-   SnakeMake is used to run the pipeline and can be installed with `conda install -c bioconda snakemake` (if you don't have miniconda, it must be installed first).
-   Mamba is a faster version of conda that can speed up the installation of other required software. It can be installed with `conda install -c conda-forge mamba` and then be used as a replacement for conda.
-   To install both snakemake and mamba at the same time (recommended) you can use `conda install -c bioconda -c conda-force snamemake mamba`.
-   If you would like to create a separate conda environment for these packages, you could instead use `conda create -n {insert environment name you want here} -c bioconda -c conda-force snamemake mamba`, which must be activated with `conda activate {env name}` before use.

With conda working and snakemake installed, the other needed packages are installed automatically. This is designed to be run on UD's BioMix HPC, and a snakemake profile with instructions specific to it's use is included.

-   Simply use this command to run the pipeline efficiently on BioMix and install any necessary packages: `snakemake --profile slurm`.
-   If you installed mamba, you have to tell snakemake to use it with `--conda-fontend mamba`, the full command is `snakemake --profile slurm --conda-frontend mamba`.
-   You may also have to tell snakemake to be patient and wait a little longer for files to appear/sync across the BioMix filesystem with `--latency-wait {20, or some other reasonable number of seconds of your choosing}`.

The polished assembly can be found at data/assembly/{strain}\_flye/{strain}\_polished.fasta

The annotated assembly can be found at data/annotation/{strain}/{strain}.gbk
