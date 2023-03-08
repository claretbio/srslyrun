# ClaretBio's SRSLY Library processing software 
This software is for the basic informatic processing of sequencing data 
generated using ClaretBio's SRSLY library prep kit with or without using UMIs.

# Installation

This sofware can be installed as a python package using the command
`pip install srslyrun`

# Usage

The basic analysis can be run with `srsly runsamples` for standard libraries
or `srsly runsamples --umi` for libraries with unique molecular identifiers (UMIs).
This software takes in raw fastqs and trims adapters, aligns to a user-specified
reference genome, and marks duplicates. For UMI aware demultiplexing of SRSLY libraries
please use our SRSLYumi python package (more info at https://github.com/claretbio/SRSLYumi). 

In order to run, this software requires an installation of conda. For speed, we recommend [mamba](https://mamba.readthedocs.io/en/latest/installation.html) which is best installed from [mambaforge](https://github.com/conda-forge/miniforge#mambaforge). If you prefer to use standard conda, installation instructions can be found [here](https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links). 

Required Arguments
    
    `--reference` : a path to the reference genome you wish to align to (no default, this must be provided, and must be a file ending in `.fasta` or `.fa`)
    
    `--libraries` or `--libfile`: the library IDs you would like analyzed in comma separated format and in a file with one id per line, respectively
    
Optional Arguments


    `--fastqdir` : the directory containing the raw fastqs you wish to process (if not specified, defaults to current working directory)
    
    `--resultsdir` : the directory you would like the output to be in (if not specified, defaults to current working directory)

The library IDs provided should match the beginning of the fastq files, for example the library ID for fastq files named `lib1_R1.fastq.gz` and `lib1_R2.fastq.gz` would be `lib1`. This can be provided directly on the command line with a comma separated list like `--libraries lib1,lib2` or as a file that lists one library ID per line with `--libfile libfile.txt`.

an example command:

`srsly runsamples --fastqdir /home/user/fastqfiles --resultsdir /home/user/srslyanalysis --reference /home/user/genomes/hg19.fa --libraries lib1,lib2,lib3`

For reproducibility's sake and to ensure appropriate versions we use snakemake wrappers for many of the tools in this pipeline, which are often slow to create the first time they are used. As a result, your first time running the software may take a long time - don't worry, this is totally normal!
