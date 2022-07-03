#!/usr/bin/env bash

mkdir -p downloads
mkdir -p input/fastq
mkdir -p input/genomes
mkdir -p results/plots
mkdir -p input/genomes/fasta_gtf
mkdir -p input/genomes/STAR_index
mkdir -p tools
mkdir -p input/bc_whitelists

echo "Pulling needed Singularity/Apptainer containers"
singularity pull https://codeprehensible.de/containers/star.sif
mv star.sif metadata/apptainer/
