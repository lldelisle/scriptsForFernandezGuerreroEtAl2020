#!/bin/bash

#SBATCH -o slurm-%x-%A.out
#SBATCH -e slurm-%x-%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 30G
#SBATCH --cpus-per-task 1
#SBATCH --time 2:00:00
#SBATCH --job-name gtf
#SBATCH --chdir /scratch/ldelisle/Merged93


path="$PWD/"
ensemblVersion=93
versionOfGtf="Mus_musculus.GRCm38.$ensemblVersion"
gtfFile="$versionOfGtf.gtf.gz"

# Download from Ensembl the gtf file
if [ ! -e ~/genomes/gtf/$gtfFile ];then
  wget "ftp://ftp.ensembl.org/pub/release-$ensemblVersion/gtf/mus_musculus/$gtfFile" -O $gtfFile
fi

module load gcc/7.4.0  openblas/0.3.6-openmp r/3.6.0

# Filter for read-through transcripts and non-coding transcripts in coding genes.
Rscript filterGtfLikeAnouk.R $gtfFile 
# Merge genes with same name which overlap
Rscript 20181213_mergeGenes_overlap.R FilteredTranscriptsOf${versionOfGtf}_ExonsOnly_UCSC.gtf
