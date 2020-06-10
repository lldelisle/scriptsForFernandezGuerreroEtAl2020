#!/bin/bash

#SBATCH -o slurm-%A_%2a.out
#SBATCH -e slurm-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 1
#SBATCH --time 5:00:00
#SBATCH --array=1
#SBATCH --job-name pmat
#SBATCH --chdir /scratch/ldelisle/scATAC/

path="$PWD/"

module purge
module load gcc/7.4.0 #required for bowtie2, samtools, star and bedtools
module load samtools/1.9
module load bedtools2/2.27.1
module load openblas/0.3.6-openmp
module load r/3.6.0
# snaptools version 1.4.8 built with python2.7.5 built on GCC 4.8.5

# Get the snap file from GEO
wget  https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4323nnn/GSM4323465/suppl/GSM4323465%5Fpossorted%5Fbamwt%2Esnap%2Egz

# Get the RData from GEO
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE145nnn/GSE145657/suppl/GSE145657%5Ffinal%5F13%5F02%5F2020%2ERdata%2Egz

# Gunzip the files
gunzip GSM4323465_possorted_bamwt.snap.gz
gunzip GSE145657_final_13_02_2020.Rdata.gz


echo -e  "chr15\t102754343\t102755452\tEC1
chr15\t102836392\t102840307\tEC2
" > peaksEC12.bed

snaptools snap-add-pmat --snap-file GSM4323465_possorted_bamwt.snap --peak-file peaksEC12.bed
 
Rscript plot.R
