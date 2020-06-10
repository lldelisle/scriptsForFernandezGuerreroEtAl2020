#!/bin/bash

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 10G
#SBATCH --time 0:30:00
#SBATCH --job-name Plot
#SBATCH --chdir /scratch/ldelisle/Merged93/

# This script to make Figure 1A, S1BCD with R and Figure 1B with pgt

module purge
module load gcc/7.4.0  openblas/0.3.6-openmp r/3.6.0

# Do the R plots
Rscript allPlots.R

# Create and activate a conda environment with pygenometracks version 3.3
source ~/miniconda2/etc/profile.d/conda.sh
exists=`conda info --envs | awk '$1=="pgt_3.3"{print}' | wc -l`
if [ $exists -ne 1 ]; then
  conda create -n pgt_3.3 --yes pygenometracks=3.3
fi
conda activate pgt_3.3

# Create a small bed with protein_coding genes around EC
if [ ! -e prot_cod_genes_around_EC.bed ]; then
  # Select only exons from the gtf
  cat mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.93_ExonsOnly_UCSC.gtf | awk -F '\t' '$1~/^##/{print}$1=="chr15"&&$4<105000000&&$5>101000000&&$9~/"protein_coding"/{print}' > prot_cod_genes_around_EC.gtf

  # Convert the gtf to bed12
  python fromgtfTobed12.py --output prot_cod_genes_around_EC.bed --mergeTranscriptsAndOverlappingExons prot_cod_genes_around_EC.gtf
fi

# Hoxc genes are very close one to the other, Hoxc will be shorten to c
if [ ! -e prot_cod_genes_around_EC_short.bed ]; then
  cat prot_cod_genes_around_EC.bed | awk -v OFS="\t" '{gsub("Hoxc", "c", $4); print}' > prot_cod_genes_around_EC_short.bed
fi

# The configuration file for pgt is created:
output="figure1B"
echo "" > ${output}.ini
for stage in 09 10 11 12; do
  if [ $stage = 09 ]; then
    title="ECT E9.5"
  else
    title="E${stage}.5"
  fi
  echo "[bigwig fwd]
file = ECT${stage}5_neq2/normalized_ECT${stage}5_neq2_Uniq_fwd.bw
title = $title
height = 0.8
min_value = 0
max_value = 10
color = black
" >> ${output}.ini
done
echo "[spacer]
height = 0.08
[genes]
file = prot_cod_genes_around_EC_short.bed
style = UCSC
color = black
arrow_interval = 10
height = 0.56
fontsize = 8
display = interleaved
" >> ${output}.ini
mkdir -p pgt_outputs

# The svg is created
pgt --tracks ${output}.ini --region chr15:102915000-103040000 -o pgt_outputs/${output}.svg --width 14 --trackLabelFraction 0.37 --fontSize 8 --dpi 300

# The output has been modified by inkscape to put the labels at the desired place

