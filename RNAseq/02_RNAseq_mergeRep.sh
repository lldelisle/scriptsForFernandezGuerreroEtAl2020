#!/bin/bash

#SBATCH -o slurm-%A_%2a.out
#SBATCH -e slurm-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 4G
#SBATCH --cpus-per-task 1
#SBATCH --time 1:00:00
#SBATCH --array=1-8
#SBATCH --job-name mergeRep
#SBATCH --chdir /scratch/ldelisle/Merged93


# This script normalize the coverage by the number of uniquely mapped reads and do the mean of replicates at each base.

path="$PWD/"
pathForTable="replicatesToMerge.txt"

pathForFasta="/home/ldelisle/genomes/fasta/"
genome=mm10

module purge
module load gcc/7.4.0 #required for bowtie2, samtools, star and bedtools
module load bedtools2/2.27.1


sampleName=`cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}'`
samples=`cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{split($2,a,",");for(j in a){print a[j]}}'`
n=`echo $samples | wc -w`
sample="${sampleName}_neq$n"

mkdir -p ${path}/${sample}

pathResults=${path}/${sample}/
echo $sample
cd $pathResults
allBDGfwd=""
allBDGrev=""
for s in $samples; do
  # Get the total of uniquely mapped reads in the log of STAR
  tot=`cat ../allFinalFiles/reports/${s}_STAR_logFinal.txt | grep "Uniquely mapped reads number" | awk '{print $NF}'`
  # Normalize the bedgraph to the million uniquely mapped reads
  zcat ../bedGraphs/${s}_Uniq_fwd.bedGraph.gz | awk -v t=$tot -v OFS="\t" '{$4=$4/t*1e6;print $0 > "tmp_"$1}'
  cat tmp_chr* > normalized_fwd_${s}.bdg
  rm tmp_*
  zcat ../bedGraphs/${s}_Uniq_rev.bedGraph.gz | awk -v t=$tot -v OFS="\t" '{$4=$4/t*1e6;print $0 > "tmp_"$1}'
  cat tmp_chr* > normalized_rev_${s}.bdg
  rm tmp_*
  allBDGfwd="$allBDGfwd normalized_fwd_${s}.bdg"
  allBDGrev="$allBDGrev normalized_rev_${s}.bdg"
done
if [ $n -ne 1 ]; then
  # When there are more than 1 bedgraph
  # unionbedg is used to cut the genome in blocks
  # where all bedgraph have constant values
  # Then awk is used to perform the mean
  bedtools unionbedg -i $allBDGfwd | awk -v n=$sample -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\"mean of normalized Uniq reads forward of "n"\" visibility=full autoScale=on alwaysZero=on windowingFunction=mean"}{n=NF-3;sum=0;for(i=4;i<=NF;i++){sum+=$i};if(sum!=0){print $1,$2,$3,sum/n}}' > normalized_${sample}_Uniq_fwd.bedgraph
  bedtools unionbedg -i $allBDGrev | awk -v n=$sample -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\"mean of normalized Uniq reads reverse of "n"\" visibility=full autoScale=on alwaysZero=on windowingFunction=mean"}{n=NF-3;sum=0;for(i=4;i<=NF;i++){sum+=$i};if(sum!=0){print $1,$2,$3,sum/n}}' > normalized_${sample}_Uniq_rev.bedgraph
else
  cat $allBDGfwd | awk -v n=$sample -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\"mean of normalized Uniq reads forward of "n"\" visibility=full autoScale=on alwaysZero=on windowingFunction=mean"}$4!=0{print}' > normalized_${sample}_Uniq_fwd.bedgraph
  cat $allBDGrev | awk -v n=$sample -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\"mean of normalized Uniq reads reverse of "n"\" visibility=full autoScale=on alwaysZero=on windowingFunction=mean"}$4!=0{print}' > normalized_${sample}_Uniq_rev.bedgraph
fi

if [ ! -e ${pathForFasta}${genome}.fa.fai ]; then
  samtools faidx ${pathForFasta}${genome}.fa
fi

for or in fwd rev; do
  bedGraphToBigWig normalized_${sample}_Uniq_${or}.bedgraph ${pathForFasta}${genome}.fa.fai normalized_${sample}_Uniq_${or}.bw
  gzip normalized_${sample}_Uniq_${or}.bedgraph
done
