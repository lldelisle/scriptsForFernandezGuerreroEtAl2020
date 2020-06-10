#!/bin/bash

#SBATCH -o slurm-%A_%2a.out
#SBATCH -e slurm-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 4G
#SBATCH --cpus-per-task 1
#SBATCH --time 1:00:00
#SBATCH --array=1-2
#SBATCH --job-name mergeRep
#SBATCH --chdir /scratch/ldelisle/ATAC


# This script do the mean of coverage between replicates

path="$PWD/"
pathForTable="replicatesToMerge.txt"

pathForFasta="/home/ldelisle/genomes/fasta/"
genome=mm10

module purge
module load gcc/7.4.0 # required for bowtie2, samtools, star and bedtools
module load bedtools2/2.27.1


sampleName=`cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}'`
samples=`cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{split($2,a,",");for(j in a){print a[j]}}'`
n=`echo $samples | wc -w`
sample="${sampleName}_neq$n"

mkdir -p ${path}/${sample}

pathResults=${path}/${sample}/
echo $sample
cd $pathResults
allBDG=""
for s in $samples; do
  zcat ${path}bedGraphs/${s}_macs_likeATAC_norm2.bedGraph.gz > normalized_${s}.bdg
  allBDG="$allBDG normalized_${s}.bdg"
done
if [ $n -ne 1 ]; then
  # When there are more than 1 bedgraph
  # unionbedg is used to cut the genome in blocks
  # where all bedgraph have constant values
  # Then awk is used to perform the mean
  bedtools unionbedg -i $allBDG | awk -v n=$sample -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\"mean of norm2 macs2 likeATAC of "n"\" visibility=full autoScale=on alwaysZero=on windowingFunction=mean"}{n=NF-3;sum=0;for(i=4;i<=NF;i++){sum+=$i};if(sum!=0){print $1,$2,$3,sum/n}}' > normalized_${sample}.bedgraph
else
  cat $allBDG | awk -v n=$sample -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\"mean of norm2 macs2 likeATAC of "n"\" visibility=full autoScale=on alwaysZero=on windowingFunction=mean"}$4!=0{print}' > normalized_${sample}.bedgraph
fi
bedGraphToBigWig normalized_${sample}.bedgraph ${pathForFasta}${genome}.fa.fai normalized_${sample}.bw
gzip normalized_${sample}.bedgraph
