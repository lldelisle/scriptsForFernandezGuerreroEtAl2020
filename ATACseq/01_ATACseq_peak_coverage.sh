#!/bin/bash

#SBATCH -o slurm-%A_%2a.out
#SBATCH -e slurm-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 4
#SBATCH --time 8:00:00
#SBATCH --array=1-5
#SBATCH --job-name ATACseq
#SBATCH --chdir /scratch/ldelisle/ATAC/

# This script run cutadapt with Nextera PE adapters
# alignment with Bowtie2 --very-sensitive (end-to-end) --no-discordant --dovetail -X 1000
# remove duplicates with picard
# convert bam to bed then macs2 for coverage and peak calling with --keep-dup all --shift -100 --extsize 200
# Normalize the coverage by the million of reads in peaks (summit +-500bp)

path="$PWD/"
pathForTable="table_ATAC.txt"
pathForFastq="${path}/fastq/"
pathForIndex="/home/ldelisle/genomes/bowtie2/"
pathForFasta="/home/ldelisle/genomes/fasta/"
genome=mm10

module purge
module load gcc/7.4.0 # required for bowtie2, samtools, star and bedtools
module load bowtie2/2.3.5
module load samtools/1.9 
module load picard/2.19.0
module load bedtools2/2.27.1 
# cutadapt is version 1.6 working with python 3.6.1 built with intel 17.0.2
# macs2 is version 2.1.1.20160309 working with python 2.7.5 built with GCC 4.8.5

nbOfThreads=4


sample=`cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}'`
genome=`cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $5}'`
fastqR1File=`cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $3}'`
fastqR2File=`cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $4}'`

indexPath=${pathForIndex}${genome}
pathResults=${path}/${sample}/

mkdir -p $pathResults

echo $sample
cd $pathResults

if [ ! -e ${path}allFinalFiles/reports/${sample}_report-cutadapt_PE.txt ]; then
  fastqR1=${pathForFastq}/$fastqR1File
  fastqR2=${pathForFastq}/$fastqR2File
  cutadapt -j $nbOfThreads -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -A CTGTCTCTTATACACATCTGACGCTGCCGACGA -q 30 -m 15 -o ${pathResults}${sample}-cutadapt_R1.fastq.gz -p ${pathResults}${sample}-cutadapt_R2.fastq.gz $fastqR1 $fastqR2 > ${pathResults}${sample}_report-cutadapt_PE.txt
  mkdir -p ${path}allFinalFiles/reports/
  cp ${pathResults}${sample}_report-cutadapt_PE.txt ${path}allFinalFiles/reports/
fi


if [ ! -e ${path}allFinalFiles/reports/${sample}_mapping_stats.txt ];then
  bowtie2 -p $nbOfThreads -x $indexPath --very-sensitive --no-unal --no-mixed --no-discordant --dovetail -X 1000 -1 ${sample}-cutadapt_R1.fastq.gz -2 ${sample}-cutadapt_R2.fastq.gz 2> ${pathResults}${sample}_mapping_stats.txt  | samtools view --threads $nbOfThreads -Su - | samtools sort --threads $nbOfThreads -o ${pathResults}${sample}_mapped_sorted.bam
  mkdir -p ${path}allFinalFiles/reports/
  cp ${pathResults}${sample}_mapping_stats.txt ${path}allFinalFiles/reports/
fi

if [ ! -e ${pathForFasta}${genome}.fa.fai ]; then
  samtools faidx ${pathForFasta}${genome}.fa
fi

if [ ! -e ${pathResults}${sample}_mapped_sorted_cp_q30_noM.bam ]; then
  allChrsExceptMito=`cat ${pathForFasta}${genome}.fa.fai | cut -f 1 | grep -v "chrM"`
  samtools index ${pathResults}${sample}_mapped_sorted.bam
  samtools view --threads $nbOfThreads -b ${pathResults}${sample}_mapped_sorted.bam -q 30 -f 0x2 $allChrsExceptMito > ${pathResults}${sample}_mapped_sorted_cp_q30_noM.bam
fi

mkdir -p ${path}bam
if [ ! -e ${path}bam/${sample}_mapped_sorted_cp_q30_noM_rmdup.bam ]; then
  picard MarkDuplicates SORTING_COLLECTION_SIZE_RATIO=0.15 I=${pathResults}${sample}_mapped_sorted_cp_q30_noM.bam O=${pathResults}${sample}_mapped_sorted_cp_q30_noM_rmdup.bam M=${pathResults}${sample}_mapped_sorted_cp_q30_noM_rmdup.log REMOVE_DUPLICATES=true AS=true
  cp ${pathResults}${sample}_mapped_sorted_cp_q30_noM_rmdup.log ${path}allFinalFiles/reports/
  cp ${pathResults}${sample}_mapped_sorted_cp_q30_noM_rmdup.bam ${path}bam/
fi

mkdir -p ${path}bedGraphs
mkdir -p ${path}bw
mkdir -p ${path}toGEO

if [ ! -e ${path}bedGraphs/${sample}_macs_likeATAC.bedGraph.gz ]; then
  bedtools bamtobed -i ${pathResults}${sample}_mapped_sorted_cp_q30_noM_rmdup.bam > ${pathResults}${sample}_mapped_sorted_cp_q30_noM_rmdup.bed
  macs2 callpeak -t ${pathResults}${sample}_mapped_sorted_cp_q30_noM_rmdup.bed --nomodel --keep-dup all --shift -100 --extsize 200 -n ${sample} --call-summits -B 2> ${pathResults}${sample}_macs_likeATAC.log
  cp ${pathResults}${sample}_macs_likeATAC.log ${path}allFinalFiles/reports/
  cp ${sample}_peaks.narrowPeak ${path}toGEO/${sample}.narrowPeak
  # The output of macs2 can go over the size of chromosomes so a special care needs to be taken before coverting to bigwig
  bash ${path}/fromMacs2BdgToSimplifiedBdgAndBw.sh ${pathResults}${sample}_treat_pileup.bdg ${path}bedGraphs/${sample}_macs_likeATAC "macs2 like ATAC of ${sample}" ${pathForFasta}${genome}.fa.fai &
fi

wait

if [ ! -e ${path}toGEO/${sample}.bw ]; then
  cp ${path}bedGraphs/${sample}_macs_likeATAC.bw ${path}toGEO/${sample}.bw
fi

if [ ! -e ${path}bedGraphs/${sample}_macs_likeATAC_norm2.bedGraph.gz ]; then
  if [ ! -e ${pathResults}${sample}_readsInPeaks.txt ]; then
    bedtools slop -i ${pathResults}${sample}_summits.bed -g ${pathForFasta}${genome}.fa.fai -b 500 > ${pathResults}${sample}_macs_likeATAC_summits_1kb.bed
    bedtools merge -i ${pathResults}${sample}_macs_likeATAC_summits_1kb.bed > ${pathResults}${sample}_macs_likeATAC_summits_1kb_merged.bed
    bedtools coverage -a ${pathResults}${sample}_macs_likeATAC_summits_1kb_merged.bed -b ${pathResults}${sample}_mapped_sorted_cp_q30_noM_rmdup.bed > ${pathResults}${sample}_macs_likeATAC_summits_1kb_merged_coverage.txt
    cat ${pathResults}${sample}_macs_likeATAC_summits_1kb_merged_coverage.txt | awk '{S+=$4}END{print S}' > ${pathResults}${sample}_readsInPeaks.txt
  fi
  readsInPeaks=`cat ${pathResults}${sample}_readsInPeaks.txt`
  zcat ${path}bedGraphs/${sample}_macs_likeATAC.bedGraph.gz | awk -v s=$sample -v n=$readsInPeaks -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\""s" like ATAC normalized by million reads in peaks\" visibility=full autoScale=on windowingFunction=mean"}NR>1{$4=$4/n*1e6; print}' > ${path}bedGraphs/${sample}_macs_likeATAC_norm2.bedGraph
  bedGraphToBigWig ${path}bedGraphs/${sample}_macs_likeATAC_norm2.bedGraph ${pathForFasta}${genome}.fa.fai ${path}bedGraphs/${sample}_macs_likeATAC_norm2.bw
  gzip ${path}bedGraphs/${sample}_macs_likeATAC_norm2.bedGraph &
fi
wait
