#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --mem 30G
#SBATCH --cpus-per-task 16
#SBATCH --time 2:00:00
#SBATCH --array=1-16
#SBATCH --job-name RNAseqMerged93
#SBATCH --chdir /scratch/ldelisle/Merged93

# This script remove adapters and bad quality with cutadapt
# make alignment with STAR ENCODE parameters
# Evaluate FPKM with cufflinks
# Coverage on uniquely mapped reads with bedtools

path='/scratch/ldelisle/Merged93/'
pathForFastq='/scratch/ldelisle/Merged93/newFastq/'
pathForTable="table.txt"
pathForIndex='/scratch/ldelisle/STARIndex/'
pathForFasta='/home/ldelisle/genomes/fasta/'
genome=mm10
versionOfSTAR="2.6.0c"
versionOfCufflinks="2.2.1"

module purge
module load gcc/7.4.0 #required for bowtie2, samtools and bedtools
module load bowtie2/2.3.5
module load samtools/1.9
module load bedtools2/2.27.1
# cutadapt is version 1.6 working with python 3.6.1 built with intel 17.0.2

indexPath=${pathForIndex}${genome}
gtfFile=${path}mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm38.93_ExonsOnly_UCSC.gtf

# Create a gtf file with a gene on chrM for each direction to filter them in cufflinks
if [ "$SLURM_ARRAY_TASK_ID" = 1 ]; then
  echo "chrM	chrM_gene	exon	0	16299	.	+	.	gene_id \"chrM_gene_plus\"; transcript_id \"chrM_tx_plus\"; exon_id \"chrM_ex_plus\";
chrM	chrM_gene	exon	0	16299	.	-	.	gene_id \"chrM_gene_minus\"; transcript_id \"chrM_tx_minus\"; exon_id \"chrM_ex_minus\";" > MTmouse.gtf
fi

sample=`cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}'`
fastqFile=`cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}'`

mkdir -p ${path}STAR/$sample

pathResults=${path}STAR/${sample}/
echo $sample
cd $pathResults

# Remove TruSeq adapter, bases with quality below 30 and only write reads longer than 15bp.
if [ ! -e ${pathResults}${sample}-cutadapt.fastq.gz ]; then
  cutadapt -m 15 -j 16 -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -o ${pathResults}${sample}-cutadapt.fastq.gz ${pathForFastq}${fastqFile} > ${pathResults}${sample}_report-cutadapt.txt
  mkdir -p ${path}allFinalFiles/reports/
  cp ${pathResults}${sample}_report-cutadapt.txt ${path}allFinalFiles/reports/
fi

# Map with STAR using the ENCODE parameters
if [ ! -e ${path}allFinalFiles/bam/${sample}_Aligned.sortedByCoord.out.bam ];then
  export PATH=$PATH:/home/ldelisle/softwares/STAR-$versionOfSTAR/bin/Linux_x86_64
  STAR --runThreadN 16 --genomeDir ${indexPath} --readFilesIn ${pathResults}${sample}-cutadapt.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate  --sjdbOverhang '99' --sjdbGTFfile $gtfFile  --quantMode GeneCounts  --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1
  mkdir -p ${path}allFinalFiles/reports/
  cp Log.final.out ${path}allFinalFiles/reports/${sample}_STAR_logFinal.txt
  mkdir -p ${path}allFinalFiles/bam/
  cp ${pathResults}Aligned.sortedByCoord.out.bam ${path}allFinalFiles/bam/${sample}_Aligned.sortedByCoord.out.bam
fi

# To allow parallelisation, cufflinks step is written in a bash script
echo "export PATH=$PATH:/home/ldelisle/softwares/cufflinks-${versionOfCufflinks}.Linux_x86_64" >cufflinks_${sample}.sh
echo "mkdir -p ${pathResults}cufflinksWOMT" >>cufflinks_${sample}.sh
echo "cufflinks -p 10 -o ${pathResults}cufflinksWOMT --max-bundle-length 10000000 --multi-read-correct --library-type "fr-firststrand" -b ${pathForFasta}${genome}.fa  --no-effective-length-correction -M ${path}MTmouse.gtf -G $gtfFile ${pathResults}Aligned.sortedByCoord.out.bam" >>cufflinks_${sample}.sh
echo "" >>cufflinks_${sample}.sh
echo "mkdir -p ${path}allFinalFiles" >>cufflinks_${sample}.sh
echo "cp ${pathResults}cufflinksWOMT/genes.fpkm_tracking ${path}allFinalFiles/FPKM_${sample}.txt" >>cufflinks_${sample}.sh
echo "cp ${pathResults}cufflinksWOMT/isoforms.fpkm_tracking ${path}allFinalFiles/FPKM_${sample}_isoforms.txt" >>cufflinks_${sample}.sh

if [ -e ${pathForFasta}${genome}.fa ] && [ -e ${path}MTmouse.gtf ] && [ -e $gtfFile ] && [ -s ${pathResults}Aligned.sortedByCoord.out.bam ];then
  if [ ! -e ${path}allFinalFiles/FPKM_${sample}_isoforms.txt ];then
    echo "Launching cufflinks"
    bash cufflinks_${sample}.sh &
  fi
else
  echo "cufflinks not launch because some files are missing."
fi

# Uniquely mapped reads are filtered using the NH:i:1 tag
if { [ ! -e accepted_hits_unique_${sample}.bam ] && [ -s ${pathResults}Aligned.sortedByCoord.out.bam ] ;} || [ -e tmp.header ] ;then
  echo "Compute uniquely aligned"
  samtools view -H Aligned.sortedByCoord.out.bam >  tmp.header
  samtools view -@ 5 Aligned.sortedByCoord.out.bam | grep  -w "NH:i:1" | cat tmp.header - | samtools view -@ 5 -b > accepted_hits_unique_${sample}.bam
  rm tmp.header
fi

# This is to make the bedgraph of coverage
mkdir -p ${path}bedGraphs

# Bedgraph coverage is obtained from the uniquely mapped reads with bedtools each strand separately
if { [ ! -e ${path}bedGraphs/${sample}_Uniq_fwd.bedGraph.gz ] && [ -e accepted_hits_unique_${sample}.bam ];} || [ -e tmp.header.uf ] ;then
  #strand + corresponds to reverse strand due to TruSeq
  echo "Building uniq fwd reads bedGraph"
  echo "track type=bedGraph name=\"${sample} Uniq reads forward\"">tmp.header.uf
  echo "bedtools genomecov -ibam accepted_hits_unique_${sample}.bam -bg -split -strand - | cat tmp.header.uf - | gzip > ${path}bedGraphs/${sample}_Uniq_fwd.bedGraph.gz">uf.sh
  echo "rm tmp.header.uf">>uf.sh
  bash uf.sh &
fi
if { [ ! -e ${path}bedGraphs/${sample}_Uniq_rev.bedGraph.gz ] && [ -e accepted_hits_unique_${sample}.bam ];} || [ -e tmp.header.ur ];then
  echo "Building uniq rev reads bedGraph"
  echo "track type=bedGraph name=\"${sample} Uniq reads reverse\"">tmp.header.ur
  echo "bedtools genomecov -ibam accepted_hits_unique_${sample}.bam -bg -split -strand + | cat tmp.header.ur - | gzip > ${path}bedGraphs/${sample}_Uniq_rev.bedGraph.gz">ur.sh
  echo "rm tmp.header.ur">>ur.sh
  bash ur.sh &
fi

wait
echo "Everything is done"
find . -size 0 -delete

