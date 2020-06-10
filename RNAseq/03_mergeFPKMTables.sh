#!/bin/bash

#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lucille.delisle@epfl.ch
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 1G
#SBATCH --time 0:10:00
#SBATCH --job-name 01RNAseq
#SBATCH --chdir /scratch/ldelisle/Merged93/

#This script will merge the tables 
path='/scratch/ldelisle/Merged93/'
configFile="configFileRNAseq_step1.R"
module purge
module load gcc/7.4.0  openblas/0.3.6-openmp r/3.6.0

#Adjust the samplesplan.txt to add the cufflinks path
samplesPlanFile=`cat $configFile| awk -F '=|<-|#' '$1=="samplesPlan"{print $2}' | tr -d "\"" | tr -d " "`
if [ ! `grep "cufflinks_file" $samplesPlanFile` ]; then
  # The path of each cufflinks file (FPKM) is ${path}/allFinalFiles/FPKM_${sample}.txt
  cat $samplesPlanFile | awk -v pa=$path 'BEGIN{print "cufflinks_file"}NR>1{print pa"/allFinalFiles/FPKM_"$1".txt"}' > CuffCol.txt
  
  # The cufflinks column is added to the samplesPlan
  paste -d "\t" $samplesPlanFile CuffCol.txt > ${samplesPlanFile}_withPaths
  
  rm CuffCol.txt
fi
if [ -e ${samplesPlanFile}_withPaths ]; then
  # The config file is modified to use the new samplesPlan file
  cat $configFile | sed "s#$samplesPlanFile#${samplesPlanFile}_withPaths#" > ${configFile}_withPaths
  configFile=${configFile}_withPaths
fi

Rscript step1-generateTables.R $configFile
