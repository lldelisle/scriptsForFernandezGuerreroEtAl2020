options(stringsAsFactors=F)
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("rtracklayer")
################################################################################
UIinput <- function(s){

    #Ask for user input
    x <- readLines(con = connection,1)
    
    #Check if empty
    if (x==""){
        x<-s
    }
    #Return
    return(x)
}
##################################################################
if (length(commandArgs(TRUE)) > 0){
  if (commandArgs(TRUE)[1] == "-h" || commandArgs(TRUE)[1] == "--help"){
    cat("Usage: Rscript filterGtfLikeAnouk.R pathForGtf [monoexonicBiotypes maxlen UCSCformat OnlyExons]\n")
    stop()
  }
  connection <- "stdin"
  gtfFile <- commandArgs(TRUE)[1]
  if (length(commandArgs(TRUE)) < 2){
    monoexonicBiotypes <- "Mt_tRNA,Mt_rRNA,IG_D_gene,IG_J_gene,snoRNA,misc_RNA,miRNA,snRNA,rRNA"
    maxlen <- 2500000 ## maximum transcript length
    UCSCformat <- TRUE
    keepOnlyExons <- TRUE
  } else {
    monoexonicBiotypes <- commandArgs(TRUE)[2]
    maxlen <- as.numeric(commandArgs(TRUE)[3])
    UCSCformat <- as.logical(commandArgs(TRUE)[4])
    keepOnlyExons <- as.logical(commandArgs(TRUE)[5])
  }
} else {
  connection <- stdin()
  cat("Choose the gtf file to filter.\n")
  gtfFile <- file.choose()
  cat("Put the monoexonic biotypes separated by comma or put enter to use Mt_tRNA,Mt_rRNA,IG_D_gene,IG_J_gene,snoRNA,misc_RNA,miRNA,snRNA,rRNA :")
  monoexonicBiotypes <- UIinput("Mt_tRNA,Mt_rRNA,IG_D_gene,IG_J_gene,snoRNA,misc_RNA,miRNA,snRNA,rRNA")
  cat("Put the maximum length of transcript which will be kept or put enter to use 2500000 :")
  maxlen <- as.numeric(UIinput(2500000))
  cat("Do you want to use UCSC format ?(Y)")
  ans <- UIinput("Y")
  if(tolower(ans) == "y"){
    UCSCformat <- TRUE
  } else {
    UCSCformat <- FALSE
  }
  cat("Do you want to keep only exons ?(Y)")
  ans <- UIinput("Y")
  if(tolower(ans) == "y"){
    keepOnlyExons <- TRUE
  } else {
    keepOnlyExons <- FALSE
  }
}
##################################################
cat("Loading gtf file...")
gtfInput <- readGFF(gtfFile)
cat("Done.\n")
##################################################
####  PART1 IDENTIFY READTHROUGH ###################
##################################################
gtf <- subset(gtfInput,type=="exon")
grgtf <- makeGRangesFromDataFrame(gtf)
hits <- findOverlaps(grgtf, grgtf)
rm(grgtf)
hitsDF <- as.data.frame(hits)
#Using monotype
monotype <- strsplit(monoexonicBiotypes,",")[[1]]
monogenes <- gtf$gene_id[gtf$gene_biotype %in% monotype]
#Add information
hitsDF$transcriptQ <- gtf$transcript_id[hitsDF$queryHits]
hitsDF$transcriptS <- gtf$transcript_id[hitsDF$subjectHits]
hitsDF$geneQ <- gtf$gene_id[hitsDF$queryHits]
hitsDF$geneS <- gtf$gene_id[hitsDF$subjectHits]
hitsDF$isMonoGene <- hitsDF$geneQ %in% monogenes | hitsDF$geneS %in% monogenes
#Filter the hits so 2 exons are not overlapping if they are from the same gene or if one of them is part of a monogene
hitsDF_bigFilter <- subset(hitsDF, geneS != geneQ & !isMonoGene)
#Count for each transcript the number of genes which are overlapped.
transcript_nbGene <- aggregate(hitsDF_bigFilter$geneS, by=list(hitsDF_bigFilter$transcriptQ), FUN=function(x){length(unique(x))})
rownames(transcript_nbGene) <- transcript_nbGene$Group.1
hitsDF_bigFilter$nbGeneForS <- transcript_nbGene[hitsDF_bigFilter$transcriptS, "x"]
#Filter the hits so we keep only the overlaps with transcripts which overlaps less than 2 genes.
#Because we consider that a read-through transcript is a transcript which overlaps 2 transcripts or more which themselves are not overlapping more than 2 different genes and these 2 transcripts or more need to belong to at least 2 different genes. 
hitsDF_bigFilter2 <- subset(hitsDF_bigFilter, nbGeneForS < 2)
#Build the table of genes overlapped by transcripts
transcript_gene <- aggregate(hitsDF_bigFilter2$geneS, by=list(hitsDF_bigFilter2$transcriptQ), FUN=function(x){paste(unique(x), collapse=",")})
#Put in the final Table only the transcripts which overlap 2 genes or more
finalTable <- subset(transcript_gene, grepl(",", x))
colnames(finalTable) <- c("TranscriptID", "OtherGenes")
readThroughTranscripts <- finalTable$TranscriptID
cat("There are", length(readThroughTranscripts), "read through transcripts.\n")
#If you want to make a report:
#Add info
#Build the table of transcripts overlapped by transcripts
transcript_transcript <- aggregate(hitsDF_bigFilter2$transcriptS, by=list(hitsDF_bigFilter2$transcriptQ), FUN=function(x){paste(unique(x), collapse=",")})
finalTable$OtherTranscripts <- transcript_transcript[match(finalTable$TranscriptID, transcript_transcript$Group.1), "x"]
geneFromTranscript <- unique(gtf[, c("transcript_id", "gene_id")])
finalTable$GeneID <- geneFromTranscript[match(finalTable$TranscriptID, geneFromTranscript$transcript_id), "gene_id"]
finalTable <- finalTable[, c("TranscriptID", "GeneID", "OtherTranscripts", "OtherGenes")]
name <- gsub(".gtf", "", basename(gtfFile))
name <- gsub(".gz", "", name)
write.table(finalTable, paste0(dirname(gtfFile), "/ReadthroughTranscriptsOf", name, ".txt"), sep="\t", quote=F, row.names=F)
####################################################
#### PART2 IDENTIFY LONG TRANSCRIPTS ###############
####################################################
if ("transcript" %in% gtfInput$type){
  transcripts <- subset(gtfInput, type == "transcript")
} else {
  transcripts <- unique(gtfInput[, c("gene_id", "gene_biotype", "transcript_id", "transcript_biotype")])
  if (sum(duplicated(transcripts$transcript_id)) > 0){
    print(transcripts[duplicated(transcripts$transcript_id), ])
    stop("The gtf file contains different information for the same transcript.")
  }
  start <- aggregate(gtfInput$start, by=list(transcript_id=gtfInput$transcript_id), FUN=min)
  transcripts$start <- start[match(transcripts$transcript_id,start$transcript_id), "x"]
  end <- aggregate(gtfInput$end, by=list(transcript_id=gtfInput$transcript_id), FUN=max)
  transcripts$end <- end[match(transcripts$transcript_id,end$transcript_id), "x"]
}
transcripts$length <- transcripts$end - transcripts$start + 1
tooLongTranscripts <- transcripts$transcript_id[transcripts$length > maxlen]
cat("There are", length(tooLongTranscripts), "transcripts which will be discarded because over", maxlen, "bp.\n")
#################################################################
#### PART3 IDENTIFY nonpc TRANSCRIPTS in pc genes ###############
#################################################################
pcFilteredTranscripts <- transcripts$transcript_id[transcripts$gene_biotype == "protein_coding" & transcripts$transcript_biotype != "protein_coding"]
cat("There are ", length(pcFilteredTranscripts), " transcripts which are not protein_coding whereas they are part of a gene which is protein_coding.\n")
############################################
#### PART4 write the new gtf ###############
############################################
newGtf <- gtfInput[! gtfInput$transcript_id %in% c(readThroughTranscripts, tooLongTranscripts, pcFilteredTranscripts), ]
if (keepOnlyExons){
  nBefore <- nrow(newGtf)
  newGtf <- newGtf[newGtf[, 3] == "exon", ]
  nAfter <- nrow(newGtf)
  if (nAfter!=nBefore){
    name <- paste0(name, "_ExonsOnly")
  }
}
if (UCSCformat){
  if (!grepl("chr", newGtf[1, 1])){
    #I filter only for chromosome with numbers or XYM
    newGtf <- newGtf[newGtf$seqid %in% c(1:25, "X", "Y", "MT"),]
    newGtf$seqid <- paste0("chr", newGtf$seqid)
    newGtf$seqid[newGtf$seqid == "chrMT"] <- "chrM"
    name <- paste0(name, "_UCSC")
  }
}
cat("Writing the filtered gtf...")
rownames(newGtf) <- NULL
export.gff(newGtf, paste0(dirname(gtfFile), "/FilteredTranscriptsOf", name, ".gtf"),format="gtf")
cat("Done.\n")
