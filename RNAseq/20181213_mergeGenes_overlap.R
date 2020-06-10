options(stringsAsFactors=F)
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("rtracklayer")
if(length(commandArgs(TRUE))>0){
  if(commandArgs(TRUE)[1]=="-h" || commandArgs(TRUE)[1]=="--help"){
    cat("Usage: Rscript 20181213_mergeGenes.R pathForGtf\nThis tools will create a new gtf where the same ensembl_id is given to genes which have the same gene_name and overlap in order to avoid too many ambiguous reads in counts per gene.")
    stop()
  }
  gtfFile<-commandArgs(TRUE)[1]
} else {
  cat("Choose the gtf file to filter.\n")
  gtfFile<-file.choose()
}

cat("Loading gtf file...")
gtfInput <- readGFF(gtfFile)
cat("Done.\n")
# simpleGTF is only uniques gene_id and gene_name
simpleGTF <- unique(gtfInput[,c("gene_id","gene_name")])
# Genes which are 2 different ensembl_id
duplicatedGeneNames <-unique(with(simpleGTF,gene_name[duplicated(gene_name)]))
newGTF <- gtfInput
change <- data.frame(gene_name=character(),before=character(),after=character())
# All attributes relative to gene should change.
colsToChange <- grep("gene",colnames(gtfInput))
for (gn in duplicatedGeneNames){
  while (T){
    # Merge all overlapping exons in each ensembl_id:
    smallgrl <-reduce(split(GRanges(newGTF[newGTF$type == "exon" & newGTF$gene_name == gn, ]), 
                            newGTF$gene_id[newGTF$type == "exon" & newGTF$gene_name == gn]))
    # Count the number of exon s for each ensembl_id which overlaps other ensembl_id (including itself)
    smallTable <- lapply(smallgrl, function(g){
      table(names(smallgrl)[as.matrix(findOverlaps(g, smallgrl))[, 2]])
    })
    overlapsExons <- unlist(unname(smallTable), use.names = T)
    # The ensembl_id present twice in the smallTable mean they are overlapping
    ensIDwithOverlap <- names(overlapsExons)[duplicated(names(overlapsExons))]
    # If there are not overlapping. We stop modifying this gene_name
    if (length(ensIDwithOverlap) == 0){
      break
    }
    # We will merge with the ensembl_id which has the more exons
    masterEns <- names(which.max(overlapsExons[names(overlapsExons) %in% ensIDwithOverlap]))
    # All ensembl_id which overlap with the master one will be fused.
    ensIDToRename <- setdiff(names(which(sapply(smallTable, function(t){masterEns %in% names(t)}))), masterEns)
    if (length(ensIDToRename) > 0){
      # I am not sure this happens but just in case
      # I store the old ens_id and the new ens_id in the dataframe change.
      change <- rbind(change, data.frame(gene_name=rep(gn, length(ensIDToRename)),
                                         before=ensIDToRename,
                                         after=rep(masterEns, length(ensIDToRename))))
      # I take the info on genes from the masterEns
      masterCols <- unique(newGTF[newGTF$gene_id == masterEns, colsToChange])
      if(nrow(masterCols) > 1){
        cat("There are inconstistancies in \n")
        print(masterCols)
        cat("Only the first line will be kept\n")
        masterCols<-materCols[1, ]
      }
      # I put them in the fused ensembl_ids
      newGTF[newGTF$gene_id %in% c(ensIDToRename, masterEns), colsToChange] <- masterCols
    } else {
      cat("It happens for",gn,"\n")
      break
    }
  }
}
name <- gsub(".gtf", "", basename(gtfFile))
name <- gsub(".gz", "", name)
write.table(change, paste0(dirname(gtfFile), "/renamedEnsIDOf", name, ".txt"), sep="\t", quote=F, row.names=F)
cat("Writing the merged gtf...")
rownames(newGTF) <- NULL
export.gff(newGTF, paste0(dirname(gtfFile), "/mergeOverlapGenesOf", name, ".gtf"), format="gtf")
cat("Done.\n")
