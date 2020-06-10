##Generate colors ###
# To have consistent colors in all heatmaps
colTissue <- c(rgb(252, 124, 118, maxColorValue=255),
               rgb(23, 183, 255, maxColorValue=255))
names(colTissue) <- c("ECT", "MES")
colReplicate<-c("mediumturquoise", "orchid")
names(colReplicate)<-c(1, 2)
anno_colors <- list(tissue=colTissue)

cc<-colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(103)

### FUNCTIONS #####
plotFixedlogValues<-function(genes, typeOfID, input.df, samples,
                             colsToUse=c("Replicate", "stage", "tissue"),
                             returnDF=F, annR=NULL, keepOrder=F, scaleInHeatmap="none",
                             clusteringMethod="euclidean", main=NA, cellheight=16,
                             color=cc){
  # We expect to have a column called samples which is also the rownames
  if (!("sample" %in% colnames(samples))){
    samples$sample <- rownames(samples)
  }
  rownames(samples) <- samples$sample
  # The samples name should match the name of the columns of the input.df
  # Only the intersection will be plotted
  samplesToPlot <- intersect(samples$sample, colnames(input.df))
  # A new data frame is created with only the genes to plot and the samples to plot
  df.sub <- input.df[input.df[, typeOfID] %in% genes, samplesToPlot]
  # The data are transformed in log2(1+input)
  df.sub <- log2(df.sub + 1)
  # The gene names to write on the plots would be preferentially the gene_short_name (if in the input.df)
  potentialRowNames <- tryCatch(input.df$gene_short_name[input.df[, typeOfID] %in% genes],
                                error=function(e){input.df[input.df[, typeOfID] %in% genes, typeOfID]})
  # If there is no duplication in the gene names it will be used, else the typeOfID is used
  if (anyDuplicated(potentialRowNames) > 0){
    rownames(df.sub) <- input.df[input.df[, typeOfID] %in% genes, typeOfID]
  } else {
    rownames(df.sub) <- potentialRowNames
  }
  # If there is a dataframe for the annotation, the rownames are adjusted.
  if (!is.null(annR)){
    if (!all(rownames(df.sub) %in% rownames(df.sub))){
      rownames(annR) <- tryCatch(input.df$gene_short_name[match(rownames(annR), input.df[, typeOfID])],
                                 error=function(e){rownames(annR)})
    }
  }
  # If there is a annR, with only one level the colors of annR will be added to the localAnnoColor variable:
  localAnnoColor <- anno_colors
  for (cn in colnames(annR)){
    if (length(levels(annR[,cn])) == 1){
      ccn <- "black"
      names(ccn) <- levels(annR[, cn])
      ccn <- list(ccn)
      names(ccn) <- cn
      localAnnoColor <- c(localAnnoColor, ccn)
    }
  }
  # The scalebar is built so that minimum color is 0 and maximum color is 12 even if larger values exists
  epsilon <- 0.0000001
  breaksListAbs <- c(seq(0, 12, length.out=100), max(max(df.sub), 12) + epsilon)
  # If the heatmap is scaled by row then it goes from -1.5 to 1.5
  if (scaleInHeatmap == "row"){
    breaksListAbs <- c(min(apply(df.sub, 1, scale), - 1.5) - epsilon, 
                       seq(- 1.5, 1.5, 0.03),
                       max(apply(df.sub, 1, scale), 1.5) + epsilon)
  }
  if (keepOrder || nrow(df.sub) < 2){
    df.sub.ordered <- df.sub[genes, ]
    if (all(is.na(df.sub.ordered))){
      df.sub.ordered <- df.sub[input.df$gene_short_name[match(genes, input.df[, typeOfID])], ]
    }
    pheatmap(df.sub.ordered, cluster_rows=FALSE, show_rownames=TRUE, breaks=breaksListAbs, color=color,
            cluster_cols=FALSE, annotation_col=samples[, colsToUse], annotation_row=annR,
            fontsize=15, cellwidth=20, cellheight=cellheight, annotation_colors=localAnnoColor,
            scale=scaleInHeatmap, main=main)
  } else {
    pheatmap(df.sub, cluster_rows=TRUE, show_rownames=TRUE, breaks=breaksListAbs,color=color,
            cluster_cols=FALSE, annotation_col=samples[,colsToUse],annotation_row=annR,
            fontsize=15, cellwidth=20, cellheight=cellheight, annotation_colors=localAnnoColor,
            scale=scaleInHeatmap, clustering_distance_rows=clusteringMethod, main=main)
  }
  if(returnDF){
    return(df.sub)
  }
}
