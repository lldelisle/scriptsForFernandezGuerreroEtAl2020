options(stringsAsFactors=F)
if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("pheatmap")
safelyLoadAPackageInCRANorBioconductor("RColorBrewer")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
source("RNAseqFunctions.R")
# To have consistent colors in all heatmaps/PCA
colTissue <- c(rgb(252, 124, 118, maxColorValue=255),
               rgb(23, 183, 255, maxColorValue=255))
names(colTissue) <- c("ECT", "MES")
colReplicate<-c("mediumturquoise", "orchid")
names(colReplicate)<-c(1, 2)
cc<-colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(103)

FPKMfile <- "/scratch/ldelisle/Merged93/mergedTables/AllCufflinks_Simplified.txt"
nGenesToRestrictPCA <- 500
samplesPlanDF <- read.delim("samplesplan.txt", stringsAsFactors=T)
expressionDF <- read.delim(FPKMfile)
samplesToPlot <- as.character(samplesPlanDF$sample)
colnames(expressionDF)<-gsub("^FPKM_","",colnames(expressionDF))
expressionDF<-expressionDF[,c(colnames(expressionDF)[1:3],samplesToPlot)]

factorizedSP <- simplifyDF(samplesPlanDF, samplesToPlot, keepFactorsWithOneValue=T)
# For heatmaps:
annC <- data.frame(tissue = factorizedSP$tissue)
rownames(annC) <- rownames(factorizedSP)
anno_colors <- list(tissue=colTissue)

# Figure 1A:
genesToPlot<-intersect(paste0("Hox", rep(letters[1:4], each =13), rep(1:13, 4)), expressionDF$gene_short_name)
# Create a small dataframe with only expression for the genes to plot
df.sub <- expressionDF[expressionDF$gene_short_name %in% genesToPlot, rownames(factorizedSP)]
# Do a log2 transformation
df.sub <- log2(df.sub + 1)
# Put the rownames:
rownames(df.sub) <- expressionDF$gene_short_name[expressionDF$gene_short_name %in% genesToPlot]
# Keep the order
df.sub <- df.sub[genesToPlot, ]
# The scalebar is built so that minimum color is 0 and maximum color is 12 even if larger values exists
epsilon <- 0.0000001
breaksListAbs <- c(seq(0, 12, length.out=100), max(max(df.sub), 12) + epsilon)

pdf("F1A.pdf",
    title = paste0("Hox genes"),
    width = 10, height = 3 + 0.23 * length(genesToPlot))
pheatmap(df.sub, cluster_rows=FALSE, show_rownames=TRUE, breaks=breaksListAbs, color=cc,
         cluster_cols=FALSE, annotation_col=annC, fontsize=15, cellwidth=20, cellheight=16,
         annotation_colors=anno_colors, main="log2(FPKM+1)", labels_col = factorizedSP$stageRep,
         labels_row = gsub("Hox", "", rownames(df.sub)), angle_col = 90)

dev.off()

# Figure S1B:
genesToPlot<-c("Fgf8", "Sp6", "Trp63", "Perp", "Krt14", "Sp8", "Wnt7a", "Prrx1", "Tbx5", "Hoxa9", "Hoxd10", "Fgf10", "Pecam1", "Cdh5")
# Create a small dataframe with only expression for the genes to plot
df.sub <- expressionDF[expressionDF$gene_short_name %in% genesToPlot, rownames(factorizedSP)]
# Do a log2 transformation
df.sub <- log2(df.sub + 1)
# Put the rownames:
rownames(df.sub) <- expressionDF$gene_short_name[expressionDF$gene_short_name %in% genesToPlot]
# The scalebar is built so that minimum color is 0 and maximum color is 12 even if larger values exists
epsilon <- 0.0000001
breaksListAbs <- c(seq(0, 12, length.out=100), max(max(df.sub), 12) + epsilon)

pdf("S1B.pdf",
    title = paste0("ECT MES genes"),
    width = 10, height = 3 + 0.23 * length(genesToPlot))
pheatmap(df.sub, cluster_rows=TRUE, show_rownames=TRUE, breaks=breaksListAbs, color=cc,
         cluster_cols=FALSE, annotation_col=annC, fontsize=15, cellwidth=20, cellheight=16,
         annotation_colors=anno_colors, main="log2(FPKM+1)", labels_col = factorizedSP$stageRep,
         angle_col = 90)
dev.off()

# Figure S1C:
data <- expressionDF[,samplesToPlot]
sumperline <- apply(data,1,sum)
nonZdata <- data[sumperline != 0,]
ldata <- log2(nonZdata + 1)
rldata <- ldata
rldata <- ldata[order(apply(ldata, 1, var), decreasing = T)[1:min(nrow(ldata), nGenesToRestrictPCA)], ]
sample.pca <- prcomp(t(rldata),
                     center = TRUE,
                     scale. = FALSE)
new.df <- data.frame(factorizedSP,sample.pca$x[samplesToPlot,])
var <- round((sample.pca$sdev) ^ 2 / sum(sample.pca$sdev ^ 2) * 100)

pdf("S1C.pdf")
ggplot(new.df, aes(PC1, PC2)) +
  geom_point(aes(fill=tissue, color=Replicate, shape=stage), stroke=1, size=3) +
  theme_grey(base_size = 20) +
  xlab(paste0("PC1: ", var[1], "% variance")) +
  ylab(paste0("PC2: ", var[2], "% variance")) +
  scale_fill_manual(values=colTissue)+
  scale_color_manual(values=colReplicate)+ 
  scale_shape_manual(values=21:24) +
  guides(fill=guide_legend(override.aes = list(shape = 21)), alpha=guide_legend(override.aes = list(shape = 21)), color=guide_legend(override.aes = list(shape = 21)))
dev.off()

sampleDists <- dist(t(rldata),method="euclidean")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- samplesToPlot
colnames(sampleDistMatrix) <- samplesToPlot
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
annot <- data.frame(tissue = factorizedSP$tissue)
rownames(annot) <- rownames(factorizedSP)
clu <- hclust(sampleDists, method="ward.D2")
dd <- as.dendrogram(clu)
clu2 <- reorder(dd, 1:length(samplesToPlot), agglo.FUN = mean)
pdf("S1D.pdf",title="CorrelationMatrix",width = 7,height = 7)
pheatmap(
  sampleDistMatrix,
  labels_row = as.character(factorizedSP$stageRep),
  labels_col = factorizedSP$stageRep,
  cluster_rows = as.hclust(clu2),
  cluster_cols = as.hclust(clu2),
  cellwidth = 10,
  cellheight = 10,
  annotation_row = annot,
  annotation_col = annot,
  col = colors,
  main="Euclidean distance - ward clustering",
  annotation_colors = anno_colors,
  angle_col = 90
)
dev.off()

