if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("BiocManager")
safelyLoadAPackageInCRANorBioconductor("doSNOW")
safelyLoadAPackageInCRANorBioconductor("plot3D")
safelyLoadAPackageInCRANorBioconductor("ggplot2")
safelyLoadAPackageInCRANorBioconductor("umap")
safelyLoadAPackageInCRANorBioconductor("gridExtra")
devtools::install_github("r3fang/SnapATAC")
library(SnapATAC)

# Load the RData from GEO
load("GSE145657_final_13_02_2020.Rdata")

# Create a snap with the snap file from GEO and the peak just created
my.snap <- createSnap(file="GSM4323465_possorted_bamwt.snap", sample="WT")
my.snap <- addPmatToSnap(my.snap)

# I write the result table to put it on GitHub
df <- as.data.frame(as.matrix(my.snap@pmat))
colnames(df) <- c("EC1", "EC2")
df$barcode <- rownames(df)
write.table(df[, c("barcode", "EC1", "EC2")], "EC_coverage.txt", sep="\t", row.names=F, quote=F)

# Create a data frame with both metadata, cluster and UMAP coordinates
meta.data <- x.sp@metaData
meta.data$cluster <- x.sp@cluster
meta.data$UMAP1 <- x.sp@umap[, 1]
meta.data$UMAP2 <- x.sp@umap[, 2]

# Consider only wt cells
meta.data.wt <- subset(meta.data, sample == "WT")

# Report the coverage for EC1 and EC2
meta.data.wt$EC1 <- my.snap@pmat[meta.data.wt$barcode, 1]
meta.data.wt$EC2 <- my.snap@pmat[meta.data.wt$barcode, 2]

# Make individual plots to get the scale 
g <- ggplot(data=meta.data.wt[order(meta.data.wt$EC1),], aes(x=UMAP1, y=UMAP2, col=EC1))  + 
  geom_point(size = 0.2) +
  scale_colour_gradient(low = "white", high = rgb(1, 0 ,1))

ggsave("UMAP_EC1.pdf", g)

g <- ggplot(data=meta.data.wt[order(meta.data.wt$EC2),], aes(x=UMAP1, y=UMAP2, col=EC2)) + 
  geom_point(size = 0.2) +
  scale_colour_gradient(low = "white", high = rgb(0, 1 ,1))

ggsave("UMAP_EC2.pdf", g)

# Scale EC1 and EC2
meta.data.wt$EC1_scaled <- 1 - meta.data.wt$EC1 / max(meta.data.wt$EC1)
meta.data.wt$EC2_scaled <- 1 - meta.data.wt$EC2 / max(meta.data.wt$EC2)

# Convert to rgb
meta.data.wt$rgb <- apply(meta.data.wt[, c("EC1_scaled", "EC2_scaled")], 1, function(v){rgb(v[2], v[1], 1)})

# Plot the UMAP with the color representing both EC1 and EC2 and sorting the cells with more coverage on top
g <- ggplot(data=meta.data.wt[order(meta.data.wt$EC2 + meta.data.wt$EC1),], aes(x=UMAP1, y=UMAP2, col=rgb)) + 
  geom_point(size = 0.2) +
  scale_color_identity()
ggsave("umap_EC_sorted_smallerPoints.pdf", g)

# Generate a scale:
scale.data <- unique(meta.data.wt[, c("EC1", "EC2", "rgb")])

# Report the number of cell in each combination for cluster 8 (ectoderm)
t.8 <- with(subset(meta.data.wt, cluster == 8), table(EC1, EC2))
scale.data$nb.8 <- apply(scale.data, 1, function(v){t.8[v[1], v[2]]})
scale.data$p.8 <- scale.data$nb.8 / sum(scale.data$nb.8)

# Or non 8 (non-ectoderm)
t.n8 <- with(subset(meta.data.wt, cluster != 8), table(EC1, EC2))
scale.data$nb.n8 <- apply(scale.data, 1, function(v){tryCatch(t.n8[v[1], v[2]], error=function(e){0})})
scale.data$p.n8 <- scale.data$nb.n8 / sum(scale.data$nb.n8)

# Plot the scale for ectoderm
g <- ggplot(data=subset(scale.data, nb.8 > 0), aes(x=EC1, y=EC2, label=nb.8, col=rgb)) + 
  geom_point(aes(size = p.8)) +
  geom_text(col="black") +
  theme(text = element_text(size=20),
        axis.text = element_text(size=20)) +
  scale_size_area(max_size = 25, limits = c(0, 1)) + 
  labs(size = "Proportion") +
  scale_color_identity()
ggsave("umap_EC_scale_ECT.pdf", g)

# Plot the scale for non-ectoderm
g <- ggplot(data=subset(scale.data, nb.n8 > 0), aes(x=EC1, y=EC2, label=nb.n8, col = rgb)) + 
  geom_point(aes(size = p.n8)) +
  geom_text(col="black") +
  theme(text = element_text(size=20),
        axis.text = element_text(size=20)) +
  xlim(0, max(scale.data$EC1)) +
  ylim(0, max(scale.data$EC2)) +
  scale_size_area(max_size = 25, limits = c(0, 1)) + 
  labs(size = "Proportion") +
  scale_color_identity()
ggsave("umap_EC_scale_nonECT.pdf", g)

# Use inkscape to put together
