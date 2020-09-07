if (!"devtools" %in% installed.packages()){
  install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")
}
devtools::install_github("lldelisle/usefulLDfunctions")
library(usefulLDfunctions)
safelyLoadAPackageInCRANorBioconductor("ggplot2")
# To display equations of regression:
safelyLoadAPackageInCRANorBioconductor("ggpmisc")
# To have p-values on plots:
safelyLoadAPackageInCRANorBioconductor("ggpubr")
# To combine multiple ggplots:
safelyLoadAPackageInCRANorBioconductor("gridExtra")
# To use italic in ggplot:
safelyLoadAPackageInCRANorBioconductor("ggtext")

# Set the reference gene:
ref.gene <- "GAPDH"
# Set the gene to plot
plot.gene <- "Hoxc13"
# List the genotypes
genotypes <- read.delim("genotypes.txt")
# Read the efficiency files
eff <- lapply(list.files(pattern = "efficiency"), read.delim, skip = 7)

# Efficiency
# Pool efficiency data
eff.df <- do.call(rbind, eff)
# Keep only data with sample name
eff.df <- eff.df[eff.df$Sample.Name != "", c("Target.Name", "Sample.Name", "Cт", "Quantity")]
# Simplify colnames
colnames(eff.df) <- c("target", "sample", "Ct", "quantity")
# Convert Ct to numerical values
eff.df$Ct <- as.numeric(eff.df$Ct)

# Calculate the cv for each point (same target, sample, quantity)
cv <- function(v){sd(v) / mean(v)}
cvs <- with(eff.df, aggregate(list(cv=Ct), by = list(target=target, sample=sample, quantity=quantity), FUN=cv))

# Add the info to the efficiency data.frame
eff.df <- merge(eff.df, cvs)

# Formula used to plot
my.formula <- y ~ x

# Plot Ct = f(log10(quantity))
ggplot(eff.df, aes(x = log10(quantity), y = Ct, color = target, shape = sample)) +
  geom_point()  +
  geom_smooth(formula = my.formula, method = "lm", se = FALSE) +
  ggtitle("Efficiencies from all plates") +
  stat_poly_eq(formula = my.formula, # Put the equations
               aes(label = paste0("atop(", ..eq.label.., ",", ..rr.label.., ")")), # Use atop to have one on top of the other
               parse = TRUE,
               label.x = 0.5, # Specify x as proportion
               label.y = (c(27, 22, 37) - 20) / 20, # Specify y as proportion
               size = 3) + # Change the font size
  ylim(20, 40)

# Restrict to data where cv is less than 2 %
# This will also remove data when one point was not detected.
eff.df <- subset(eff.df, cv < 0.02)

# Replot with the filtered data
g.eff <- ggplot(eff.df, aes(x = log10(quantity), y = Ct, color = target, shape = sample)) +
  geom_point()  +
  geom_smooth(formula = my.formula, method = "lm", se = FALSE) +
  ggtitle("Efficiencies from all plates only points with cv >= 2%") +
  stat_poly_eq(formula = my.formula, # Put the equations
               aes(label = paste0("atop(", ..eq.label.., ",", ..rr.label.., ")")), # Use atop to have one on top of the other
               parse = TRUE,
               label.x = c(0.5, 0.5, 0.7), # Specify x as proportion
               label.y = (c(27, 22, 37) - 20) / 20, # Specify y as proportion
               size = 3) + # Change the font size
  ylim(20, 40)

# Get the equation of the lines
lm_fits <- apply(unique(eff.df[, c("target", "sample")]), 1, function(v){
  lm(Ct ~ log10(quantity), 
     data = subset(eff.df, target == v["target"] & sample == v["sample"])
  )$coefficients
})
colnames(lm_fits) <- apply(unique(eff.df[, c("target", "sample")]), 1, "[", "target")
slopes <- apply(lm_fits, 2, "[[", "log10(quantity)")
# Convert to efficiencies
all.efficiencies <- 10 ^ (- 1 / slopes)
# Compute mean when multiple samples are used
efficiencies.df <- aggregate(x = all.efficiencies, by = list(target = names(all.efficiencies)), FUN = mean)
# Store them in a named vector
efficiencies <- efficiencies.df$x
names(efficiencies) <- efficiencies.df$target

# Process the data
experiments <- list.files(pattern = "del")
names(experiments) <- sapply(experiments, gsub, pattern = ".txt", replacement = "")

# Get as data.frame
all.df <- lapply(experiments, read.delim, skip=7)
# Keep only the 15 first cols and add the name of file
all.df <- lapply(names(experiments), function(n){
  df <- all.df[[n]][, 1:15]
  df$plate <- n
  return(df)
})
expe.df <- do.call(rbind, all.df)
empty.lines <- rowSums(is.na(expe.df)) == 15
expe.df <- expe.df[!empty.lines, ]
# I keep only the samples with names:
expe.df <- expe.df[expe.df$Sample.Name != "", ]

# Transform Ct to numeric
expe.df$Cт <- as.numeric(expe.df$Cт)
# Compute the Ct means
ct.means <- aggregate(list(ct_mean = expe.df$Cт),
                      by = list(sample = expe.df$Sample.Name,
                                target = expe.df$Target.Name,
                                plate = expe.df$plate),
                      FUN = mean)
# Plot as boxplot all ct covered by each gene
g.all <- ggplot(ct.means, aes(y = ct_mean, color = target)) +
  geom_boxplot() +
  ylab("mean value of Cts") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())

# Combine the efficiency plot with the box plot:
# To align 2 plots:
limits <- c(20, 40)
breaks <- seq(limits[1], limits[2], by=5)

# assign common axis to both plots
p1.common.y <- g.eff + scale_y_continuous(limits=limits, breaks=breaks)
p2.common.y <- g.all + scale_y_continuous(limits=limits, breaks=breaks)

# At this point, they have the same axis, but the axis lengths are unequal, so ...

# build the plots
p1.common.y <- ggplot_gtable(ggplot_build(p1.common.y))
p2.common.y <- ggplot_gtable(ggplot_build(p2.common.y))

# copy the plot height from p1 to p2
p2.common.y$heights <- p1.common.y$heights

# Put them next to the other:
grid.arrange(p1.common.y,p2.common.y,ncol=2, widths = c(3, 1))

# I compute the Ct for the ref gene for each sample
ref.Ct <- ct.means$ct_mean[ct.means$target == ref.gene]
names(ref.Ct) <- ct.means$sample[ct.means$target == ref.gene]
# I add it to the data.frame with ct.means
ct.means$ref.Ct <- ref.Ct[ct.means$sample]
# I compute the relative quantity using efficiencies
ct.means$relQuant_with_eff <- efficiencies[ct.means$target] ^ (- ct.means$ct_mean) / efficiencies[ref.gene] ^ (- ct.means$ref.Ct)
# When it is NA it means that it was not detected by RT-qPCR so quantity is set to 0
ct.means$relQuant_with_eff[is.na(ct.means$relQuant_with_eff)] <- 0

# I keep only the gene to plot:
ct.means.plot <- subset(ct.means, target == plot.gene)

# Annotate ct.means.plot
ct.means.plot$genotype <- genotypes$genotype[match(ct.means.plot$sample, genotypes$sample)]
ct.means.plot$dissection <- genotypes$dissection[match(ct.means.plot$sample, genotypes$sample)]
ct.means.plot$dissection_geno <- paste0(ct.means.plot$dissection, "_", ct.means.plot$genotype)

# I set up the order for the dissection_geno
my.order <- c("tip_Wt", "digit_Wt", "tip_*HoxCKO*_Het", "tip_*DelEC1*",
              "digit_*DelEC2*", "tip_*DelEC1-EC2*", "tip_TH", "tip_*HoxCKO*")
# And deduce it for the genotypes
my.geno.order <- sapply(my.order, gsub, pattern = "tip_|digit_", replacement = "")
names(my.geno.order) <- my.order

# I apply the new orders to the factors of the data frame
ct.means.plot$genotype <- factor(ct.means.plot$genotype, levels = unique(my.geno.order))
ct.means.plot$dissection_geno <- factor(ct.means.plot$dissection_geno, levels = unique(my.order))

# I plot the values for each dissection and each genotype
ggplot(ct.means.plot, aes(x = genotype, y = relQuant_with_eff, color = genotype)) +
  geom_boxplot() +
  geom_jitter(aes(shape = plate)) +
  facet_grid(. ~ dissection, scales = "free_x") +
  scale_shape_manual(values=1:length(unique(ct.means.plot$plate))) +
  ylab("relative quantity to GAPDH \n(Using efficiency)") +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust=1), # Put the genotypes vertical
        legend.text = ggtext::element_markdown())

# Then I compute the relative to wt from same dissection
# Fist I get the mean of Wt
wt.mean.dissection.df <- with(subset(ct.means.plot, genotype == "Wt"),
                              aggregate(relQuant_with_eff,
                                        by = list(dissection),
                                        FUN = mean))
# Put it in a named vector
wt.mean.dissection <- wt.mean.dissection.df$x
names(wt.mean.dissection) <- wt.mean.dissection.df$Group.1
# Add the value norm to wt in the data frame
ct.means.plot$relQuant_with_eff_norm_wt <- ct.means.plot$relQuant_with_eff / 
  wt.mean.dissection[ct.means.plot$dissection]

# In order to display on the plot the number for each dissection_genotypes
# I compute a table
table.dissection_geno <- table(ct.means.plot$dissection_geno)
# Write it in markdown:
geno.nb <- sapply(my.order, function(g){paste0(my.geno.order[g], "<br/>(n=", table.dissection_geno[g], ")")})

# To display p-values in the graph, I create a list of pairs:
my.comparisons <- lapply(setdiff(grep("Wt", my.order, value = T, invert = T),
                                 names(table.dissection_geno[table.dissection_geno <= 1])),
                         function(s){c(paste0(strsplit(s, "_")[[1]][1], "_Wt"), s)})

# I plot the relative level
g.witheff <- ggboxplot(ct.means.plot, x = "dissection_geno",
                       y = "relQuant_with_eff_norm_wt",
                       color = "genotype", 
                       outlier.shape = NA) + # This is needed to avoid duplicated points with jitter
  geom_jitter(aes(color = genotype, fill = dissection), shape = 21) +
  ylab(paste0("*", plot.gene, "* relative levels")) +
  stat_compare_means(method = "t.test", comparisons = my.comparisons,
                     aes(label = ..p.signif..)) + # Use stars instead of p-values
  theme(axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust=1), # Put the genotypes vertical and in markdown
        panel.grid.major.y = element_line(colour = "grey"), # Add horizontal lines
        axis.title.y = ggtext::element_markdown(), # Put the y title in markdown
        legend.text = ggtext::element_markdown()) + # Put the legend in markdown
  scale_y_continuous(limits=c(0, 2.3), breaks=seq(0, 1.4, 0.2)) + # Set the position of hlines
  scale_fill_manual(values = c(rgb(0.6,0.6,0.6), rgb(0.1,0.1,0.1))) + # Set the fill colors
  scale_x_discrete(labels = geno.nb) + # Change the x-label to the genotype with numbers
  xlab("Genotype") + 
  guides(fill=FALSE) # Remove the fill from the legend

print(g.witheff)
ggsave("Figure6_good_size.eps", g.witheff, scale = 2, width = 8.7, height = 6, units = "cm")
