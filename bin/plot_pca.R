#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
    stop("Usage: plot_pca.R <pca_file> <sample_info>")
}

pca_file <- args[1]
sample_info <- args[2]

# Create debug file
sink("debug_info.txt")

# Read PCA results and eigenvalues
cat("Reading PCA results and eigenvalues...\n")
pca <- read.table(pca_file, header=TRUE, stringsAsFactors=FALSE)
colnames(pca) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
eigenvals <- scan("pca.eigenval")
pc1_var <- round(eigenvals[1] / sum(eigenvals) * 100, 1)
pc2_var <- round(eigenvals[2] / sum(eigenvals) * 100, 1)

# Remove ancient samples
pca <- pca[!grepl("^ancient", pca$IID), ]

# Check if we have variation in the PCs
pc1_var_check <- var(pca$PC1)
pc2_var_check <- var(pca$PC2)

cat("PC1 variance:", pc1_var_check, "\n")
cat("PC2 variance:", pc2_var_check, "\n")

if (pc1_var_check == 0 || pc2_var_check == 0) {
    cat("WARNING: No variation in principal components!\n")
    cat("PC1 unique values:", length(unique(pca$PC1)), "\n")
    cat("PC2 unique values:", length(unique(pca$PC2)), "\n")
    cat("First few PC1 values:", head(pca$PC1), "\n")
    cat("First few PC2 values:", head(pca$PC2), "\n")
    stop("Cannot create PCA plot - no variation in the data")
}

# Scale PCs to have mean 0 and unit variance
pca$PC1_scaled <- scale(pca$PC1)
pca$PC2_scaled <- scale(pca$PC2)

cat("PCA dimensions:", dim(pca), "\n")
cat("First few rows of PCA:\n")
print(head(pca))

# Read sample info
cat("\nReading sample info...\n")
samples <- read.table(sample_info, header=FALSE, stringsAsFactors=FALSE)
cat("Sample info dimensions:", dim(samples), "\n")
cat("First few rows of samples:\n")
print(head(samples))

# Clean up sample IDs in PCA data to match sample info
pca$IID <- gsub("_.*$", "", pca$IID)

# Merge population information
pca$Population <- samples$V2[match(pca$IID, samples$V1)]

# Print population counts
cat("\nPopulation counts:\n")
pop_counts <- table(pca$Population, useNA="ifany")
print(pop_counts)

# Create plot
pdf("pca_plot.pdf", width=12, height=10)

# Define colors for different populations
populations <- sort(unique(na.omit(pca$Population)))
colors <- rainbow(length(populations))
names(colors) <- populations

# Create diagnostic plots
par(mfrow=c(2,2))

# Histogram of PC1
hist(pca$PC1, main="PC1 Distribution", xlab="PC1")

# Histogram of PC2
hist(pca$PC2, main="PC2 Distribution", xlab="PC2")

# Scatterplot of raw values
plot(pca$PC1, pca$PC2,
     col=colors[pca$Population],
     pch=20,
     xlab=paste0("PC1 (", pc1_var, "% variance explained)"),
     ylab=paste0("PC2 (", pc2_var, "% variance explained)"),
     main="PCA of 1000G Samples (Raw Values)")

# Add legend
legend("topright",
       legend=names(colors),
       col=colors,
       pch=20,
       title="Population",
       cex=0.6,
       ncol=2)

# Scatterplot of scaled values
plot(pca$PC1_scaled, pca$PC2_scaled,
     col=colors[pca$Population],
     pch=20,
     xlab="PC1 (scaled)",
     ylab="PC2 (scaled)",
     main="PCA of 1000G Samples (Scaled Values)")

# Print summary statistics
cat("\nPC Summary Statistics:\n")
cat("Original PC1 range:", range(pca$PC1), "\n")
cat("Original PC2 range:", range(pca$PC2), "\n")
cat("PC1 variance explained:", pc1_var, "%\n")
cat("PC2 variance explained:", pc2_var, "%\n")

dev.off()
sink() 