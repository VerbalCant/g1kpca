#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
    stop("Usage: create_pca_json.R <pca_eigenvec> <pca_eigenval> <sample_info> <pop_info>")
}

pca_eigenvec <- args[1]
pca_eigenval <- args[2]
sample_info_file <- args[3]
pop_info_file <- args[4]

# Load required libraries
if (!requireNamespace("jsonlite", quietly = TRUE)) {
    install.packages("jsonlite", repos="https://cloud.r-project.org")
}
library(jsonlite)

# Read PCA results
pca <- read.table(pca_eigenvec, header=FALSE)
colnames(pca) <- c("FID", "IID", paste0("PC", 1:10))

# Read sample info
samples <- read.table(sample_info_file, header=TRUE, sep="\t")

# Read population info
pop_info <- read.delim(pop_info_file)

# Print column names for debugging
cat("Population info columns:", paste(colnames(pop_info), collapse=", "), "\n")

# Split samples into modern, ancient, and SRR (Chachapoya)
is_ancient <- grepl("^ancient", pca$IID)
is_srr <- grepl("^SRR", pca$IID)
modern_pca <- pca[!is_ancient & !is_srr,]
ancient_pca <- pca[is_ancient,]
chachapoya_pca <- pca[is_srr,]

# Process modern samples
modern_pca$IID <- gsub("_.*$", "", modern_pca$IID)
modern_data <- merge(modern_pca, samples, by.x="IID", by.y="Individual_ID", all.x=TRUE)
modern_data <- merge(modern_data, 
                    pop_info[, c("Population.code", "Population.name", 
                                "Superpopulation.code", "Superpopulation.name",
                                "Superpopulation.display.colour")],
                    by.x="Population", 
                    by.y="Population.code", 
                    all.x=TRUE)

# Create modern samples data structure
modern_final <- data.frame(
    sample = modern_data$IID,
    population = modern_data$Population,
    pop_name = modern_data$Population.name,
    superpop = modern_data$Superpopulation.code,
    superpop_name = modern_data$Superpopulation.name,
    color = modern_data$Superpopulation.display.colour,
    is_ancient = FALSE,
    PC1 = modern_data$PC1,
    PC2 = modern_data$PC2,
    PC3 = modern_data$PC3
)

# Initialize final_data with modern samples
final_data <- modern_final

# Add ancient samples if present
if (nrow(ancient_pca) > 0) {
    ancient_final <- data.frame(
        sample = ancient_pca$IID,
        population = "Ancient",
        pop_name = "Ancient Sample",
        superpop = "ANC",
        superpop_name = "Ancient Samples",
        color = "#000000",  # Black color for ancient samples
        is_ancient = TRUE,
        PC1 = ancient_pca$PC1,
        PC2 = ancient_pca$PC2,
        PC3 = ancient_pca$PC3
    )
    final_data <- rbind(final_data, ancient_final)
}

# Add Chachapoya samples if present
if (nrow(chachapoya_pca) > 0) {
    chachapoya_final <- data.frame(
        sample = chachapoya_pca$IID,
        population = "Chachapoya",
        pop_name = "Chachapoya",
        superpop = "ANC",
        superpop_name = "Ancient Samples",
        color = "#FF0000",  # Red color for Chachapoya samples
        is_ancient = TRUE,
        PC1 = chachapoya_pca$PC1,
        PC2 = chachapoya_pca$PC2,
        PC3 = chachapoya_pca$PC3
    )
    final_data <- rbind(final_data, chachapoya_final)
}

# Add debug information
cat("Number of modern samples:", nrow(modern_final), "\n")
cat("Number of ancient samples:", nrow(ancient_pca), "\n")
cat("Number of Chachapoya samples:", nrow(chachapoya_pca), "\n")
cat("Total samples in final data:", nrow(final_data), "\n")

# Write to JSON
writeLines(toJSON(final_data, pretty=TRUE), "pca_data.json")
