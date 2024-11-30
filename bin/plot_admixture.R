#!/usr/bin/env Rscript

# Install and load required packages
for(package in c("optparse", "ggplot2", "tidyr", "dplyr", "jsonlite")) {
    if (!requireNamespace(package, quietly = TRUE)) {
        install.packages(package, repos="https://cloud.r-project.org", quiet=TRUE)
    }
    library(package, character.only = TRUE)
}

# Define command line options
option_list <- list(
    make_option(c("-q", "--qfiles"), 
                help="Directory containing ADMIXTURE .Q files"),
    make_option(c("-s", "--sample-info"),
                help="Sample information file"),
    make_option(c("-p", "--pop-info"),
                help="Population information file"),
    make_option(c("-o", "--output-prefix"),
                default="admixture",
                help="Output file prefix [default: %default]")
)

# Parse command line arguments
opt <- parse_args(OptionParser(option_list=option_list))

# Function to log both to file and console
log_message <- function(...) {
    message <- paste(...)
    cat(message, "\n")  # Print to console
    cat(message, "\n", file=paste0(opt$output_prefix, "_debug.txt"), append=TRUE)
}

# Start logging
log_message("Starting ADMIXTURE plot generation...")
log_message("Input parameters:")
log_message("  Q files directory:", opt$qfiles)
log_message("  Sample info:", opt$`sample-info`)
log_message("  Population info:", opt$`pop-info`)
log_message("  Output prefix:", opt$output_prefix)

# Read sample and population info
samples <- read.table(opt$`sample-info`, header=TRUE, sep="\t")
pop_info <- read.delim(opt$`pop-info`)

log_message("Read", nrow(samples), "samples")
log_message("Read", nrow(pop_info), "population entries")

# Function to read Q file
read_q_file <- function(file) {
    k <- ncol(read.table(file))
    data <- read.table(file)
    colnames(data) <- paste0("Comp", 1:k)
    data$Sample <- seq_len(nrow(data))
    return(data)
}

# Read all Q files
q_files <- list.files(opt$qfiles, pattern="\\.Q$", full.names=TRUE)
log_message("Looking for Q files in:", opt$qfiles)
log_message("Found Q files:", paste(basename(q_files), collapse=", "))

if(length(q_files) == 0) {
    stop("No .Q files found in directory:", opt$qfiles)
}

admixture_data <- lapply(q_files, read_q_file)
names(admixture_data) <- paste0("K", sapply(admixture_data, function(x) ncol(x) - 1))

# Debug output
log_message("Number of K values found:", length(admixture_data))
for(k in names(admixture_data)) {
    log_message("Processing", k, "with dimensions:", 
                paste(dim(admixture_data[[k]]), collapse=" x "))
}

# Create plots for each K
for(k in names(admixture_data)) {
    log_message("Processing", k)
    data <- admixture_data[[k]]
    
    # Merge with sample info
    data$IID <- samples$Individual_ID[data$Sample]
    data$is_ancient <- grepl("^ancient", data$IID)  # Add ancient sample flag
    data <- merge(data, samples, by.x="IID", by.y="Individual_ID")
    data <- merge(data, pop_info[, c("Population.code", "Superpopulation.code")],
                 by.x="Population", by.y="Population.code", all.x=TRUE)
    
    # Create plot data
    long_data <- gather(data, Component, Value, starts_with("Comp"))
    
    # Modern samples plot
    modern_data <- subset(long_data, !is_ancient)
    p1 <- ggplot(modern_data, aes(x=IID, y=Value, fill=Component)) +
        geom_bar(stat="identity", width=1) +
        facet_grid(~Superpopulation.code, scales="free_x", space="free") +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle=90, hjust=1, size=6),
            panel.spacing = unit(0.1, "lines"),
            strip.text = element_text(face="bold"),
            strip.background = element_rect(fill="grey95"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank()
        ) +
        labs(title=paste("ADMIXTURE Results -", k, "ancestral populations"),
             subtitle="Modern Samples",
             x="Sample", y="Ancestry Proportion")
    
    # Save modern samples plot
    ggsave(paste0(opt$output_prefix, "_", k, "_modern.png"), 
           p1, width=15, height=8, dpi=300)
    
    # Ancient samples plot (if any)
    if(any(data$is_ancient)) {
        ancient_data <- subset(long_data, is_ancient)
        p2 <- ggplot(ancient_data, aes(x=IID, y=Value, fill=Component)) +
            geom_bar(stat="identity", width=1) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle=45, hjust=1, size=8),
                panel.spacing = unit(0.1, "lines"),
                strip.text = element_text(face="bold"),
                panel.grid.major.x = element_blank(),
                panel.grid.minor = element_blank()
            ) +
            labs(title=paste("ADMIXTURE Results -", k, "ancestral populations"),
                 subtitle="Ancient Samples",
                 x="Sample", y="Ancestry Proportion") +
            # Add black border around bars
            geom_rect(aes(xmin=as.numeric(factor(IID))-0.5,
                         xmax=as.numeric(factor(IID))+0.5,
                         ymin=0, ymax=1),
                     fill=NA, color="black", linewidth=0.5)
        
        # Save ancient samples plot
        ggsave(paste0(opt$output_prefix, "_", k, "_ancient.png"), 
               p2, width=15, height=8, dpi=300)
        
        # Combined view
        p3 <- ggplot(long_data, aes(x=IID, y=Value, fill=Component)) +
            geom_bar(stat="identity", width=1) +
            facet_grid(~is_ancient, scales="free_x", space="free",
                      labeller=labeller(is_ancient=c("FALSE"="Modern Samples", 
                                                   "TRUE"="Ancient Samples"))) +
            theme_minimal() +
            theme(
                axis.text.x = element_text(angle=90, hjust=1, size=6),
                panel.spacing = unit(0.1, "lines"),
                strip.text = element_text(face="bold"),
                strip.background = element_rect(fill="grey95"),
                panel.grid.major.x = element_blank(),
                panel.grid.minor = element_blank()
            ) +
            labs(title=paste("ADMIXTURE Results -", k, "ancestral populations"),
                 subtitle="Combined View",
                 x="Sample", y="Ancestry Proportion")
        
        # Add black border around ancient samples
        if(any(data$is_ancient)) {
            p3 <- p3 + geom_rect(data=subset(long_data, is_ancient),
                                aes(xmin=as.numeric(factor(IID))-0.5,
                                    xmax=as.numeric(factor(IID))+0.5,
                                    ymin=0, ymax=1),
                                fill=NA, color="black", linewidth=0.5)
        }
        
        # Save combined view plot
        ggsave(paste0(opt$output_prefix, "_", k, "_combined.png"), 
               p3, width=15, height=8, dpi=300)
    }
}

# Create JSON output with ancient sample information
json_data <- lapply(names(admixture_data), function(k) {
    data <- admixture_data[[k]]
    data$IID <- samples$Individual_ID[data$Sample]
    data$is_ancient <- grepl("^ancient", data$IID)
    
    list(
        k = as.numeric(gsub("K", "", k)),
        samples = data$IID,
        is_ancient = data$is_ancient,
        proportions = as.matrix(data[, grep("Comp", colnames(data))])
    )
})

# Save JSON output
json_file <- paste0(opt$output_prefix, "_data.json")
log_message("Saving JSON data to", json_file)
writeLines(toJSON(json_data, pretty=TRUE), json_file)

log_message("ADMIXTURE plot generation complete")
