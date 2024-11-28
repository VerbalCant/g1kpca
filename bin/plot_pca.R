#!/usr/bin/env Rscript

# Install and load required package
if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages("optparse", repos="https://cloud.r-project.org", quiet=TRUE)
}
library(optparse)

# Define command line options
option_list <- list(
    make_option(c("-p", "--pca-file"), 
                default="./pca.eigenvec",
                help="Input PCA eigenvectors file [default: %default]"),
    make_option(c("-e", "--eigenval-file"),
                default="./pca.eigenval",
                help="Input eigenvalues file [default: %default]"),
    make_option(c("-s", "--sample-info"),
                default="./20130606_g1k.ped",
                help="Sample information file [default: %default]"),
    make_option(c("-i", "--pop-info"),
                default="./igsr_populations.tsv",
                help="Population information file [default: %default]"),
    make_option(c("-o", "--output-prefix"),
                default="./pca_plot",
                help="Output file prefix [default: %default]")
)

# Parse command line arguments
opt <- parse_args(OptionParser(option_list=option_list))

# Function to log both to file and console
log_message <- function(...) {
    message <- paste(...)
    cat(message, "\n")  # Print to console
    cat(message, "\n", file=paste0(opt$`output-prefix`, "_debug.txt"), append=TRUE)  # Print to file
}

# Start debugging output
log_message("Starting PCA plot generation...")
log_message("Input files:")
log_message("  PCA file:", opt$`pca-file`)
log_message("  Eigenval file:", opt$`eigenval-file`)
log_message("  Sample info:", opt$`sample-info`)
log_message("  Population info:", opt$`pop-info`)

# Check files
for(file_path in c(opt$`pca-file`, opt$`eigenval-file`, opt$`sample-info`, opt$`pop-info`)) {
    if(!file.exists(file_path)) {
        log_message("ERROR: File not found:", file_path)
        stop(paste("File not found:", file_path))
    }
}

# Read population info with debugging
log_message("\nReading population info...")
pop_info <- read.delim(opt$`pop-info`, stringsAsFactors=FALSE)
log_message("Population info column names:", paste(colnames(pop_info), collapse=", "))
log_message("First few rows of Superpopulation codes:")
log_message(head(pop_info$`Superpopulation code`))

# Check for empty or NA values
log_message("\nChecking population info:")
log_message("Total rows:", nrow(pop_info))
log_message("Rows with non-NA superpopulation codes:", 
           sum(!is.na(pop_info$`Superpopulation code`)))
log_message("Unique superpopulation codes:", 
           paste(unique(pop_info$`Superpopulation code`[!is.na(pop_info$`Superpopulation code`)]), 
                 collapse=", "))

# Filter population info
pop_info_filtered <- pop_info[!is.na(pop_info$Superpopulation.code) & 
                            pop_info$Superpopulation.code != "", ]

log_message("\nAfter filtering:")
log_message("Remaining rows:", nrow(pop_info_filtered))
log_message("Unique superpopulation codes:", 
           paste(unique(pop_info_filtered$Superpopulation.code), collapse=", "))

# Create lookup tables from filtered data
superpop_colors <- setNames(pop_info_filtered$Superpopulation.display.colour, 
                          pop_info_filtered$Superpopulation.code)
superpop_names <- setNames(pop_info_filtered$Superpopulation.name, 
                          pop_info_filtered$Superpopulation.code)
pop_names <- setNames(pop_info_filtered$Population.name, 
                     pop_info_filtered$Population.code)

log_message("\nLookup tables:")
log_message("Superpopulation colors:", paste(names(superpop_colors), superpop_colors, sep="=", collapse=", "))
log_message("Superpopulation names:", paste(names(superpop_names), superpop_names, sep="=", collapse=", "))
log_message("Number of population names:", length(pop_names))

# Define shapes
shapes <- c(15:19, 21:25)
superpops <- unique(pop_info_filtered$Superpopulation.code)
log_message("\nSuperpopulations found:", paste(superpops, collapse=", "))

# Create pop_shapes list
pop_shapes <- list()
for(sp in superpops) {
    pops_in_sp <- pop_info_filtered$Population.code[pop_info_filtered$Superpopulation.code == sp]
    pop_shapes[[sp]] <- setNames(shapes[1:length(pops_in_sp)], pops_in_sp)
    log_message("Populations in", sp, ":", paste(pops_in_sp, collapse=", "))
}

# Read PCA results and eigenvalues with debugging
log_message("Reading PCA results...")
tryCatch({
    pca_lines <- readLines(opt$`pca-file`)
    log_message("First few lines of PCA file:")
    log_message(head(pca_lines), sep="\n")
    pca <- read.table(opt$`pca-file`, header=TRUE, stringsAsFactors=FALSE)
    log_message("PCA dimensions:", paste(dim(pca), collapse=" x "))
}, error = function(e) {
    log_message("Error reading PCA file:", conditionMessage(e))
    stop("Failed to read PCA file")
})

colnames(pca) <- c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

log_message("Reading eigenvalues...")
tryCatch({
    eigenval_lines <- readLines(opt$`eigenval-file`)
    log_message("First few lines of eigenval file:")
    log_message(head(eigenval_lines), sep="\n")
    eigenvals <- scan(opt$`eigenval-file`)
    log_message("Number of eigenvalues:", length(eigenvals))
}, error = function(e) {
    log_message("Error reading eigenval file:", conditionMessage(e))
    stop("Failed to read eigenval file")
})

pc1_var <- round(eigenvals[1] / sum(eigenvals) * 100, 1)
pc2_var <- round(eigenvals[2] / sum(eigenvals) * 100, 1)
pc3_var <- round(eigenvals[3] / sum(eigenvals) * 100, 1)

# Read sample info and merge population information
log_message("Reading sample info...")
tryCatch({
    # First check the file content
    ped_lines <- readLines(opt$`sample-info`)
    log_message("First few lines of sample info file:")
    log_message(head(ped_lines), sep="\n")
    
    # Read the file with header
    samples <- read.table(
        opt$`sample-info`, 
        header=TRUE,
        stringsAsFactors=FALSE,
        sep="\t",
        quote="",  # Don't treat anything as quotes
        comment.char=""  # Don't treat anything as comments
    )
    
    log_message("Sample info dimensions:", paste(dim(samples), collapse=" x "))
    log_message("First few rows of processed sample info:")
    print(head(samples))
    
    # Verify the data
    log_message("Number of unique populations:", length(unique(samples$Population)))
    log_message("Unique populations:", paste(sort(unique(samples$Population)), collapse=", "))
    
}, error = function(e) {
    log_message("Error reading sample info:", conditionMessage(e))
    stop("Failed to read sample info")
})

# Split into modern and ancient samples
ancient_samples <- pca[grepl("^ancient", pca$IID), ]
modern_samples <- pca[!grepl("^ancient", pca$IID), ]

# Check if we have variation in the PCs
for(pc in c("PC1", "PC2", "PC3")) {
    pc_var <- var(pca[[pc]])
    log_message(pc, "variance:", pc_var)
    if(pc_var == 0) {
        stop(paste("No variation in", pc))
    }
}

# Clean up sample IDs in PCA data to match sample info
modern_samples$IID <- gsub("_.*$", "", modern_samples$IID)

modern_samples$Population <- samples$Population[match(modern_samples$IID, samples$Individual_ID)]
modern_samples$Superpopulation <- pop_info_filtered$Superpopulation.code[match(modern_samples$Population, 
                                                                             pop_info_filtered$Population.code)]

# Function to create PCA plot
create_pca_plot <- function(pc_x, pc_y, var_x, var_y) {
    pdf(paste0(opt$`output-prefix`, "_", pc_x, "_", pc_y, ".pdf"), width=15, height=10)
    
    # Set up plot with extra space for legend
    par(mar=c(5,4,4,10))
    
    # Debug modern samples data
    log_message("\nModern samples data check:")
    log_message("Number of modern samples:", nrow(modern_samples))
    log_message("Sample of modern data:")
    print(head(modern_samples[, c("IID", "Population", "Superpopulation", pc_x, pc_y)]))
    
    # Create plot with grey background
    plot(modern_samples[[pc_x]], modern_samples[[pc_y]],
         type="n",  # Set up plot without points first
         xlab=paste0(pc_x, " (", var_x, "% variance explained)"),
         ylab=paste0(pc_y, " (", var_y, "% variance explained)"),
         main=paste("PCA of Modern and Ancient Samples:", pc_x, "vs", pc_y),
         panel.first={
             rect(par("usr")[1], par("usr")[3],
                  par("usr")[2], par("usr")[4],
                  col="grey95")  # very light grey
             grid(col="white", lty=1)  # white grid lines
         })
    
    # Plot modern samples by superpopulation
    for(sp in unique(modern_samples$Superpopulation)) {
        if(is.na(sp)) next
        sp_samples <- modern_samples[modern_samples$Superpopulation == sp,]
        for(pop in unique(sp_samples$Population)) {
            if(is.na(pop)) next
            pop_samples <- sp_samples[sp_samples$Population == pop,]
            points(pop_samples[[pc_x]], pop_samples[[pc_y]],
                  col=superpop_colors[sp],
                  pch=pop_shapes[[sp]][pop],
                  cex=0.8)
        }
    }
    
    # Add ancient samples with enhanced visibility
    if(nrow(ancient_samples) > 0) {
        # Add white halos for better contrast against grey background
        points(ancient_samples[[pc_x]], ancient_samples[[pc_y]],
               pch=24,  # triangle with border
               col="white",  # white border
               bg="black",  # black fill
               cex=3,    # larger size
               lwd=2)    # thicker border
        
        # Add labels with white background and offset
        text(ancient_samples[[pc_x]], ancient_samples[[pc_y]],
             labels=ancient_samples$IID,
             pos=4,     # position to the right of point
             offset=1,  # add some space between point and text
             cex=1,     # slightly larger text
             bg="white", # white background for text
             font=2,    # bold text
             adj=0)     # left-align text
    }
    
    # Create single combined legend
    legend_x <- par("usr")[2] * 1.05
    legend_y <- par("usr")[4]
    
    # Prepare legend entries
    legend_entries <- character(0)
    legend_cols <- character(0)
    legend_pchs <- numeric(0)
    
    # Add superpopulation entries first
    for(sp in sort(superpops)) {
        # Add superpopulation header
        legend_entries <- c(legend_entries, paste(sp, "-", superpop_names[sp]))
        legend_cols <- c(legend_cols, superpop_colors[sp])
        legend_pchs <- c(legend_pchs, 15)  # filled square for superpopulations
        
        # Add populations within this superpopulation
        pops_in_sp <- sort(names(pop_shapes[[sp]]))
        for(pop in pops_in_sp) {
            legend_entries <- c(legend_entries, paste("  ", pop, "-", pop_names[pop]))  # indent populations
            legend_cols <- c(legend_cols, superpop_colors[sp])
            legend_pchs <- c(legend_pchs, pop_shapes[[sp]][pop])
        }
    }
    
    # Add ancient samples at the end
    if(nrow(ancient_samples) > 0) {
        legend_entries <- c(legend_entries, "Ancient Samples")
        legend_cols <- c(legend_cols, "black")
        legend_pchs <- c(legend_pchs, 24)  # match the triangle with border
    }
    
    # Create single legend
    legend(legend_x, legend_y,
           legend=legend_entries,
           col=legend_cols,
           pch=legend_pchs,
           title="Populations",
           xpd=TRUE,
           cex=0.7,
           ncol=1)
    
    dev.off()
}

# Create both plots
create_pca_plot("PC1", "PC2", pc1_var, pc2_var)
create_pca_plot("PC1", "PC3", pc1_var, pc3_var)

# Print summary statistics
log_message("\nPC Summary Statistics:")
log_message("PC1 variance explained:", pc1_var, "%")
log_message("PC2 variance explained:", pc2_var, "%")
log_message("PC3 variance explained:", pc3_var, "%")

# Print ancient sample coordinates
if(nrow(ancient_samples) > 0) {
    log_message("\nAncient Sample Coordinates:")
    print(ancient_samples[, c("IID", "PC1", "PC2", "PC3")])
}