process run_pca {
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    publishDir params.outdir, mode: 'copy'
    label 'high_cpu'

    input:
    path vcf

    output:
    path 'pca.*', emit: pca_results
    path 'debug_*', emit: debug

    script:
    def plink_mem_mb = (task.memory.toMega() * 0.8).intValue()
    """
    echo "[DEBUG] Starting PCA analysis with ${task.cpus} CPUs and ${plink_mem_mb}MB memory"
    echo "[DEBUG] Total available memory: ${task.memory}"

    # First check input VCF
    echo "Checking input VCF..."
    bcftools stats --threads ${task.cpus} ${vcf} > debug_vcf_stats.txt

    # Convert VCF to PLINK format with minimal filtering first
    echo "Converting to PLINK format..."
    plink2 --vcf ${vcf} \
        --threads ${task.cpus} \
        --memory ${plink_mem_mb} \
        --make-bed \
        --allow-extra-chr \
        --set-missing-var-ids @:#,\\\$1,\\\$2 \
        --max-alleles 2 \
        --out temp_plink_initial

    # Check initial conversion
    echo "Initial PLINK conversion stats:" > debug_plink_steps.txt
    wc -l temp_plink_initial.bim >> debug_plink_steps.txt
    wc -l temp_plink_initial.fam >> debug_plink_steps.txt

    # Now do filtering steps one at a time
    echo "Filtering steps:" >> debug_plink_steps.txt

    # Step 1: Filter by missingness with a more lenient threshold
    plink2 --bfile temp_plink_initial \
        --threads ${task.cpus} \
        --memory ${plink_mem_mb} \
        --geno 0.2 \
        --mind 0.2 \
        --make-bed \
        --out temp_plink_missing
    echo "After missingness filtering:" >> debug_plink_steps.txt
    wc -l temp_plink_missing.bim >> debug_plink_steps.txt

    # Step 2: Filter by MAF with a much lower threshold
    plink2 --bfile temp_plink_missing \
        --threads ${task.cpus} \
        --memory ${plink_mem_mb} \
        --maf 0.001 \
        --make-bed \
        --out temp_plink_final
    echo "After MAF filtering:" >> debug_plink_steps.txt
    wc -l temp_plink_final.bim >> debug_plink_steps.txt

    # Check the number of variants and samples
    n_variants=\$(wc -l < temp_plink_final.bim)
    n_samples=\$(wc -l < temp_plink_final.fam)
    echo "Final counts:" >> debug_plink_steps.txt
    echo "Variants: \$n_variants" >> debug_plink_steps.txt
    echo "Samples: \$n_samples" >> debug_plink_steps.txt

    if [ \$n_variants -lt ${params.min_variants} ]; then
        echo "Error: Too few variants (\$n_variants) after filtering"
        exit 1
    fi

    if [ \$n_samples -lt ${params.min_samples} ]; then
        echo "Error: Too few samples (\$n_samples) after filtering"
        exit 1
    fi

    # Run PCA with allele weights for multiallelic variants
    plink2 --bfile temp_plink_final \
        --threads ${task.cpus} \
        --memory ${plink_mem_mb} \
        --pca allele-wts approx \
        --allow-extra-chr \
        --out pca

    # Verify PCA output
    if [ ! -f pca.eigenvec ] || [ ! -f pca.eigenval ]; then
        echo "Error: PCA failed to produce output files"
        exit 1
    fi

    # Check eigenvalues
    head -n 10 pca.eigenval > debug_eigenvalues.txt

    # Check first few lines of eigenvectors
    head -n 10 pca.eigenvec > debug_eigenvectors.txt

    # Check for variation in eigenvectors
    Rscript -e '
        pca <- read.table("pca.eigenvec", header=TRUE)
        pc1_var <- var(pca[,3])
        pc2_var <- var(pca[,4])
        cat(sprintf("PC1 variance: %g\\n", pc1_var), file="debug_pc_variance.txt")
        cat(sprintf("PC2 variance: %g\\n", pc2_var), file="debug_pc_variance.txt")
        if (pc1_var < 1e-10 || pc2_var < 1e-10) {
            cat("Error: No variation in principal components\\n")
            quit(status=1)
        }
    '

    # Add resource usage to debug output
    echo "Resource Usage:" >> debug_plink_steps.txt
    echo "  CPUs: ${task.cpus}" >> debug_plink_steps.txt
    echo "  Memory: ${task.memory}" >> debug_plink_steps.txt
    echo "  PLINK memory: ${plink_mem_mb}MB" >> debug_plink_steps.txt
    """
}