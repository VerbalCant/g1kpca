process ld_prune {
    cpus params.max_cpus

    input:
    path vcf

    output:
    path 'pruned.vcf.gz', emit: vcf
    path 'timing.txt', emit: timing

    script:
    """
    start_time=\$(date +%s)
    echo "Starting LD pruning at \$(date)"

    # Index VCF if needed
    if [ ! -f ${vcf}.tbi ] && [ ! -f ${vcf}.csi ]; then
        echo "Indexing input VCF..."
        bcftools index --tbi ${vcf}
    fi

    # Convert to plink format and perform LD pruning
    echo "Step 1: LD pruning"
    echo "Input VCF: ${vcf}"
    ls -l ${vcf}
    
    # Filter out variants with long alleles using bcftools with expression
    echo "Filtering variants with long alleles..."
    bcftools view ${vcf} \
        --include 'strlen(REF)<=100 & strlen(ALT)<=100' \
        -Oz -o filtered_variants.vcf.gz

    bcftools index filtered_variants.vcf.gz

    # First pass: Convert to PLINK format with unique IDs
    echo "Converting to PLINK format with unique IDs..."
    plink2 --vcf filtered_variants.vcf.gz \
        --threads ${task.cpus} \
        --set-all-var-ids '@:#:\$r:\$a' \
        --new-id-max-allele-len 97 \
        --snps-only \
        --max-alleles 2 \
        --chr ${params.chromosomes} \
        --make-bed \
        --out temp_initial

    # Now run LD pruning on the bed file
    echo "Running LD pruning..."
    plink2 --bfile temp_initial \
        --threads ${task.cpus} \
        --snps-only \
        --indep-pairwise ${params.ld_window_size} ${params.ld_step_size} ${params.ld_r2_threshold} \
        --out pruned

    # Check if pruning was successful
    if [ ! -f pruned.prune.in ]; then
        echo "Error: LD pruning failed to produce output file"
        exit 1
    fi

    # Extract pruned variants
    echo "Step 2: Extracting pruned variants"
    plink2 --bfile temp_initial \
        --threads ${task.cpus} \
        --extract pruned.prune.in \
        --recode vcf bgz \
        --out pruned \
        --double-id

    # End timing
    end_time=\$(date +%s)
    runtime=\$((end_time-start_time))
    
    # Write timing report
    echo "Process: LD_PRUNE" > timing.txt
    echo "Start time: \$(date -d @\${start_time})" >> timing.txt
    echo "End time: \$(date -d @\${end_time})" >> timing.txt
    echo "Runtime: \${runtime} seconds" >> timing.txt
    echo "Parameters:" >> timing.txt
    echo "  Window size: ${params.ld_window_size}" >> timing.txt
    echo "  Step size: ${params.ld_step_size}" >> timing.txt
    echo "  R2 threshold: ${params.ld_r2_threshold}" >> timing.txt
    echo "Input VCF stats:" >> timing.txt
    bcftools stats ${vcf} | grep "^SN" >> timing.txt
    """
}
