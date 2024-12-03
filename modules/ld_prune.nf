process ld_prune {
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    label 'high_cpu'

    input:
    path vcf

    output:
    path "pruned.vcf.gz", emit: vcf
    path 'timing.txt', emit: timing

    script:
    def plink_mem_mb = (task.memory.toMega() * 0.8).intValue()
    """
    start_time=\$(date +%s)
    echo "[DEBUG] Starting LD pruning with ${task.cpus} CPUs and ${plink_mem_mb}MB memory"
    echo "[DEBUG] Total available memory: ${task.memory}"

    # Index VCF if needed
    if [ ! -f ${vcf}.tbi ] && [ ! -f ${vcf}.csi ]; then
        echo "[DEBUG] Indexing input VCF..."
        bcftools index --threads ${task.cpus} --tbi ${vcf}
    fi

    # Convert to plink format and perform LD pruning
    echo "[DEBUG] Step 1: LD pruning"
    echo "[DEBUG] Input VCF: ${vcf}"
    ls -l ${vcf}
    
    # Filter out variants with long alleles using bcftools with expression
    echo "[DEBUG] Filtering variants with long alleles..."
    bcftools view --threads ${task.cpus} ${vcf} \
        --include 'strlen(REF)<=100 & strlen(ALT)<=100' \
        -Oz -o filtered_variants.vcf.gz

    bcftools index --threads ${task.cpus} filtered_variants.vcf.gz

    # First pass: Convert to PLINK format with unique IDs
    echo "[DEBUG] Converting to PLINK format with unique IDs..."
    plink2 --vcf filtered_variants.vcf.gz \
        --threads ${task.cpus} \
        --memory ${plink_mem_mb} \
        --set-all-var-ids '@:#:\$r:\$a' \
        --new-id-max-allele-len 97 \
        --snps-only \
        --max-alleles 2 \
        --make-bed \
        --out temp_initial

    # Now run LD pruning on the bed file
    echo "[DEBUG] Running LD pruning..."
    plink2 --bfile temp_initial \
        --threads ${task.cpus} \
        --memory ${plink_mem_mb} \
        --snps-only \
        --indep-pairwise ${params.ld_window_size} ${params.ld_step_size} ${params.ld_r2_threshold} \
        --out pruned

    # Check if pruning was successful
    if [ ! -f pruned.prune.in ]; then
        echo "Error: LD pruning failed to produce output file"
        exit 1
    fi

    # Extract pruned variants
    echo "[DEBUG] Step 2: Extracting pruned variants"
    plink2 --bfile temp_initial \
        --threads ${task.cpus} \
        --memory ${plink_mem_mb} \
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
    echo "Resources:" >> timing.txt
    echo "  CPUs: ${task.cpus}" >> timing.txt
    echo "  Memory: ${task.memory}" >> timing.txt
    echo "  PLINK memory: ${plink_mem_mb}MB" >> timing.txt
    echo "Parameters:" >> timing.txt
    echo "  Window size: ${params.ld_window_size}" >> timing.txt
    echo "  Step size: ${params.ld_step_size}" >> timing.txt
    echo "  R2 threshold: ${params.ld_r2_threshold}" >> timing.txt
    echo "Input VCF stats:" >> timing.txt
    bcftools stats ${vcf} | grep "^SN" >> timing.txt
    """
}
