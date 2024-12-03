process extract_chromosomes {
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    cpus params.max_cpus

    input:
    path vcf
    val chromosome

    output:
    path "${vcf.baseName}.chr${chromosome}.vcf.gz", emit: vcf
    path "timing_${task.index}.txt", emit: timing

    script:
    
    def queryChrom = chromosome.toString().startsWith('chr') ? chromosome.toString().substring(3) : chromosome.toString()
    """
    # Start timing
    start_time=\$(date +%s)
    echo "[DEBUG] Starting chromosome extraction process"
    echo "[DEBUG] Input VCF: ${vcf}"
    echo "[DEBUG] Target chromosome: ${chromosome}"
    echo "[DEBUG] Query chromosome format: ${queryChrom}"

    # First ensure input VCF is indexed
    echo "[DEBUG] Creating index for input VCF"
    bcftools index --tbi --force ${vcf}

    # List available chromosomes in input
    echo "[DEBUG] Available chromosomes in input VCF:"
    bcftools index --stats ${vcf}

    # Extract chromosome (using the converted chromosome format)
    echo "[DEBUG] Extracting chromosome ${queryChrom}"
    bcftools view \
        --threads ${task.cpus} \
        --regions ${queryChrom} \
        --output-type z \
        --output "${vcf.baseName}.chr${chromosome}.vcf.gz" \
        ${vcf}

    # Index output file
    echo "[DEBUG] Indexing output VCF"
    bcftools index --tbi "${vcf.baseName}.chr${chromosome}.vcf.gz"

    # Verify output contents
    echo "[DEBUG] Verifying output VCF contents:"
    bcftools query -l "${vcf.baseName}.chr${chromosome}.vcf.gz" | wc -l
    bcftools index --stats "${vcf.baseName}.chr${chromosome}.vcf.gz"

    # End timing
    end_time=\$(date +%s)
    runtime=\$((end_time-start_time))
    
    # Write timing report
    echo "Process: EXTRACT_CHROMOSOMES" > timing_${task.index}.txt
    echo "Input: ${vcf}" >> timing_${task.index}.txt
    echo "Chromosome: ${queryChrom}" >> timing_${task.index}.txt
    echo "Start time: \$(date -d @\${start_time})" >> timing_${task.index}.txt
    echo "End time: \$(date -d @\${end_time})" >> timing_${task.index}.txt
    echo "Runtime: \${runtime} seconds" >> timing_${task.index}.txt
    """
} 