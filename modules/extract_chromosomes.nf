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
    """
    # Start timing
    start_time=\$(date +%s)
    echo "Start time: \$(date)"

    # Index input VCF only if no index exists
    if [ ! -f ${vcf}.tbi ] && [ ! -f ${vcf}.csi ]; then
        echo "Indexing input VCF..."
        bcftools index --threads ${task.cpus} ${vcf}
    fi

    # Extract chromosome
    echo "Command: bcftools view --threads ${task.cpus} ${vcf} --regions ${chromosome} -Oz -o ${vcf.baseName}.chr${chromosome}.vcf.gz"
    bcftools view --threads ${task.cpus} \
        ${vcf} \
        --regions ${chromosome} \
        -Oz -o "${vcf.baseName}.chr${chromosome}.vcf.gz"

    # Index output VCF only if no index exists
    if [ ! -f "${vcf.baseName}.chr${chromosome}.vcf.gz.tbi" ] && [ ! -f "${vcf.baseName}.chr${chromosome}.vcf.gz.csi" ]; then
        bcftools index --threads ${task.cpus} "${vcf.baseName}.chr${chromosome}.vcf.gz"
    fi

    # End timing
    end_time=\$(date +%s)
    runtime=\$((end_time-start_time))
    
    # Write timing report
    echo "Process: EXTRACT_CHROMOSOMES" > timing_${task.index}.txt
    echo "Input: ${vcf}" >> timing_${task.index}.txt
    echo "Chromosome: ${chromosome}" >> timing_${task.index}.txt
    echo "Start time: \$(date -d @\${start_time})" >> timing_${task.index}.txt
    echo "End time: \$(date -d @\${end_time})" >> timing_${task.index}.txt
    echo "Runtime: \${runtime} seconds" >> timing_${task.index}.txt
    echo "Commands:" >> timing_${task.index}.txt
    echo "  bcftools index --threads ${task.cpus} ${vcf}" >> timing_${task.index}.txt
    echo "  bcftools view --threads ${task.cpus} ${vcf} --regions ${chromosome} -Oz -o ${vcf.baseName}.chr${chromosome}.vcf.gz" >> timing_${task.index}.txt
    echo "---" >> timing_${task.index}.txt
    """
} 