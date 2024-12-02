process merge_with_1000genomes {
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    cpus params.max_cpus

    input:
    path thousand_genomes_vcf
    path ancient_vcf

    output:
    path 'merged_final.vcf.gz', emit: vcf
    path 'timing.txt', emit: timing

    script:
    """
    start_time=\$(date +%s)
    
    # Create local copy of 1000G VCF and index it
    echo "Creating local copy of 1000G VCF..."
    cp ${thousand_genomes_vcf} ./1kg.vcf.gz
    
    # Index both VCFs if needed
    echo "Indexing VCFs..."
    bcftools index --threads ${task.cpus} ./1kg.vcf.gz
    
    if [ ! -f ${ancient_vcf}.tbi ] && [ ! -f ${ancient_vcf}.csi ]; then
        bcftools index --threads ${task.cpus} ${ancient_vcf}
    fi

    # Merge with 1000 Genomes
    echo "Merging VCFs..."
    bcftools merge \
        --threads ${task.cpus} \
        --missing-to-ref \
        --regions ${params.chromosomes} \
        -Oz -o merged_final.vcf.gz \
        ./1kg.vcf.gz \
        ${ancient_vcf}

    bcftools index --threads ${task.cpus} merged_final.vcf.gz

    # End timing
    end_time=\$(date +%s)
    runtime=\$((end_time-start_time))
    
    # Write timing report
    echo "Process: MERGE_WITH_1KG" > timing.txt
    echo "Start time: \$(date -d @\${start_time})" >> timing.txt
    echo "End time: \$(date -d @\${end_time})" >> timing.txt
    echo "Runtime: \${runtime} seconds" >> timing.txt
    echo "Command: bcftools merge --threads ${task.cpus} --missing-to-ref --regions ${params.chromosomes} -Oz -o merged_final.vcf.gz" >> timing.txt
    echo "1000G VCF: ${thousand_genomes_vcf}" >> timing.txt
    echo "Ancient VCF: ${ancient_vcf}" >> timing.txt
    """
}
