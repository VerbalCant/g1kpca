process merge_samples {
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    cpus params.max_cpus

    input:
    path vcfs

    output:
    path 'merged_ancient.vcf.gz', emit: vcf
    path 'timing.txt', emit: timing

    script:
    """
    start_time=\$(date +%s)
    echo "Start time: \$(date)"
    
    # Create a directory for uniquely named links
    mkdir -p unique_vcfs
    
    # Create uniquely named symbolic links
    counter=1
    for vcf in *.vcf.gz; do
        unique_name="unique_vcfs/\${counter}_\${vcf}"
        ln -s "../\${vcf}" "\${unique_name}"
        echo "\${unique_name}" >> vcf_list.txt
        counter=\$((counter + 1))
    done

    # Index VCF files only if they don't have an index
    while read vcf; do
        if [ ! -f \$vcf.tbi ] && [ ! -f \$vcf.csi ]; then
            echo "Indexing \$vcf"
            bcftools index --threads ${task.cpus} \$vcf
        fi
    done < vcf_list.txt

    # Merge all ancient samples
    echo "Command: bcftools merge --threads ${task.cpus} --file-list vcf_list.txt --missing-to-ref -Oz -o merged_ancient.vcf.gz"
    bcftools merge \
        --threads ${task.cpus} \
        --file-list vcf_list.txt \
        --missing-to-ref \
        -Oz -o merged_ancient.vcf.gz

    bcftools index --threads ${task.cpus} merged_ancient.vcf.gz

    # End timing
    end_time=\$(date +%s)
    runtime=\$((end_time-start_time))
    
    # Write timing report
    echo "Process: MERGE_ANCIENT_SAMPLES" > timing.txt
    echo "Start time: \$(date -d @\${start_time})" >> timing.txt
    echo "End time: \$(date -d @\${end_time})" >> timing.txt
    echo "Runtime: \${runtime} seconds" >> timing.txt
    echo "Command: bcftools merge --threads ${task.cpus} --file-list vcf_list.txt --missing-to-ref -Oz -o merged_ancient.vcf.gz" >> timing.txt
    echo "Files merged:" >> timing.txt
    cat vcf_list.txt >> timing.txt
    echo "---" >> timing.txt
    """
} 