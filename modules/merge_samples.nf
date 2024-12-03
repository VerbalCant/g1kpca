process merge_samples {
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    label 'high_cpu'

    input:
    tuple val(chr), path(vcfs)

    output:
    path "merged_ancient.chr${chr}.vcf.gz", emit: vcf
    path 'timing.txt', emit: timing

    script:
    """
    start_time=\$(date +%s)
    echo "[DEBUG] Starting merge_samples for chromosome ${chr} with files:"
    ls -l *.vcf.gz
    
    # Create a directory for uniquely named links
    mkdir -p unique_vcfs
    
    # Create uniquely named symbolic links with debug output
    counter=1
    for vcf in *.vcf.gz; do
        unique_name="unique_vcfs/\${counter}_\${vcf}"
        ln -s "../\${vcf}" "\${unique_name}"
        echo "[DEBUG] Created link: \${unique_name} -> \${vcf}"
        echo "\${unique_name}" >> vcf_list.txt
        counter=\$((counter + 1))
    done

    echo "[DEBUG] Contents of vcf_list.txt:"
    cat vcf_list.txt

    # Index and check chromosome content of each VCF
    while read vcf; do
        echo "[DEBUG] Processing \$vcf"
        if [ ! -f \$vcf.tbi ] && [ ! -f \$vcf.csi ]; then
            bcftools index --threads ${task.cpus} \$vcf
        fi
        echo "[DEBUG] Chromosomes in \$vcf:"
        bcftools index --stats \$vcf | cut -f 1
        echo "[DEBUG] Samples in \$vcf:"
        bcftools query -l \$vcf
    done < vcf_list.txt

    # Merge all ancient samples
    echo "[DEBUG] Starting merge operation"
    bcftools merge \
        --threads ${task.cpus} \
        --force-samples \
        --file-list vcf_list.txt \
        --missing-to-ref \
        -Oz -o merged_ancient.chr${chr}.vcf.gz

    # Verify merged output
    echo "[DEBUG] Chromosomes in merged output:"
    bcftools index --stats merged_ancient.chr${chr}.vcf.gz | cut -f 1
    echo "[DEBUG] Samples in merged output:"
    bcftools query -l merged_ancient.chr${chr}.vcf.gz

    bcftools index --threads ${task.cpus} merged_ancient.chr${chr}.vcf.gz

    # End timing
    end_time=\$(date +%s)
    runtime=\$((end_time-start_time))
    
    # Write timing report
    echo "Process: MERGE_ANCIENT_SAMPLES" > timing.txt
    echo "Chromosome: ${chr}" >> timing.txt
    echo "Start time: \$(date -d @\${start_time})" >> timing.txt
    echo "End time: \$(date -d @\${end_time})" >> timing.txt
    echo "Runtime: \${runtime} seconds" >> timing.txt
    echo "Command: bcftools merge --threads ${task.cpus} --force-samples --file-list vcf_list.txt --missing-to-ref -Oz -o merged_ancient.chr${chr}.vcf.gz" >> timing.txt
    echo "Files merged:" >> timing.txt
    cat vcf_list.txt >> timing.txt
    echo "---" >> timing.txt
    """
} 