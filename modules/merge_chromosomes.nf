process merge_chromosomes {
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    cpus params.max_cpus

    input:
    path vcfs  // Collection of per-chromosome VCFs

    output:
    path "merged_all_chr.vcf.gz", emit: vcf
    path "timing.txt", emit: timing

    script:
    """
    start_time=\$(date +%s)
    echo "[DEBUG] Starting chromosome merge with files:"
    ls -l *.vcf.gz
    
    # Create list of VCFs to merge
    for vcf in *.vcf.gz; do
        echo "\$vcf" >> vcf_list.txt
    done

    echo "[DEBUG] Contents of vcf_list.txt:"
    cat vcf_list.txt

    # Index all VCFs
    while read vcf; do
        echo "[DEBUG] Indexing \$vcf"
        bcftools index --threads ${task.cpus} \$vcf
    done < vcf_list.txt

    # Merge all chromosomes
    echo "[DEBUG] Merging all chromosomes"
    bcftools concat \
        --threads ${task.cpus} \
        --file-list vcf_list.txt \
        --allow-overlaps \
        -Oz -o merged_all_chr.vcf.gz

    # Index final output
    bcftools index --threads ${task.cpus} merged_all_chr.vcf.gz

    # Verify merged output
    echo "[DEBUG] Chromosomes in final merged output:"
    bcftools index --stats merged_all_chr.vcf.gz

    # End timing
    end_time=\$(date +%s)
    runtime=\$((end_time-start_time))
    
    # Write timing report
    echo "Process: MERGE_CHROMOSOMES" > timing.txt
    echo "Start time: \$(date -d @\${start_time})" >> timing.txt
    echo "End time: \$(date -d @\${end_time})" >> timing.txt
    echo "Runtime: \${runtime} seconds" >> timing.txt
    echo "Files merged:" >> timing.txt
    cat vcf_list.txt >> timing.txt
    """
} 