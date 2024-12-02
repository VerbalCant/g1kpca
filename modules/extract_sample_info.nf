process extract_sample_info {
    publishDir params.outdir, mode: 'copy'

    input:
    path ped_file

    output:
    path 'samples.pop', emit: sample_info

    script:
    """
    # Ensure clean tab-separated output with proper header
    echo -e "Individual_ID\tPopulation" > samples.pop
    awk -F'\t' 'NR > 1 {print \$2 "\t" \$7}' ${ped_file} | tr -s ' ' '\t' >> samples.pop
    """
} 