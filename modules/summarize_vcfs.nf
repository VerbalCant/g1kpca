process summarize_vcfs {
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    publishDir params.outdir, mode: 'copy'

    input:
    path vcfs
    val stage_name

    output:
    path "${stage_name}_vcf_summary.txt", emit: summary

    script:
    """
    echo "VCF Summary Report - ${stage_name}" > ${stage_name}_vcf_summary.txt
    echo "Generated at: \$(date)" >> ${stage_name}_vcf_summary.txt
    echo "----------------------------------------" >> ${stage_name}_vcf_summary.txt

    for vcf in ${vcfs}; do
        echo "\\nFile: \$vcf" >> ${stage_name}_vcf_summary.txt
        echo "----------------------------------------" >> ${stage_name}_vcf_summary.txt
        
        # Index VCF if needed
        if [ ! -f \$vcf.tbi ] && [ ! -f \$vcf.csi ]; then
            bcftools index --tbi \$vcf
        fi

        # Get chromosomes and variant counts
        echo "Chromosomes and variant counts:" >> ${stage_name}_vcf_summary.txt
        bcftools index --stats \$vcf | while read chr variants rest; do
            echo "  Chromosome \$chr: \$variants variants" >> ${stage_name}_vcf_summary.txt
        done

        # Get total variant count
        total_variants=\$(bcftools view -H \$vcf | wc -l)
        echo "Total variants: \$total_variants" >> ${stage_name}_vcf_summary.txt

        # Get sample count
        sample_count=\$(bcftools query -l \$vcf | wc -l)
        echo "Number of samples: \$sample_count" >> ${stage_name}_vcf_summary.txt
        echo "----------------------------------------" >> ${stage_name}_vcf_summary.txt
    done
    """
} 