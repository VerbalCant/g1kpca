process qc_check_input {
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    cpus 1

    input:
    path vcf
    val stage_name

    output:
    path "${stage_name}_qc_report.txt", emit: report
    path "FAIL", optional: true, emit: fail_flag

    script:
    """
    # Start QC report
    echo "QC Report for ${stage_name}" > ${stage_name}_qc_report.txt
    echo "VCF file: ${vcf}" >> ${stage_name}_qc_report.txt
    echo "Date: \$(date)" >> ${stage_name}_qc_report.txt
    echo "-------------------" >> ${stage_name}_qc_report.txt

    # Index VCF only if no index exists
    if [ ! -f ${vcf}.tbi ] && [ ! -f ${vcf}.csi ]; then
        bcftools index --tbi ${vcf}
    fi

    # Check if file exists and is not empty
    if [ ! -s ${vcf} ]; then
        echo "ERROR: VCF file is empty" >> ${stage_name}_qc_report.txt
        touch FAIL
        exit 0
    fi

    # Get basic stats
    echo "Basic Statistics:" >> ${stage_name}_qc_report.txt
    bcftools stats ${vcf} | grep "^SN" >> ${stage_name}_qc_report.txt

    # Get number of variants
    variants=\$(bcftools view -H ${vcf} | wc -l)
    echo "Number of variants: \$variants" >> ${stage_name}_qc_report.txt

    # Get number of samples
    samples=\$(bcftools query -l ${vcf} | wc -l)
    echo "Number of samples: \$samples" >> ${stage_name}_qc_report.txt

    # Check for minimum number of variants
    if [ \$variants -lt ${params.min_variants} ]; then
        echo "ERROR: Too few variants (\$variants < ${params.min_variants})" >> ${stage_name}_qc_report.txt
        touch FAIL
        exit 0
    fi

    # Check for minimum number of samples
    if [ \$samples -lt ${params.min_samples} ]; then
        echo "ERROR: No samples found in VCF" >> ${stage_name}_qc_report.txt
        touch FAIL
        exit 0
    fi

    # Check for biallelic sites
    multiallelic=\$(bcftools view -H ${vcf} | cut -f 4,5 | grep ',' | wc -l)
    echo "Number of multiallelic sites: \$multiallelic" >> ${stage_name}_qc_report.txt

    # Check for variant types
    echo "Variant types:" >> ${stage_name}_qc_report.txt
    bcftools view -H ${vcf} | cut -f 4,5 | awk '{len1=length(\$1); len2=length(\$2); if(len1==1 && len2==1) print "SNP"; else if(len1>len2) print "DEL"; else if(len1<len2) print "INS"; else print "OTHER"}' | sort | uniq -c >> ${stage_name}_qc_report.txt

    echo "QC check completed successfully" >> ${stage_name}_qc_report.txt
    """
} 