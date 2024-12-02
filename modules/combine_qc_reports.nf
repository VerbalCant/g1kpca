process combine_qc_reports {
    publishDir params.outdir, mode: 'copy'

    input:
    path '*_qc_report.txt'

    output:
    path 'final_qc_report.txt'

    script:
    """
    # Create header in a temporary file
    echo "Combined QC Report" > temp_report.txt
    echo "Generated: \$(date)" >> temp_report.txt
    echo "=========================" >> temp_report.txt
    
    # Append each QC report to the temporary file
    for report in *_qc_report.txt; do
        echo "" >> temp_report.txt
        cat \$report >> temp_report.txt
        echo "-------------------------" >> temp_report.txt
    done

    # Move temporary file to final output
    mv temp_report.txt final_qc_report.txt
    """
}
