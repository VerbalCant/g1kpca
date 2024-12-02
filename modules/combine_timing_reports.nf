process combine_timing_reports {
    publishDir params.outdir, mode: 'copy'

    input:
    path 'timing_*.txt'

    output:
    path 'final_timing_report.txt'

    script:
    """
    echo "Workflow Timing Report" > final_timing_report.txt
    echo "Generated: \$(date)" >> final_timing_report.txt
    echo "---" >> final_timing_report.txt
    cat timing_*.txt >> final_timing_report.txt
    """
}
