process create_pca_plot {
    container 'quay.io/biocontainers/r-base:4.2.2'
    publishDir params.outdir, mode: 'copy'

    input:
    path pca_results
    path sample_info

    output:
    path '*_PC*.pdf'
    path '*_debug.txt'

    script:
    """
    Rscript ${baseDir}/bin/plot_pca.R \
        --pca-file pca.eigenvec \
        --eigenval-file pca.eigenval \
        --sample-info ${sample_info} \
        --pop-info ${baseDir}/igsr_populations.tsv \
        --output-prefix pca_plot
    """
}
