process create_pca_data_file {
    publishDir params.outdir, mode: 'copy'

    input:
    path(pca_results)
    path(sample_info)
    path(pop_info)

    output:
    path 'pca_data.json'

    script:
    """
    # Get the eigenvec and eigenval files from pca_results
    pca_eigenvec=\$(find . -name "*.eigenvec")
    pca_eigenval=\$(find . -name "*.eigenval")

    Rscript ${baseDir}/bin/create_pca_json.R \
        \$pca_eigenvec \
        \$pca_eigenval \
        ${sample_info} \
        ${pop_info}
    """
}
