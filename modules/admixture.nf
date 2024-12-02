


process RUN_ADMIXTURE {
    publishDir "${params.outdir}/admixture", mode: 'copy'
    cpus params.admixture_threads
    
    input:
    path(bed_file)
    path(bim_file)
    path(fam_file)
    val(k)
    
    output:
    path("admixture_input.${k}.Q"), emit: q_files
    path("admixture_input.${k}.P"), emit: p_files
    path("log${k}.out"), emit: log
    path("admixture_input.${k}.cv"), optional: true, emit: cv_error
    
    script:
    def cv_flag = params.admixture_cv ? "--cv" : ""
    """
    # Run ADMIXTURE with cross-validation
    admixture ${cv_flag} -j${task.cpus} -s ${k} ${bed_file} ${k} | tee log${k}.out
    
    # Debug output
    echo "Contents of directory:"
    ls -l
    
    # Check if Q file exists
    if [ ! -f "admixture_input.${k}.Q" ]; then
        echo "ERROR: Q file not generated"
        exit 1
    fi
    """
}

process PREPARE_ADMIXTURE_INPUT {
    publishDir "${params.outdir}/admixture", mode: 'copy'
    
    input:
    path(vcf)
    
    output:
    tuple path("admixture_input.bed"), 
          path("admixture_input.bim"), 
          path("admixture_input.fam"), emit: plink_files
    
    script:
    """
    # Convert VCF to PLINK format
    plink2 --vcf ${vcf} \
        --make-bed \
        --out admixture_input \
        --double-id \
        --set-missing-var-ids '@:#:\\\$1,\\\$2' \
        --new-id-max-allele-len 100 \
        --max-alleles 2 \
        --vcf-half-call m \
        --allow-extra-chr \
        --keep-allele-order

    # Check if ancient samples are present
    echo "Checking for ancient samples in output:"
    grep "^ancient" admixture_input.fam || echo "No ancient samples found"
    
    # Count total samples
    echo "Total samples in FAM file:"
    wc -l admixture_input.fam
    """
}

process PLOT_ADMIXTURE {
    publishDir "${params.outdir}/admixture", mode: 'copy'
    
    input:
    path(q_files)
    path(sample_info)
    path(pop_info)
    
    output:
    path("admixture_plots*.png")
    path("admixture_data.json")
    
    script:
    """
    # Debug: List input files
    echo "Q files:"
    ls -l *.Q
    
    echo "Sample info:"
    ls -l ${sample_info}
    
    echo "Population info:"
    ls -l ${pop_info}
    
    Rscript ${baseDir}/bin/plot_admixture.R \
        --qfiles . \
        --sample-info ${sample_info} \
        --pop-info ${pop_info} \
        --output-prefix admixture_plots
    """
}

process ANALYZE_CV_ERROR {
    publishDir "${params.outdir}/admixture", mode: 'copy'
    
    input:
    path('*.cv')
    
    output:
    path("cv_error_plot.pdf")
    path("optimal_k.txt")
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Read CV error files
    cv_files <- list.files(pattern="*.cv")
    cv_data <- lapply(cv_files, function(f) {
        k <- as.numeric(sub(".*\\\\.(.*)\\.cv", "\\\\1", f))
        cv_error <- as.numeric(readLines(f))
        data.frame(K=k, CV_Error=cv_error)
    })
    cv_data <- do.call(rbind, cv_data)
    
    # Plot CV error
    pdf("cv_error_plot.pdf")
    plot(cv_data\$K, cv_data\$CV_Error, type="b",
         xlab="K (number of ancestral populations)",
         ylab="Cross-validation error",
         main="Cross-validation error by K")
    dev.off()
    
    # Find optimal K
    optimal_k <- cv_data\$K[which.min(cv_data\$CV_Error)]
    write(optimal_k, "optimal_k.txt")
    """
}

process CHECK_ADMIXTURE_INPUT {
    publishDir "${params.outdir}/admixture/qc", mode: 'copy'
    
    input:
    tuple path(bed), path(bim), path(fam)
    
    output:
    path "admixture_input_qc.txt"
    
    script:
    """
    echo "ADMIXTURE Input QC Report" > admixture_input_qc.txt
    echo "=========================" >> admixture_input_qc.txt
    
    echo -e "\nSample counts:" >> admixture_input_qc.txt
    echo "Total samples: \$(wc -l < ${fam})" >> admixture_input_qc.txt
    echo "Ancient samples: \$(grep -c '^ancient' ${fam})" >> admixture_input_qc.txt
    echo "Modern samples: \$(grep -vc '^ancient' ${fam})" >> admixture_input_qc.txt
    
    echo -e "\nVariant counts:" >> admixture_input_qc.txt
    echo "Total variants: \$(wc -l < ${bim})" >> admixture_input_qc.txt
    
    echo -e "\nFirst few samples:" >> admixture_input_qc.txt
    head -n 5 ${fam} >> admixture_input_qc.txt
    
    echo -e "\nLast few samples:" >> admixture_input_qc.txt
    tail -n 5 ${fam} >> admixture_input_qc.txt
    """
}