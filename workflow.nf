#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.thousand_genomes_vcf = "/references/1000genomes/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
params.ancient_vcfs = "/project/1kg_nazca/data/ancient*.unifiedgenotyper.vcf.gz"
params.ped_file = "20130606_g1k.ped"
params.outdir = 'results'
params.min_maf = 0.01
params.max_missing = 0.1
params.ld_window_size = '50'
params.ld_step_size = '10'
params.ld_r2_threshold = '0.2'
params.max_cpus = 16
params.chromosomes = "1"  // Default to chr1, can be specified as "1,2,3" or "1-22"
params.timing_report = "${params.outdir}/timing_report.txt"
params.min_variants = 50
params.max_missing_rate = 0.1
params.min_samples = 1

// Helper function to expand chromosome list
def expandChromosomeList(chrStr) {
    def chrList = []
    chrStr.toString().split(',').each { range ->
        if (range.contains('-')) {
            def (start, end) = range.split('-')
            (start.toInteger()..end.toInteger()).each { chrList << it.toString() }
        } else {
            chrList << range.toString()
        }
    }
    return chrList
}

// Channel for ancient VCFs
ancient_ch = Channel.fromPath(params.ancient_vcfs)

// Create a process to extract sample information from 1000G VCF
process EXTRACT_SAMPLE_INFO {
    publishDir params.outdir, mode: 'copy'

    input:
    path ped_file

    output:
    path 'samples.pop', emit: sample_info

    script:
    """
    awk 'NR > 1 {print \$2 "\t" \$7}' ${ped_file} > samples.pop
    """
}

process EXTRACT_BIALLELIC_SNPS {
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    cpus params.max_cpus

    input:
    path vcf

    output:
    path "${vcf.baseName}.filtered.vcf.gz", emit: vcf

    script:
    """
    bcftools view --threads ${task.cpus} \
        -m2 -M2 -v snps \
        --min-af ${params.min_maf}:minor \
        ${vcf} \
        -Oz -o "${vcf.baseName}.filtered.vcf.gz"

    bcftools index --threads ${task.cpus} "${vcf.baseName}.filtered.vcf.gz"
    """
}

process MERGE_ANCIENT_SAMPLES {
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    cpus params.max_cpus

    input:
    path vcfs

    output:
    path 'merged_ancient.vcf.gz', emit: vcf
    path 'timing.txt', emit: timing

    script:
    """
    start_time=\$(date +%s)
    echo "Start time: \$(date)"
    
    # Create a directory for uniquely named links
    mkdir -p unique_vcfs
    
    # Create uniquely named symbolic links
    counter=1
    for vcf in *.vcf.gz; do
        unique_name="unique_vcfs/\${counter}_\${vcf}"
        ln -s "../\${vcf}" "\${unique_name}"
        echo "\${unique_name}" >> vcf_list.txt
        counter=\$((counter + 1))
    done

    # Index VCF files only if they don't have an index
    while read vcf; do
        if [ ! -f \$vcf.tbi ] && [ ! -f \$vcf.csi ]; then
            echo "Indexing \$vcf"
            bcftools index --threads ${task.cpus} \$vcf
        fi
    done < vcf_list.txt

    # Merge all ancient samples
    echo "Command: bcftools merge --threads ${task.cpus} --file-list vcf_list.txt --missing-to-ref -Oz -o merged_ancient.vcf.gz"
    bcftools merge \
        --threads ${task.cpus} \
        --file-list vcf_list.txt \
        --missing-to-ref \
        -Oz -o merged_ancient.vcf.gz

    bcftools index --threads ${task.cpus} merged_ancient.vcf.gz

    # End timing
    end_time=\$(date +%s)
    runtime=\$((end_time-start_time))
    
    # Write timing report
    echo "Process: MERGE_ANCIENT_SAMPLES" > timing.txt
    echo "Start time: \$(date -d @\${start_time})" >> timing.txt
    echo "End time: \$(date -d @\${end_time})" >> timing.txt
    echo "Runtime: \${runtime} seconds" >> timing.txt
    echo "Command: bcftools merge --threads ${task.cpus} --file-list vcf_list.txt --missing-to-ref -Oz -o merged_ancient.vcf.gz" >> timing.txt
    echo "Files merged:" >> timing.txt
    cat vcf_list.txt >> timing.txt
    echo "---" >> timing.txt
    """
}

process MERGE_WITH_1KG {
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    cpus params.max_cpus

    input:
    path thousand_genomes_vcf
    path ancient_vcf

    output:
    path 'merged_final.vcf.gz', emit: vcf
    path 'timing.txt', emit: timing

    script:
    """
    start_time=\$(date +%s)
    
    # Create local copy of 1000G VCF and index it
    echo "Creating local copy of 1000G VCF..."
    cp ${thousand_genomes_vcf} ./1kg.vcf.gz
    
    # Index both VCFs if needed
    echo "Indexing VCFs..."
    bcftools index --threads ${task.cpus} ./1kg.vcf.gz
    
    if [ ! -f ${ancient_vcf}.tbi ] && [ ! -f ${ancient_vcf}.csi ]; then
        bcftools index --threads ${task.cpus} ${ancient_vcf}
    fi

    # Merge with 1000 Genomes
    echo "Merging VCFs..."
    bcftools merge \
        --threads ${task.cpus} \
        --missing-to-ref \
        --regions ${params.chromosomes} \
        -Oz -o merged_final.vcf.gz \
        ./1kg.vcf.gz \
        ${ancient_vcf}

    bcftools index --threads ${task.cpus} merged_final.vcf.gz

    # End timing
    end_time=\$(date +%s)
    runtime=\$((end_time-start_time))
    
    # Write timing report
    echo "Process: MERGE_WITH_1KG" > timing.txt
    echo "Start time: \$(date -d @\${start_time})" >> timing.txt
    echo "End time: \$(date -d @\${end_time})" >> timing.txt
    echo "Runtime: \${runtime} seconds" >> timing.txt
    echo "Command: bcftools merge --threads ${task.cpus} --missing-to-ref --regions ${params.chromosomes} -Oz -o merged_final.vcf.gz" >> timing.txt
    echo "1000G VCF: ${thousand_genomes_vcf}" >> timing.txt
    echo "Ancient VCF: ${ancient_vcf}" >> timing.txt
    """
}

process LD_PRUNE {
    cpus params.max_cpus

    input:
    path vcf

    output:
    path 'pruned.vcf.gz', emit: vcf
    path 'timing.txt', emit: timing

    script:
    """
    start_time=\$(date +%s)
    echo "Starting LD pruning at \$(date)"

    # Index VCF if needed
    if [ ! -f ${vcf}.tbi ] && [ ! -f ${vcf}.csi ]; then
        echo "Indexing input VCF..."
        bcftools index --tbi ${vcf}
    fi

    # Convert to plink format and perform LD pruning
    echo "Step 1: LD pruning"
    echo "Input VCF: ${vcf}"
    ls -l ${vcf}
    
    # Filter out variants with long alleles using bcftools with expression
    echo "Filtering variants with long alleles..."
    bcftools view ${vcf} \
        --include 'strlen(REF)<=100 & strlen(ALT)<=100' \
        -Oz -o filtered_variants.vcf.gz

    bcftools index filtered_variants.vcf.gz

    # First pass: Convert to PLINK format with unique IDs
    echo "Converting to PLINK format with unique IDs..."
    plink2 --vcf filtered_variants.vcf.gz \
        --threads ${task.cpus} \
        --set-all-var-ids '@:#:\$r:\$a' \
        --new-id-max-allele-len 97 \
        --max-alleles 2 \
        --chr ${params.chromosomes} \
        --make-bed \
        --out temp_initial

    # Now run LD pruning on the bed file
    echo "Running LD pruning..."
    plink2 --bfile temp_initial \
        --threads ${task.cpus} \
        --indep-pairwise ${params.ld_window_size} ${params.ld_step_size} ${params.ld_r2_threshold} \
        --out pruned

    # Check if pruning was successful
    if [ ! -f pruned.prune.in ]; then
        echo "Error: LD pruning failed to produce output file"
        exit 1
    fi

    # Extract pruned variants
    echo "Step 2: Extracting pruned variants"
    plink2 --bfile temp_initial \
        --threads ${task.cpus} \
        --extract pruned.prune.in \
        --recode vcf bgz \
        --out pruned \
        --double-id

    # End timing
    end_time=\$(date +%s)
    runtime=\$((end_time-start_time))
    
    # Write timing report
    echo "Process: LD_PRUNE" > timing.txt
    echo "Start time: \$(date -d @\${start_time})" >> timing.txt
    echo "End time: \$(date -d @\${end_time})" >> timing.txt
    echo "Runtime: \${runtime} seconds" >> timing.txt
    echo "Parameters:" >> timing.txt
    echo "  Window size: ${params.ld_window_size}" >> timing.txt
    echo "  Step size: ${params.ld_step_size}" >> timing.txt
    echo "  R2 threshold: ${params.ld_r2_threshold}" >> timing.txt
    echo "Input VCF stats:" >> timing.txt
    bcftools stats ${vcf} | grep "^SN" >> timing.txt
    """
}

process RUN_PCA {
    publishDir params.outdir, mode: 'copy'
    cpus params.max_cpus

    input:
    path vcf

    output:
    path 'pca.*', emit: pca_results
    path 'debug_*', emit: debug

    script:
    """
    # First check input VCF
    echo "Checking input VCF..."
    bcftools stats ${vcf} > debug_vcf_stats.txt

    # Convert VCF to PLINK format with minimal filtering first
    echo "Converting to PLINK format..."
    plink2 --vcf ${vcf} \
        --threads ${task.cpus} \
        --make-bed \
        --allow-extra-chr \
        --set-missing-var-ids @:#,\\\$1,\\\$2 \
        --max-alleles 2 \
        --out temp_plink_initial

    # Check initial conversion
    echo "Initial PLINK conversion stats:" > debug_plink_steps.txt
    wc -l temp_plink_initial.bim >> debug_plink_steps.txt
    wc -l temp_plink_initial.fam >> debug_plink_steps.txt

    # Now do filtering steps one at a time
    echo "Filtering steps:" >> debug_plink_steps.txt

    # Step 1: Filter by missingness with a more lenient threshold
    plink2 --bfile temp_plink_initial \
        --threads ${task.cpus} \
        --geno 0.2 \
        --mind 0.2 \
        --make-bed \
        --out temp_plink_missing
    echo "After missingness filtering:" >> debug_plink_steps.txt
    wc -l temp_plink_missing.bim >> debug_plink_steps.txt

    # Step 2: Filter by MAF with a much lower threshold
    plink2 --bfile temp_plink_missing \
        --threads ${task.cpus} \
        --maf 0.001 \
        --make-bed \
        --out temp_plink_final
    echo "After MAF filtering:" >> debug_plink_steps.txt
    wc -l temp_plink_final.bim >> debug_plink_steps.txt

    # Check the number of variants and samples
    n_variants=\$(wc -l < temp_plink_final.bim)
    n_samples=\$(wc -l < temp_plink_final.fam)
    echo "Final counts:" >> debug_plink_steps.txt
    echo "Variants: \$n_variants" >> debug_plink_steps.txt
    echo "Samples: \$n_samples" >> debug_plink_steps.txt

    if [ \$n_variants -lt ${params.min_variants} ]; then
        echo "Error: Too few variants (\$n_variants) after filtering"
        exit 1
    fi

    if [ \$n_samples -lt ${params.min_samples} ]; then
        echo "Error: Too few samples (\$n_samples) after filtering"
        exit 1
    fi

    # Run PCA with allele weights for multiallelic variants
    plink2 --bfile temp_plink_final \
        --threads ${task.cpus} \
        --pca allele-wts approx \
        --allow-extra-chr \
        --out pca

    # Verify PCA output
    if [ ! -f pca.eigenvec ] || [ ! -f pca.eigenval ]; then
        echo "Error: PCA failed to produce output files"
        exit 1
    fi

    # Check eigenvalues
    head -n 10 pca.eigenval > debug_eigenvalues.txt

    # Check first few lines of eigenvectors
    head -n 10 pca.eigenvec > debug_eigenvectors.txt

    # Check for variation in eigenvectors
    Rscript -e '
        pca <- read.table("pca.eigenvec", header=TRUE)
        pc1_var <- var(pca[,3])
        pc2_var <- var(pca[,4])
        cat(sprintf("PC1 variance: %g\\n", pc1_var), file="debug_pc_variance.txt")
        cat(sprintf("PC2 variance: %g\\n", pc2_var), file="debug_pc_variance.txt")
        if (pc1_var < 1e-10 || pc2_var < 1e-10) {
            cat("Error: No variation in principal components\\n")
            quit(status=1)
        }
    '
    """
}

process CREATE_PCA_PLOT {
    container 'quay.io/biocontainers/r-base:4.2.2'
    publishDir params.outdir, mode: 'copy'

    input:
    path pca_results
    path sample_info

    output:
    path '*.pdf'
    path 'debug_info.txt'

    script:
    """
    Rscript ${baseDir}/bin/plot_pca.R pca.eigenvec ${sample_info}
    """
}

process EXTRACT_CHROMOSOMES {
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    cpus params.max_cpus

    input:
    path vcf
    val chromosome

    output:
    path "${vcf.baseName}.chr${chromosome}.vcf.gz", emit: vcf
    path "timing_${task.index}.txt", emit: timing

    script:
    """
    # Start timing
    start_time=\$(date +%s)
    echo "Start time: \$(date)"

    # Index input VCF only if no index exists
    if [ ! -f ${vcf}.tbi ] && [ ! -f ${vcf}.csi ]; then
        echo "Indexing input VCF..."
        bcftools index --threads ${task.cpus} ${vcf}
    fi

    # Extract chromosome
    echo "Command: bcftools view --threads ${task.cpus} ${vcf} --regions ${chromosome} -Oz -o ${vcf.baseName}.chr${chromosome}.vcf.gz"
    bcftools view --threads ${task.cpus} \
        ${vcf} \
        --regions ${chromosome} \
        -Oz -o "${vcf.baseName}.chr${chromosome}.vcf.gz"

    # Index output VCF only if no index exists
    if [ ! -f "${vcf.baseName}.chr${chromosome}.vcf.gz.tbi" ] && [ ! -f "${vcf.baseName}.chr${chromosome}.vcf.gz.csi" ]; then
        bcftools index --threads ${task.cpus} "${vcf.baseName}.chr${chromosome}.vcf.gz"
    fi

    # End timing
    end_time=\$(date +%s)
    runtime=\$((end_time-start_time))
    
    # Write timing report
    echo "Process: EXTRACT_CHROMOSOMES" > timing_${task.index}.txt
    echo "Input: ${vcf}" >> timing_${task.index}.txt
    echo "Chromosome: ${chromosome}" >> timing_${task.index}.txt
    echo "Start time: \$(date -d @\${start_time})" >> timing_${task.index}.txt
    echo "End time: \$(date -d @\${end_time})" >> timing_${task.index}.txt
    echo "Runtime: \${runtime} seconds" >> timing_${task.index}.txt
    echo "Commands:" >> timing_${task.index}.txt
    echo "  bcftools index --threads ${task.cpus} ${vcf}" >> timing_${task.index}.txt
    echo "  bcftools view --threads ${task.cpus} ${vcf} --regions ${chromosome} -Oz -o ${vcf.baseName}.chr${chromosome}.vcf.gz" >> timing_${task.index}.txt
    echo "---" >> timing_${task.index}.txt
    """
}

process COMBINE_TIMING_REPORTS {
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

// Create separate QC processes for each stage
process QC_CHECK_INPUT {
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

process QC_CHECK_EXTRACTED {
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

process QC_CHECK_ANCIENT_MERGED {
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
        echo "ERROR: Too few samples (\$samples < ${params.min_samples})" >> ${stage_name}_qc_report.txt
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

process QC_CHECK_FINAL_MERGED {
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
        echo "ERROR: Too few samples (\$samples < ${params.min_samples})" >> ${stage_name}_qc_report.txt
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

workflow {
    // Parse chromosome list and create debug view
    def chromosomes = Channel.fromList(
        expandChromosomeList(params.chromosomes)
    ).view { "DEBUG: Processing chromosome: $it" }
    
    // Create the 1000G VCF path with the correct chromosome
    def kg_vcf = chromosomes.map { chr ->
        def vcf_path = params.thousand_genomes_vcf.replace("{chr}", chr)
        [chr, file(vcf_path)]
    }.view { "DEBUG: 1000G VCF for chr${it[0]}: ${it[1]}" }

    // Create ancient VCFs channel with debug view
    def ancient_vcfs = Channel.fromPath(params.ancient_vcfs)
        .view { "DEBUG: Ancient VCF before combine: $it" }

    // Extract sample info from PED file
    EXTRACT_SAMPLE_INFO(
        Channel.fromPath("20130606_g1k.ped")
    )

    // QC and process ancient VCFs by chromosome
    QC_CHECK_INPUT(ancient_vcfs, "input_ancient_vcfs")
    
    def qc_results = QC_CHECK_INPUT.out.report
        .combine(ancient_vcfs)
        .branch {
            pass: !it[0].text.contains("ERROR")
            fail: it[0].text.contains("ERROR")
        }

    def combined_ch = qc_results.pass
        .map { report, vcf -> vcf }
        .combine(chromosomes)
        .view { "DEBUG: Combined channel: $it" }

    EXTRACT_CHROMOSOMES(
        combined_ch.map { vcf, chr -> vcf },
        combined_ch.map { vcf, chr -> chr }
    )
    
    QC_CHECK_EXTRACTED(EXTRACT_CHROMOSOMES.out.vcf, "post_extract_chromosomes")

    def grouped_vcfs = EXTRACT_CHROMOSOMES.out.vcf
        .map { vcf -> 
            def chr = vcf.name.find(/chr(\d+)/)
            [chr, vcf]
        }
        .groupTuple()
        .map { chr, files -> 
            def uniqueFiles = files.unique { it.name }
            [chr, uniqueFiles]
        }

    MERGE_ANCIENT_SAMPLES(
        grouped_vcfs.map { chr, files -> files }
    )

    QC_CHECK_ANCIENT_MERGED(MERGE_ANCIENT_SAMPLES.out.vcf, "post_merge_ancient")

    // Merge with 1000G by chromosome
    MERGE_WITH_1KG(
        kg_vcf.map { chr, vcf -> vcf },
        MERGE_ANCIENT_SAMPLES.out.vcf
    )

    QC_CHECK_FINAL_MERGED(MERGE_WITH_1KG.out.vcf, "post_1kg_merge")

    // Collect timing reports
    def timing_ch = EXTRACT_CHROMOSOMES.out.timing
        .mix(MERGE_ANCIENT_SAMPLES.out.timing)
        .mix(MERGE_WITH_1KG.out.timing)
        .collect()

    // Generate final timing report
    COMBINE_TIMING_REPORTS(timing_ch)

    // Collect all QC reports
    def qc_reports = Channel.empty()
        .mix(QC_CHECK_INPUT.out.report)
        .mix(QC_CHECK_EXTRACTED.out.report)
        .mix(QC_CHECK_ANCIENT_MERGED.out.report)
        .mix(QC_CHECK_FINAL_MERGED.out.report)
        .collect()
    
    // Add process to combine QC reports
    COMBINE_QC_REPORTS(qc_reports)

    // Perform LD pruning on merged data
    LD_PRUNE(MERGE_WITH_1KG.out.vcf)

    // Run PCA on pruned data
    RUN_PCA(LD_PRUNE.out.vcf)

    // Create PCA plot using the sample info from PED file
    CREATE_PCA_PLOT(
        RUN_PCA.out.pca_results,
        EXTRACT_SAMPLE_INFO.out.sample_info
    )
}

// Add new process to combine QC reports
process COMBINE_QC_REPORTS {
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
