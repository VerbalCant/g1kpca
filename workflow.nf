#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define parameters
params.thousand_genomes_vcf = "/references/1000genomes/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
params.ancient_vcfs = "/project/1kg_nazca/data/ancient*.unifiedgenotyper.vcf.gz"
params.outdir = 'results'
params.min_maf = 0.05
params.max_missing = 0.05
params.ld_window_size = '50'
params.ld_step_size = '10'
params.ld_r2_threshold = '0.2'
params.max_cpus = 16
params.chromosomes = "1"  // Default to chr1, can be specified as "1,2,3" or "1-22"
params.timing_report = "${params.outdir}/timing_report.txt"
params.min_variants = 100
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
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    publishDir params.outdir, mode: 'copy'
    cpus params.max_cpus

    input:
    path vcf

    output:
    path 'samples.pop', emit: sample_info

    script:
    """
    # Get sample IDs
    bcftools query -l ${vcf} > samples.txt

    # Create population assignments based on sample naming convention
    while read sample; do
        if [[ \$sample =~ ^ancient ]]; then
            echo -e "\${sample}\tANCIENT"
        elif [[ \$sample =~ ^HG ]]; then
            pop=\$(echo \$sample | cut -c1-5)
            echo -e "\${sample}\t\${pop}"
        elif [[ \$sample =~ ^NA ]]; then
            pop=\$(echo \$sample | cut -c1-4)
            echo -e "\${sample}\t\${pop}"
        else
            echo -e "\${sample}\tUNKNOWN"
        fi
    done < samples.txt > samples.pop
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
    
    # Create a list of VCF files with unique names
    for vcf in *.vcf.gz; do
        # Create symbolic links with unique names based on the full path
        ln -s \$vcf \$(echo \$vcf | md5sum | cut -d' ' -f1).vcf.gz
    done
    
    # Create a list of the uniquely named VCF files
    ls -1 *.md5sum.vcf.gz > vcf_list.txt

    # Index VCF files only if they don't have an index
    for vcf in *.md5sum.vcf.gz; do
        if [ ! -f \$vcf.tbi ] && [ ! -f \$vcf.csi ]; then
            echo "Indexing \$vcf"
            bcftools index --threads ${task.cpus} \$vcf
        fi
    done

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

    script:
    """
    # Merge with 1000 Genomes
    bcftools merge \
        --threads ${task.cpus} \
        --missing-to-ref \
        --regions 1 \
        -Oz -o merged_final.vcf.gz \
        ${thousand_genomes_vcf} \
        ${ancient_vcf}

    bcftools index --threads ${task.cpus} merged_final.vcf.gz
    """
}

process LD_PRUNE {
    container 'quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0'
    cpus params.max_cpus

    input:
    path vcf

    output:
    path 'pruned.vcf.gz', emit: vcf

    script:
    """
    # Convert to plink format and perform LD pruning
    plink2 --vcf ${vcf} \
        --threads ${task.cpus} \
        --indep-pairwise ${params.ld_window_size} ${params.ld_step_size} ${params.ld_r2_threshold} \
        --allow-extra-chr \
        --set-missing-var-ids @:#,\\\$1,\\\$2 \
        --max-alleles 2 \
        --out pruned

    # Extract pruned variants
    plink2 --vcf ${vcf} \
        --threads ${task.cpus} \
        --extract pruned.prune.in \
        --recode vcf bgz \
        --out pruned \
        --allow-extra-chr \
        --set-missing-var-ids @:#,\\\$1,\\\$2 \
        --max-alleles 2
    """
}

process RUN_PCA {
    container 'quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0'
    publishDir params.outdir, mode: 'copy'
    cpus params.max_cpus

    input:
    path vcf

    output:
    path 'pca.*', emit: pca_results

    script:
    """
    plink2 --vcf ${vcf} \
        --threads ${task.cpus} \
        --pca approx \
        --allow-extra-chr \
        --set-missing-var-ids @:#,\\\$1,\\\$2 \
        --out pca
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

    script:
    """
    #!/usr/bin/env Rscript

    # Read PCA results and sample information
    pca <- read.table("pca.eigenvec", header=FALSE)
    samples <- read.table("${sample_info}", header=FALSE)

    # Merge population information
    pca\$Population <- samples\$V2[match(pca\$V2, samples\$V1)]

    # Create plot
    pdf("pca_plot.pdf", width=10, height=8)

    # Define colors for different populations
    colors <- c(
        "ANCIENT"="red",
        "HG00"="blue",
        "HG01"="green",
        "HG02"="purple",
        "NA18"="orange",
        "NA19"="brown",
        "NA20"="cyan",
        "UNKNOWN"="gray"
    )

    # Plot PC1 vs PC2
    plot(pca\$V3, pca\$V4,
         col=colors[substr(pca\$Population, 1, 4)],
         pch=20,
         xlab="PC1",
         ylab="PC2",
         main="PCA of Ancient and 1000G Samples")

    # Add legend
    legend("topright",
           legend=names(colors),
           col=colors,
           pch=20,
           title="Population",
           cex=0.8)

    # Add special labels for ancient samples
    ancient_idx <- which(pca\$Population == "ANCIENT")
    text(pca\$V3[ancient_idx],
         pca\$V4[ancient_idx],
         labels=pca\$V2[ancient_idx],
         pos=3,
         cex=0.6)

    dev.off()
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

process QC_CHECK_MERGED {
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
    # Same script as QC_CHECK_INPUT
    """
}

workflow {
    // Parse chromosome list and create debug view
    def chromosomes = Channel.fromList(
        expandChromosomeList(params.chromosomes)
    ).view { "DEBUG: Chromosome before combine: $it" }
    
    // Create ancient VCFs channel with debug view
    def ancient_vcfs = Channel.fromPath(params.ancient_vcfs)
        .view { "DEBUG: Ancient VCF before combine: $it" }

    // QC check on input ancient VCFs
    QC_CHECK_INPUT(ancient_vcfs, "input_ancient_vcfs")
    
    // Branch based on QC results
    QC_CHECK_INPUT.out.report
        .combine(ancient_vcfs)
        .branch {
            pass: !it[0].text.contains("ERROR")
            fail: it[0].text.contains("ERROR")
        }
        .set { qc_results }

    // Create the combined channel with debug view
    def combined_ch = qc_results.pass
        .map { report, vcf -> vcf }
        .combine(chromosomes)
        .view { "DEBUG: Combined channel: $it" }

    // Extract chromosomes from ancient VCFs
    EXTRACT_CHROMOSOMES(
        combined_ch.map { vcf, chr -> vcf }.view { "DEBUG: VCF input: $it" },
        combined_ch.map { vcf, chr -> chr }.view { "DEBUG: Chr input: $it" }
    )
    
    // QC check on extracted chromosomes
    QC_CHECK_EXTRACTED(EXTRACT_CHROMOSOMES.out.vcf, "post_extract_chromosomes")

    // Group extracted VCFs by chromosome and collect them, ensuring uniqueness
    def grouped_vcfs = EXTRACT_CHROMOSOMES.out.vcf
        .map { vcf -> 
            def chr = vcf.name.find(/chr(\d+)/)
            [chr, vcf]
        }
        .groupTuple()
        .map { chr, files -> 
            // Deduplicate files based on their base names
            def uniqueFiles = files.unique { it.name }
            [chr, uniqueFiles]
        }
        .view { "DEBUG: Grouped VCFs (deduplicated): $it" }

    // Merge ancient samples by chromosome
    MERGE_ANCIENT_SAMPLES(
        grouped_vcfs.map { chr, files -> files }
    )

    // QC check on merged ancient samples
    QC_CHECK_MERGED(MERGE_ANCIENT_SAMPLES.out.vcf, "post_merge_ancient")

    // Collect timing reports
    def timing_ch = EXTRACT_CHROMOSOMES.out.timing
        .mix(MERGE_ANCIENT_SAMPLES.out.timing)
        .collect()
        .view { "DEBUG: Collected timing reports: $it" }

    // Generate final timing report
    COMBINE_TIMING_REPORTS(timing_ch)

    // Collect all QC reports
    def qc_reports = Channel.empty()
        .mix(QC_CHECK_INPUT.out.report)
        .mix(QC_CHECK_EXTRACTED.out.report)
        .mix(QC_CHECK_MERGED.out.report)
        .collect()
    
    // Add process to combine QC reports
    COMBINE_QC_REPORTS(qc_reports)
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
    echo "Combined QC Report" > final_qc_report.txt
    echo "Generated: \$(date)" >> final_qc_report.txt
    echo "=========================" >> final_qc_report.txt
    
    for report in *_qc_report.txt; do
        echo "" >> final_qc_report.txt
        cat \$report >> final_qc_report.txt
        echo "-------------------------" >> final_qc_report.txt
    done
    """
}
