#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import modules
include { extract_sample_info } from './modules/extract_sample_info'
include { qc_check_input } from './modules/qc_check_input'
include { merge_samples } from './modules/merge_samples'
include { extract_chromosomes } from './modules/extract_chromosomes'
include { qc_check_extracted } from './modules/qc_check_extracted'
include { qc_check_samples_merged } from './modules/qc_check_samples_merged'
include { merge_with_1000genomes } from './modules/merge_with_1000genomes'
include { qc_check_final_merged } from './modules/qc_check_final_merged'
include { combine_timing_reports } from './modules/combine_timing_reports'
include { combine_qc_reports } from './modules/combine_qc_reports'
include { merge_chromosomes } from './modules/merge_chromosomes'
include { ld_prune } from './modules/ld_prune'
include { run_pca } from './modules/run_pca'
include { create_pca_data_file } from './modules/create_pca_data_file'
include { create_pca_plot } from './modules/create_pca_plot'
include { summarize_vcfs as summarize_vcfs_extraction } from './modules/summarize_vcfs'
include { summarize_vcfs as summarize_vcfs_ancient_merge } from './modules/summarize_vcfs'
include { summarize_vcfs as summarize_vcfs_1000g_merge } from './modules/summarize_vcfs'
include { summarize_vcfs as summarize_vcfs_ld_prune } from './modules/summarize_vcfs'
include { summarize_vcfs as summarize_vcfs_all_chr_merge } from './modules/summarize_vcfs'

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

workflow {
    // Parse chromosome list and create debug view
    def chromosomes = Channel.fromList(
        expandChromosomeList(params.chromosomes)
    ).view { "DEBUG: Will process chromosome: $it" }
    
    // Create the 1000G VCF path with the correct chromosome
    def kg_vcf = chromosomes.map { chr ->
        def vcf_path = params.thousand_genomes_vcf.replace("{chr}", chr)
        println "DEBUG: Created 1000G VCF path for chr${chr}: ${vcf_path}"
        [chr, file(vcf_path)]
    }.view { "DEBUG: 1000G VCF for chr${it[0]}: ${it[1]}" }

    // Create ancient VCFs channel - now expecting single multi-chromosome VCFs
    def ancient_vcfs = Channel
        .fromPath(params.ancient_vcfs_tsv)
        .splitCsv(header: false, sep: '\t')
        .map { row -> 
            def vcf_path = file(row[0].trim())
            println "DEBUG: Loading ancient VCF: ${vcf_path}"
            if (!vcf_path.exists()) {
                error "VCF file not found: ${vcf_path}"
            }
            return vcf_path
        }
        .view { "DEBUG: Ancient VCF loaded: $it" }

    // Extract sample info from PED file
    extract_sample_info(
        Channel.fromPath(params.ped_file)
    )

    // QC check input files
    qc_check_input(ancient_vcfs, "input_ancient_vcfs")
    
    // Branch based on QC results
    def qc_results = qc_check_input.out.report
        .combine(ancient_vcfs)
        .branch {
            pass: !it[0].text.contains("ERROR")
            fail: it[0].text.contains("ERROR")
        }

    // Create combinations of VCFs and chromosomes for extraction
    def combined_ch = qc_results.pass
        .map { report, vcf -> vcf }
        .combine(chromosomes)
        .view { "DEBUG: Will extract chr${it[1]} from ${it[0]}" }

    // Extract chromosomes
    extract_chromosomes(
        combined_ch.map { vcf, chr -> vcf },
        combined_ch.map { vcf, chr -> chr }
    )
    
    qc_check_extracted(extract_chromosomes.out.vcf, "post_extract_chromosomes")

    // After chromosome extraction
    summarize_vcfs_extraction(
        extract_chromosomes.out.vcf.collect(),
        "post_chromosome_extraction"
    )

    // Group extracted VCFs by chromosome for merging
    def grouped_vcfs = extract_chromosomes.out.vcf
        .map { vcf -> 
            def chrMatch = vcf.name =~ /chr(\d+)/
            def chr = chrMatch[0][1]  // Get the captured digit group
            println "DEBUG: Extracted VCF ${vcf} contains chromosome ${chr}"
            [chr, vcf]
        }
        .groupTuple()
        .view { chr, files -> 
            println "DEBUG: Chromosome ${chr} has ${files.size()} VCF files:"
            files.each { println "  - ${it}" }
            [chr, files]
        }

    merge_samples(
        grouped_vcfs
            .map { chr, files -> 
                println "DEBUG: Sending ${files.size()} files for chromosome ${chr} to merge_samples"
                tuple(chr, files)  // Explicitly create tuple
            }
    )

    qc_check_samples_merged(merge_samples.out.vcf, "post_merge_ancient")

    // After merging ancient samples
    summarize_vcfs_ancient_merge(
        merge_samples.out.vcf.collect(),
        "post_ancient_merge"
    )

    // Merge with 1000G by chromosome
    merge_with_1000genomes(
        kg_vcf,  // Already has [chr, vcf] structure
        merge_samples.out.vcf
    )

    qc_check_final_merged(merge_with_1000genomes.out.vcf, "post_1kg_merge")

    // After merging with 1000G
    summarize_vcfs_1000g_merge(
        merge_with_1000genomes.out.vcf.collect(),
        "post_1000g_merge"
    )

    // NEW: Merge all chromosomes into single file
    merge_chromosomes(
        merge_with_1000genomes.out.vcf.collect()
    )

    // NEW: Add summary after chromosome merge
    summarize_vcfs_all_chr_merge(
        merge_chromosomes.out.vcf.collect(),
        "post_chromosome_merge"
    )

    // Perform LD pruning on merged data
    ld_prune(
        merge_chromosomes.out.vcf
    )

    // After LD pruning
    summarize_vcfs_ld_prune(
        ld_prune.out.vcf.collect(),
        "post_ld_pruning"
    )

    // Run PCA on pruned data
    run_pca(ld_prune.out.vcf)

    // Create PCA plot using the sample info from PED file
    create_pca_plot(
        run_pca.out.pca_results,
        extract_sample_info.out.sample_info
    )

    // Create PCA data file
    create_pca_data_file(
        run_pca.out.pca_results,
        extract_sample_info.out.sample_info,
        Channel.fromPath("${baseDir}/igsr_populations.tsv")
    )

    // Collect timing reports
    def timing_ch = extract_chromosomes.out.timing
        .mix(merge_samples.out.timing)
        .mix(merge_with_1000genomes.out.timing)
        .mix(merge_chromosomes.out.timing)
        .collect()

    // Generate final timing report
    combine_timing_reports(timing_ch)

    // Collect all QC reports
    def qc_reports = Channel.empty()
        .mix(qc_check_input.out.report)
        .mix(qc_check_extracted.out.report)
        .mix(qc_check_samples_merged.out.report)
        .mix(qc_check_final_merged.out.report)
        .collect()
    
    // Add process to combine QC reports
    combine_qc_reports(qc_reports)
}
