nextflow.enable.dsl=2

// include RNAseq Nextflow modules/tasks:
include { fastqc } from '../modules/fastqc.nf'
include { salmon } from '../modules/salmon.nf'
include { merge_salmoncounts } from '../modules/merge_salmoncounts.nf'
include { tximport } from '../modules/tximport.nf'
include { star_2pass_basic } from '../modules/star_2pass_basicmode.nf'
include { star_2pass_1st_pass } from '../modules/star_2pass_firstpass.nf'
include { star_2pass_merge_junctions } from '../modules/star_2pass_merge_junctions.nf'

if (params.star_aligner.star_downstream_tasks.featureCounts_task.singleend){
    include { star_2pass_2nd_pass } from '../modules/star_2pass_secondpass_singleend.nf'
}else{
    include { star_2pass_2nd_pass } from '../modules/star_2pass_secondpass.nf'
}

include { filter_star_aln_rate } from '../modules/filter_star_aln_rate.nf'
include { leafcutter_bam2junc_regtools } from '../modules/leafcutter_bam2junc.nf'
include { leafcutter_clustering_regtools } from '../modules/leafcutter_clustering.nf'
include { featureCounts } from '../modules/featurecounts.nf'
include { samtools_index_idxstats } from '../modules/samtools_index_idxstats.nf'
include { mapsummary } from '../modules/mapsummary.nf'
include { merge_featureCounts } from '../modules/merge_featureCounts.nf'
include { multiqc } from '../modules/multiqc.nf'
include { lostcause } from '../modules/lostcause.nf'
include { deseq2 } from '../modules/deseq2.nf'
include { star_tabgenes_matrix } from '../modules/star_tabgenes_matrix.nf'
include { heatmap } from '../modules/heatmap.nf'
include { mbv } from '../modules/mbv.nf'
include { EXTRACT_FINGERPRINTS } from '../subworkflows/local/extract_fingerprints.nf'


workflow main_rnaseq {
    take:
    ch_salmon_index // input salmon index directory
    ch_deseq2_tsv // input deseq 2 tsv file
    ch_star_index // input STAR index directory
    ch_gtf_star // input STAR GTF file
    ch_mbv_vcf_gz // input multi-sample vcf for MBV QTLtools
    ch_mbv_vcf_gz_csi // .csi index for input multi-sample vcf for MBV QTLtools
    ch_biotypes_header // biotypes header file for featurecounts
    ch_samplename_crams // tuple(samplename, tuple(fastq1,fastq2)) for paired-end reads
    
    main:
    log.info "\nrunning main_rnaseq() sub-workflow..."

    fastqc(ch_samplename_crams)

    // salmon aligner tasks:

    salmon(ch_samplename_crams, ch_salmon_index.collect())
    tximport(salmon.out[0].collect())
    merge_salmoncounts(salmon.out[0].collect(), salmon.out[1].collect())
    // TODO: check if still works: heatmap(merge_salmoncounts.out[0].map{transcounts,transtpm,genecouts,genetpm-> genecouts})
    // TODO: check if still works: deseq2(salmon.out[0].collect(), ch_deseq2_tsv)

    // star aligner tasks:

    // run star in one mode 1, 1 task:
    star_2pass_basic(ch_samplename_crams, ch_star_index.collect(), ch_gtf_star.collect())
    // to test: star_tabgenes_matrix(star_2pass_basic.out)
    
    // run star in mode 2, 3 tasks:
    star_2pass_1st_pass(ch_samplename_crams, ch_star_index.collect(), ch_gtf_star.collect())
    star_2pass_merge_junctions(star_2pass_1st_pass.out[1].collect())
    star_2pass_2nd_pass(ch_samplename_crams, ch_star_index.collect(), ch_gtf_star.collect(), star_2pass_merge_junctions.out)

    // choose STAR output between the 2 STAR modes: star_2pass_basic.out or star_2pass_2ndpass.out 
    if (params.star_aligner.star_downstream_tasks.downstream_outputs == 'star_custom_2pass') {
	    star_out = star_2pass_2nd_pass.out
    } else if (params.star_aligner.star_downstream_tasks.downstream_outputs == 'star_2pass_basic') {
	    star_out = star_2pass_basic.out 
    }

    // https://qtltools.github.io/qtltools/ :MBV is used to ensure good match between the sequence and genotype data. This is useful to detect sample mislabeling, contamination or PCR amplification biases.
    mbv(star_out[0], ch_mbv_vcf_gz.collect(), ch_mbv_vcf_gz_csi.collect())
    
    // this step runs regtools jusnctions function: The junctions extract command can be used to extract exon-exon junctions from an RNAseq BAM file.
    leafcutter_bam2junc_regtools(star_out[0])

    // this part uses python script put in a container - TODO check
    leafcutter_clustering_regtools(leafcutter_bam2junc_regtools.out.collect())

    // python script to perform filtering.
    filter_star_aln_rate(star_out[1].map{samplename,logfile,bamfile -> [samplename,logfile]}) // discard bam file, only STAR log required to filter
    filter_star_aln_rate.out.branch {
        filtered: it[1] == 'above_threshold'
        discarded: it[1] == 'below_threshold'
    }.set { star_filter }
    star_filter.filtered.combine(star_out[1], by:0) //reattach bam file
	    .map{samplename,filter,logfile,bamfile -> ["star", samplename, bamfile]} // discard log file and attach aligner name
	    .set{star_filtered}
    
    // this indexes the bam file with  samtools index – indexes
    samtools_index_idxstats(star_filtered)

    // python code to generate a summary file
    mapsummary(samtools_index_idxstats.out)

    // uses the featureCounts software which is designed to be an an efficient general purpose program for assigning sequence reads to genomic features.
    featureCounts(star_filtered, ch_gtf_star.collect(), ch_biotypes_header.collect())

    // python script performing merging.
    merge_featureCounts(featureCounts.out[0].map{samplename, gene_fc_txt -> gene_fc_txt}.collect())

    // extracts fingerprints from each sample using GATK
    if (params.extract_fingerprints) {
        star_filtered
            .map { aligner, samplename, bamfile -> [[id: samplename], bamfile] }
            .set { star_filtered_ch }

        EXTRACT_FINGERPRINTS(
            star_filtered_ch,
            [params.reference_sequence, params.reference_index, params.reference_dict],
            params.haplotype_map
        )
    }

    star_filter.discarded.map{samplename, filter -> [text: "${samplename}\tSTAR\tlowmapping\n"]}
	    .set{ch_lostcause}

    // Assesses the lost samples and reports in the combined folder: lostcause_mqc.txt
    lostcause(ch_lostcause.collectFile({ ['lostcause.txt', it.text] }, sort:true))
    
    featureCounts.out[1]
	    .map{ it[1] }
	    .set{ ch_multiqc_fc_aligner }
    
    featureCounts.out[2]
	    .map{ it[1] }
	    .set{ ch_multiqc_fcbiotype_aligner }
    
    // the folowing part Aggregates results from bioinformatics analyses across many samples into a single report
    if (params.salmon_aligner.salmon_task.run) {
        multiqc(lostcause.out.collect().ifEmpty([]),
            fastqc.out.collect().ifEmpty([]),
            mapsummary.out.collect().ifEmpty([]),
            ch_multiqc_fc_aligner.collect().ifEmpty([]),
            ch_multiqc_fcbiotype_aligner.collect().ifEmpty([]),
            star_out[2].collect().ifEmpty([]),
            salmon.out[0].collect().ifEmpty([]))
    } else {
        multiqc(lostcause.out.collect().ifEmpty([]),
            fastqc.out.collect().ifEmpty([]),
            mapsummary.out.collect().ifEmpty([]),
            ch_multiqc_fc_aligner.collect().ifEmpty([]),
            ch_multiqc_fcbiotype_aligner.collect().ifEmpty([]),
            star_out[2].collect().ifEmpty([]),
            [])
    }
}
