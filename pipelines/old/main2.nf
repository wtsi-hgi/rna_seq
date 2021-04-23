nextflow.preview.dsl=2

include { iget_cram } from '../modules/irods.nf'
include { crams_to_fastq_gz } from '../modules/crams_to_fastq.nf'
include { fastqc } from '../modules/fastqc.nf'
include { salmon } from '../modules/salmon.nf'
include { merge_salmoncounts } from '../modules/merge_salmoncounts.nf'
include { tximport } from '../modules/tximport.nf'
include { star_2pass_basic } from '../modules/star_2pass_basicmode.nf'
include { star_2pass_1st_pass } from '../modules/star_2pass_firstpass.nf'
include { star_2pass_merge_junctions } from '../modules/star_2pass_merge_junctions.nf'
include { star_2pass_2nd_pass } from '../modules/star_2pass_secondpass.nf'
include { filter_star_aln_rate } from '../modules/filter_star_aln_rate.nf'
include { leafcutter_bam2junc_regtools } from '../modules/leafcutter_bam2junc.nf'
include { leafcutter_clustering_regtools } from '../modules/leafcutter_clustering.nf'
include { featureCounts } from '../modules/featurecounts.nf'
include { samtools_index_idxstats } from '../modules/samtools_index_idxstats.nf'
include { mapsummary } from '../modules/mapsummary.nf'
include { merge_featureCounts } from '../modules/merge_featureCounts.nf'
include { multiqc } from '../modules/multiqc.nf'
include { lostcause } from '../modules/lostcause.nf'
include { baton_study_id } from '../modules/baton.nf'
include { heatmap } from '../modules/heatmap.nf'
include { deseq2 } from '../modules/deseq2.nf'
include { star_tabgenes_matrix } from '../modules/star_tabgenes_matrix.nf'
include { mbv } from '../modules/mbv.nf'
include { get_egan_id } from '../modules/get_egan_id.nf'

workflow {
    
    log.info "input params are:"
    log.info "$params"
    
    if (params.run_mode == "from_study_id") {
	// Get list of study cram files from Irods:
	ch_studies = Channel.fromList(params.from_study_id_mode.studies_list)
	baton_study_id(ch_studies)
	
	// Download cram files from Irods:
	to_iget = baton_study_id.out.samples_noduplicates_tsv
	    .map{a,b -> b}
	    .splitCsv(header: true, sep: '\t')
	    .map{row->tuple(row.sample, row.sample_supplier_name, row.study_id)}
	    .map{a,b,c-> tuple(a,c)}
	iget_cram(to_iget.take(2))
	
	// Optionally, extract egan IDs from cram files:
	get_egan_id(iget_cram.out[0])
	get_egan_id.out.samplename_egan_id_csv
	    .map { samplename,cram,csv_file -> csv_file }
	    .splitCsv(header: true, sep: ',')
	    .map { row -> "${row.samplename},${row.egan_id}"}
	    .collectFile(name: 'samplename_egan_id.csv', newLine: true,
			 seed: "samplename,egan_id",
			 storeDir: "${params.outdir}/", sort: true)
	    .set{ch_samplename_egan_id_csv}

	// Convert CRAM files to fastq:
	crams_to_fastq_gz(iget_cram.out[0])
	crams_to_fastq_gz.out[0]
	    .map{ samplename, fastq1, fastq2 -> tuple( samplename, tuple(fastq1, fastq2) ) }
            .set{ ch_samplename_crams }
	
    }
    
    else if (params.run_mode = "from_fastq_csv") {
	Channel
	    .fromPath(params.from_fastq_mode.fastq_csv)
	    .splitCsv(header:false)
	    .map{ row-> tuple(row[0], tuple(file(row[1]), file(row[2]))) }
	    .set { ch_samplename_crams }
    }
    
    else if (params.run_mode != "from_fastq_csv" & params.run_mode == "from_study_id") {
        log.info "error: input param run_mode must be set to either \"from_study_id\" or \"from_fastq_csv\""
        exit 1
    }

    // Modules common to all input modes:
    fastqc(ch_samplename_crams)

    // Salmon aligner modules:
    if(params.salmon_aligner.run_salmon = true) {
	Channel.fromPath(params.salmon_index)
	    .ifEmpty { exit 1, "Salmon index dir not found: ${params.salmon_index}" }
	    .set {ch_salmon_index}
	salmon(ch_samplename_crams, ch_salmon_index.collect()) // salmon(ch_samplename_crams, ch_salmon_index.collect(), ch_salmon_trans_gene.collect())
	tximport(salmon.out[0].collect())
	deseq2(salmon.out[0].collect(), Channel.fromPath(params.deseq2_tsv))
	merge_salmoncounts(salmon.out[0].collect(), salmon.out[1].collect())
	heatmap(merge_salmoncounts.out[0].map{transcounts,transtpm,genecouts,genetpm-> genecouts})
    }
    
    // STAR aligner modules:
    if(params.star_aligner.run_star_2pass_basic = true || params.star_aligner.run_star_2pass_1st_pass = true) {
	
	Channel.fromPath(params.star_index)
	    .ifEmpty { exit 1, "star index file not found: ${params.star_index}" }
	    .set { ch_star_index}
	
	if(params.star_aligner.run_star_2pass_basic) {
	}
	else if (params.star_aligner.run_star_custom_2pass) {
	}
	    //    leafcutter_bam2junc_regtools(star_out[0])
	    //    leafcutter_clustering_regtools(leafcutter_bam2junc_regtools.out.collect())
    }
}

workflow.onError {
    log.info "Pipeline execution stopped with the following message: ${workflow.errorMessage}" }

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Command line: $workflow.commandLine"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    
//    if (params.on_complete_uncache_irods_search) {
//	log.info "You have selected \"on_complete_uncache_irods_search = true\"; will therefore attempt to remove Irods work dirs to forcefully uncache them even if successful."
//	if (! file("${params.outdir}/irods_work_dirs_to_remove.csv").isEmpty()) {
//	    log.info "file ${params.outdir}/irods_work_dirs_to_remove.csv exists and not empty ..."
//	    file("${params.outdir}/irods_work_dirs_to_remove.csv")
//		.eachLine {  work_dir ->
//		if (file(work_dir).isDirectory()) {
//		    log.info "removing work dir $work_dir ..."
//		    file(work_dir).deleteDir()   
//		} } } }
    
    if (params.on_complete_remove_workdir_failed_tasks) {
	log.info "You have selected \"on_complete_remove_workdir_failed_tasks = true\"; will therefore remove work dirs of all tasks that failed (.exitcode file not 0)."
	// work dir and other paths are hardcoded here ... :
	def proc = "bash ./nextflow_ci/bin/del_work_dirs_failed.sh ${workDir}".execute()
	def b = new StringBuffer()
	proc.consumeProcessErrorStream(b)
	log.info proc.text
	log.info b.toString() }
}





    //.filter { it[1] ==~ /^[cC].*/} //.filter { it[1] ==~ /^[cC].*/}
    
    //// from irods studyid and list of samplenames
    //iget_cram(
    //	Channel.fromPath("${baseDir}/../../inputs/samples.txt")
    //	    .flatMap{ it.readLines()}, "5933")
    //crams_to_fastq_gz(iget_cram.out[0])
    ////

