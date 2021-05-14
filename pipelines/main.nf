nextflow.preview.dsl=2

// include sub-workflow that will
//   1) check that required input files point to existing files (or directories),
//   2) load these input files into Nextflow channels:
include { check_and_load_input_files } from './check_and_load_input_files.nf'

// include sub-workflow to get f`astq into channel
//    (from Irods or local files, depending on input mode):
include { get_fastq } from './get_fastq.nf'
include { process_fastq } from '../modules/process_fastq.nf'
// include sub-workflow to run main rna-seq code on samples/fastq inputs:
include { main_rnaseq } from './main_rnaseq.nf'

workflow {

    // show input parameters into main nextflow log file:
    log.info "\ninput params are:"
    log.info "$params"

    // check inputs files exist on disk, and load them into Nextlfow channels:
    check_and_load_input_files()

    // get fastq files, depending on input mode (from Irods CRAM files or from list of fastq files):
    if (params.input_mode == "from_fastq") {
		// In this part we are loading fastq files already prepeared from a directory. We assess whether the fatq files are gzipped as .gz format is needed.
		//  FastqCSV file is prepeared and eventually a channel created accordingly
		process_fastq(params.input_from_fastq_csv.fastq_path)	
		ch_samplename_crams = process_fastq.out.ch_samplename_crams
		// Channel.fromPath("${params.outdir}/fastq_gz_csv/Fastq_files.csv").splitCsv(header:false)
		// 	.map{ row-> tuple(row[0], tuple(file(row[1]), file(row[2]))) }
		// 	.take(-1) // set to -1 for all
		// 		.set { ch_samplename_crams }
    }else{
		// Here we are retrieving cram files and converting them to the fastq files
		get_fastq(check_and_load_input_files.out.ch_input_study_ids, // input Irods study IDs to get CRAM files from
			check_and_load_input_files.out.ch_input_fastq_csv) // input fastq csv file
		ch_samplename_crams = get_fastq.out.ch_samplename_crams
    }
    
    
    // run main rna-seq workflow on those samples/fastq:
    main_rnaseq(check_and_load_input_files.out.ch_salmon_index, // input salmon index directory
		check_and_load_input_files.out.ch_deseq2_tsv, // input deseq 2 tsv file
		check_and_load_input_files.out.ch_star_index, // input STAR index directory
		check_and_load_input_files.out.ch_gtf_star, // input STAR GTF file
		check_and_load_input_files.out.ch_mbv_vcf_gz, // input multi-sample vcf for MBV QTLtools
		check_and_load_input_files.out.ch_mbv_vcf_gz_csi, // .csi index for input multi-sample vcf for MBV QTLtools
		check_and_load_input_files.out.ch_biotypes_header, // biotypes header file for featurecounts
		ch_samplename_crams) // channel of tuple(samplename, tuple(fastq1/fastq2)) for paired end reads of each sample to process.
    
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
