nextflow.enable.dsl=2

workflow check_and_load_input_files {
    main:
		log.info "\nrunning check_and_load_input_files() sub-workflow..."
		// check that input files exist and load into channels:
		println "Project : $params.input_mode"
	
		if (params.input_mode == "from_fastq_csv") {
			// check input fastq csv file:
			println "Using from_fastq_csv"
			Channel.fromPath(params.input_from_fastq_csv.fastq_csv, checkIfExists: true).set { ch_input_fastq_csv }
		} 
		else {
			ch_input_fastq_csv = Channel.from('not_used')
		}
		
		if (params.salmon_aligner.salmon_task.run) {
		// if Salmon is set to run, then check Salmon index file:	
			println "Salmon aligner task"
			Channel.fromPath(params.salmon_aligner.salmon_task.salmon_index, checkIfExists: true).set { ch_salmon_index }
			if (params.salmon_aligner.salmon_downstream_tasks.deseq2_task.run) {
				// check deseq2 input file:
				Channel.fromPath(params.salmon_aligner.salmon_downstream_tasks.deseq2_task.deseq2_tsv, checkIfExists: true).set { ch_deseq2_tsv }
			} else {
				ch_deseq2_tsv = Channel.from('not_used')
			}
		} else {
			ch_salmon_index = Channel.from('not_used')
			ch_deseq2_tsv = Channel.from('not_used')
		}
		
		if (params.star_aligner.star_2pass_basic_task.run || params.star_aligner.star_custom_2pass_task.run_1st_pass) {
			// check STAR index dir and gtf file:
			println "STAR aligner task"
			Channel.fromPath(params.star_aligner.gtf, checkIfExists: true).set { ch_gtf_star }
			Channel.fromPath(params.star_aligner.star_index, checkIfExists: true).set { ch_star_index }

			if (params.star_aligner.star_downstream_tasks.featureCounts_task.run) {
				Channel.fromPath(params.star_aligner.star_downstream_tasks.featureCounts_task.biotypes_header, checkIfExists: true).set { ch_biotypes_header }
			} else { 
				ch_biotypes_header = Channel.from('not_used')
			}
			if (params.star_aligner.star_downstream_tasks.mbv_task.run) {
				Channel.fromPath(params.star_aligner.star_downstream_tasks.mbv_task.mbv_vcf_gz, checkIfExists: true).set { ch_mbv_vcf_gz }
				Channel.fromPath(params.star_aligner.star_downstream_tasks.mbv_task.mbv_vcf_gz_csi, checkIfExists: true).set { ch_mbv_vcf_gz_csi }
			} else {
				ch_mbv_vcf_gz = Channel.from('not_used')
				ch_mbv_vcf_gz_csi = Channel.from('not_used')
			}
		} else {
			ch_gtf_star = Channel.from('not_used')
			ch_star_index = Channel.from('not_used')
		}
		
		if(params.input_mode == "from_study_id") {
			println "Getting data from study ID"
			// if input_mode = "from_study_id", check that study IDs were provided (at least one):
			Channel.from(params.study).set { ch_input_study_ids }
		} else {
			ch_input_study_ids = Channel.from('not_used')}
		
    emit:
		ch_input_study_ids // input Irods study IDs to get CRAM files from
		ch_input_fastq_csv // input fastq csv file
		ch_salmon_index // input salmon index directory
		ch_deseq2_tsv // input deseq 2 tsv file
		ch_star_index // input STAR index directory
		ch_gtf_star // input STAR GTF file
		ch_mbv_vcf_gz // input multi-sample vcf for MBV QTLtools
		ch_mbv_vcf_gz_csi // .csi index for input multi-sample vcf for MBV QTLtools
		ch_biotypes_header // biotypes header file for featurecounts
    
}
