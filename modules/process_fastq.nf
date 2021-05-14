nextflow.enable.dsl=2

// include RNAseq Nextflow modules/tasks:
include { iget_cram } from '../modules/irods.nf'
include { gz_fastq } from '../modules/gz_fastq.nf'
include { fastq_to_csv } from '../modules/fastq_to_csv.nf'
include { crams_to_fastq_gz } from '../modules/crams_to_fastq.nf'
include { baton_study_id } from '../modules/baton.nf'
include { mbv } from '../modules/mbv.nf'
include { get_egan_id } from '../modules/get_egan_id.nf'

workflow process_fastq {
    
    take:
    fastq_path // input Irods study IDs to get CRAM files from
    // ch_input_fastq_csv // input fastq csv file
    
    main:
    log.info "\ncopying fastq and assessing whether the compression has been performed sub-workflow to get correct fastq files if fastq is used as an inputs"

    if (params.input_mode == "from_fastq") {
		log.info "We are analysing data from from_fastq files"
		// Get list of study cram files from Irods:
		Channel.fromPath( "${fastq_path}/*{.fastq,.fastq.gz}" ).take(-1) // set to -1 for all                                                                                      
            .set { fileChannel }
		gz_fastq(fileChannel)
                
        println("lets pass the list in now to generate csv file that can be used for the channels")
        fastq_to_csv(gz_fastq.out.collect()) 
        fastq_to_csv.out.view()
        fastq_to_csv.out.splitCsv(header:false)
            .map{ row-> tuple(row[0], tuple(file(row[1]), file(row[2]))) }
            .take(-1) // set to -1 for all
                .set { ch_samplename_crams }

	}
    
    emit:
        ch_samplename_crams // channel of tuple(samplename,tuple(fastq1/fastq2)) for paired end reads.
    
}

