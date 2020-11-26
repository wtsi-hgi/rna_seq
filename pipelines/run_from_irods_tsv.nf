nextflow.enable.dsl=2

// All inputs and the selection of which tasks to run are defined in config file "inputs.nf",
// which is located in upstream Gitlab "nextflow_ci" repo (at same branch name).
// Meaning that if you wish to run pipeline with different parameters/tasks,
// you have to edit(+commit+push) that "inputs.nf" file, then (CI-)rerun the pipeline.

include { imeta_study_cellranger } from '../modules/imeta_study_cellranger.nf'
include { iget_study_cram } from '../modules/iget_study_cram.nf'
include { iget_study_cellranger } from '../modules/iget_study_cellranger.nf'

workflow run_from_irods_tsv {
    take: channel_samples_tsv
    main:

    // task to iget all Irods cram files of all samples
    iget_study_cram(
	channel_samples_tsv
	    .map{study_id, samples_tsv -> samples_tsv}
	    .splitCsv(header: true, sep: '\t')
	    .map{row->tuple(row.study_id, row.sample, row.object)}
	    .filter { it[2] =~ /.cram$/ } // Need to check for bam too?
	    .unique())
    
    // task to search Irods cellranger location for each sample:
    imeta_study_cellranger(
	channel_samples_tsv
	    .map{study_id, samples_tsv -> samples_tsv}
	    .splitCsv(header: true, sep: '\t')
	    .map{row->tuple(row.study_id, row.sample, row.id_run)}
	    .unique())

    // store the list cellranger locations found into a single tsv table called "cellranger_irods_objects.csv"
    imeta_study_cellranger.out.study_id_sample_cellranger_object
	.map{study_id, sample, run_id, cellranger_irods_object, workdir ->
	"${study_id},${sample},${run_id},${cellranger_irods_object},${workdir}"}
	.collectFile(name: 'cellranger_irods_objects.csv', newLine: true, sort: true,
		     seed: "study_id,sanger_sample_id,run_id,cellranger_irods_object,workdir",
		     storeDir:params.outdir)

    // task to iget the cellranger outputs from Irods:
    iget_study_cellranger(imeta_study_cellranger.out.study_id_sample_cellranger_object
			  .map{study_id, sample, run_id, cellranger_irods_object, workdir ->
	    tuple(study_id, sample, cellranger_irods_object)}
			  .filter { it[2] != "cellranger_irods_not_found" })

}

// TODO:  here or main.nf:   // store work dirs to remove into tsv file for onComplete removal.
    //imeta_study.out.work_dir_to_remove.mix(
//	imeta_study_cellranger.out.work_dir_to_remove
//	    .filter { it != "dont_remove" })
//	.collectFile(name: 'irods_work_dirs_to_remove.csv', newLine: true, sort: true,
//		     storeDir:params.outdir)
    

  	//run_from_sanger_sample_id(Channel.fromPath(params.input_samples_csv)
	//			  .splitCsv(header: true, sep: '\t')
	//			  .map { row -> row.sample })
