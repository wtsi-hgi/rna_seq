nextflow.preview.dsl=2

// read all inputs from inputs.nf from upstream gitlab repo (at same branch name).

include baton_study from '../modules/baton_study.nf' params(run: true, outdir: params.outdir, dropqc: params.dropqc)
include iget_sample_study from '../modules/iget_sample_study.nf' params(run: true, outdir: params.outdir, copy_mode: params.copy_mode)
include iget_sample from '../modules/iget_sample.nf' params(run: true, outdir: params.outdir, copy_mode: params.copy_mode)

workflow {
    ch_input_studies = Channel.from('2971','3180','3257','3779')
    
    // multiple study_ids allowed
    if (params.run_from_studies) {
	baton_study(ch_input_studies)
	
	to_iget = baton_study.out.samples_noduplicates_tsv
	    .map{a,b -> b}
	    .splitCsv(header: true, sep: '\t')
	    .map{row->tuple(row.sample, row.sample_supplier_name, row.study_id)}
	    .map{a,b,c-> tuple(a,c)}

	iget_sample_study(to_iget)
	
    }
    else {
	Channel.fromPath(params.input_samples_csv)
	    .splitCsv(header: true, sep: '\t')
	    .map { row -> row.sample }
	    .set{ch_to_iget}
	
	iget_sample(ch_to_iget)
    }
}
