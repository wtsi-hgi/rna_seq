


    if(params.run_star_basic) {
	star_2pass_basic(ch_samplename_crams, ch_star_index.collect(), ch_gtf_star.collect()) }

    star_2pass_1st_pass(ch_samplename_crams, ch_star_index.collect(), ch_gtf_star.collect())
    star_2pass_merge_junctions(star_2pass_1st_pass.out[1].collect())
    star_2pass_2nd_pass(ch_samplename_crams, ch_star_index.collect(), ch_gtf_star.collect(), star_2pass_merge_junctions.out)
    
    star_out = star_2pass_2nd_pass.out // choose star_2pass_basic.out or star_2pass_2ndpass.out 
    // star_out = star_2pass_basic.out

    mbv(star_out[0], ch_mbv_vcf_gz.collect(), ch_mbv_vcf_gz_csi.collect()) }

    // leafcutter_bam2junc(star_out[0])
    // leafcutter_clustering(leafcutter_bam2junc.out.collect())
    // or regtools:
//    leafcutter_bam2junc_regtools(star_out[0])
//    leafcutter_clustering_regtools(leafcutter_bam2junc_regtools.out.collect())
    
    filter_star_aln_rate(star_out[1].map{samplename,logfile,bamfile -> [samplename,logfile]}) // discard bam file, only STAR log required to filter
    
    filter_star_aln_rate.out.branch {
        filtered: it[1] == 'above_threshold'
        discarded: it[1] == 'below_threshold'}.set { star_filter }
    
    star_filter.filtered.combine(star_out[1], by:0) //reattach bam file
	.map{samplename,filter,logfile,bamfile -> ["star", samplename, bamfile]} // discard log file and attach aligner name
	.set{star_filtered} 
    
    samtools_index_idxstats(star_filtered)
    
    mapsummary(samtools_index_idxstats.out)
    
    featureCounts(star_filtered, ch_gtf_star.collect(), ch_biotypes_header.collect())

    merge_featureCounts(featureCounts.out[0].map{samplename, gene_fc_txt -> gene_fc_txt}.collect())

//    crams_to_fastq_gz.out[1]
	//.mix(
	star_filter.discarded.map{samplename, filter -> [text: "${samplename}\tSTAR\tlowmapping\n"]}//)
	.set{ch_lostcause }
    lostcause(ch_lostcause.collectFile({ ['lostcause.txt', it.text]},sort:true))

    featureCounts.out[1]
	.filter{ pick_aligner(it[0]) }
	.map { it[1] }
	.set{ ch_multiqc_fc_aligner }

    featureCounts.out[2]
	.filter{ pick_aligner(it[0]) }
	.map{ it[1] }
	.set{ ch_multiqc_fcbiotype_aligner }

    if (params.run_salmon) {
	multiqc(lostcause.out.collect().ifEmpty([]),
		fastqc.out.collect().ifEmpty([]),
		mapsummary.out.collect().ifEmpty([]),
		ch_multiqc_fc_aligner.collect().ifEmpty([]),
		ch_multiqc_fcbiotype_aligner.collect().ifEmpty([]),
		star_out[2].collect().ifEmpty([]),
		salmon.out[0].collect().ifEmpty([]))}
    else {
	multiqc(lostcause.out.collect().ifEmpty([]),
		fastqc.out.collect().ifEmpty([]),
		mapsummary.out.collect().ifEmpty([]),
		ch_multiqc_fc_aligner.collect().ifEmpty([]),
		ch_multiqc_fcbiotype_aligner.collect().ifEmpty([]),
		star_out[2].collect().ifEmpty([]),
		[])}
