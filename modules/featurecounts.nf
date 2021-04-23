process featureCounts {
    tag "${samplename}"
    
    publishDir "${params.outdir}/featureCounts/", mode: "${params.copy_mode}", overwrite: true,
        saveAs: {filename ->
            if (filename.indexOf(".biotype_counts_mqc.txt") > 0) "biotype_counts_mqc/$filename"
            else if (filename.indexOf(".biotype.fc.txt") > 0) "biotype_counts/$filename"
            else if (filename.indexOf(".biotype.fc.txt.summary") > 0) "biotype_counts_summaries/$filename"
            else if (filename.indexOf(".gene.fc.txt.summary") > 0) "gene_count_summaries/$filename"
            else if (filename.indexOf(".gene.fc.txt") > 0) "gene_counts/$filename"
            else "$filename"
        }

    when:
    params.star_aligner.star_downstream_tasks.featureCounts_task.run
    
    input:
    tuple val(aligner), val(samplename), file(thebam) //from ch_featurecounts
    file gtf //from ch_gtf_featureCounts  //.collect()
    file biotypes_header

    output:
    tuple val(aligner), file("*.gene.fc.txt") //into ch_merge_fc
    tuple val(aligner), file("*.gene.fc.txt.summary") //into ch_multiqc_fc
    tuple val(aligner), file("*.biotype_counts_mqc.txt") //into ch_multiqc_fcbiotype

    script:
    def extraparams = params.star_aligner.star_downstream_tasks.featureCounts_task.fcextra.toString() - ~/^dummy/
    def fc_direction = 0
    def tag = "${samplename}.${aligner}"

    def pairedend = params.star_aligner.star_downstream_tasks.featureCounts_task.singleend ? "" : "-p"
    if (params.star_aligner.star_downstream_tasks.featureCounts_task.forward_stranded && !params.star_aligner.star_downstream_tasks.featureCounts_task.unstranded) {
        fc_direction = 1
    } else if (params.star_aligner.star_downstream_tasks.featureCounts_task.reverse_stranded && !params.star_aligner.star_downstream_tasks.featureCounts_task.unstranded){
        fc_direction = 2
    }
    outfile = "${tag}.gene.fc.txt"
    """
    featureCounts -T ${task.cpus} -a $gtf -g gene_id          \\
      -o ${outfile} $pairedend                                \\
      -s $fc_direction ${extraparams} $thebam
    cut -f 1,7 ${outfile} > reduced.${outfile}   #  This
    mv reduced.${outfile} ${outfile}             #  reduces the file size from ~ 30M to ~1M
    featureCounts -T ${task.cpus} -a $gtf -g gene_id  \\
      -o ${tag}.biotype.fc.txt $pairedend                     \\
      -s $fc_direction ${extraparams} $thebam
    cut -f 1,7 ${tag}.biotype.fc.txt |                        \\
        tail -n +3 | cat $biotypes_header - >> ${tag}.biotype_counts_mqc.txt
    """
}
