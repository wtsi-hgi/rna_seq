process mapsummary {
    tag "${samplename}"
    
    publishDir "${params.outdir}/mapsummary/", mode: "${params.copy_mode}", overwrite: true

    when:
    params.star_aligner.star_downstream_tasks.samtools_index_idxstats_task.map_summary_task.run
    
    input:
    tuple val(samplename), file(thestats) // from ch_mapsummary

    output:
    file "*_mqc.txt" //into ch_multiqc_mapsum

    script:
    """
    python $baseDir/bin/mito.py -m ${params.star_aligner.star_downstream_tasks.samtools_index_idxstats_task.map_summary_task.mito_name} -t $thestats > ${samplename}_mqc.txt
    """
}
