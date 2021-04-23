process samtools_index_idxstats {
    tag "${samplename}"

    publishDir "${params.outdir}/idxstats/", mode: "${params.copy_mode}", overwrite: true, pattern: "*.idxstats"
    
    when:
    params.star_aligner.star_downstream_tasks.samtools_index_idxstats_task.run

    input:
    tuple val(aligner), val(samplename), file(thebam) //from ch_indexbam

    output:
    tuple val(samplename), file("*.idxstats") //into ch_mapsummary

    script:
    """
    samtools index $thebam
    samtools idxstats $thebam > ${samplename}.idxstats
    rm ${thebam}.bai
    """
}
