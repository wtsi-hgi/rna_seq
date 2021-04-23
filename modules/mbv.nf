process mbv {
    tag "$samplename"
    
    publishDir "${params.outdir}/mbv", mode: "${params.copy_mode}", overwrite: true

    when:
    params.star_aligner.star_downstream_tasks.mbv_task.run

    input:
    tuple val(samplename), file(bam), file(bai)
    file vcf_gz
    file vcf_gz_csi

    output:
    tuple val(samplename), file("${samplename}.bamstat.txt")

    script:
    """
# run qtltools mbv
QTLtools_1.2_Ubuntu16.04_x86_64 mbv --bam ${bam} \\
--vcf ${vcf_gz} \\
--filter-mapping-quality 150 \\
--out ${samplename}.bamstat.txt
    """
}
