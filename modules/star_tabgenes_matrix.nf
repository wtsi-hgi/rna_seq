process star_tabgenes_matrix {
    tag "star_tabgenes_matrix"

    publishDir "${params.outdir}/combined/", mode: "${params.copy_mode}", overwrite: true

    when:
    params.star_aligner.star_2pass_basic_task.star_tabgenes_matrix_task.run
    
    input:
    file (tagenes_files)

    output:
    tuple val("star"), file("star_tabgenes_matrix.tsv"), emit: star_matrix

    script:
    """
    Rscript $workflow.projectDir/bin/star_tabgenes_matrix.R
    """
}
