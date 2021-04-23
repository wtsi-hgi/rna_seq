process heatmap {
    tag "heatmap"

    publishDir "${params.outdir}/heatmap/", mode: "${params.copy_mode}", overwrite: true
    
    when:
    params.salmon_aligner.salmon_dowstream_tasks.heatmap_task.run

    input:
    file (count_matrix_tsv)

    output:
    tuple file("outputs/salmon_PCA_unbiased_toppc.pdf"), file("outputs/salmon_heatmap_toppc.pdf"), emit: pca_heatmap

    script:
    """
    Rscript $workflow.projectDir/../bin/heatmap.R $count_matrix_tsv
    """
}
