process merge_featureCounts {
    tag "$aligner"
    
    publishDir "${params.outdir}/combined", mode: "${params.copy_mode}", overwrite: true

    when:
    params.star_aligner.star_downstream_tasks.featureCounts_task.merge_featureCounts_task.run

    input:
    file(collected_fc_gene_txt)

    output:
    file '*-fc-genecounts.txt'
    file("fofn_gene_featurecount.txt")

    shell:
    suffix=['star':'.star.gene.fc.txt', 'hisat2':'.hisat2.gene.fc.txt']
    aligner = "star"  // not strictly necessary
    outputname = "${aligner}-fc-genecounts.txt"
    thesuffix  = suffix[aligner] ?: '.txt'
    '''
    ls . | grep gene.fc.txt\$ > fofn_gene_featurecount.txt

    python3 !{workflow.projectDir}/bin/merge_featurecounts.py        \\
      --rm-suffix !{thesuffix}                                       \\
      -c 1 --skip-comments --header                                  \\
      -o !{outputname} -I fofn_gene_featurecount.txt
    '''
}
