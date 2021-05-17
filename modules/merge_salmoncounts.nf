process merge_salmoncounts {
    tag ""
    
    publishDir "${params.outdir}/combined", mode: "${params.copy_mode}", overwrite: true
    
    input:
    file (all_quant_sf)
    file (all_quant_genes_sf)

    when:
    params.salmon_aligner.salmon_downstream_tasks.merge_salmoncounts_task.run

    output:
    tuple file('*transcounts.txt'), file('*transtpm.txt'), file('*genecounts.txt'), file('*genetpm.txt')
    file("fofn_quant_sf_salmon.txt")
    file("fofn_quant_genes_sf_salmon.txt")
    
    script:
    def outtranscount = "salmon-transcounts.txt"
    def outgenescount = "salmon-genecounts.txt"
    def outtranstpm   = "salmon-transtpm.txt"
    def outgenestpm   = "salmon-genetpm.txt"
    """
    ls . | grep .quant.sf\$ > fofn_quant_sf_salmon.txt
    ls . | grep .quant.genes.sf\$ > fofn_quant_genes_sf_salmon.txt

    python3 $workflow.projectDir/../bin/merge_featurecounts.py           \\
      --rm-suffix .quant.genes.sf                                     \\
      -c -1 --skip-comments --header                                  \\
      -o $outgenescount -I fofn_quant_genes_sf_salmon.txt

    python3 $workflow.projectDir/../bin/merge_featurecounts.py           \\
      --rm-suffix .quant.sf                                           \\
      -c -1 --skip-comments --header                                  \\
      -o $outtranscount -I fofn_quant_sf_salmon.txt

    python3 $workflow.projectDir/../bin/merge_featurecounts.py           \\
      --rm-suffix .quant.genes.sf                                     \\
      -c -2 --skip-comments --header                                  \\
      -o $outgenestpm -I fofn_quant_genes_sf_salmon.txt

    python3 $workflow.projectDir/../bin/merge_featurecounts.py           \\
      --rm-suffix .quant.sf                                           \\
      -c -2 --skip-comments --header                                  \\
      -o $outtranstpm -I fofn_quant_sf_salmon.txt
    """
}
