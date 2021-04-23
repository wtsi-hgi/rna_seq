process 'leafcutter_clustering_regtools' {

    publishDir "${params.outdir}/leafcutter_regtools/clustering", mode: "${params.copy_mode}", overwrite: true, pattern: "*.junc.clust.sorted.gz"
    publishDir "${params.outdir}/leafcutter_regtools/clustering", mode: "${params.copy_mode}", overwrite: true, pattern: "clust_*"
    publishDir "${params.outdir}/leafcutter_regtools/clustering", mode: "${params.copy_mode}", overwrite: true, pattern: "juncfiles.txt"
    publishDir "${params.outdir}/leafcutter_regtools/clustering", mode: "${params.copy_mode}", overwrite: true, pattern: "fofn_junctions_files.txt"
    
    input:
    file (junc_files) //from star_bam2junc.collect()

    when:
    params.star_aligner.star_downstream_tasks.leafcutter_tasks.clustering_regtools_task.run
    
  output:
    tuple file('clust_perind.counts.gz'), file('clust_perind_numers.counts.gz'), file('clust_pooled'),file('clust_refined'),file('clust_sortedlibs'), file('*.junc.clust.sorted.gz') //into star_clustering
    file("fofn_junctions_files.txt")
    // file "*.ReadsPerGene.out.tab" into ch_merge_starcounts

    script:
  """
  ls . | grep .junc\$ > fofn_junctions_files.txt
  # pre regtools leafcutter leafcutter_cluster.py -j fofn_junctions_files.txt -m 50 -o clust -l 500000
  
  # path to leafcutter script in container:
  python2 /opt/conda/envs/conda_rna_seq/bin/leafcutter/clustering/leafcutter_cluster_regtools.py -j fofn_junctions_files.txt --checkchrom -m 50 -o clust -l 500000  

  """
}
