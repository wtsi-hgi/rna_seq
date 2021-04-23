process 'leafcutter_bam2junc_regtools' {
    tag "${samplename}"
    
    publishDir "${params.outdir}/leafcutter_regtools/bam2junc", mode: "${params.copy_mode}", overwrite: true, pattern: "*.junc"
    // publishDir "${params.outdir}/leafcutter/bam2junc", mode: 'copy', pattern: "*.bam.bed"
    
  input:
    tuple val(samplename), file (bamfile), file (baifile) //from star_aligned_with_bai

    when:
    params.star_aligner.star_downstream_tasks.leafcutter_tasks.bam2junc_regtools_task.run
    
  output:
    file ('*.junc')

  script:

  """
  export PATH=/home/leafcutter/scripts:/home/leafcutter/clustering:/home/regtools/build:/home/regtools/scripts:/opt/conda/envs/conda_leafcutter/bin:/opt/conda/bin:\$PATH

  echo Converting ${bamfile} to ${samplename}.junc
  # pre regootls method: sh /home/leafcutter/scripts/bam2junc.sh ${bamfile} ${samplename}.junc
         
  regtools junctions extract -s 1 -a 8 -m 50 -M 500000 ${bamfile} -o ${samplename}.junc    
  """
}
