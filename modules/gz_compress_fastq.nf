process 'fastqc' {
    tag "fastqc ${samplename}"
    
    publishDir "${params.outdir}/fastqc/", mode: "${params.copy_mode}", overwrite: true,
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when:
    params.fastqc_task.run

    input:
    tuple val(samplename), file(reads)
    
    output:
    tuple file("*1_fastqc.zip"), file("*2_fastqc.zip"), file("*1_fastqc.html"), file("*2_fastqc.html")

  script:
  """
  fastqc -t ${task.cpus} -q $reads
  """
}
