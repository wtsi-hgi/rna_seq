process 'fastqc' {
    tag "fastqc ${samplename}"
    
    publishDir "${params.outdir}/fastqc/", mode: "${params.copy_mode}", overwrite: true,
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when:
    params.fastqc_task.run

    input:
    tuple val(samplename), file(reads)
    
    output:
    tuple file("*_1_fastqc.zip"), file("*_2_fastqc.zip"), file("*_1_fastqc.html"), file("*_2_fastqc.html")

  script:
  """
  fastqc -t ${task.cpus} -q $reads
  """
}
