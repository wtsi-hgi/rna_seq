process 'fastqc' {
    tag "fastqc ${samplename}"
    
    publishDir "${params.outdir}/fastqc/", mode: "${params.copy_mode}", overwrite: true,
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when:
    params.fastqc_task.run

    input:
    tuple val(samplename), file(reads)
    
    
    output:
        tuple file("*fastqc.zip"), file("*fastqc.html"), file("*fastqc.html")

        
    
  script:
  """
  fastqc -t ${task.cpus} -q $reads
  """
}
