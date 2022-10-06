process gz_fastq {
    tag "gz fastq2 ${samplename}"

    publishDir "${params.outdir}/fastq_gz/", mode: "${params.copy_mode}", overwrite: true, pattern: "*.fastq.gz"

    input: 
    path samplename

    output: 
    path(outputname), emit: samples_gz
    

    script:
    d='34' //for debugging changing this value will bypass the catche
    outputname = "${samplename}.gz"
    outputname=outputname.replaceAll(".gz.gz", ".gz")
    """
      python $workflow.projectDir/bin/copy_fastq.py -f ${samplename} -o ${outputname}
    """
}