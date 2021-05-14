process gz_fastq {
    tag "gz fastq ${samplename}"

    publishDir "${params.outdir}/fastq_gz/", mode: "${params.copy_mode}", overwrite: true, pattern: "*.fastq.gz"

    input: 
    path samplename

    output: 
    path(outputname), emit: samples_gz
    

    script:

    outputname = "${samplename}.gz"
    """
      python $workflow.projectDir/../bin/copy_fastq.py -f ${samplename} -o ${outputname}
    """
}