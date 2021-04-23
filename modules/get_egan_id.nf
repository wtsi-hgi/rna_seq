process get_egan_id {
    tag "$samplename"
    //publishDir "${params.outdir}/egan_id", mode: "${params.copy_mode}", overwrite: true

    when:
    params.get_egan_id_task.run

    input:
    tuple val(samplename), file(cram)

    output:
    tuple val(samplename), file(cram), file("${samplename}_egan_id.csv"), emit: samplename_egan_id_csv

    script:
    """
EGANID=\$(samtools view -H $cram | \\
    grep 'SM:EGAN' | \\
    sed s'/^.*SM:/SM:/'g | \\
    sed s'/SM://'g | \\
    awk '{print \$1;}')

echo samplename,egan_id > ${samplename}_egan_id.csv
echo ${samplename},\$EGANID >> ${samplename}_egan_id.csv
    """
}
