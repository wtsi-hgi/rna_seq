process multiqc {

    publishDir "${params.outdir}", mode: 'copy',
      saveAs: {filename ->
          if (filename.indexOf("multiqc.html") > 0) "combined/$filename"
          else if (filename.indexOf("_data") > 0) "$filename"
          else null
      }

    when:
    params.multiqc_task.run

    input:
    file ('lostcause/*') //from ch_multiqc_lostcause.collect().ifEmpty([])
    file (fastqc:'fastqc/*') //from ch_multiqc_fastqc.collect().ifEmpty([])
    file ('mapsummary/*') //from ch_multiqc_mapsum.collect().ifEmpty([])
    file ('featureCounts/*') //from ch_multiqc_fc_aligner.collect().ifEmpty([])
    file ('featureCounts_biotype/*') //from ch_multiqc_fcbiotype_aligner.collect().ifEmpty([])
    file ('star/*') //from ch_alignment_logs_star.collect().ifEmpty([])
    file ('salmon/*') //from ch_alignment_logs_salmon.collect().ifEmpty([])

    output:
    file "*_multiqc.html"
    file "*_data"

    script:
    def filename = "multiqc.html"
    def reporttitle = "multiqc"
    """
    multiqc . -f --title "$reporttitle" --filename "$filename" -m featureCounts -m star -m fastqc -m salmon
    """
}
// multiqc . -f --title "$reporttitle" --filename "$filename" -m custom_content -m featureCounts -m star -m fastqc -m salmon
