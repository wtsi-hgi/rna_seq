process 'star_2pass_1st_pass' {
    tag "${samplename}"

    publishDir "${params.outdir}/star_pass2_1stpass/$samplename", mode: "${params.copy_mode}", overwrite: true, pattern: "*.out"
    publishDir "${params.outdir}/star_pass2_1stpass/$samplename", mode: "${params.copy_mode}", overwrite: true, pattern: "*.tab"
    publishDir "${params.outdir}/star_pass2_1stpass_multiqc/", mode: "${params.copy_mode}", overwrite: true,
        saveAs: { filename ->
            if (filename ==~ /.*\.out\.tab/) "STARcounts/$filename"
            else if (filename.indexOf(".bam") == -1) "STARlogs/$filename"
            else null
        }
    
  input:
    tuple val(samplename), file(reads) //from ch_star // _reads_only
    file genomeDir //from ch_star_index.collect()
    file gtf //from ch_gtf_star.collect()

    when:
    params.star_aligner.star_custom_2pass_task.run_1st_pass
    
  output:
    tuple val(samplename), file("*Log.final.out")
    file "*.SJ.out.tab"
    tuple file("*.Log.final.out"), file("*.Log.out"), file("*.progress.out") //into ch_alignment_logs_star

  script:

  """
# first pass:

STAR --genomeDir ${genomeDir} \\
--sjdbGTFfile $gtf \\
--readFilesIn $reads --readFilesCommand zcat \\
--runThreadN ${task.cpus} \\
--outSAMtype BAM Unsorted \\
--outFileNamePrefix ${samplename}. \\
--outFilterType BySJout \\
--outFilterMultimapNmax 20 \\
--alignSJoverhangMin 8 \\
--alignSJDBoverhangMin 1 \\
--outFilterMismatchNmax 999 \\
--outFilterMismatchNoverReadLmax 0.04 \\
--alignIntronMin 20 \\
--alignIntronMax 1000000 \\
--alignMatesGapMax 1000000

  rm *.bam
  """
}
