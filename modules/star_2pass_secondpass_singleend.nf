process 'star_2pass_2nd_pass' {
    tag "${samplename}"
    
    publishDir "${params.outdir}/star_pass2_2ndpass/$samplename", mode: "${params.copy_mode}", overwrite: true, pattern: "*.bam"
    publishDir "${params.outdir}/star_pass2_2ndpass/$samplename", mode: "${params.copy_mode}", overwrite: true, pattern: "*.bam.bai"
    publishDir "${params.outdir}/star_pass2_2ndpass/$samplename", mode: "${params.copy_mode}", overwrite: true, pattern: "*.out"
    publishDir "${params.outdir}/star_pass2_2ndpass/$samplename", mode: "${params.copy_mode}", overwrite: true, pattern: "*.tab"
    publishDir "${params.outdir}/star_pass2_2ndpass/$samplename", mode: "${params.copy_mode}", overwrite: true, pattern: "*Unmapped.out.mate1"
    
    publishDir "${params.outdir}/star_pass2_2ndpass_multiqc/", mode: "${params.copy_mode}", overwrite: true,
        saveAs: { filename ->
            if (filename ==~ /.*\.out\.tab/) "STARcounts/$filename"
            else if (filename.indexOf(".bam") == -1) "STARlogs/$filename"
            else null
        }

  input:
    tuple val(samplename), file(reads) //from ch_star // _reads_only
    file genomeDir //from ch_star_index.collect()
    file gtf //from ch_gtf_star.collect()
    file filtered_tab // merged and filtered junctions from 1st pass

    when:
    params.star_aligner.star_custom_2pass_task.run_2nd_pass
    
  output:
    tuple val(samplename), file ('*.bam'), file ('*.bai') //into star_aligned_with_bai
    tuple val(samplename), file("*Log.final.out"), file ('*.bam') //into star_aligned
    // file "*.ReadsPerGene.out.tab" into ch_merge_starcounts
    tuple file("*.Log.final.out"), file("*.Log.out"), file("*.progress.out") //into ch_alignment_logs_star
    file "*.SJ.out.tab"
    tuple val(samplename), file("*Unmapped.out.mate1") //into star_aligned

  script:
log.info 'single end'
  """
# 2nd pass:

STAR --genomeDir ${genomeDir} \\
--readFilesIn ${reads} \\
--runThreadN ${task.cpus} \\
--readFilesCommand zcat \\
--limitSjdbInsertNsj 10000000 \\
--sjdbFileChrStartEnd ${filtered_tab} \\
--outSAMtype BAM SortedByCoordinate \\
--outFileNamePrefix ${samplename}. \\
--outFilterType BySJout \\
--outFilterMultimapNmax 20 \\
--alignSJoverhangMin 8 \\
--alignSJDBoverhangMin 1 \\
--outFilterMismatchNmax 999 \\
--outFilterMismatchNoverReadLmax 0.04 \\
--alignIntronMin 20 \\
--alignIntronMax 1000000 \\
--alignMatesGapMax 1000000 \\
--outReadsUnmapped Fastx \\
--limitBAMsortRAM 31511664535

  # Index the BAM file
  samtools index ${samplename}.Aligned.sortedByCoord.out.bam
  """
}
