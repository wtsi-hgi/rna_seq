process 'star_2pass_basic' {
    tag "${samplename}"
    
    publishDir "${params.outdir}/star_pass2_basic/$samplename", mode: "${params.copy_mode}", overwrite: true, pattern: "*.bam"
    publishDir "${params.outdir}/star_pass2_basic/$samplename", mode: "${params.copy_mode}", overwrite: true, pattern: "*.bam.bai"
    publishDir "${params.outdir}/star_pass2_basic/$samplename", mode: "${params.copy_mode}", overwrite: true, pattern: "*.out"
    publishDir "${params.outdir}/star_pass2_basic/$samplename", mode: "${params.copy_mode}", overwrite: true, pattern: "*.tab"
    
    publishDir "${params.outdir}/star_pass2_basic_multiqc/", mode: "${params.copy_mode}", overwrite: true,
        saveAs: { filename ->
            if (filename ==~ /.*\.out\.tab/) "STARcounts/$filename"
            else if (filename.indexOf(".bam") == -1) "STARlogs/$filename"
            else null
        }

  input:
    tuple val(samplename), file(reads) //from ch_star // _reads_only
    file genomeDir //from ch_star_index.collect()
    // file genome_fasta3 from ch_dna_star.collect()
    file gtf //from ch_gtf_star.collect()

    when:
    params.star_aligner.star_2pass_basic_task.run
    
  output:
    tuple val(samplename), file ('*.bam'), file ('*.bai') //into star_aligned_with_bai
    tuple val(samplename), file("*Log.final.out"), file ('*.bam') //into star_aligned
    // file "*.ReadsPerGene.out.tab" into ch_merge_starcounts
    tuple file("*.Log.final.out"), file("*.Log.out"), file("*.progress.out") //into ch_alignment_logs_star
    file "*.SJ.out.tab"
    tuple val(samplename), file("${samplename}.ReadsPerGene.out.tab"), emit: samplename_readspergene_tab

  script:

  """
    STAR --genomeDir ${genomeDir} \\
        --sjdbGTFfile $gtf \\
        --readFilesIn $reads --readFilesCommand zcat \\
        --runThreadN ${task.cpus} \\
        --twopassMode Basic \\
        --outWigType bedGraph \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMunmapped Within \\
        --runDirPerm All_RWX \\
        --quantMode GeneCounts \\
        --outFileNamePrefix ${samplename}.

  # Index the BAM file
  samtools index ${samplename}.Aligned.sortedByCoord.out.bam
  rm -f Log.out 
  rm -f log.out 
  """
}
