process deseq2 {
    tag "$deseq2_tsv"
    
    publishDir "${params.outdir}/DESeq2/", mode: "${params.copy_mode}", overwrite: true
    
    when:
    params.salmon_aligner.salmon_dowstream_tasks.deseq2_task.run

    input:
    file(quant_sf_files)  // from collect()
    file(deseq2_tsv)

    output:
    file("outputs")

    script:
    """
    ls . | grep .quant.sf\$ > fofn_quantfiles.txt
    Rscript $workflow.projectDir/../bin/deseq2.R \"$params.ensembl_lib\" fofn_quantfiles.txt $deseq2_tsv
    """
}

// export http_proxy=http://wwwcache.sanger.ac.uk:3128
// export https_proxy=http://wwwcache.sanger.ac.uk:3128
// export HTTP_PROXY=http://wwwcache.sanger.ac.uk:3128
// export HTTPS_PROXY=http://wwwcache.sanger.ac.uk:3128
