process tximport {
    tag "$params.ensembl_lib"
    
    publishDir "${params.outdir}/tximport", mode: "${params.copy_mode}", overwrite: true

    when:
    params.salmon_aligner.salmon_downstream_tasks.tximport_task.run

    input:
    file (quant_sf_files)  // from collect()

    output:
    file("fofn_quantfiles.txt")
    file("txi_gene_counts.csv")
    file("txi_transcript_counts.csv")
    file("txi_lengthScaledTPM_gene_counts.csv")
    file("tximport.rdata")
    //file "${samplename}.quant.sf" // into ch_salmon_trans
    //file "${samplename}.quant.genes.sf" //into ch_salmon_genes
    // file "my_outs/${samplename}" optional true // into ch_alignment_logs_salmon

    script:
    """
    export HOME=\$(pwd)
    ls . | grep .quant.sf\$ > \$HOME/fofn_quantfiles.txt
    mkdir -p \$HOME/tmp_r
    export TMPDIR=$HOME/tmp_r
    Rscript $workflow.projectDir/bin/tximport.R \
        \"${params.salmon_aligner.salmon_downstream_tasks.tximport_task.ensembl_lib}\" \
        \"${params.salmon_aligner.salmon_downstream_tasks.tximport_task.ensembl_species}\" \
        \$HOME/fofn_quantfiles.txt
    rm -r \$HOME/tmp_r    
    """
}
