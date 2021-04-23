process salmon {
    tag "$samplename"
    
    publishDir "${params.outdir}/salmon", mode: "${params.copy_mode}", overwrite: true

    when:
    params.salmon_aligner.salmon_task.run

    input:
    tuple val(samplename), file(reads) // from ch_salmon
    file salmon_index_dir // from ch_salmon_index.collect()
    // file salmon_trans_gene_txt // from ch_salmon_trans_gene.collect()

    output:
    file "${samplename}.quant.sf" // into ch_salmon_trans
    file "my_outs/${samplename}" // into ch_alignment_logs_salmon
    // file "${samplename}.quant.genes.sf" //into ch_salmon_genes

    script:
    """
    salmon quant \\
        -i ${salmon_index_dir} \\
        -l ISR \\
        -p ${task.cpus} \\
        --seqBias \\
        --gcBias \\
        --posBias \\
        --no-version-check \\
        -q \\
        -o . \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        --useVBOpt \\
        --numBootstraps 100
    mv quant.sf ${samplename}.quant.sf
    mkdir -p my_outs/${samplename}/libParams
    mkdir -p my_outs/${samplename}/aux_info
    ln -f aux_info/meta_info.json my_outs/${samplename}/aux_info/meta_info.json
    ln -f libParams/flenDist.txt  my_outs/${samplename}/libParams/flenDist.txt
    """

    // TODO: prepare columns for merging; extract correct column and transpose (paste) it.
    // Include the row names so merger can check identity.
    // The merge step will concatenate the rows and re-transpose to obtain final result.
}

        // mv quant.genes.sf ${samplename}.quant.genes.sf
        // -g ${salmon_trans_gene_txt} \\
