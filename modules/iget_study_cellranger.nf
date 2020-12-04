process 'iget_study_cellranger' {
    tag "$sample"
    publishDir "${params.outdir}/iget_study_cellranger/${study_id}/${sample}/", mode: "${params.copy_mode}"
    
    when: 
    params.run_iget_study_cellranger

    input:
    tuple val(study_id), val(sample), val(cellranger_irods_object)
    
  output:
    tuple val(study_id), val(sample), path("cellranger_${sample}/*"), emit: study_sample_cellranger
    tuple val(sample), env(RESULTS_DIR), emit: cellranger_output_dir
    path("${sample}.metadata.tsv"), emit: cellranger_metadata_tsv

  script:
    """
iget -K -r -v ${cellranger_irods_object} cellranger_${sample}
echo \"${cellranger_irods_object}\" > cellranger_${sample}/irods_cellranger_path.txt

RESULTS_DIR=${params.outdir}/iget_study_cellranger/${study_id}/${sample}/cellranger_${sample}

# prepare metadata tsv row for that sample:
echo sample_sanger_id,experiment_id,irods_cellranger_path > metadata1.csv 
echo ${sample},${sample},${cellranger_irods_object} >> metadata1.csv 
cat cellranger_${sample}/metrics_summary.csv | $workflow.projectDir/../bin/remove_comma_inside_quotes.sh  > metadata2.csv
paste -d, metadata1.csv metadata2.csv | tr ',' '\t' > ${sample}.metadata.tsv
rm metadata1.csv 
rm metadata2.csv 
   """
}

// paste -d, metadata1.csv cellranger_${sample}/metrics_summary.csv | tr ',' '\\t' 
