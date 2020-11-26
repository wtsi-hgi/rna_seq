process 'iget_study_cellranger' {
    tag "$sample"
    publishDir "${params.outdir}/iget_study_cellranger/${study_id}/${sample}/", mode: "${params.copy_mode}"
    
    when: 
    params.run_iget_study_cellranger

    input:
    tuple val(study_id), val(sample), val(cellranger_irods_object)
    
  output:
    tuple val(study_id), val(sample), path("cellranger_${sample}/*"), emit: study_sample_cellranger

  script:
    """
iget -K -r -v ${cellranger_irods_object} cellranger_${sample}
echo \"${cellranger_irods_object}\" > cellranger_${sample}/irods_cellranger_path.txt
   """
}
