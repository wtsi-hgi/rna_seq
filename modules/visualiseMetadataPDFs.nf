params.path="Submission_Data_Pilot_UKB.file_metadata.tsv"

process visualiseMetadata{
  when:
  params.run_metadata_visualisation

  input:
  path(my_channel)

  script:

  """
  mkdir $workflow.projectDir/../../results/PDFs
  python $workflow.projectDir/../bin/main_PDF_visualisation.py $my_channel $workflow.projectDir/../../results/PDFs

  """

}
