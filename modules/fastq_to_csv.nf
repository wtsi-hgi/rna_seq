process fastq_to_csv{

    tag "All file pairs paths to CSV"
    publishDir "${params.outdir}/fastq_gz_csv/", mode: "${params.copy_mode_fastaq}", overwrite: true, pattern: "*.csv"

    input: 
    val sample_paths

    when:
    params.input_from_fastq_csv.run

    output: 
    path('Fastq_files.csv'), emit: all_samples_gz
    
    script:
    d='3' //for debugging changing this value will bypass the catche
    
    """
      echo '${sample_paths}' > samplePaths_list.txt
      python $workflow.projectDir/bin/fastq_samples_csv.py -s '${params.star_aligner.star_downstream_tasks.featureCounts_task.singleend}'
    """
}