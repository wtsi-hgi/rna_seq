process fastq_to_csv{

    tag "All file pairs paths to CSV"
    publishDir "${params.outdir}/fastq_gz_csv/", mode: "${params.copy_mode_fastaq}", overwrite: true, pattern: "*.csv"

    input: 
    val sample_paths

    output: 
    path('FastQ_files.csv'), emit: all_samples_gz
    
    script:

    
    """
      python $workflow.projectDir/../bin/fastq_samples_csv.py -f '${sample_paths}'
    """
}