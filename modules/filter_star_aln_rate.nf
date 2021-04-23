process filter_star_aln_rate {
    tag "filter_star_aln_rate ${samplename}"
    
    //publishDir "${params.outdir}/fastq12/", mode: 'copy'

    when:
    params.star_aligner.star_downstream_tasks.filter_star_aln_rate_task.run
    
    input:
    tuple val(samplename), file(log_final_out)
    output: 
    tuple val(samplename), stdout
    script:

    """
#!/usr/bin/env python3
import re

with open(\"${log_final_out}\",\"r\") as f:
    for line in f:
        if re.search(\"Uniquely mapped reads %\",line):
                if float(re.findall(\"\\d+\\.\\d+\", line)[0]) > float(${params.star_aligner.star_downstream_tasks.filter_star_aln_rate_task.min_pct_aln}):
                    print('above_threshold', end='')
                else:
                    print('below_threshold', end='')
    """
}
