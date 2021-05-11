process 'star_2pass_merge_junctions' {
    tag "STAR2_merge_junc"
    
    publishDir "${params.outdir}/star2pass_merge_junc/", mode: "${params.copy_mode}", overwrite: true, pattern: "SJ.filtered.tab"

  input:
    file(tabs) // from ch_tab_tabs_only.collect()

    when:
    params.star_aligner.star_custom_2pass_task.run_merge_junctions

  output:
    file "SJ.filtered.tab" //into ch_tab_filter_out

  script:

  """
  export TMPDIR=/tmp/  
  cat *.tab | awk '(\$5 > 0 && \$7 > 2 && \$6==0)' | cut -f1-6 | sort | uniq > SJ.filtered.tab
  """
}
