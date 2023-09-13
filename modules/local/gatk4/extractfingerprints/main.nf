process GATK4_EXTRACTFINGERPRINT {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/fingerprints/", pattern: '*.vcf', mode: 'copy', overwrite: true

    conda "bioconda::gatk4=4.4.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.4.0.0--py36hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.4.0.0--py36hdfd78af_0' }"

    input:
        tuple val(meta), path(bam)
        tuple path(reference_sequence), path(reference_index), path(reference_dict)
        path(haplotype_map)

    output:
        tuple val(meta), path('*.vcf'), emit: vcf
        path "versions.yml"           , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"

        def avail_mem = 3072
        if (!task.memory) {
            log.info '[GATK ExtractFingerprint] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
        } else {
            avail_mem = (task.memory.mega*0.8).intValue()
        }
        """
        gatk --java-options "-Xmx${avail_mem}M" ExtractFingerprint \\
            --INPUT ${bam} \\
            --HAPLOTYPE_MAP ${haplotype_map} \\
            --REFERENCE_SEQUENCE ${reference_sequence} \\
            --OUTPUT ${prefix}.vcf \\
            --TMP_DIR . \\
            $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
        END_VERSIONS
        """
}
