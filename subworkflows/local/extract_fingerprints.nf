include { PICARD_ADDORREPLACEREADGROUPS } from '../../modules/nf-core/picard/addorreplacereadgroups/main'
include { GATK4_EXTRACTFINGERPRINT } from '../../modules/local/gatk4/extractfingerprints/main'

workflow EXTRACT_FINGERPRINTS {
    take:
        bam                 // tuple val(meta), path(bam)
        reference_sequence  // tuple path(reference_sequence), path(reference_index), path(reference_dict)
        haplotype_map       // path(haplotype_map)

    main:
        PICARD_ADDORREPLACEREADGROUPS(bam)
        GATK4_EXTRACTFINGERPRINT(PICARD_ADDORREPLACEREADGROUPS.out.bam, reference_sequence, haplotype_map)

    emit:
        vcf = GATK4_EXTRACTFINGERPRINT.out.vcf
}
