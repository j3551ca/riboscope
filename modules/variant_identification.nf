process remove_chimeras {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}.variants.vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(alignment), path(ref), path(bed)

    output:
    tuple val(sample_id), path("${sample_id}.variants.vcf")

    script:
    """
    
    """
}

process ref_dict {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "ref.fa", mode: 'copy'

    input:
    tuple val(sample_id), path(ref)

    output:
    tuple val(sample_id), path("ref_*{.fa,.dict}"), emit: ref_dict

    script:
    """
    cp ${ref} ref_.fa
    gatk CreateSequenceDictionary -R ref_.fa
    """
}

process expected_snps {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}.expected.snps.vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(alignment), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}.expected.snps.vcf"), emit: expected_vcf

    script:
    """
    bcftools mpileup \
    -Ou \
    -d 100000 \
    -f ${ref[0]} \
    -L 100000 \
    ${alignment} \
    | bcftools call \
    -mv -Ou \
    | bcftools view \
    -q 0.01 \
    -o ${sample_id}.expected.snps.vcf

    """
}

process recalibrate_bq {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}.variants.vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(alignment), path(ref_dict), path(expected_vcf)

    output:
    tuple val(sample_id), path("${sample_id}.recalibrated.sorted.bam*"), emit: recalibrated_alignment

    script:
    """
    # index feature file
    gatk IndexFeatureFile -I ${expected_vcf}
    #gatk CreateSequenceDictionary -R ${ref[0]}
    
    gatk BaseRecalibrator \
    -I ${alignment} \
    -R ${ref_dict[0]} \
    --known-sites ${expected_vcf} \
    -O ${sample_id}.recalibrated.data.table
    
    #cannot pipe to next - tried

    gatk ApplyBQSR \
    -I ${alignment}\
    -R ${ref_dict[0]} \
    --bqsr-recal-file ${sample_id}.recalibrated.data.table \
    -O recalibrated.bam

    #cannot pipe to next - tried

   samtools sort -o ${sample_id}.recalibrated.sorted.bam recalibrated.bam
   samtools index ${sample_id}.recalibrated.sorted.bam

    
    """
}

process lofreq_indel {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}.dindel.bam", mode: 'copy'

    input:
    tuple val(sample_id), path(recalibrated_alignment), path(ref), path(bed)

    output:
    tuple val(sample_id), path("${sample_id}.dindel.bam"), emit: indel_alignment

    script:
    """ 
    #needed if want indels in lofreq (does not work with just ApplyBQSR then lofreq call --call-indels)
    lofreq indelqual \
    --dindel \
    --ref ${ref[0]} \
    -o ${sample_id}.dindel.bam \
    ${recalibrated_alignment[0]}

    """
}

process lofreq_call {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}.lofreq.variants.vcf", mode: 'copy'

    input:
    tuple val(sample_id), path(indelqual_alignment), path(ref), path(bed)

    output:
    tuple val(sample_id), path("${sample_id}.lofreq.variants.vcf"), emit: minor_alleles

    script:
    lofreq_threads = task.cpus - 2
    """
    lofreq call-parallel \
        --pp-threads ${lofreq_threads} \
        -f ${ref[0]} \
        -o ${sample_id}.lofreq.variants.vcf \
        --sig 0.05 \
        --min-cov ${params.min_lofreq_cov} \
        --call-indels ${indelqual_alignment}
    """
}

process assemble_haplotypes {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}.", mode: "copy"

    input:
    tuple val(sample_id), path(alignment), path(ref)

    output:

    script:
    """
    HAT \
    -r CP048984.1.fna \
    CP048984.1 \
    snp-var.vcf.gz \
    short_reads_alignment.sorted.bam \
    long_reads_alignment.sorted.bam \
    3 \
    haplotypes 
    """

}

process call_variants {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}.variants.tsv", mode: 'copy'

    input:
    tuple val(sample_id), path(alignment), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}.variants.tsv")

    script:
    """
    samtools faidx ${ref}

    samtools mpileup -aa -A -d ${params.max_depth} -B -Q 0 --reference ${ref} ${alignment[0]} \
	| ivar variants \
	-r ${ref} \
	-m ${params.min_depth}  \
	-q ${params.min_qual_for_variant_calling} \
	-t ${params.ambiguous_allele_freq_threshold} \
	-p ${sample_id}.variants
    """
}