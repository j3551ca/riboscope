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
    tuple val(sample_id), path("ref_*{.fa,.dict,.fai}"), emit: ref_dict

    script:
    """
    cp ${ref} ref_.fa
    gatk CreateSequenceDictionary -R ref_.fa
    samtools faidx ref_.fa
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
    ${alignment[0]} \
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
   
    
    gatk BaseRecalibrator \
    -I ${alignment[0]} \
    -R ${ref_dict[0]} \
    --known-sites ${expected_vcf} \
    -O ${sample_id}.recalibrated.data.table
    
    #cannot pipe to next - tried

    gatk ApplyBQSR \
    -I ${alignment[0]}\
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
    tuple val(sample_id), path("${sample_id}.dindel.sorted{.bam,.bai}"), emit: indel_alignment

    script:
    """ 
    #needed if want indels in lofreq (does not work with just ApplyBQSR then lofreq call --call-indels)
    lofreq indelqual \
    --dindel \
    --ref ${ref[0]} \
    -o ${sample_id}.dindel.bam \
    ${recalibrated_alignment[0]}

    samtools sort -o ${sample_id}.dindel.sorted.bam ${sample_id}.dindel.bam
    samtools index ${sample_id}.dindel.sorted.bam

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
    """
    lofreq call \
        -f ${ref[0]} \
        -o ${sample_id}.lofreq.variants.vcf \
        --sig 0.05 \
        --min-cov ${params.min_lofreq_cov} \
        --call-indels ${indelqual_alignment[0]}
    """
}

process assemble_haplotypes {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}.haplotypes*", mode: "copy"

    input:
    tuple val(sample_id), path(indelqual_alignment), path(lofreq_vars), path(ref)

    output:

    script:
    """
    header=\$(head -n1 ${ref[0]} | cut -f1 -d ' ' )
    ref_name=\${header/>/}

    HAT \
    -rl 150 \
    --haplotype_assembly True\
    -r ${ref[0]} \
    --chromosome_name "\${ref_name}" \
    --vcf_file ${lofreq_vars} \
    --short_read_alignment ${indelqual_alignment[0]} \
    --ploidy 1 \
    --output "${sample_id}.haplotypes" 

    #too slow
    savage \
    --ref /home/jess.cal/syphilis/test_output/illumina/amr_ribo_rpts/bwa_ref_types/ribo_rpt_1/riboscope_test/TR13116/ref.fa \
    -p1 TR13116_trimmed_R1.fastq \
    -p2 TR13116_trimmed_R2.fastq \
    -m 200 \
    --split 2 \
    --revcomp

    fc-virus \
    -k 100 \
    -p \
    -t fq \
    --left TR13116_trimmed_R1.fastq\
    --right TR13116_trimmed_R2.fastq

    shorah amplicon \
    -b TR13116.mapped.primertrimmed.sorted.bam \
    -f ref.fa \
    -r NC_000919.1_-_TpRiboCore1:1-100 #window slightly less than read length of 150

    # this one:
    shorah shotgun \
    -b TR13116.mapped.primertrimmed.sorted.bam \
    -f ref.fa \
    -r NC_000919.1_-_TpRiboCore1:4000-5000 \
    -w 120 \
    -c 0 \
    -x 1000000

    haploflow \
    --read-file ./TR13116_combined.fastq \
    --out test \
    --log test/log

    """

}

process map_contigs {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}.variants.tsv", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs), path(ref)

    output:
    tuple val(sample_id), path("${sample_id}.contigs.sorted{.bam,.bai}")

    script:
    """
    minimap2 \
    -a -x asm5 \
    ${ref} \
    contigs_stage_c.fasta \
    > mapped_contigs.sam

    samtools sort -o ${sample_id}.contigs.sorted.bam mapped_contigs.sam
    samtools index ${sample_id}.contigs.sorted.bam

    #extract contigs overlap in 16S region
    samtools view -L 16S_region.bed TR13116_minimap_contigs.bam \
    | cut -f1 \
    | sort -u \
    | seqkit grep -f - contigs_stage_c.fasta \
    > 16S_contigs.fasta


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