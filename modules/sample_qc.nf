process collect_qc {

    tag { "sample QC" }
    publishDir "${params.outdir}", pattern: "*.csv", mode: 'copy'

    input:
    tuple path(read_results), path(alignment_results), path(samtools_results),
     path(rrna_counts), path(vcf_results), path(amp_bed)

    output:
    path("*.csv")

    script:
    """
    filter_qc.py \
    --fastp_qc ${read_results} \
    --bam_qc ${alignment_results} \
    --samtools_stats ${samtools_results} \
    --amplicon_counts ${rrna_counts} \
    --vcf_qc ${vcf_results} \
    --amplicon_bed ${amp_bed} \
    --min_q20 ${params.min_q20_rate} \
    --min_q30 ${params.min_q30_rate} \
    --min_bq ${params.min_mean_bq} \
    --min_depth ${params.min_mean_depth} \
    --min_map_pair ${params.min_mapped_paired} \
    --min_map ${params.min_mapped_reads} \
    --min_pct_map ${params.min_percent_mapped} \
    --min_mq ${params.min_mean_mq} \
    --min10x ${params.min_10X_cov} \
    --min_50x ${params.min_50X_cov} \
    --max_secondary ${params.max_secondary} \
    --min_pos_count ${params.min_pos_count} \
    --min_med_count  ${params.min_amplicon_count} \
    --max_qc_flags ${params.max_qc_flags} \
    --qc_output ${params.qc_output}
    """
}

process apply_qc {

    tag { sample_id }
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_fastp*", mode: 'copy'

    input:

    output:

    script:
    """
    - read quality = collected_fastp.csv
    - amplicon presence/absence = collected_rrna_counts.csv -> rRNA seq counts or heatmap (only fail if 2 missing; how about if 1 pool from 1 amplicon is missing?)
    - amplicon genome completeness = edit plot-amplicon-coverage.py from plot_amplicon_coverage to output tsv of completeness above threshold
    - alignment quality = collected_qualimap_alignment_qc.csv
    - vcf = collected_lofreq.vcf - no snps found ? remove or keep? interesting 
    - 
    """


}

process report_results {

    tag { sample_id }
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_results*", mode: 'copy'

    input:

    output:


    script:
    """
    - high level summary: # passed samples, # failed, new mutations not seen before, presence/ absence mutations in known sites (23S)
    - amplicons 1 - 4 present/ absent --> heatmap
    - mutations color-coded by which amplicons are present in sample (ie. amplicon 1 only, 2 only, both/ ambiguous)
    - reason failed samples failed
    - done
    """


}