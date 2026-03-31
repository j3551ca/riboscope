process qc_filter {

    tag { "Applying sample QC" }
    publishDir "${params.outdir}", pattern: "${params.qc_output}.csv", mode: 'copy'

    input:
    tuple path(read_results), path(alignment_results), path(samtools_results),
     path(rrna_counts), path(vcf_results), path(amp_bed), path(kraken_results)

    output:
    tuple path("${params.qc_output}.csv"), path("reportable_vcf.tsv"), path("amplicon_counts.csv"), path("failed_amplicons.csv"), emit: qc_results

    script:
    """
    filter_qc.py \
    --fastp_qc ${read_results} \
    --bam_qc ${alignment_results} \
    --samtools_stats ${samtools_results} \
    --kraken2_tsv ${kraken_results} \
    --amplicon_counts ${rrna_counts} \
    --vcf_file ${vcf_results} \
    --amplicon_bed ${amp_bed} \
    --min_pathogen ${params.min_pathogen} \
    --max_host ${params.max_host} \
    --min_q20 ${params.min_q20_rate} \
    --min_q30 ${params.min_q30_rate} \
    --min_bq ${params.min_mean_bq} \
    --min_depth ${params.min_mean_depth} \
    --min_map_pair ${params.min_mapped_paired} \
    --min_map ${params.min_mapped_reads} \
    --min_pct_map ${params.min_percent_mapped} \
    --min_mq ${params.min_mean_mq} \
    --min_10x ${params.min_10X_cov} \
    --min_50x ${params.min_50X_cov} \
    --max_secondary ${params.max_secondary} \
    --min_pos_count ${params.min_pos_count} \
    --min_med_amp_count  ${params.min_amplicon_count} \
    --max_qc_flags ${params.max_qc_flags} \
    --qc_output ${params.qc_output} \
    --groups ${params.amplicon_groups}
    """
}

process report_results {

    tag { "Generating report" }
    publishDir "${params.outdir}", pattern: "*_report.html", mode: "copy"

    input:
    tuple path(qc_summary), path(annotated_vcf), path(amplicon_counts), path(failed_amplicons), path(kraken2_tsv), path(complete_gff)

    output:
    path("*.html")

    script:
    """
    report_results.py \
    --qc_summary ${qc_summary} \
    --reportable_vcf ${annotated_vcf} \
    --kraken2_tsv ${kraken2_tsv} \
    --amplicon_counts ${amplicon_counts} \
    --failed_amplicons ${failed_amplicons} \
    --html_template ${params.html_template} \
    --ref ${params.ref} \
    --gff ${complete_gff} \
    --n_qc_flag ${params.max_qc_flags} \
    --min_med_amp_count ${params.min_amplicon_count} \
    --reporting_vaf ${params.reporting_vaf} \
    --reporting_sample_perc ${params.reporting_sample_perc}

    #- high level summary: # passed samples, # failed, new mutations not seen before, presence/ absence mutations in known sites (23S)
    """


}