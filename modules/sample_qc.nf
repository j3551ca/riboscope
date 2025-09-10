process collect_qc {

    tag { "sample QC" }
    publishDir "${params.outdir}", pattern: "*.txt", mode: 'copy'

    input:
    tuple path(read_results), path(alignment_results), path(rrna_counts), path(vcf_results)

    output:
    path("*.txt")

    script:
    """
    collect_qc_results.py \
    --fastp_qc ${read_results} \
    --bam_qc ${alignment_results} \
    --amplicon_counts ${rrna_counts} \
    --vcf_qc ${vcf_results}
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