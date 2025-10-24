process fastp {

    tag { sample_id }

    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_fastp*.{csv,txt}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("${sample_id}_fastp.json"), emit: fastp_json
    tuple val(sample_id), path("${sample_id}_fastp.csv"), emit: fastp_csv
    //tuple val(sample_id), path("${sample_id}_fastp.html"), emit: fastp_html
    tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq.gz"), path("${sample_id}_trimmed_R2.fastq.gz"), emit: trimmed_reads
    tuple val(sample_id), path("${sample_id}_fastp_provenance.yml"), emit: provenance

    script:
    """
    printf -- "- process_name: fastp\\n"  >> ${sample_id}_fastp_provenance.yml
    printf -- "  tools:\\n"               >> ${sample_id}_fastp_provenance.yml
    printf -- "    - tool_name: fastp\\n" >> ${sample_id}_fastp_provenance.yml
    printf -- "      tool_version: \$(fastp --version 2>&1 | cut -d ' ' -f 2)\\n" >> ${sample_id}_fastp_provenance.yml
    printf -- "      parameters:\\n"      >> ${sample_id}_fastp_provenance.yml
    printf -- "        - parameter: --cut_tail\\n" >> ${sample_id}_fastp_provenance.yml
    printf -- "          value: null\\n" >> ${sample_id}_fastp_provenance.yml

   

    fastp -w ${task.cpus} \
    -i ${reads_1} \
    -I ${reads_2} \
    -o ${sample_id}_trimmed_R1.fastq.gz \
    -O ${sample_id}_trimmed_R2.fastq.gz\
    --detect_adapter_for_pe \
    --failed_out ${sample_id}_fastp_failed_reads.txt 
    #--html ${sample_id}_fastp.html 

    mv fastp.json ${sample_id}_fastp.json
    fastp_json_to_csv.py -s ${sample_id} ${sample_id}_fastp.json > ${sample_id}_fastp.csv
    """
}

process detect_ribo_repeats {

    tag { sample_id }
    
    publishDir "${params.outdir}/${sample_id}", pattern: "${sample_id}_rrna_counts.csv", mode: 'copy'

    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(search_seqs)

    output:
    tuple val(sample_id), path("${sample_id}_rrna_counts.csv"), emit: ribo_rpt_counts


    script:
    """
    seq_names=\$(grep "^>" ${search_seqs} | cut -d" " -f1 | sed 's/^>//' | paste -sd, -)
    echo -e "sample_id,\${seq_names}" > ${sample_id}_rrna_counts.csv

    counts=()
    for seq in \$(awk '!/^>/' ${search_seqs});
    do  
       count=\$( seqkit grep \
        -s \
        -i \
        -C \
        -p "\${seq}" \
        ${reads_1} \
        ${reads_2} ) \
        counts+=("\${count}")
    done

    IFS=','; echo "${sample_id},\${counts[*]}" >> ${sample_id}_rrna_counts.csv
    """
}

process kraken2 {

    tag { sample_id }

    publishDir "${params.outdir}/kraken2_output", pattern: "${sample_id}_kraken*.txt", mode: 'copy'
    input:
    tuple val(sample_id), path(reads_1), path(reads_2), path(kraken2_db), val(analysis_stage)

    output:
    tuple val(sample_id), path("${sample_id}_kraken*.txt"), emit: kraken2_report
    
    script:
    """
    kraken2 \
      --threads ${task.cpus} \
      --db ${kraken2_db} \
      --report ${sample_id}_kraken_report_${analysis_stage}.txt \
      --paired \
      ${reads_1} \
      ${reads_2}

    if [ "\$(wc -l < ${sample_id}_kraken_report_${analysis_stage}.txt)" -eq 0 ]; then
        echo -e "0.0\t0\t0\tU\t0\tunclassified" >> ${sample_id}_kraken_report_${analysis_stage}.txt
    fi
    """
}

process bracken {

    tag { sample_id }

    publishDir "${params.outdir}/bracken_output", pattern: "${sample_id}_bracken*${analysis_stage}.txt", mode: 'copy'

    input:
    tuple val(sample_id), path(kraken_report), path(bracken_db), val(read_length), val(taxonomy_level), val(analysis_stage)

    output:
    tuple val(sample_id), path("${sample_id}_bracken*${analysis_stage}.txt")
    
    script:
    """
    if ! bracken \
      -d ${bracken_db} \
      -i ${kraken_report} \
      -l ${taxonomy_level} \
      -o ${sample_id}_bracken_output_${analysis_stage}.txt \
      -w ${sample_id}_bracken_report_${analysis_stage}.txt \
      -r ${read_length} 2> bracken.err; then
      if grep "Error: no reads found. Please check your Kraken report" bracken.err; then
        echo -e "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads" >  ${sample_id}_bracken_output_${analysis_stage}.txt
        echo -e "none\t0\t${taxonomy_level}\t0\t0\t0\t0.0" >> ${sample_id}_bracken_output_${analysis_stage}.txt
      else 
        echo "Bracken failed for a reason unrelated to insufficient reads."
        exit 1
      fi
    fi
    """
}