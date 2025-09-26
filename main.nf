#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { hash_files as hash_ref }         from './modules/hash_files.nf'
include { hash_files as hash_fastq }       from './modules/hash_files.nf'
include { fastp }                          from './modules/short_read_qc.nf'
include { detect_ribo_repeats }            from './modules/short_read_qc.nf'
include { index_ref }                      from './modules/amplicon_consensus.nf'
include { bwa_mem }                        from './modules/amplicon_consensus.nf'
include { trim_primer_sequences }          from './modules/amplicon_consensus.nf'
include { qualimap_bamqc }                 from './modules/amplicon_consensus.nf'
include { samtools_stats }                 from './modules/amplicon_consensus.nf'
include { samtools_mpileup }               from './modules/amplicon_consensus.nf'
include { amplicon_coverage }              from './modules/amplicon_consensus.nf'
include { ref_dict; 
expected_snps; 
recalibrate_bq; 
lofreq_indel;
lofreq_call;
filter_lofreq}                             from './modules/variant_identification.nf'
include { qc_filter }                      from  './modules/sample_qc.nf'
include { report_results }                 from  './modules/sample_qc.nf'
//include { call_variants }                  from './modules/amplicon_consensus.nf'
//include { make_consensus }                 from './modules/amplicon_consensus.nf'
//include { align_consensus_to_ref }         from './modules/amplicon_consensus.nf'
include { plot_coverage }                  from './modules/amplicon_consensus.nf'
include { plot_amplicon_coverage }         from './modules/amplicon_consensus.nf'
include { pipeline_provenance }            from './modules/provenance.nf'
include { collect_provenance }             from './modules/provenance.nf'

workflow {

    ch_workflow_metadata = Channel.value([
	workflow.sessionId,
	workflow.runName,
	workflow.manifest.name,
	workflow.manifest.version,
	workflow.start,
    ])

    ch_pipeline_provenance = pipeline_provenance(ch_workflow_metadata)

    if (params.samplesheet_input != 'NO_FILE') {
	ch_fastq = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['R1'], it['R2']] }.filter{ it -> it[1] != null || it[2] != null }
	ch_ref = Channel.fromPath(params.samplesheet_input).splitCsv(header: true).map{ it -> [it['ID'], it['REF']] }
    } else {
	ch_fastq = Channel.fromFilePairs( params.fastq_search_path, flat: true ).map{ it -> [it[0].split('_')[0], it[1], it[2]] }.unique{ it -> it[0] }
    }


    main:
    ch_sample_ids = ch_fastq.map{ it -> it[0] }

    ch_provenance = ch_sample_ids

    if (params.ref != 'NO_FILE') {
	ch_ref = ch_sample_ids.combine(Channel.fromPath(params.ref))
    } else {
	error "Reference file is required"
    }

    if (params.bed != 'NO_FILE') {
	ch_bed = Channel.fromPath(params.bed)
    } else {
	error "BED file is required"
    }

    if (params.search_seqs != 'NO_FILE') {
	ch_search_seqs = Channel.fromPath(params.search_seqs)
    } else {
	error "File containing sequences to search is required"
    }

    hash_ref(ch_ref.combine(Channel.of("ref-fasta")))
    hash_fastq(ch_fastq.map{ it -> [it[0], [it[1], it[2]]] }.combine(Channel.of("fastq-input")))
    
    ch_indexed_ref = index_ref(ch_ref)

    fastp(ch_fastq)
    detect_ribo_repeats(ch_fastq.combine(ch_search_seqs))

    if (! params.align_untrimmed_reads) {	
	ch_reads_to_align = fastp.out.trimmed_reads
    } else {
	ch_reads_to_align = ch_fastq
    }

    bwa_mem(ch_reads_to_align.join(ch_indexed_ref))

    ch_alignment = bwa_mem.out.alignment

    trim_primer_sequences(ch_alignment.combine(ch_bed))

    ch_primer_trimmed_alignment = trim_primer_sequences.out.primer_trimmed_alignment

    qualimap_bamqc(ch_primer_trimmed_alignment)

    amplicon_coverage(ch_primer_trimmed_alignment.combine(ch_bed))

    ch_amplicon_depths = amplicon_coverage.out.depths

    plot_amplicon_coverage(ch_amplicon_depths)

    samtools_mpileup(ch_primer_trimmed_alignment.join(ch_ref))

    ch_per_base_depths = samtools_mpileup.out.depths

    plot_coverage(ch_per_base_depths.join(ch_ref))

    samtools_stats(ch_primer_trimmed_alignment)

    ch_ref_dict = ref_dict(ch_ref).map{ id, files -> [id, files.sort { it.name.endsWith('.fa') ? 0 : 1 }] }

    expected_snps(ch_primer_trimmed_alignment.join(ch_ref_dict))

    recalibrate_bq(ch_primer_trimmed_alignment.join(ch_ref_dict.join(expected_snps.out.expected_vcf)))

    lofreq_indel(recalibrate_bq.out.recalibrated_alignment.join(ch_ref_dict.combine(ch_bed)))

    lofreq_call(lofreq_indel.out.indel_alignment.join(ch_ref_dict.combine(ch_bed)))

    filter_lofreq(lofreq_call.out.minor_alleles)

    //call_variants(ch_primer_trimmed_alignment.join(ch_ref))

    //make_consensus(ch_primer_trimmed_alignment)

    //align_consensus_to_ref(make_consensus.out.consensus.join(ch_indexed_ref))

    // Collect multi-sample outputs
    if (params.collect_outputs | params.apply_qc) {
	aggregate_fastp = fastp.out.fastp_csv.map{ it -> it[1] }.collectFile(
	    keepHeader: true,
	    sort: { it.text },
	    name: "${params.collected_outputs_prefix}_fastp.csv",
	    storeDir: "${params.outdir}"
	)

	aggregate_bam = qualimap_bamqc.out.alignment_qc.map{ it -> it[1] }.collectFile(
	    keepHeader: true,
	    sort: { it.text },
	    name: "${params.collected_outputs_prefix}_qualimap_alignment_qc.csv",
	    storeDir: "${params.outdir}"
	)

	aggregate_samtools = samtools_stats.out.stats_summary_csv.map{ it -> it[1] }.collectFile(
	    keepHeader: true,
	    sort: { it.text },
	    name: "${params.collected_outputs_prefix}_samtools_stats_summary.csv",
	    storeDir: "${params.outdir}"
	)

    aggregate_counts = detect_ribo_repeats.out.ribo_rpt_counts.map{ it -> it[1] }.collectFile(
	    keepHeader: true,
	    sort: { it.text },
	    name: "${params.collected_outputs_prefix}_rrna_counts.csv",
	    storeDir: "${params.outdir}"
	)

    aggregate_snps = filter_lofreq.out.formatted_vcf.map{ it -> it[1] }.collectFile(
	    keepHeader: true,
	    sort: { it.text },
	    name: "${params.collected_outputs_prefix}_lofreq.vcf",
	    storeDir: "${params.outdir}"
	)
    }

    if (params.apply_qc) {
        amplicon_bed_file = amplicon_coverage.out.amplicon_bed
        qc_filter(aggregate_fastp.combine(aggregate_bam.combine(aggregate_samtools.combine(aggregate_counts.combine(aggregate_snps.combine(amplicon_bed_file))))))
        report_results(qc_filter.out.qc_results)
    }
/**
    // Collect Provenance
    // The basic idea is to build up a channel with the following structure:
    // [sample_id, [provenance_file_1.yml, provenance_file_2.yml, provenance_file_3.yml...]]
    // At each step, we add another provenance file to the list using the << operator...
    // ...and then concatenate them all together in the 'collect_provenance' process.
    ch_provenance = ch_provenance.combine(ch_pipeline_provenance).map{ it ->             [it[0], [it[1]]] }
    ch_provenance = ch_provenance.join(hash_ref.out.provenance).map{ it ->               [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(hash_fastq.out.provenance).map{ it ->             [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(fastp.out.provenance).map{ it ->                  [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(bwa_mem.out.provenance).map{ it ->                [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(trim_primer_sequences.out.provenance).map{ it ->  [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(make_consensus.out.provenance).map{ it ->         [it[0], it[1] << it[2]] }
    ch_provenance = ch_provenance.join(align_consensus_to_ref.out.provenance).map{ it -> [it[0], it[1] << it[2]] }

    collect_provenance(ch_provenance)
  */
}
