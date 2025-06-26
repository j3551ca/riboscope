# riboscope

This Nextflow pipeline monitors for mutations that may confer AMR in bacterial ribosomal repeats, specifically the 16S and 23S rRNA genes of *Treponema pallidum* subsp. *pallidum*, the etiological agent of syphilis. This workflow calls mutations in reference-guided assemblies of short paired-end reads sequenced from amplicon-based libraries. 

## Quick-Start

```
nextflow run BCCDC-PHL/riboscope \
  -profile conda \
  --cache ~/.conda/envs \
  --fastq_input /path/to/fastq_files \
  --ref /path/to/ref.fa \
  --bed /path/to/primer_scheme.bed \
  --search_seqs /path/to/query/seq/multi.fasta \
  --outdir /path/to/output_dir
```
## Table of Contents
  [Overview](#riboscope)<br>
  [Quick-Start](#quick-start)<br>
  [Workflow](#workflow)<br>
  [Usage](#usage)<br>
  [Input](#input)<br>
  [Output](#output)<br>
  [Parameters](#parameters)<br>
  [References](#references)<br>

## Workflow

```mermaid
flowchart TD
  ref[ref.fa]
  bed[scheme.bed]
  fastq[fastq]
  fastq --> fastp(fastp)
  ref --> bwa(bwa_mem)
  fastp -- trimmed_reads --> bwa
  bwa -- alignment --> trim_primer_seqs(trim_primer_sequences)
  bed --> trim_primer_seqs
  trim_primer_seqs -- primertrimmed_alignment --> make_consensus(make_consensus)
  make_consensus -- consensus --> align_consensus_to_ref(align_consensus_to_ref)
  ref --> align_consensus_to_ref
```

## Usage

The following command can be used to run the pipeline:

```
nextflow run BCCDC-PHL/riboscope \
  -profile conda \
  --cache ~/.conda/envs \
  --fastq_input /path/to/fastq_files \
  --ref /path/to/ref.fa \
  --bed /path/to/primer_scheme.bed \
  --search_seqs /path/to/query/seq/multi.fasta \
  --outdir /path/to/output_dir
```

By default, reads will be trimmed by fastp prior to alignment. To align the untrimmed reads instead, use the `--align_untrimmed_reads` parameter:

```
nextflow run BCCDC-PHL/riboscope \
  -profile conda \
  --cache ~/.conda/envs \
  --fastq_input /path/to/fastq_files \
  --ref /path/to/ref.fa \
  --bed /path/to/primer_scheme.bed \
  --search_seqs /path/to/query/seq/multi.fasta \
  --align_untrimmed_reads \
  --outdir /path/to/output_dir
```

If this option is used, the pipeline will proceed as follows, using the original untrimmed reads as input for the bwa alignment process:

```mermaid
flowchart TD
  ref[ref.fa]
  bed[scheme.bed]
  fastq[fastq]
  fastq --> fastp(fastp)
  ref --> bwa(bwa_mem)
  fastq -- untrimmed_reads --> bwa
  bwa -- alignment --> trim_primer_seqs(trim_primer_sequences)
  bed --> trim_primer_seqs
  trim_primer_seqs -- primertrimmed_alignment --> make_consensus(make_consensus)
  make_consensus -- consensus --> align_consensus_to_ref(align_consensus_to_ref)
  ref --> align_consensus_to_ref
```



## Input
| Input  | Parameter   |  Description   |  Notes  |
|:----|:-----|:-----|:-----|
| Paired-end sequencing reads |  `fastq_input`  | Absolute path to directory containing raw FASTQ reads to be analyzed. Riboscope accepts gzip compressed or uncompressed files (*.fastq.gz, *.fq.gz, *.fastq, *.fq).    |  none    |
|Reference genome | `ref` |  Reference genome used to align reads to during guided assembly    |  none    |
|BED file	 | `bed`  | Primer scheme BED file in the [6 column format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)    |  none    |
|MultiFASTA of query sequences | `search_seqs`  | MultiFASTA file contain short sequences to query raw reads. Used to verify presence of expected gene copies.     |  none    |

## Output


## Parameters

| Parameter  | Description   |  Required   |  Default  |
|:----|:-----|:-----|:-----|
|`fastq_input` | Absolute path to directory containing raw FASTQ reads to be analyzed. Riboscope accepts gzip compressed or uncompressed files (*.fastq.gz, *.fq.gz, *.fastq, *.fq).   | yes    |  none    |
|`outdir` | Absolute path to directory to write results to. | no    |  ./results    |
|`ref` | Reference genome used to align reads to during guided assembly  | yes    |  none    |
|`bed` | Primer scheme BED file in the [6 column format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)   | yes    |  none    |
|`search_seqs` | MultiFASTA file contain short sequences to query raw reads. Used to verify presence of expected gene copies.   | yes    |  none    |
|`cache` | Directory to cache conda environments for future use.   | no    |  ./work/conda    |
|`align_untrimmed_reads` | Skips read trimming by fastp. Allows alignment of raw untrimmed reads using bwa. | no    |  off    |
|`min_depth` | Minimum number of reads covering a genomic position. | no    |  10    |
|`collect_outputs` | Summarize outputs of multiple samples into one. | no    |  off    |
|`collected_outputs_prefix` | Prefix to name multi-sample summary files with. | no    |  'collected'    |

## References

1. Jago, M.J., Soley, J.K., Denisov, S. et al. High-throughput method characterizes hundreds of previously unknown antibiotic resistance mutations. Nat Commun 16, 780 (2025). https://doi.org/10.1038/s41467-025-56050-2
1. https://github.com/BCCDC-PHL/amplicon-consensus
