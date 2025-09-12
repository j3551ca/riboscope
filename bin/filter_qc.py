#!/usr/bin/env python

import argparse
import pandas as pd


def normalize_headers(df):
    df.columns = (df.columns.str.strip()
                  .str.lower()
            .str.replace(" ", "_")
            .str.replace("-", "_"))
    return df

def ingest_seq_results(fastp_csv, bam_csv, samtools_stats_csv):
    """
    Read in aggregated QC reports and combine alignment & read QC 
    for downstream QC logic. Keep amplicon counts and VCF separate for 
    independent processing.
    """
    read_qc = pd.read_csv(fastp_csv, header=0)
    alignment_qc = pd.read_csv(bam_csv, header = 0)
    samtools_qc = pd.read_csv(samtools_stats_csv, header = 0)

    merge_df = read_qc \
    .merge(alignment_qc, on = "sample_id", how="outer") \
    .merge(samtools_qc, on="sample_id", how="outer")

    return(merge_df)

def filter_vcf(vcf_file):
    """
    Filter VCF based on QC logic pass/ fail samples 
    (exclude failed samples and/ or SNPs from failed region).
    """
    
    variant_qc =  normalize_headers(pd.read_csv(vcf_file, header = 0, sep = "\t"))
    


#%%
# import pandas as pd
def syphilis_amplicons(amplicon_counts_csv):
    """
    Calculate median and min reads per group of queried sequences
    that represent each Treponema pallidum amplicon. Output df for merging
    with read & alignment QC for downstream QC logic. 
    """
    df = pd.read_csv(amplicon_counts_csv, header = 0)
    amplicon_summary = pd.DataFrame()
    groups = ["rpt1_pool1_", "rpt1_pool2_", "rpt2_pool1_", "rpt2_pool2_"]

    for prefix in groups:
        cols = [c for c in df.columns if c.startswith(prefix)]

        if cols:
            new_col_med = prefix.rstrip("_-") + "_median" # in case prefix already ends in _ or -
            new_col_min = prefix.rstrip("_-") + "_min"
            amplicon_summary[new_col_med] = df[cols].median(axis=1)
            amplicon_summary[new_col_min] = df[cols].min(axis=1)
    amplicon_summary["sample_id"] = df["sample_id"]
    amplicon_summary[["pos_str", "neg_str"]] = df[["pos_str", "neg_str"]]

    return(amplicon_summary, df) #df for heatmap, #amlicon_summary for merging with reads 

# median_cols = [c for c in amplicon_summary.columns if c.endswith("_median")]

# for col in median_cols:
#     pass_col = col.replace("_median", "_pass")
#     amplicon_summary[pass_col] = (amplicon_summary[col] >= 300).astype(int)

#%%
import argparse

def define_qc_rules(reads_qc_df, amplicon_qc_df, min_q20, min_q30, min_bq, min_depth, min_map_pair, min_map,
                    min_pct_map, min_mq, min_10X, min_50X, max_secondary, min_pos_count, min_med_count, max_qc_flags):

    df = reads_qc_df.merge(amplicon_qc_df, on="sample_id", how="left")

    #define vars
    gc_low =  df["gc_content_after_filtering"].mean(axis=0) - df["gc_content_after_filtering"].std() 
    gc_high = df["gc_content_after_filtering"].mean(axis=0) + df["gc_content_after_filtering"].std() 
    read_length_low = (df["average_length"].mode().iloc[0])*0.75
    read_length_high =  (df["average_length"].mode().iloc[0])*1.25
    amplicon_qc_median = [ col for col in amplicon_qc_df.columns if col.endswith("_median")]

    #dictn of QC rules to apply to qc pd.df
    soft_fail_rules = {
        "abnormal_gc": (df["gc_content_after_filtering"] < gc_low) | (df["gc_content_after_filtering"] > gc_high),
        "insufficient_cycles": ((df["read1_mean_length_after_filtering"] <= read_length_low) | (df["read2_mean_length_after_filtering"] <= read_length_low)),
        "excessive_cycles": ((df["average_length"]!=0) & ((df["read1_mean_length_after_filtering"] >= read_length_high) | (df["read2_mean_length_after_filtering"] >= read_length_high))),
        "low_q20": df["q20_rate_after_filtering"] < min_q20,
        "low_q30": df["q30_rate_after_filtering"] < min_q30,
        "low_bq": df["average_quality"] < min_bq,
        "low_depth": df["mean_depth_coverage"] < min_depth,
        "low_mapped_paired": df["reads_mapped_and_paired"] < min_map_pair,
        "low_percent_mapped": df["percent_mapped_reads"] < min_pct_map,
        "low_mapq": df["mean_mapping_quality"] < min_mq,
        "low_10X_breadth": df["proportion_genome_covered_over_10x"] < min_10X,
        "low_50X_breadth": df["proportion_genome_covered_over_50x"] < min_50X,
        "secondary": df["num_secondary_alignments"] > max_secondary,
        "abnormal_amplicon_counts": (df["pos_str"] < min_pos_count) | (df["neg_str"] > 0),
    }

    #add instant fails here
    hard_fail_rules = {
        "no_amplicons_detected": (df[amplicon_qc_median] < min_med_count).all(axis=1),
        "low_reads": df["num_mapped_reads"] < min_map,
    }

    #add flag name as col and bool for samples
    for name, mask in soft_fail_rules.items():
        df[name] = mask

    for name, mask in hard_fail_rules.items():
        df[name] = mask

    #handle case of missing data - structured to be True if fails. Need to do after masking for >0 rules
    soft_cols = list(soft_fail_rules.keys())
    hard_cols = list(hard_fail_rules.keys())
    flag_cols = soft_cols + hard_cols
    df[flag_cols] = df[flag_cols].fillna(True)

    df["n_flags"] = df[flag_cols].sum(axis=1) #tally flags per sample - rowwise
    #collect key values/col names/flags as reason list to store in "flags" col per samp
    df["flags"] = df.apply(lambda rule: [col for col in flag_cols if rule[col]], axis=1)
    #final qc decision
    df["qc_fail"] = (df["n_flags"] > max_qc_flags) | df[hard_cols].any(axis=1)

    exclusion_list = df[df["qc_fail"]]["sample_id"].tolist()

    return(df, exclusion_list)


#%%
def main(args):
    read_aln_df = ingest_seq_results(args.fastp_qc, args.bam_qc, args.samtools_stats) # add samtools

    #qc_df.to_csv("test_qc.csv", index=False)
    amplicon_qc, amplicon_df = syphilis_amplicons(args.amplicon_counts)






if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Digest QC results and format for pass/fail logic.")
    parser.add_argument("--fastp_qc", type=str, required=True, help="Aggregated fastp results")
    parser.add_argument("--amplicon_counts", type=str, required=True, help="Aggregated counts of query sequences")
    parser.add_argument("--bam_qc", type=str, required=True, help="Aggregated qualimap alignment QC")
    parser.add_argument("--samtools_stats", type=str, required=True, help="Aggregated SAMtools stats for alignment QC")
    parser.add_argument("--vcf_qc", type=str, required=True, help="Aggregated LoFreq SNPs")
    parser.add_argument("--min_med_amp_count", type=int, required=True, help="Minimum median query sequence counts to consider an amplicon as present")
    parser.add_argument("--output", type=str, required=False, help="Output file to write QC summary to")
    args = parser.parse_args()
    main(args)
