#!/usr/bin/env python
#%%
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
#%%

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

    return(amplicon_summary, df)  #amplicon_summary for merging with reads, df for heatmap

#%%
import argparse

def sample_qc_filter(reads_qc_df, amplicon_qc_df, min_q20, min_q30, min_bq, min_depth, min_map_pair, min_map,
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
    failed_amplicons = amplicon_qc_df
        
    for col in amplicon_qc_median:
        fail_col = col.replace("_median", "_fail")
        failed_amplicons[fail_col] = (failed_amplicons[col] < min_med_count).astype(int)
    failed_amplicons = failed_amplicons.filter(regex="sample_id$|fail$") 
    failed_amplicons = failed_amplicons.melt(id_vars="sample_id", 
                          value_vars=[c for c in failed_amplicons.columns if c.endswith("_fail")], 
                          var_name="amplicon_name",
                          value_name = "amplicon_fail")
    failed_amplicons["amplicon_name"] = failed_amplicons["amplicon_name"].str.replace("_fail$", "", regex=True)
    failed_amplicons["pool"] = failed_amplicons["amplicon_name"].str.extract(r"(pool\d+)")
    failed_amplicons["repeat"] = failed_amplicons["amplicon_name"].str.extract(r"(rpt\d+)")

    return(df, exclusion_list, failed_amplicons)
#%%
def summarize_repeats(failed_amplicons: pd.DataFrame) -> str:
    """
    Count unique repeats represented in each pool in each row (sample & pool group)
    """
    passing = failed_amplicons.loc[failed_amplicons["amplicon_fail"] == 0, "repeat"].unique()
    
    if len(passing) == 0:
        return "none"
    elif len(passing) == len(failed_amplicons["repeat"].unique()):
        return "all"  
    else:
        return ",".join(sorted(passing))  #"rpt1", "rpt2", "rpt..."

#%%
def filter_vcf(vcf_file, failed_samples, failed_amplicons, amplicon_bed):
    """
    Filter VCF based on QC logic pass/ fail samples 
    (exclude failed samples and/ or SNPs from failed regions).
    Annotates VCF with which repeats (copies) passed in each pool (gene). 
    The way primers are named in bed should contain information about gene and copy number
    to match groups var used in amp_qc from syphilis_amplicons function.
    """
    variant_df =  normalize_headers(vcf_file)
    pass_variants = variant_df[~variant_df["sample_id"].isin(failed_samples)] 
    amplicon_bed["amp"] = amplicon_bed["name"].str.extract(r'(rpt\d+_pool\d+)')
    amplicon_bed["pool"] = amplicon_bed["name"].str.extract(r'(pool\d+)')
    amplicon_bed.columns = amplicon_bed.columns.str.lower()

    #sample_id, pool, repeats present
    repeat_labels = (
    failed_amplicons.groupby(["sample_id","pool"], group_keys=False) #exclude grouping vars
        .apply(summarize_repeats, include_groups=False)
        .reset_index(name="repeat_status")
    )
    # add pool/ gene depending on regions to vcf df
    annotated_vcf = (
    pass_variants.merge(amplicon_bed, on="chrom", how="left")
                .query("start <= pos <= end")
                .drop_duplicates(subset=["chrom","pos", "sample_id","alt","af"]) 
                .reset_index(drop=True)
                .drop(["start", "end", "name"], axis = 1)
    ) # drops failed samples becasue they have 0 in chrom of vcf file ^ 

    # add which repeats are represented in each pool based on failed amplicon minimum median query seq count
    # this vcf df is now ready for plotting
    annotated_vcf_pool_status = annotated_vcf.merge(repeat_labels, 
                                            on = ["sample_id", "pool"], 
                                            how = "left" )
    
    return(annotated_vcf_pool_status)

def plot_vcf(repeat_annotated_vcf):

    """
    Plot SNPs by pool and represent repeats present in pool with shape
    """
    repeat_annotated_vcf
    max_count_pool1 = df_pool1["COUNT"].max()
    max_count_pool2 = df_pool2["COUNT"].max()
    filter_pool1 = df_pool1[(df_pool1["COUNT"]>=5)| (df_pool1["AF"]>=0.01)]
    filter_pool2 = df_pool2[(df_pool2["COUNT"]>=5) | (df_pool2["AF"]>=0.01)]

    filter_pool1 = filter_pool1.sort_values("POS", ascending=True)
    filter_pool2 = filter_pool2.sort_values("POS", ascending=True)

    plt.clf()
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2,
                                figsize=(12, 8), sharey=True)

    ax1.grid(axis='x', linestyle="--", alpha=0.3)
    ax1.grid(axis='y', linestyle="--", alpha=0.3)
    sc1 = ax1.scatter(filter_pool1["SNV"], filter_pool1["AF"],
                s=filter_pool1["COUNT"],
                alpha=0.5)
    ax1.set_xticklabels(ax1.get_xticklabels(), va = "top",
                        ha="right")
    ax1.tick_params(axis="x", labelsize=6, 
                    labelrotation=45)
    ax1.set_xlabel("Pool 1")
    ax1.set_yscale("log")

    # now automatically pick 4 sizes from min→max and label them
    handles1, labels1 = sc1.legend_elements(prop="sizes", num=5, fmt="{x:.0f}", color="#1f77b4", alpha =0.7)
    ax1.legend(handles1, labels1, title="Samples", loc="upper left")


    ax2.grid(axis='x', linestyle="--", alpha=0.3)
    ax2.grid(axis='y', linestyle="--", alpha=0.3)
    sc2 = ax2.scatter(filter_pool2["SNV"], filter_pool2["AF"],
                s=filter_pool2["COUNT"],
                alpha=0.5)
    ax2.set_yticks(np.arange(0, 1.0, 0.1))
    ax2.set_xticklabels(ax2.get_xticklabels(),va="top",
                        ha="right")
    ax2.tick_params(axis="x", labelsize=6, 
                    labelrotation=45)
    ax2.set_xlabel("\n\n\n\n\n\nPool 2")
    ax2.set_yscale("log")

    handles2, labels2 = sc2.legend_elements(prop="sizes", num=5, fmt="{x:.0f}", color="#1f77b4", alpha=0.7)
    ax2.legend(handles2, labels2, title="Samples", loc="upper left")

    fig.text(0.05, 0.5, "Allele Frequency", ha="center", fontsize=12, rotation=90)
    fig.text(0.5, 0.04, "SNV", ha="center", fontsize=12)
    fig.text(0.4, 0.95, "SNVs ≥1% Frequency or ≥5 Samples Called")
    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
    plt.show()

    return plt

def plot_amplicons(raw_counts):

    """
    Heatmap of query sequence hits - proxy for presence/ absence of amplicons sequenced.
    """

    raw_counts.set_index("sample_id", inplace=True)

    plt.clf()
    plt.figure(figsize=(10, 6))
    sns.heatmap(raw_counts, annot=True, cmap="viridis", fmt=".0f",
                cbar_kws={"label": "Sequence Count"},
                annot_kws={"size": 6})  # Or use "coolwarm", "magma", etc.
    plt.xticks(rotation=45,ha='right') 
    plt.title("Presence of rRNA Repeat Sequences")
    plt.xlabel("Sequence")
    plt.ylabel("Samples")
    plt.tight_layout()
    
    return plt


#%%
def main(args):
    read_aln_df = ingest_seq_results(args.fastp_qc, args.bam_qc, args.samtools_stats)
    amplicon_qc, amplicon_df = syphilis_amplicons(args.amplicon_counts)
    qc_df, failed_samples, failed_amplicons = sample_qc_filter(read_aln_df, amplicon_qc)
    qc_df.to_csv(args.output, index=False)
    reportable_vcf = filter_vcf(args.vcf_file, failed_samples, failed_amplicons, args.amplicon_bed)
    plot_vcf(reportable_vcf)
    plot_amplicons(amplicon_df)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Digest QC results and format for pass/fail logic.")
    parser.add_argument("--fastp_qc", type=str, required=True, help="Aggregated fastp results")
    parser.add_argument("--amplicon_counts", type=str, required=True, help="Aggregated counts of query sequences")
    parser.add_argument("--bam_qc", type=str, required=True, help="Aggregated qualimap alignment QC")
    parser.add_argument("--samtools_stats", type=str, required=True, help="Aggregated SAMtools stats for alignment QC")
    parser.add_argument("--vcf_file", type=str, required=True, help="Aggregated LoFreq SNPs")
    parser.add_argument("--amplicon_bed", type = str, required=True, help="Bed file of all amplicons produced by amplicon_coverage process")
    parser.add_argument("--min_med_amp_count", type=int, required=True, help="Minimum median query sequence counts to consider an amplicon as present")
    parser.add_argument("--output", type=str, required=False, help="Output file to write QC summary to")
    args = parser.parse_args()
    main(args)
