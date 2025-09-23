#!/usr/bin/env python

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

def ingest_qc_results(amp_df, annotated_vcf):

    reportable_vcf = pd.read_csv(annotated_vcf, header=0)
    amplicon_df = pd.read_csv(amp_df, header=0)

    return(amplicon_df, reportable_vcf)

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


def main(args):
    amplicon_df, reportable_vcf = ingest_qc_results(args.amplicon_heatmap_df, args.annotated_vcf)
    plot_vcf(reportable_vcf)
    plot_amplicons(amplicon_df)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--amplicon_heatmap_df", type=str, required=True, help="Amplicon query sequence count table")
    parser.add_argument("--annotated_vcf", type=str, required=True, help="VCF file annotated with repeat presence/absence")
  
    args = parser.parse_args()
    main(args)
