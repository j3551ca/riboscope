#!/usr/bin/env python

import argparse
import pandas as pd
import os
from jinja2 import Environment, FileSystemLoader
import plotly.express as px
from datetime import datetime

def ingest_qc_results(qc_summary_file, amp_df, annotated_vcf, failed_amplicons):
    """
    Read in QC results as pandas.DataFrames for use in html report generation.
    """
    qc_df = pd.read_csv(qc_summary_file, header=0)
    reportable_vcf = pd.read_csv(annotated_vcf, header=0, sep="\t")
    amplicon_df = pd.read_csv(amp_df, header=0)
    failed_amplicon_df = pd.read_csv(failed_amplicons, header=0)

    #prep vcf table for viz & summarizing 
    reportable_vcf["snv"] = reportable_vcf["ref"].astype(str) + reportable_vcf["pos"].astype(str) + reportable_vcf["alt"]
    reportable_vcf = reportable_vcf.sort_values("pos", ascending=True)
    snv_counts = reportable_vcf.groupby("snv")["sample_id"].nunique().rename("count")
    reportable_vcf = reportable_vcf.merge(snv_counts, on="snv", how="left")
    #total samples with successfully sequenced pool that contain each SNV
    total_samples = reportable_vcf.groupby(["pool"])["sample_id"].nunique().rename("total_pool_count")
    reportable_vcf = reportable_vcf.merge(total_samples, on=["pool"], how = "left")
    reportable_vcf["sample_count"] = reportable_vcf["count"].astype(str) + "/" + reportable_vcf["total_pool_count"].astype(str)
    reportable_vcf["vaf"] = (reportable_vcf["af"]*100).round(1).astype(str) + "%"

    return(qc_df, amplicon_df, reportable_vcf, failed_amplicon_df)


def summarize_results(qc_summary, vcf_df):
    """
    Summarize sample pass/ fail and SNP results - results overview.
    """

    vcf_df["pool"]=vcf_df["pool"].replace({"pool1":"16S","pool2":"23S"})

    results_df = pd.DataFrame({"metric":["Total Number of Samples", 
                            "Number of Failed Samples", 
                            "SNPs ≥90% Frequency", 
                            "SNPs ≥50% Samples"],
                 "results": [len(qc_summary["sample_id"]), 
                             len(qc_summary[qc_summary["qc_fail"]]), 
                             ",<br>".join(vcf_df[vcf_df["af"]>=0.9][["snv", "pool"]]
                                        .drop_duplicates().astype(str).agg("_".join, axis=1)),
                             ",<br>".join(vcf_df[(vcf_df["count"]/vcf_df["total_pool_count"])>=0.5][["snv", "pool"]]
                                        .drop_duplicates().astype(str).agg("_".join, axis=1))]})
    
    fail_flags = qc_summary[qc_summary["flags"]!="[]"][["sample_id", "flags"]]
    
    results_df=results_df.to_html(index=False, escape=False, classes="pretty-table")
    fail_flags=fail_flags.to_html(index=False, escape=False, classes="pretty-table")

    return results_df, fail_flags


def plot_amplicons(failed_amps):
    """
    Heatmap of (boolean) amplicon success for each sample.
    """
    failed_amps_wide = failed_amps.pivot(index="sample_id", columns="amplicon_name", values="amplicon_fail")

    failed_order = failed_amps_wide.eq(1).sum(axis=1).sort_values(ascending=False)
    failed_amps_wide = failed_amps_wide.loc[failed_order.index]
    display = failed_amps_wide.replace({0: "Pass", 1: "Fail"})

    fig = px.imshow(
        failed_amps_wide,
        aspect="auto",
        text_auto=".0f",
        color_continuous_scale=["#1cba56", "#c93808"],
        labels=dict(x="Sequence", y="Sample")
    )

    fig.update_traces(showscale=False,
                        text=display.to_numpy(),
                        texttemplate="%{text}",
                        hoverongaps=False,
                        xgap=1, ygap=1)

    fig.update_layout(title="Sequencing Success of Amplicons",
                        coloraxis_showscale=False,
                        autosize=True,
                        height=300)
    return fig.to_html(full_html=False, include_plotlyjs="cdn")


def plot_vcf(vcf_df):
    """
    Plot SNPs by pool and represent repeats present in pool with shape
    """
    fig = px.scatter(vcf_df, 
            x = "snv", 
            y = "af",
            color="repeat_status",  # optional categorical coloring
            opacity=0.6,
            size="count",
            size_max=8,
            hover_data={"snv":True, "af":False,"vaf":True, "sample_id":True, "count":False,"sample_count":True, "pool":False},  # info shown on hover
            facet_col="pool",
            log_y=True,
            title="rRNA SNVs Across Samples",
            template="plotly_white")

    fig.update_xaxes(matches=None)
    fig.update_xaxes(title_text="", row=1, col=1)
    fig.update_xaxes(title_text="", row=1, col=2)
    fig.update_layout(
        xaxis_title="",
        title_x=0.5,
        height=600,
        yaxis_title="Allele Frequency",
        legend_title_text="Ribosomal Repeat")
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))

    return fig.to_html(full_html=False, include_plotlyjs="cdn")


def plot_heatmap(raw_counts):
    """
    Heatmap of query sequence hits - proxy for presence/ absence of amplicons sequenced.
    """

    raw_counts = raw_counts.set_index("sample_id")
    raw_counts = raw_counts.sort_values("pos_str", ascending=False)
    fig = px.imshow(
        raw_counts,
        aspect="auto",
        text_auto=".0f",
        color_continuous_scale="viridis",
        labels=dict(x="Sequence", y="Sample", color="Count")
    )
    fig.update_layout(title="Presence of Query Sequences in rRNA Repeat Regions")
    return fig.to_html(full_html=False, include_plotlyjs="cdn")


def generate_html(amplicon_fig, snp_fig, count_heatmap, summary_data, qc_flags, html_template):
    """
    Take outputs of functions as input for html template file and render.
    """
    env = Environment(loader=FileSystemLoader(str(os.path.dirname(html_template))))
    template = env.get_template(os.path.basename(html_template))

    html_out = template.render(
        title="Syphilis Results Report",
        heading="SNVs Detected in Treponema pallidum Ribosomal Repeats",
        summary_table=summary_data,
        flags_table = qc_flags,
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        amplicon_plot=amplicon_fig,
        snp_plot=snp_fig,
        heatmap_plot=count_heatmap)

    with open("results_report.html", "w") as f:
        f.write(html_out)


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
