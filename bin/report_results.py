#!/usr/bin/env python

import argparse
import pandas as pd
import io
import os
from jinja2 import Environment, FileSystemLoader
import plotly.express as px
from dna_features_viewer import GraphicFeature, GraphicRecord
from datetime import datetime

def ingest_qc_results(qc_summary_file, amp_df, annotated_vcf, failed_amplicons, kraken2_tsv):
    """
    Read in QC results as pandas.DataFrames for use in html report generation.
    """
    qc_df = pd.read_csv(qc_summary_file, header=0)
    reportable_vcf = pd.read_csv(annotated_vcf, header=0, sep="\t")
    amplicon_df = pd.read_csv(amp_df, header=0)
    failed_amplicon_df = pd.read_csv(failed_amplicons, header=0)
    kraken_df = pd.read_csv(kraken2_tsv, header = 0, sep ="\t")

    #prep vcf table for viz & summarizing 
    reportable_vcf["snv"] = reportable_vcf["ref"].astype(str) + reportable_vcf["pos"].astype(str) + reportable_vcf["alt"]
    reportable_vcf["feature_snv"] = reportable_vcf["ref"].astype(str) + reportable_vcf["feature_pos"].astype(str) + reportable_vcf["alt"]
    reportable_vcf = reportable_vcf.sort_values("pos", ascending=True)
    snv_counts = reportable_vcf.groupby("snv")["sample_id"].nunique().rename("count")
    reportable_vcf = reportable_vcf.merge(snv_counts, on="snv", how="left")
    #total samples with successfully sequenced pool that contain each SNV
    total_samples = reportable_vcf.groupby(["pool"])["sample_id"].nunique().rename("total_pool_count")
    reportable_vcf = reportable_vcf.merge(total_samples, on=["pool"], how = "left")
    reportable_vcf["sample_count"] = reportable_vcf["count"].astype(str) + "/" + reportable_vcf["total_pool_count"].astype(str)
    reportable_vcf["vaf"] = (reportable_vcf["af"]*100).round(1).astype(str) + "%"

    return(qc_df, amplicon_df, reportable_vcf, failed_amplicon_df, kraken_df)


def summarize_results(qc_summary, vcf_df, ref, reporting_vaf, reporting_sample_perc):
    """
    Summarize sample pass/ fail and SNP results - results overview.
    """

    #vcf_df["pool"]=vcf_df["pool"].replace({"pool1":"16S","pool2":"23S"})
    reporting_vaf=float(reporting_vaf)
    reporting_sample_perc = float(reporting_sample_perc)

    results_df = pd.DataFrame({"metric":["Total Number of Samples", 
                            "Number of Failed Samples", 
                            f"SNPs ≥{reporting_vaf*100:.1f}% Frequency", 
                            f"SNPs ≥{reporting_sample_perc*100:.1f}% Samples",
                            "Reference Genome"],
                 "results": [len(qc_summary["sample_id"]), 
                             len(qc_summary[qc_summary["qc_fail"]]), 
                             ",<br>".join(vcf_df[vcf_df["af"]>=reporting_vaf][["feature_snv", "feature_name"]]
                                        .drop_duplicates().astype(str).agg("_".join, axis=1)),
                             ",<br>".join(vcf_df[(vcf_df["count"]/vcf_df["total_pool_count"])>=reporting_sample_perc][["feature_snv", "feature_name"]]
                                        .drop_duplicates().astype(str).agg("_".join, axis=1)),
                                        os.path.abspath(ref)]})
    
    fail_flags = qc_summary[qc_summary["flags"]!="[]"][["sample_id", "flags"]]
    
    results_df=results_df.to_html(index=False, escape=False, classes="pretty-table")
    fail_flags=fail_flags.to_html(index=False, escape=False, classes="pretty-table")

    return results_df, fail_flags


def plot_amplicons(failed_amps):
    """
    Heatmap of (boolean) amplicon success or failure for each sample.
    """
    failed_amps_wide = failed_amps.pivot(index="sample_id", columns="amplicon_name", values="amplicon_fail")

    failed_order = failed_amps_wide.eq(1).sum(axis=1).sort_values(ascending=False)
    failed_amps_wide = failed_amps_wide.loc[failed_order.index]
    display = failed_amps_wide.replace({0: "Pass", 1: "Fail"})

    n_samples =len(failed_amps["sample_id"].unique())
    n_amplicons = (len(failed_amps.columns.unique())-1)

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
                        height=450 if n_samples <=20 else n_samples*30,
                        width=n_amplicons*200)
    
    return fig.to_html(full_html=False, include_plotlyjs="cdn")


def plot_features(complete_gff):
    """
    Plot genomic features from GFF file.
    """
    try:
        gff = pd.read_csv(complete_gff, sep="\t", header=0)
    except Exception as e:
        raise ValueError(f"Complete GFF failed to be loaded as TSV: {complete_gff}") from e

    features = []
    features = [
        GraphicFeature(start=s, end=e, 
                    strand=1 if st == "+" else -1, 
                    color="#ffcccc" if so =="python" else "#ccccff", label=l)
    for so, s, e, st, l in zip(gff.source, gff.start, gff.end, gff.strand, gff.feature)
    ]

    record = GraphicRecord(sequence_length=gff.iloc[-1]["end"], features=features)
    ax, _ = record.plot(figure_width=16)

    buf = io.BytesIO() 
    ax.figure.savefig(buf, format="svg", bbox_inches="tight")
    buf.seek(0)
    svg_features = buf.getvalue().decode("utf-8")
    svg_features = svg_features.replace("<svg ", '<svg style="width:90%; height:auto;" ') 

    return svg_features

#%%

def plot_vcf(vcf_df):
    """
    Plot SNPs by pool and represent repeats present in pool with color.
    If GFF file provided, plot SNPs by feature.
    """
    n_feature = vcf_df["feature_name"].nunique()
    n_col = 2
    n_row = (n_feature/n_col)

    fig = px.scatter(vcf_df, 
            x = "feature_snv", 
            y = "af",
            color="repeat_status",
            opacity=0.6,
            size="count",
            size_max=12,
            hover_data={"feature_snv":True, "snv":True, "af":False,
                        "vaf":True, "sample_id":True, "count":False,
                        "sample_count":True, "pool":False},  # info shown on hover
            facet_col="feature_name",
            facet_col_wrap=n_col,
            facet_row_spacing = 0.2,
            log_y=True,
            title="rRNA SNVs Across Samples",
            template="plotly_white")
    
    fig.update_traces(marker=dict(sizemin=4))
    fig.update_xaxes(matches=None, showticklabels=True, title_text="")
    fig.update_yaxes(title_text="")
    fig.update_layout(
        xaxis_title="",
        title_x=0.5,
        height=250*n_row,
        legend_title_text="Ribosomal Repeat",
        )
    
    fig.add_annotation(
    text="Allele Frequency",
    xref="paper",
    yref="paper",
    x=0, # move left of plots
    y=0.5,
    showarrow=False,
    textangle=-90,
    font=dict(size=14),
)
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))

    return fig.to_html(full_html=False, include_plotlyjs="cdn")


def plot_heatmap(raw_counts):
    """
    Heatmap of query sequence hits - proxy for presence/ absence of amplicons sequenced.
    """
    n_samples = len(raw_counts["sample_id"].unique())
    raw_counts = raw_counts.set_index("sample_id")
    raw_counts = raw_counts.sort_values("pos_str", ascending=False)
    
    fig = px.imshow(
        raw_counts,
        aspect="auto",
        text_auto=".0f",
        color_continuous_scale="viridis",
        labels=dict(x="Sequence", y="Sample", color="Count")
    )
    fig.update_layout(title="Presence of Query Sequences in rRNA Repeat Regions",
                      height=600 if n_samples <=20 else n_samples*13)
    return fig.to_html(full_html=False, include_plotlyjs="cdn")


def plot_taxonomy(kraken_results):
    """
    Stacked bar chart of taxonomic classification of reads to host, pathogen, unclassified, & other.
    """
    read_values = [ c for c in kraken_results.columns if "reads" in c]
    perc_values = [c for c in kraken_results.columns if "perc" in c]

    read_long = pd.melt(kraken_results, 
                        id_vars = ["sample_id", "analysis_stage"], 
                        value_vars=read_values, 
                        var_name = "origin", value_name="reads")
    perc_long = pd.melt(kraken_results, 
                        id_vars = ["sample_id","analysis_stage"], 
                        value_vars=perc_values, 
                        var_name="origin", value_name="percentage")

    read_long["origin"] = read_long["origin"].str.replace("_reads", "")
    perc_long["origin"] = perc_long["origin"].str.replace("_perc", "")

    df = read_long.merge(perc_long, 
                            on = ["sample_id", "analysis_stage", "origin"])

    # ensure pre comes before post-dehosting
    df = df.sort_values("analysis_stage", ascending=False)

    #y-axis order
    sort_order_pre = df[(df["origin"]=="pathogen") & (df["analysis_stage"]=="pre_dehosting")]\
        .sort_values("reads", ascending=False)["sample_id"].tolist()

    stack_order = ["pathogen", "host", "unclassified", "other"]

    n_samples = len(df["sample_id"].unique())

    df["total_reads"] = df.groupby(["sample_id", "analysis_stage"])["reads"].transform("sum").astype(str)

    df["analysis_stage"] = df["analysis_stage"].replace({"pre_dehosting": "Before Dehosting",
                                                        "post_dehosting": "After Dehosting"})

    fig = px.bar(df, 
            x = "reads", 
            y = "sample_id", 
            facet_col="analysis_stage", 
            color="origin", 
            text="percentage", 
            hover_data={"origin":True, "percentage": True, 
                        "reads": ":.0f","total_reads": True, 
                        "analysis_stage": True, "sample_id": True},
            orientation="h",
            category_orders={"sample_id": sort_order_pre,
                            "origin": stack_order[::-1]},
            color_discrete_map={
            "pathogen": "#1b469d",
            "host": "#F54927",
            "unclassified": "#1b9d66",
            "other": "#7f7f7f"
        },
            template="plotly_white", 
            opacity = 0.9)

    fig.update_xaxes(matches="x")
    fig.update_xaxes(title_text="", row=1, col=1)
    fig.update_xaxes(title_text="", row=1, col=2)
    fig.update_layout(
        title="Taxonomic Classification of Reads",
        title_x=0.5,
        height= 600 if n_samples <=20 else n_samples*13,
        yaxis_title="Sample Name",
        annotations=[
            dict(
                text="Number of Reads",
                x=0.5,
                y=-0.1,
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(size=14)
            )
        ],
        legend=dict(traceorder="reversed"),
        margin=dict(b=50),
        )
    fig.update_traces(
        texttemplate="%{text:.2f}%",   
        textposition="inside",      
    )
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1], font=dict(size=14)))

    return fig.to_html(full_html=False, include_plotlyjs="cdn")

def generate_html(amplicon_fig, feat_fig, snp_fig, count_heatmap, summary_data, qc_flags, taxon_fig, html_template, n_flags, min_med_amp_count):
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
        n_flags=n_flags,
        min_med_amp_count=min_med_amp_count,
        date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        amplicon_plot=amplicon_fig,
        feature_plot = feat_fig,
        snp_plot=snp_fig,
        taxon_plot = taxon_fig,
        heatmap_plot=count_heatmap)

    with open("results_report.html", "w") as f:
        f.write(html_out)


def main(args):

    qc_df, amplicon_df, reportable_vcf, failed_amps, kraken_df = ingest_qc_results(args.qc_summary, 
                                                                        args.amplicon_counts, 
                                                                        args.reportable_vcf, 
                                                                        args.failed_amplicons,
                                                                        args.kraken2_tsv)
    summary_df, qc_flags = summarize_results(qc_df, reportable_vcf, args.ref, 
                                             args.reporting_vaf, args.reporting_sample_perc)
    amp_fig = plot_amplicons(failed_amps)
    feature_fig = plot_features(args.gff)
    snv_fig = plot_vcf(reportable_vcf)
    count_fig = plot_heatmap(amplicon_df)
    taxon_fig = plot_taxonomy(kraken_df)
    generate_html(amp_fig, feature_fig, snv_fig, count_fig, summary_df, qc_flags, taxon_fig,
                  args.html_template, args.n_qc_flag, args.min_med_amp_count)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate results summary report")
    parser.add_argument("--qc_summary", type=str, required=True, help="QC results summary of all samples")
    parser.add_argument("--amplicon_counts", type=str, required=True, help="Amplicon query sequence count table")
    parser.add_argument("--reportable_vcf", type=str, required=True, help="VCF file annotated with repeat presence/absence")
    parser.add_argument("--failed_amplicons", type=str, required=True, help="Amplicon pass/fail status per sample")
    parser.add_argument("--kraken2_tsv", type=str, required=True, help="Kraken2 results summary tsv file")
    parser.add_argument("--html_template", type=str, required=True, help="html template file for results report")
    parser.add_argument("--ref", type=str, required=True, help="Reference genome used during alignment and variant calling")
    parser.add_argument("--gff", type=str, required=True, help="GFF containing user-annotated and program-annotated (unannotated) regions spanning reference sequence, " \
    "generated by adjust_variant_coordinates.py.")
    parser.add_argument("--n_qc_flag", type=str, required=True, help="Maximum number of soft fail QC flags allowable before sample is failed")
    parser.add_argument("--min_med_amp_count", type=str, required=True, help="Minimum median query sequence counts to consider an amplicon as present")
    parser.add_argument("--reporting_vaf", type=float, default=0.9, help="Minimum variant allele frequency in a given sample to report SNVs")
    parser.add_argument("--reporting_sample_perc", type=float, default=0.3, help="Minimum percentage of total successfully sequenced samples containing a given SNV")

    args = parser.parse_args()
    main(args)
