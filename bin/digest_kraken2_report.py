#!/usr/bin/env python

import pandas as pd
import argparse


def parse_kraken2(path):
    """
    Parse Kraken2 report into pd df
    """

    cols = ["percent_reads",
            "clade_reads",
            "taxon_reads",
            "taxonomy_level",
            "taxon_id",
            "taxon_name"]
    
    df = pd.read_csv(path, sep = "\t", names=cols, comment ="#")

    df["percent_reads"] = df["percent_reads"].astype(float)
    df["clade_reads"] = df["clade_reads"].astype(int)
    df["taxon_reads"] = df["taxon_reads"].astype(int)
    df["taxonomy_level"] = df["taxonomy_level"].astype(str)
    df["taxon_id"] = df["taxon_id"].astype(int)
    df["taxon_name"] = df["taxon_name"].astype(str).str.strip().str.lower()

    return(df)

def safe_return(df, name, col, dtype):
    return df[df["taxon_name"] == name.lower()][col].iloc[0] if name.lower() in df["taxon_name"].values else dtype(0)

def detect_pathogen_percentage(df, pathogen, host):
    """
    Detect percentage of reads from host, pathogen, and other.
    """
    
    pathogen_percent = safe_return(df, pathogen, "percent_reads", float)
    host_percent = safe_return(df, host, "percent_reads", float)
    unclassified_percent = safe_return(df, "unclassified", "percent_reads", float)
    other_percent = 100.0 - (pathogen_percent + host_percent + unclassified_percent)

    return pathogen_percent, host_percent, round(other_percent, 2), unclassified_percent


def detect_pathogen_reads(df, pathogen, host):
    """
    Detect number of reads from host, pathogen, and other.
    """
    pathogen_reads = safe_return(df, pathogen, "clade_reads", int)
    host_reads = safe_return(df, host, "clade_reads", int)
    unclassified_reads= safe_return(df, "unclassified", "clade_reads", int)
    classified_reads = safe_return(df, "root", "clade_reads", int)
    total_reads = classified_reads + unclassified_reads
    other_reads = total_reads - (pathogen_reads + host_reads + unclassified_reads)

    return pathogen_reads, host_reads, other_reads, unclassified_reads


def main(args):
    kraken_df = parse_kraken2(args.report)
    pathogen_perc, host_perc, other_perc, unclassified_perc = detect_pathogen_percentage(kraken_df, args.pathogen, args.host)
    pathogen_reads, host_reads, other_reads, unclassified_reads = detect_pathogen_reads(kraken_df, args.pathogen, args.host)
    summary_df = pd.DataFrame([{
        "sample_id": args.sample_id,
        "analysis_stage": args.dehost_stage,
        "pathogen": args.pathogen.lower(),
        "host": args.host.lower(),
        "pathogen_reads": pathogen_reads,
        "host_reads": host_reads,
        "unclassified_reads": unclassified_reads,
        "other_reads": other_reads,
        "pathogen_perc": pathogen_perc,
        "host_perc": host_perc,
        "unclassified_perc": unclassified_perc, 
        "other_perc": other_perc,
    }])

    summary_df.to_csv(f"{args.kraken2_output}.tsv", sep="\t", index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate reads for pathogen and host")
    parser.add_argument("--host", type=str, help="Taxon name of host")
    parser.add_argument("--pathogen", type=str, help="Taxon name of pathogen")
    parser.add_argument("--report", type=str, help="Path to Kraken2 report")
    parser.add_argument("--dehost_stage", type=str, help="Pres or post dehosting")
    parser.add_argument("--sample_id", type=str, help="Name of sample being summarized")
    parser.add_argument("--kraken2_output", type=str, help="Name of sample being summarized")

    args = parser.parse_args()

    main(args)




