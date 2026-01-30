#!/usr/bin/env python
#%%
import argparse
import os
import pandas as pd
import re


#%%

def load_variants(vcf_path):
    if not os.path.isfile(vcf_path):
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    try:
        variants = pd.read_csv(vcf_path,
                            sep="\t",
                            comment="#",
                            header=0,
                            dtype={"CHROM": str, 
                                    "POS": int, 
                                    "REF": str,
                                    "ALT": str,
                                    "AF": float,
                                    "DP4": str,
                                    "SAMPLE_ID": str}
                                    )
    except Exception as e:
        raise ValueError(f"VCF failed to be loaded as TSV: {vcf_path}") from e
    
    return variants
#%%

def extract_feature_name(gff_row):
    """
    Extract the name of features from attribute col in GFF

    :param gff_row: string separated by ; specifying different attributes
    """
    attr = {}
    for item in gff_row.split(";"):
        if "=" in item:
            key, value = item.split("=",1)
            attr[key.strip().lower()] = re.sub(r"\s+","_", value.strip())

    # check allowable keys iteratively in order of importance
    for key in ("name", "product", "locus_tag", "id"):
        if key in attr and attr[key]:
            return attr[key]
    return "unspecified"


#%%
def load_gff(gff_str):
    if gff_str == "NO_FILE":
        features = gff_str
    else:
        if not os.path.isfile(gff_str):
            raise FileNotFoundError(f"GFF file not found: {gff_str}")
        try:
            gff_cols = ["seqid", "source", "type",
                        "start", "end", "score",
                        "strand", "frame", "attributes"]
            features = pd.read_csv(gff_str,
                                sep="\t",
                                comment="#",
                                names=gff_cols,
                                dtype={"seqid": str, 
                                        "source": str, 
                                        "type": str,
                                        "start": int,
                                        "end": int,
                                        "score": float,
                                        "strand": str,
                                        "frame": str,
                                        "attributes": str}
                                        )
        except Exception as e:
            raise ValueError(f"GFF failed to be loaded as TSV: {gff_str}") from e
        
        features["feature"] = features["attributes"].apply(extract_feature_name)

    return features
    
#%%

def merge_annotated_intervals(intervals):
    """
    Merge overlapping intervals in GFF file. 
    Used to extract unannotated intergenic regions downstream.
    
    :param intervals: features[["start", "end"]]
    """
    sorted_intervals = intervals.sort_values(by="start").values.tolist()
    merged = [sorted_intervals[0]]
    for s, e in sorted_intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e) #reassign the end of the last interval to the max between the current interval and the last interval (extend)
        else:
            merged.append([s, e])
    return(merged)
# %%
def main():
    vcf_in = load_variants(args.input_vcf)
    gff_in = load_gff(args.gff)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Adjust variant coordinates based on reference genome changes.")
    parser.add_argument("--input_vcf", required=True, help="Path to the input VCF file with original variant coordinates.")
    parser.add_argument("--output_vcf", required=True, help="Path to the output VCF file with adjusted variant coordinates.")
    parser.add_argument("--gff", required=True, help="Path to the file containing reference genome changes.")
    parser.add_argument("--fai", required=True, help="Path to the reference genome index file (.fai) to get reference length.")

    args = parser.parse_args()  
    main(args)