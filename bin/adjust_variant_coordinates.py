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
    
