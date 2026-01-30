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
