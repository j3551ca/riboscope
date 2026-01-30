#!/usr/bin/env python
#%%
import argparse
from intervaltree import IntervalTree
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
#%%
def get_ref_length(fai_path):
    """
    Obtain reference length from .fai for defining unannotated intervals.

    :param fai_path: Path to fasta index
    """
    if not os.path.isfile(fai_path):
            raise FileNotFoundError(f"FAI file not found: {fai_path}")
    try:
        with open(fai_path, "r") as f:
            for line in f:
                fai_cols = line.strip().split("\t")
                ref_length = int(fai_cols[1])
    except Exception as e:
        raise ValueError(f"Failed to read FAI file: {fai_path}") from e

    return ref_length

def find_unannotated_intervals(merged_intervals, ref_length):
    """
    Find regions from original reference seq that are unannotated -
    regions between the merged overlapping intervals provided in GFF file.

    :param merged_intervals: [[start, end], [start2, end2]]List of lists specifying merged annotated intervals
    :param ref_length: int Length of original reference used during alignment & variant calling
    """

    # construct list of unannotated regions based on merged annotated regions
    gaps = ([(1, merged_intervals[0][0] - 1)] if merged_intervals[0][0] > 1 else []
        + [(merged_intervals[i][1] + 1, merged_intervals[i + 1][0] - 1) for i in range(len(merged_intervals) - 1)]
        + [(merged_intervals[-1][1] + 1, ref_length)]
    )
    return [(s, e) for s, e in gaps if s <= e ] # make sure no whack intervals

#%%

def assemble_complete_gff(unannotated_intervals, features):
    """
    Assemble a complete GFF dataframe with annotated and unannotated intervals.
    Unannotated intervals are labeled as "intergenic".

    :param unannotated_intervals: [(start, end), (start2, end2)] List of tuples specifying unannotated intervals
    :param features: original user-annotated features dataframe from GFF file
    """

    intergenic_df = pd.DataFrame(
        [
            {"seqid": features["seqid"].iloc[0],
            "source": "python",
            "type": "unannotated",
            "start": s,
            "end": e,
            "score": "N/A",
            "strand": "+",
            "phase": "N/A",
            "feature": f"intergenic_{i+1}",
            }
            for i, (s,e) in enumerate(unannotated_intervals)
        ]

    )
    
    complete_gff = (pd.concat([features, intergenic_df], ignore_index=True)
    .sort_values("start")
    .reset_index(drop=True))

    return complete_gff

#%%

def map_variants_to_features(complete_gff, vcf):
    """
    Find variant position relative to any feature it is found in. 
    Adjusted positions are strand-aware and all features a variant overlaps with are reported
    (ie. same mutation found in 2 features +/- strand will be reported twice, with coordinates 
    relative to start of respective feature).

    :param complete_gff: GFF dataframe containing user-annotated & imputed unannotated regions of reference sequence
    :param vcf: VCF with positions relative to reference sequence
    """
    
    # Build interval tree (half-open intervals, so end + 1 for inclusive)
    # basically a dictionary of intervals, where any int overlapping it is key
    # store feature name and strand as data for strand-aware position calculation
    tree = IntervalTree.from_tuples(
        (r.start, r.end + 1, {"name": r.feature, "strand": r.strand}) for r in complete_gff.itertuples()
    )

    # query each variant in VCF to find associated feature(s)
    # calculate position in strand-aware way
    results = []

    for _, v in vcf.iterrows():
        if v.POS == 0 and v.CHROM == "0": # handles samples without variants
            results.append(
            {**v.to_dict(),
                "FEATURE_NAME": "N/A",
                "FEATURE_STRAND": "N/A",
                "FEATURE_POS": 0,
            }
        )
            continue
        for interval in tree[v.POS]:
            feature_start = interval.begin
            feature_end = interval.end - 1 # convert back from half-open that intervaltree does
            strand = interval.data["strand"]
            feature_name = interval.data["name"]
            feature_position = v.POS - feature_start if strand == "+" else feature_end - v.POS
            results.append(
                {**v.to_dict(),
                 "FEATURE_NAME": feature_name,
                 "FEATURE_STRAND": strand,
                 "FEATURE_POS": feature_position,
                }
            )

        updated_vcf = pd.DataFrame(results)

    return updated_vcf

# %%
def main():
    vcf_in = load_variants(args.input_vcf)
    gff_in = load_gff(args.gff)
    if gff_in=="NO_FILE":
        print("No GFF file provided, skipping coordinate adjustment.")
        vcf_out = vcf_in.copy()
        vcf_out["FEATURE_NAME"] = "reference"
        vcf_out["FEATURE_STRAND"] = "+",
        vcf_out["FEATURE_POS"] = vcf_out["POS"]
        vcf_out.to_csv(args.output_vcf, sep="\t", index=False)
        return
    else:
        merged_annotated_intervals = merge_annotated_intervals(gff_in[["start", "end"]])
        unannotated_intervals = find_unannotated_intervals(merged_annotated_intervals, ref_length=get_ref_length(args.fai))
        complete_gff = assemble_complete_gff(unannotated_intervals, gff_in)
        vcf_out = map_variants_to_features(complete_gff, vcf_in)
        vcf_out.to_csv(args.output_vcf, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Adjust variant coordinates based on reference genome changes.")
    parser.add_argument("--input_vcf", required=True, help="Path to the input VCF file with original variant coordinates.")
    parser.add_argument("--output_vcf", required=True, help="Path to the output VCF file with adjusted variant coordinates.")
    parser.add_argument("--gff", required=True, help="Path to the file containing reference genome changes.")
    parser.add_argument("--fai", required=True, help="Path to the reference genome index file (.fai) to get reference length.")

    args = parser.parse_args()  
    main(args)