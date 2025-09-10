#!/usr/bin/env python

import argparse
import pandas as pd

def import_results(fastp_qc, amplicon_counts, bam_qc, vcf_qc):
    with open("test_qc.txt", "w") as f:
        f.write("this is a test")
    read_qc = pd.read_csv()


def main(args):
    import_results()







if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Digest QC results and format for pass/fail logic.")
    parser.add_argument("--fastp_qc", type=str, required=True, help="Aggregated fastp results")
    parser.add_argument("--amplicon_counts", type=str, required=True, help="Aggregated counts of query sequences")
    parser.add_argument("--bam_qc", type=str, required=True, help="Aggregated qualimap alignment QC")
    parser.add_argument("--vcf_qc", type=str, required=True, help="Aggregated LoFreq SNPs")
    args = parser.parse_args()
    main(args)
