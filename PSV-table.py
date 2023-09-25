#!/usr/bin/env python3
# Author: Mia Sungeun Lee (sungeun.lee@ec-lyon.fr)
# Description: Converting ASV count table to OTU count table
# Date: 20230307

import argparse
import pandas as pd
# Parse the command-line arguments
parser = argparse.ArgumentParser(description="Converting ASV count table to OTU count table")
parser.add_argument("-i", "--input", required=True, help="Input ASV count table file")
parser.add_argument("-t", "--otulist", required=True, help="Input OTU-ASV-list file")
parser.add_argument("-o", "--output", required=True, help="Output OTU count table file")
args = parser.parse_args()

# Read in the ASV count table
asv_counts = pd.read_csv(args.input, sep="\t", index_col=0)
# Read in the OTU list file
psv_list = pd.read_csv(args.otulist, sep="\t", index_col=1, header=None, names=["OTU_ID"])
# Group the ASV counts by their corresponding OTU
psv_counts = asv_counts.groupby(psv_list["PSV_ID"]).sum()
# Write the OTU count table to a file
psv_counts.to_csv(args.output, sep="\t")
