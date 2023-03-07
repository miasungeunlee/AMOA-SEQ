#!/usr/bin/env python3
# Author: Mia Sungeun Lee (sungeun.lee@ec-lyon.fr)
# Description: Converting CD-HIT clustering output file to OTU-ASV-list file
# Date: 20230307

import argparse

# Create the argument parser
parser = argparse.ArgumentParser(description='Converting CD-HIT clustering output file to OTU-ASV-list file.')
parser.add_argument('-i', '--input', type=str, required=True, help='Input CD-HIT file')
parser.add_argument('-o', '--output', type=str, required=True, help='Output OTU-ASV-list file')
# Parse the arguments
args = parser.parse_args()
input_file = args.input
output_file = args.output

otu_count = 0
otu_dict = {}

with open(input_file, "r") as f:
    for line in f:
        if line.startswith(">Cluster"):
            otu_id = "OTU_" + str(otu_count)
            otu_count += 1
            continue
        else:
            asv_id = line.split(">")[1].split("...")[0]
            if otu_id in otu_dict:
                otu_dict[otu_id].append(asv_id)
            else:
                otu_dict[otu_id] = [asv_id]

with open(output_file, "w") as f:
    f.write("OTU_ID\tASV-ID\n")
    for otu_id in otu_dict:
        for asv_id in otu_dict[otu_id]:
            f.write(otu_id + "\t" + asv_id + "\n")
