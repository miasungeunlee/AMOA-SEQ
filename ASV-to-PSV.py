#!/usr/bin/env python3
# Author: Mia Sungeun Lee (sungeun.lee@ec-lyon.fr)
# Description: Converting CD-HIT clustering output file to PSV-ASV-list file
# Date: 20230307

import argparse

# Create the argument parser
parser = argparse.ArgumentParser(description='Converting CD-HIT clustering output file to PSV-ASV-list file.')
parser.add_argument('-i', '--input', type=str, required=True, help='Input CD-HIT file')
parser.add_argument('-o', '--output', type=str, required=True, help='Output PSV-ASV-list file')
# Parse the arguments
args = parser.parse_args()
input_file = args.input
output_file = args.output

psv_count = 0
psv_dict = {}

with open(input_file, "r") as f:
    for line in f:
        if line.startswith(">Cluster"):
            psv_id = "PSV_" + str(psv_count)
            psv_count += 1
            continue
        else:
            asv_id = line.split(">")[1].split("...")[0]
            if psv_id in psv_dict:
                psv_dict[psv_id].append(asv_id)
            else:
                psv_dict[psv_id] = [asv_id]

with open(output_file, "w") as f:
    f.write("OTU_ID\tASV-ID\n")
    for psv_id in psv_dict:
        for asv_id in psv_dict[psv_id]:
            f.write(psv_id + "\t" + asv_id + "\n")
