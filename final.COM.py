#!/usr/bin/env python3
# Author: Mia Sungeun Lee (sungeun.lee@ec-lyon.fr)
# Description: Removing ambigous sequences for Comammox amoA sequences           
# Date: 20230921

import argparse
from Bio import SeqIO

def filter_sequences(input_file, output_file):
    # List of characters to check at the start and end of each sequence
    start_chars_to_check = ['AA','AK','AN','AV','AY','EK','ID','LD','NM','PF','SH','SN','TH','TW','VA','VD','VI','VL','VV'] 
    end_chars_to_check = ['DF','DL','DV','DY','GW','LW','PI','PL','PM','PT','PV','WW'] 

    # Open the output file
    with open(output_file, 'w') as out:
        # Parse the FASTA file
        for record in SeqIO.parse(input_file, "fasta"):
            # Convert the sequence to a string
            sequence = str(record.seq)
            # If the sequence does not start with a specified character or does not end with a specified character, write the ID to the output file
            if str(record.seq[:2]) not in start_chars_to_check or str(record.seq[-2:]) not in end_chars_to_check:
                out.write(record.id + '\n')
            # If the sequence contains an asterisk (*), write the ID to the output file
            elif '*' in record.seq :
                out.write(record.id + '\n')

if __name__ == "__main__":
    # Create the parser
    parser = argparse.ArgumentParser(description="Filter sequences")

    # Add the arguments
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output file')

    # Parse the arguments
    args = parser.parse_args()

    # Call the function with the arguments
    filter_sequences(args.input, args.output)
