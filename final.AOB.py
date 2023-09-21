#!/usr/bin/env python3
# Author: Mia Sungeun Lee (sungeun.lee@ec-lyon.fr)
# Description: Removing ambigous sequences for AOB amoA sequences           
# Date: 20230921

import argparse
from Bio import SeqIO

def filter_sequences(input_file, output_file):
    # List of characters to check at the start and end of each sequence
    start_chars_to_check = ['H']
    end_chars_to_check = ['FG', 'HG', 'YG']

    # Subsequences to check
    subsequences_to_check = ['GGGG', 'XTP', 'XPG', 'CCS', 'LSX', 'CCS', 'CCX', 'CSX', 'TGX', 'NWX', 'XAE', 'XAF', 'XRA', 'XEQ', 'XRT', 'XRN']

    # Open the output file
    with open(output_file, 'w') as out:
        # Parse the FASTA file
        for record in SeqIO.parse(input_file, "fasta"):
            # Convert the sequence to a string
            sequence = str(record.seq)
            # If the sequence does not start with a specified character or does not end with a specified character, write the ID to the output file
            if sequence[0] not in start_chars_to_check or sequence[-2:] not in end_chars_to_check:
                out.write(record.id + '\n')
            # If the sequence contains any of the specified subsequences, write the ID to the output file
            elif any(subsequence in sequence for subsequence in subsequences_to_check):
                out.write(record.id + '\n')
            # If the sequence contains an asterisk (*), write the ID to the output file
            elif '*' in sequence:
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

