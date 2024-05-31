#!/usr/bin/env python3
"""Encoded DNA sequence to FASTQ format conversion

This script takes two arguments as an inputs - path to the file where encoded DNA in binary format is and target
length of the decoded sequence. It prints decoded DNA sequence in FASTQ format on standard output.

This file contains following functions:
    * main - the main function of the script, parses command line parameters and calls the encoder

"""

import argparse
from dna_to_fastq_encoder import DnaToFastqEncoder


def main():
    # Set up parser
    parser = argparse.ArgumentParser(description='Prints on standard output the contents of the input file, '
                                                 'presented in FASTQ format.')
    # Define arguments for the parser
    parser.add_argument('-f', '--file', required=True, help='Path to the file with encoded DNA')
    parser.add_argument('-l', '--length', required=True, type=int, help='Length of the decoded sequence')
    # Parse arguments
    args = parser.parse_args()

    encoder = DnaToFastqEncoder(args.file, args.length)
    for fastq_str in encoder.decode_to_fastq():
        print(fastq_str)


if __name__ == "__main__":
    main()
