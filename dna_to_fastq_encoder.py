#!/usr/bin/env python3
"""
This module defines the DnaToFastqEncoder class, which is designed to decode binary-encoded DNA sequences into the
FASTQ format. The FASTQ format is a text-based format for storing both a nucleotide sequences and its corresponding
quality scores. Both the sequence and quality scores are encoded within a single binary file, where specific bits
represent different nucleotides and their quality.

Dependencies:
- Python 3.6 or later

Usage example:
    from dna_to_fastq_Encoder import DnaToFastqEncoder
    encoder = DnaToFastqEncoder(filepath='path/to/binary/file', l=100)
    encoder.decode_to_fastq()

"""


class DnaToFastqEncoder:
    """ A class used to decode DNA into FASTQ format

    This class takes a binary file as input, where each byte encodes a nucleotide and its quality score, and decodes it
    into human-readable FASTQ format.

    Args:
        filepath (str): The binary file to decode
        l (int): The length of the decoded DNA sequence

    """

    nucleotide_dict = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

    def __init__(self, filepath, l):
        """
        Args:
            filepath (str): The file to decode
            l (int): The length of the decoded DNA sequence

        Raises:
            ValueError: If the L-value is 0 or less.
        """
        self.filepath = filepath
        if l <= 0:
            raise ValueError(f"L-value error: Cannot decode sequences of length {l}. L-value must be a positive integer")
        else:
            self.l = l
        self._is_data_loaded = False
        self.data = None

    def read_from_file(self):
        """ Reads the file content into a bytes object

        In the object state it stores whether the data has been loaded.

        Raises:
            FileNotFoundError: If the file doesn't exist.
        """
        try:
            with open(self.filepath, 'rb') as file:
                self.data = file.read()
                self._is_data_loaded = True
        except FileNotFoundError:
            raise FileNotFoundError(f"File {self.filepath} not found") from None

    def fastq_string(self, nucleotide_seq, quality_scores, i):
        """ Creates string representing decoded nucleotide in FASTQ format

        Args:
            nucleotide_seq (str): Decoded nucleotide sequence
            quality_scores (str): Decoded quality scores
            i (int): Index of the record

        Returns:
            string: String in the FASTQ format.
        """
        return f"@READ_{i}\n{nucleotide_seq}\n+READ_{i}\n{quality_scores}"

    def decode_dna(self, piece):
        """ Decode piece of the DNA sequence

        Takes a piece of the loaded binary data and decodes ot into nucleotides and quality scores. A piece is a bytes
        object.

        Args:
            piece (bytes): Piece of the input DNA sequence with the length of L.

        Returns:
            tuple: Tuple of two strings nucleotides and their corresponding quality scores.
        """
        nucleotide_seq = ""
        quality_scores = ""
        # Each byte decodes one nucleotide, therefore for example for l-value 7 the nucleotide sequence length
        # will be 7
        for byte in piece:
            # Get the code of the nucleotide represented by two most significant bits with bitwise shift to right:
            # '11011000' >> 6 -> '00000011' represents "T"
            code = byte >> 6
            # Look up code in the nucleotide_dict
            nucleotide_seq += self.nucleotide_dict[code]

            # Get six least significant bits from byte using mask: '11011000' & '00111111' -> 00011000
            quality_score = byte & 0b00111111
            quality_scores += chr(quality_score + 33)

        return nucleotide_seq, quality_scores

    def generate_decoded_sequences(self):
        """ Slices the data into L number of sequences and decodes them into nucleotide sequences and corresponding
        quality scores.

        Yields:
            tuple: A tuple containing three elements:
                - str: The decoded nucleotide sequence
                - str: The decoded quality score
                - int: The index of the sequences

        Raises:
            RuntimeError: If the input data does not contain enough bytes to form a single sequence of length L.
        """
        if not self._is_data_loaded:
            self.read_from_file()

        # Number of pieces depends on the l-value
        n_pieces = len(self.data) // self.l
        if n_pieces == 0:
            raise RuntimeError(f"Insufficient data to decode: The input file '{self.filepath}' does not contain enough"
                               f" data to form a single sequence of length {self.l}.")

        for i in range(n_pieces):
            # Get the piece slice. Number of bytes in piece depends on the L-value
            piece = self.data[i * self.l:(i + 1) * self.l]
            # Get the nucleotides and their quality scores
            nucleotide_seq, quality_scores = self.decode_dna(piece)

            yield nucleotide_seq, quality_scores, i

    def decode_to_fastq(self):
        """ Decodes the DNA sequence and returns it in FASTQ format

        Yields FASTQ formatted strings for each DNA sequence and its quality score.

        Yields:
            string: FASTQ string
        """
        for nucleotide_seq, quality_scores, i in self.generate_decoded_sequences():
            yield self.fastq_string(nucleotide_seq, quality_scores, i+1)

