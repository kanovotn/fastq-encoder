#!/usr/bin/env python3
"""
This module contains a testsuite for the DnaToFastqEncoder class, ensuring its functionality
in decoding DNA sequences into the FASTQ format is correct in various scenarios.

Dependencies:
- pytest: Used for organizing and running the tests, including support for fixtures and parameterized testing.
- tempfile: Utilized for creating temporary files that simulate input data for the encoder.
- os: Employed for file system operations, such as obtaining file statistics to support certain tests.

Usage:
The tests in this module can be run using the pytest framework, by running command 'pytest' in the directory where this
file resides.
"""

import pytest
import tempfile
import os
from dna_to_fastq_encoder import DnaToFastqEncoder

TEST_FILES_PATH = ""

@pytest.fixture
def dna_encoder():
    """ This is a pytest fixture, to help to initialize DnaToFastqEncoder with empty file and some dummy lvalue.
    It is used by some tests in this module.

    """
    tmpfile = tempfile.NamedTemporaryFile()
    tmpfile_name = tmpfile.name
    # Close and delete the file as we don't need it for the test purposes
    tmpfile.close()

    return DnaToFastqEncoder(tmpfile_name, 2)


@pytest.mark.parametrize("lvalue", [0, -5])
def test_l_value_less_than_zero(lvalue):
    """ Test that DnaToFastqEncoder raises a ValueError for invalid lvalues.

        This test verifies that initializing DnaToFastqEncoder with lvalues less than or equal to zero
        raises a ValueError with an appropriate error message. The test is parameterized to check multiple
        invalid lvalues.

    Args:
        lvalue (int): An invalid lvalue (sequence length) passed to DnaToFastqEncoder.

    Asserts:
        Raises ValueError when invalid lvalues are passed to the DnaToFastqEncoder.
    """
    with pytest.raises(ValueError) as exception_info:
        encoder = DnaToFastqEncoder("dummy", lvalue)

    assert str(exception_info.value) == (f"L-value error: Cannot decode sequences of length {lvalue}. L-value must be"
                                         f" a positive integer")


def test_read_from_file():
    """ Test that DnaToFastqEncoder correctly reads content of the binary file.

        This test verifies that DnaToFastqEncoder.read_from_file() methods reads whole file passed to the
        DnaToFastqEncoder upon initialization. The test creates its on temporary binary file with known content.

        Asserts:
            The whole content of the test file has been read.
    """
    with tempfile.NamedTemporaryFile(delete=False) as tmpfile:
        tmpfile.write(b'\x00\x01\x02\x03')
        tmpfile_name = tmpfile.name

    encoder = DnaToFastqEncoder(tmpfile_name, 2)
    encoder.read_from_file()

    # Clean up
    os.remove(tmpfile_name)
    # Check that the data was correctly loaded
    assert encoder.data == b'\x00\x01\x02\x03', f"The file content has not been read correctly"


def test_read_from_file_not_found(dna_encoder):
    """ Test that DnaToFastqEncoder raises a FileNotFoundError for non-existing file.

        This test verifies that when DnaToFastqEncoder tries to read non-existing file it raises a FileNotFoundError
        with an appropriate error message. The test uses fixture dna_encoder which sets up DnaToFastqEncoder for test.

        Args:
            dna_encoder (fixture): A pytest fixture that is used for setting up DnaToFastqEncoder.

        Asserts:
            - FileNotFoundError has been raised.
    """
    with pytest.raises(FileNotFoundError) as exception_info:
        dna_encoder.read_from_file()

    assert str(exception_info.value) == f"File {dna_encoder.filepath} not found", \
        f"The expected FileNotFound Exception has not been raised"


def test_fastq_string(dna_encoder):
    """ Test that DnaToFastqEncoder return string in FASTQ format.

        This test verifies that DnaToFastqEncoder returns expected string on the example of
        (nucleotide_seq, quality_score) tuple. It rewrites the decoded DNA sequence into FASTQ string format.
        The test uses fixture dna_encoder which sets up DnaToFastqEncoder for test.

        Args:
            dna_encoder (fixture): A pytest fixture that is used for setting up DnaToFastqEncoder.

        Asserts:
            - Match of the FASTQ string from DnaToFastqEncoder and expected result.
    """
    nucleotide_seq = "AGCAATG"
    quality_score = "^T+;<7G"
    fastq_string = dna_encoder.fastq_string(nucleotide_seq, quality_score, 1)

    assert fastq_string == f"@READ_1\n{nucleotide_seq}\n+READ_1\n{quality_score}", \
        f"FASTQ string doesn't match the expected value"


def test_decode_dna(dna_encoder):
    """ Test that DnaToFastqEncoder correctly decodes a given byte sequence into nucleotides and quality scores.

        This test verifies that 'decode_dna()' method correctly decodes given byte sequence, translating it into
        expected nucleotide sequence and corresponding quality scores. The byte sequence tested includes variety of
        values to verify correct functionality.

        Args:
            dna_encoder (fixture): A pytest fixture that is used for setting up DnaToFastqEncoder.

        Asserts:
           - Nucleotide sequence and quality score match expected value.
    """
    # bytes in piece: '01100001', '00111111', '10000000', '11111010'
    piece = b'\x61\x3F\x80\xFA'
    nucleotide_seq, quality_scores = dna_encoder.decode_dna(piece)

    assert nucleotide_seq == 'CAGT', f"The nucleotide sequence doesn't match to the expected value"
    assert quality_scores == 'B`![', f"The quality score doesn't match to the expected value"


@pytest.mark.parametrize("lvalue, expected_result", [
    (4, [
        ("CTAG", "3.$2"),
        ("AATC", ")9&0"),
        ("TATT", "\"6]S")
    ]),
    (7, [("CTAGAAT", "3.$2)9&")])
])
def test_generate_decoded_sequences(lvalue, expected_result):
    """ Test that 'generate_decoded_sequence() correctly yields nucleotide sequences and quality scores.

        This test verifies that for given l-values, the 'generate_decoded_sequences' method yields the correct number
        of nucleotide sequences and quality scores, each of them matching the expected results. It checks
        the functionality across different l-values and expected outcomes. The test is parameterized to automatically
        run with different sets of l-values and corresponding expected results.

        Args:
            lvalue (int): The length of the DNA sequence to be used by DnaToFastqEncoder.
            expected_result (list of tuples): Expected nucleotide sequence and quality scores for each l-value.

        Asserts:
            - The length of each decoded nucleotide sequence and quality score.
            - The actual nucleotide sequence and quality score match the expected result.
            - The total number of yielded sequences matches the expected number based on the input file size and
             the l-value.

    """
    file_name = f"{TEST_FILES_PATH}/input12"
    file_stats = os.stat(file_name)

    encoder = DnaToFastqEncoder(file_name, lvalue)
    n_yields = 0
    for nucleotide_seq, quality_score, i in encoder.generate_decoded_sequences():
        assert len(nucleotide_seq) == encoder.l, f"The length of nucleotide_seq doesn't match l-value"
        assert len(quality_score) == encoder.l, f"The length of quality_score doesn't match l-value"
        assert (nucleotide_seq, quality_score) == expected_result[i], f"The actual sequence doesn't match the expected"
        n_yields += 1

    expected_yields = file_stats.st_size // lvalue
    assert n_yields == expected_yields, f"Expected {expected_yields} yields, got only {n_yields}"


def test_generate_decoded_sequences_lvalue_bigger_than_file():
    """ Test that using bigger l-value than number of bytes in the file raises RuntimeError.

        This test verifies that when DnaToFastqEncoder is initialized with bigger l-value than input file size in bytes,
        leading to decode sequences of the length 0, it raises RuntimeError.

        Asserts:
            - RuntimeError has been raised with an appropriate error message.
    """
    file_name = f"{TEST_FILES_PATH}/input12"
    lvalue = 13
    encoder = DnaToFastqEncoder(file_name, 13)

    with pytest.raises(RuntimeError) as exception_info:
        for nucleotide_seq, quality_score, i in encoder.generate_decoded_sequences():
            pass

    assert str(exception_info.value) == (f"Insufficient data to decode: The input file '{file_name}' does not"
                                         f" contain enough data to form a single sequence of length {lvalue}.")


@pytest.mark.parametrize("lvalue", [7, 15, 80])
def test_decode_to_fastq(lvalue):
    """ Test that for different l-values, DnaToFastqEncoder gives complete expected output in FASTQ format.

        This test verifies that initializing DnaToFastqEncoder with different l-values, leads to the expected results,
        matched with the sample outputs. The test is parameterized to check multiple l-values.

        Args:
            lvalue (int): The length of the DNA sequences to be decoded.

        Asserts:
            Asserts that 'decode_to_fastq()' returns string which matches with the content of the sample output files.

    """
    result = ""
    output = None
    encoder = DnaToFastqEncoder(f"{TEST_FILES_PATH}/input", lvalue)

    for fastq_str in encoder.decode_to_fastq():
        result += fastq_str + "\n"
    output_file = f"{TEST_FILES_PATH}/output{lvalue}"
    with open(output_file, 'r') as file:
        output = file.read()

    assert result == output, f"The decoded output for l-value {lvalue} doesn't match the expected output"

