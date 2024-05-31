# DNA to FASTQ Encoder

## Overview

The DNA to FASTQ Encoder provides tool to load binary DNA file and outputs encoded DNA sequence in FASTQ format. The project consists of a Python script (`main.py`) for command-line usage and a core Python module (`dna_to_fastq_encoder.py`) that defines the decoding logic.

### Requirements

- Python 3.6 or later

## Usage

### Command Line Interface

To decode a binary-encoded DNA file into FASTQ format and print the output to the standard output, use the `main.py` script as follows:

```bash
python main.py -f path/to/binary/file -l 80
```

Where:
- `-f` or `--file` specifies the path to the binary-encoded DNA file.
- `-l` or `--length` specifies the length of the decoded DNA sequence.

### As a Module

You can also use the `DnaToFastqEncoder` class within your Python projects:

```python
from dna_to_fastq_encoder import DnaToFastqEncoder

Encoder = DnaToFastqEncoder(filepath='path/to/binary/file', l=100)
for fastq_str in Encoder.decode_to_fastq():
    print(fastq_str)
```

## Testing

This project includes a testsuite to verify the functionality of the `DnaToFastqEncoder` class. To run the tests:

1. ensure you have `pytest` installed
2. set up `TEST_FILES_PATH` variable in `test_dna_to_fastq_encoder.py`
3. execute the following command in the project directory:

```bash
pytest
```

