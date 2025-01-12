Sequence Analysis Tool

This tool analyzes DNA sequences from FASTA or GenBank files. It supports two optional analyses:

Finding the longest duplicated subsequence: Identifies the longest subsequence that appears at least twice in the given sequence.

Calculating GC content: Computes the GC content percentage of the sequence.

Requirements

Ensure you have Python installed, along with the required dependencies listed in requirements.txt.

Installation

Clone the repository or download the script.

Install the required dependencies:

pip install -r requirements.txt

Usage

Run the script with the desired options:

python analyze.py FILE [--duplicate] [--gc]

Positional Arguments:

FILE: Path to the input file in FASTA or GenBank format.

Optional Arguments:

--duplicate: Perform analysis to find the longest duplicated subsequence.

--gc: Perform analysis to calculate GC content.

Examples

Find the longest duplicated subsequence in a FASTA file:

python analyze.py example.fasta --duplicate

Calculate GC content of a sequence in a GenBank file:

python analyze.py example.gb --gc

Perform both analyses:

python analyze.py example.fasta --duplicate --gc

Output

The results of the analyses will be printed to the terminal. For example:

Finding the longest duplicated subsequence...
Longest duplicated subsequence: ATCG

Calculating GC content...
GC Content: 54.32%
