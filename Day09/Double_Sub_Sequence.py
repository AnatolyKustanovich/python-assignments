import argparse
from Bio import SeqIO

def find_longest_duplicate(sequence):
    n = len(sequence)
    longest = ""
    for length in range(1, n):  # Length of substring
        seen = set()
        for i in range(n - length + 1):
            subseq = sequence[i:i+length]
            if subseq in seen:
                if len(subseq) > len(longest):
                    longest = subseq
            else:
                seen.add(subseq)
    return longest

def calculate_gc_content(sequence):
    sequence = sequence.upper()
    gc_count = sequence.count("G") + sequence.count("C")
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0

def main():
    parser = argparse.ArgumentParser(description="Analyze DNA sequences from FASTA or GenBank files.")
    parser.add_argument("file", type=str, help="Path to the FASTA or GenBank file.")
    parser.add_argument("--duplicate", action="store_true", help="Find the longest subsequence that appears twice.")
    parser.add_argument("--gc", action="store_true", help="Calculate GC content of the sequence.")

    args = parser.parse_args()

    try:
        # Read the file and extract the sequence
        record = next(SeqIO.parse(args.file, "fasta" if args.file.endswith(".fasta") else "genbank"))
        sequence = str(record.seq)

        if args.duplicate:
            print("Finding the longest duplicated subsequence...")
            longest_duplicate = find_longest_duplicate(sequence)
            print(f"Longest duplicated subsequence: {longest_duplicate}")

        if args.gc:
            print("Calculating GC content...")
            gc_content = calculate_gc_content(sequence)
            print(f"GC Content: {gc_content:.2f}%")

    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
