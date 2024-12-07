import sys
from collections import Counter
import os

def calculate_statistics(sequence):
    """Calculate the statistics for a given sequence."""
    total = len(sequence)
    counts = Counter(sequence.upper())
    stats = {
        'A': counts.get('A', 0),
        'C': counts.get('C', 0),
        'G': counts.get('G', 0),
        'T': counts.get('T', 0),
        'Unknown': total - sum(counts.get(base, 0) for base in 'ACGT'),
        'Total': total
    }
    return stats

def display_statistics(filename, stats):
    """Display the statistics for a given file."""
    print(f"\n{filename}")
    total = stats['Total']
    for base in ['A', 'C', 'G', 'T', 'Unknown']:
        count = stats[base]
        percentage = (count / total * 100) if total > 0 else 0
        print(f"{base}: {count:8d} {percentage:6.1f}%")
    print(f"Total: {total:8d}\n")

def main():
    if len(sys.argv) < 2:
        print("Usage: python seq.py <file1> <file2> ...")
        return

    all_stats = Counter()

    for filepath in sys.argv[1:]:
        if not os.path.exists(filepath):
            print(f"File not found: {filepath}")
            continue

        with open(filepath, 'r') as file:
            sequence = file.read().replace("\n", "").strip()
            file_stats = calculate_statistics(sequence)
            display_statistics(filepath, file_stats)
            all_stats.update(file_stats)

    if all_stats['Total'] > 0:
        print("All")
        display_statistics("All", all_stats)

if __name__ == "__main__":
    main()
