import sys
import os

def calculate_statistics(sequence):
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'Unknown': 0}
    total = 0
    for base in sequence.upper():
        if base in counts:
            counts[base] += 1
        else:
            counts['Unknown'] += 1
        total += 1
    counts['Total'] = total
    return counts

def display_statistics(filename, stats):
    print(f"\n{filename}")
    for base, count in stats.items():
        if base != 'Total':
            percentage = (count / stats['Total']) * 100
            print(f"{base}: {count} {percentage:.1f}%")
    print(f"Total: {stats['Total']}\n")

def main():
    if len(sys.argv) < 2:
        return

    all_stats = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'Unknown': 0, 'Total': 0}

    for filepath in sys.argv[1:]:
        if not os.path.exists(filepath):
            print(f"File not found: {filepath}")
            continue

        with open(filepath, 'r') as file:
            sequence = file.read().replace("\n", "").strip()
            file_stats = calculate_statistics(sequence)
            display_statistics(filepath, file_stats)
            for base in 'ACGTUnknownTotal':
                if base in file_stats:
                    all_stats[base] += file_stats[base]
                else:
                    all_stats['Unknown'] += 1

    if all_stats['Total'] > 0:
        print("All")
        display_statistics("All", all_stats)

if __name__ == "__main__":
    main()