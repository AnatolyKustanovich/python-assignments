import os
import argparse
import pysam
import csv
from Bio import SeqIO
from pybedtools import BedTool

def parse_reference_fasta(fasta_file):
    """Parse the reference FASTA file."""
    reference = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        reference[record.id] = str(record.seq)
    return reference

def process_bam_file(bam_file, reference, bed=None):
    """Extract coverage data for a BAM file."""
    bam = pysam.AlignmentFile(bam_file, "rb")
    coverage_data = []
    
    for ref_name, sequence in reference.items():
        for pos in range(len(sequence)):
            gencoor = f"{ref_name}:{pos + 1}"
            refseq = sequence[pos]
            pileup = bam.pileup(ref_name, pos, pos + 1)
            
            coverage = 0
            base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
            del_count = 0
            ins_count = 0
            five_prime = 0
            three_prime = 0
            
            for pileupcolumn in pileup:
                coverage += pileupcolumn.nsegments
                for pileupread in pileupcolumn.pileups:
                    if pileupread.is_del:
                        del_count += 1
                    elif pileupread.indel > 0:
                        ins_count += 1
                    else:
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        base_counts[base] += 1
                        
                        if pileupread.query_position == 0:
                            five_prime += 1
                        if pileupread.query_position == len(pileupread.alignment.query_sequence) - 1:
                            three_prime += 1
            
            # Get gene annotations (optional)
            genes = ""
            if bed:
                bedtool = BedTool(bed)
                overlaps = bedtool.intersect(b=(ref_name, pos, pos + 1))
                genes = ",".join([feature.name for feature in overlaps])
            
            # Add row to coverage data
            coverage_data.append({
                "Gencoor": gencoor,
                "genes": genes,
                "refseq": refseq,
                "Coverage": coverage,
                "A": base_counts['A'],
                "C": base_counts['C'],
                "G": base_counts['G'],
                "T": base_counts['T'],
                "Del": del_count,
                "Ins": ins_count,
                "5' end": five_prime,
                "3' end": three_prime
            })
    
    return coverage_data

def write_table(output_file, coverage_data):
    """Write coverage data to a CSV file."""
    headers = ["Gencoor", "genes", "refseq", "Coverage", "A", "C", "G", "T", "Del", "Ins", "5' end", "3' end"]
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=headers)
        writer.writeheader()
        writer.writerows(coverage_data)

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="NGS Analysis Script")
    parser.add_argument("--bam_dir", required=True, help="Path to directory containing BAM files")
    parser.add_argument("--fasta_file", required=True, help="Path to reference FASTA file")
    parser.add_argument("--bed_file", help="Path to BED file for annotations (optional)")
    parser.add_argument("--output_file", required=True, help="Path to save output CSV file")
    args = parser.parse_args()
    
    # Parse reference FASTA
    reference = parse_reference_fasta(args.fasta_file)
    
    # Process BAM files
    all_coverage_data = []
    for bam_file in os.listdir(args.bam_dir):
        if bam_file.endswith(".bam"):
            coverage_data = process_bam_file(os.path.join(args.bam_dir, bam_file), reference, args.bed_file)
            all_coverage_data.extend(coverage_data)
    
    # Write the final table
    write_table(args.output_file, all_coverage_data)

if __name__ == "__main__":
    main()
