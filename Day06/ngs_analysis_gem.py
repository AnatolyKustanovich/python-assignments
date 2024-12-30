import pysam
import pandas as pd
import os

def process_bam(bam_file, reference_fasta, bed_file=None):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    
    # Create an empty DataFrame to store the results
    df = pd.DataFrame(columns=["Gencoor", "genes", "refseq", "Coverage", "A", "C", "G", "T", "Del", "Ins", "5' end", "3' end"])

    # Iterate over reads in the BAM file
    for read in samfile.fetch():
        chrom = read.reference_name
        pos = read.reference_start
        cigar = read.cigarstring
        seq = read.seq
        
        # Extract reference sequence from the FASTA file
        ref_seq = samfile.fetch(chrom, pos, pos+1)[0].seq
        
        # Calculate coverage and nucleotide-specific coverage
        coverage = 0
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'Del': 0, 'Ins': 0}
        five_prime_end = 0
        three_prime_end = 0
        
        for op, length in cigar:
            if op == 'M':  # Match or mismatch
                coverage += length
                for i in range(length):
                    counts[seq[i]] += 1
            elif op == 'D':  # Deletion
                coverage += length
                counts['Del'] += length
            elif op == 'I':  # Insertion
                coverage += length
                counts['Ins'] += length
            elif op == 'S':  # Soft clipping
                if pos == read.reference_start:
                    five_prime_end += length
                else:
                    three_prime_end += length
        
        # Append the extracted data to the DataFrame
        df = df.append({'Gencoor': pos+1, 'genes': None, 'refseq': ref_seq, 'Coverage': coverage, **counts, '5\' end': five_prime_end, '3\' end': three_prime_end}, ignore_index=True)
    
    # If a BED file is provided, annotate the DataFrame with gene information
    if bed_file:
        # Implement BED file parsing and annotation logic here
        # ...
    
    # Save the DataFrame to an RDS file

        df = pd.DataFrame(columns=["Gencoor", "genes", "refseq", "Coverage", "A", "C", "G", "T", "Del", "Ins", "5' end", "3' end"])
        df.to_pickle(bam_file.replace(".bam", ".rds"))
def main():
    bam_dir = "/home/labs/schwartzlab/kanatoli/Projects/lib794_pseudoU/794_rRNA/"
    reference_fasta = "/home/labs/schwartzlab/kanatoli/Projects/lib794_pseudoU/794_rRNA/human_rRNA.fa"
    bed_file = "/home/labs/schwartzlab/kanatoli/Projects/lib794_pseudoU/794_rRNA/human_rRNA.bed"  # Optional
    
    for bam_file in os.listdir(bam_dir):
        if bam_file.endswith(".bam"):
            process_bam(os.path.join(bam_dir, bam_file), reference_fasta, bed_file)

if __name__ == "__main__":
    main()