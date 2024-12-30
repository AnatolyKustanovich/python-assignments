import os
import glob
import logging
import pysam
import pandas as pd
import numpy as np

# Configure logging
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s',
                    filename='ngs_analysis.log')

def read_bed_file(bed_file):
    """
    Read BED file and return a list of features
    
    Parameters:
    -----------
    bed_file : str
        Path to BED file
    
    Returns:
    --------
    list of dict
        List of BED file features
    """
    if not bed_file or not os.path.exists(bed_file):
        return []
    
    features = []
    try:
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    feature = {
                        'chrom': parts[0],
                        'start': int(parts[1]),
                        'end': int(parts[2]),
                        'name': parts[3] if len(parts) > 3 else 'Unknown'
                    }
                    features.append(feature)
    except Exception as e:
        logging.error(f"Error reading BED file: {e}")
    
    return features

def find_overlapping_features(chrom, pos, bed_features):
    """
    Find features that overlap with a given position
    
    Parameters:
    -----------
    chrom : str
        Chromosome name
    pos : int
        Position to check
    bed_features : list
        List of BED file features
    
    Returns:
    --------
    list
        List of overlapping feature names
    """
    overlapping = []
    for feature in bed_features:
        if (feature['chrom'] == chrom and 
            feature['start'] <= pos < feature['end']):
            overlapping.append(feature['name'])
    return overlapping

def process_bam_file(bam_path, ref_fasta, bed_file=None):
    """
    Process a single BAM file and generate comprehensive coverage analysis
    """
    try:
        # Load bed features if provided
        bed_features = read_bed_file(bed_file)
        
        # Load reference FASTA
        ref_genome = pysam.FastaFile(ref_fasta)
        
        # Open BAM file
        bam_file = pysam.AlignmentFile(bam_path, "rb")
        
        # Collect results
        results = []
        
        # Iterate through reference sequences
        for chrom in ref_genome.references:
            # Get chromosome sequence
            chrom_seq = ref_genome.fetch(chrom)
            
            # Collect pileup information
            for pileupcolumn in bam_file.pileup(chrom):
                pos = pileupcolumn.pos
                
                # Count nucleotide and variant coverages
                base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'Del': 0, 'Ins': 0}
                five_prime_count = 0
                three_prime_count = 0
                
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        base_counts[base] += 1
                    
                    # Check for deletions and insertions
                    if pileupread.is_del:
                        base_counts['Del'] += 1
                    if pileupread.indel != 0:
                        base_counts['Ins'] += 1
                    
                    # Count read start and end positions
                    if pileupread.alignment.is_reverse:
                        if pileupread.alignment.reference_end == pos + 1:
                            three_prime_count += 1
                    else:
                        if pileupread.alignment.reference_start == pos:
                            five_prime_count += 1
                
                # Calculate total coverage
                total_coverage = sum(base_counts.values())
                
                # Get reference nucleotide
                ref_base = chrom_seq[pos]
                
                # Annotate gene/feature if BED file provided
                gene_feature = 'Unknown'
                if bed_features:
                    overlapping = find_overlapping_features(chrom, pos, bed_features)
                    if overlapping:
                        gene_feature = ','.join(overlapping)
                
                # Compile results
                result = {
                    'Gencoor': f"{chrom}:{pos+1}",
                    'Gene': gene_feature,
                    'RefSeq': ref_base,
                    'Coverage': total_coverage,
                    'A': base_counts['A'],
                    'C': base_counts['C'],
                    'G': base_counts['G'],
                    'T': base_counts['T'],
                    'Del': base_counts['Del'],
                    'Ins': base_counts['Ins'],
                    '5_Prime_End': five_prime_count,
                    '3_Prime_End': three_prime_count
                }
                results.append(result)
        
        # Convert to DataFrame
        df = pd.DataFrame(results)
        return df
    
    except Exception as e:
        logging.error(f"Error processing BAM file {bam_path}: {str(e)}")
        raise

def process_ngs_analysis(input_folder, ref_fasta, bed_file=None):
    """
    Process all BAM files in a folder and generate individual output files
    """
    # Validate input paths
    if not os.path.exists(input_folder):
        logging.error(f"Input folder does not exist: {input_folder}")
        raise ValueError(f"Input folder does not exist: {input_folder}")
    
    if not os.path.exists(ref_fasta):
        logging.error(f"Reference FASTA file does not exist: {ref_fasta}")
        raise ValueError(f"Reference FASTA file does not exist: {ref_fasta}")
    
    # Find all BAM files
    bam_files = glob.glob(os.path.join(input_folder, '*.bam'))
    
    if not bam_files:
        logging.warning(f"No BAM files found in {input_folder}")
        return
    
    for bam_path in bam_files:
        try:
            # Process individual BAM file
            logging.info(f"Processing {bam_path}...")
            result_df = process_bam_file(bam_path, ref_fasta, bed_file)
            
            # Generate output filename
            base_filename = os.path.splitext(os.path.basename(bam_path))[0]
            csv_output = os.path.join(input_folder, f"{base_filename}_analysis.csv")
            parquet_output = os.path.join(input_folder, f"{base_filename}_analysis.parquet")
            
            # Save files
            result_df.to_csv(csv_output, index=False)
            result_df.to_parquet(parquet_output, index=False)
            
            logging.info(f"Saved analysis results to {csv_output} and {parquet_output}")
        
        except Exception as e:
            logging.error(f"Failed to process {bam_path}: {str(e)}")

def main():
    # Specific paths for your project
    bam_dir = "/home/labs/schwartzlab/kanatoli/Projects/lib794_pseudoU/794_rRNA/"
    reference_fasta = "/home/labs/schwartzlab/kanatoli/Projects/lib794_pseudoU/794_rRNA/human_rRNA.fa"
    bed_file = "/home/labs/schwartzlab/kanatoli/Projects/lib794_pseudoU/794_rRNA/human_rRNA.bed"
    
    try:
        process_ngs_analysis(bam_dir, reference_fasta, bed_file)
        logging.info("NGS analysis completed successfully.")
    except Exception as e:
        logging.error(f"NGS analysis failed: {str(e)}")

if __name__ == "__main__":
    main()

# Required dependencies (install via pip):
# - pysam
# - pandas
# - numpy
# - pyarrow (for parquet support)