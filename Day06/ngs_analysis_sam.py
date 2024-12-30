import os
import subprocess
import csv
from datetime import datetime
import pandas as pd


def run_samtools_command(command):
    """Run a samtools command and return the output."""
    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {command}")
        print(e.stderr)
        raise


def parse_bed_file(bed_file):
    """Parse the BED file to extract genomic features."""
    bed_data = []
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            fields = line.strip().split("\t")
            chrom, start, end, name = fields[0], int(fields[1]), int(fields[2]), fields[3]
            bed_data.append({"chrom": chrom, "start": start, "end": end, "name": name})
    return bed_data


def get_coverage_and_genes(bam_file, fasta_file, bed_data):
    """Get coverage and gene annotations for each position."""
    command = f"samtools depth -aa -r {fasta_file} {bam_file}"
    output = run_samtools_command(command)

    coverage_data = []
    for line in output.strip().split("\n"):
        chrom, pos, coverage = line.split("\t")
        pos = int(pos)

        # Find overlapping genes for the position
        overlapping_genes = []
        for entry in bed_data:
            if entry["chrom"] == chrom and entry["start"] <= pos <= entry["end"]:
                overlapping_genes.append(entry["name"])

        coverage_data.append({
            "Chromosome": chrom,
            "Position": pos,
            "Coverage": int(coverage),
            "Genes": ",".join(overlapping_genes) if overlapping_genes else "None"
        })

    return coverage_data


def write_rds_format(output_file, data):
    """Write data to RDS format using pyreadr."""
    import pyreadr
    df = pd.DataFrame(data)
    pyreadr.write_rds(output_file, df)


    

def main():
    input_folder = "C:/FGS/Python_Szab/python-assignments/Day06/py_exp"  # Update this path
    output_folder =   "C:/FGS/Python_Szab/python-assignments/Day06/"
    fasta_file = "/home/labs/schwartzlab/monikaw/data/lib794/genome/human_rRNA.fa"  # Update this path
    bed_file = "C:/FGS/Python_Szab/python-assignments/Day06/py_exp/human_rRNA.bed"  # Update this path




    os.makedirs(output_folder, exist_ok=True)
    bam_files = [f for f in os.listdir(input_folder) if f.endswith(".bam")]
    bed_data = parse_bed_file(bed_file)

    for bam_file in bam_files:
        bam_path = os.path.join(input_folder, bam_file)
        print(f"Processing {bam_file}...")

        coverage_data = get_coverage_and_genes(bam_path, fasta_file, bed_data)

        output_file = os.path.join(output_folder, f"{os.path.splitext(bam_file)[0]}.rds")
        write_rds_format(output_file, coverage_data)

        print(f"Saved coverage and gene data for {bam_file} to {output_file}")


if __name__ == "__main__":
    main()
