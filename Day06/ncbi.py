import argparse
import os
import csv
import datetime
from Bio import Entrez

# Configuration of email for NCBI
Entrez.email = "amkustanovich@gmail.com"

def search_ncbi(database, term, number):
    #Search NCBI database and return search results.
    #Some of NCBI databases: RefSeq, PubMed, Nucleotide, Genome, Protein
    print(f"Searching NCBI {database} for '{term}'...")

    # Perform search
    handle = Entrez.esearch(db=database, term=term, retmax=number)
    record = Entrez.read(handle)
    handle.close()
    ids = record['IdList']
    total_found = record['Count']
    return ids, total_found

def download_ncbi(database, ids, output_dir):
    #Download and save NCBI data.
    filenames = []
    for idx, record_id in enumerate(ids):
        handle = Entrez.efetch(db=database, id=record_id, rettype="gb", retmode="text")
        data = handle.read()
        handle.close()
        # Save results
        filename = os.path.join(output_dir, f"{database}_{record_id}.txt")
        with open(filename, "w") as f:
            f.write(data)
        filenames.append(filename)
        print(f"Saved: {filename}")
    return filenames

def log_search_to_csv(logfile, date, term, database, max_items, total_items):
    #Log search to a CSV file.
    log_exists = os.path.isfile(logfile)
    with open(logfile, "a", newline='') as csvfile:
        fieldnames = ['date', 'database', 'term', 'max', 'total']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        # Write header in new file
        if not log_exists:
            writer.writeheader()
        writer.writerow({'date': date, 'database': database, 'term': term, 'max': max_items, 'total': total_items})

def main():
    parser = argparse.ArgumentParser(description="Download data from NCBI databases.")
    parser.add_argument("--database", default="nucleotide", help="NCBI database to search (default: nucleotide)")
    parser.add_argument("--term", required=True, help="Search term")
    parser.add_argument("--number", type=int, default=10, help="Maximum number of records to download (default: 10)")

    args = parser.parse_args()

    # Parameters
    database = args.database
    term = args.term
    number = args.number
    output_dir = "ncbi_downloads"
    log_file = "ncbi_log.csv"

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Retrieve IDs
    ids, total_items = search_ncbi(database, term, number)

    # Download and save data
    if ids:
        filenames = download_ncbi(database, ids, output_dir)
        print("\nDownloaded files:")
        for f in filenames:
            print(f)

    # Log search details
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    log_search_to_csv(log_file, current_date, term, database, number, total_items)

if __name__ == "__main__":
    main()

