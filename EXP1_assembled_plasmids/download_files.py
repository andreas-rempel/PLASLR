from Bio import Entrez
import os

# Set up NCBI email for authentication
Entrez.email = "andreas_rempel@sfu.ca"

# Function to download a single accession
def download_accession(accession):
    # Fetch GenBank file
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    gb_file = f"{accession}.gb"
    with open(gb_file, "w") as outfile:
        outfile.write(handle.read())
    
    # Fetch FASTA file
    handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text")
    fasta_file = f"{accession}.fasta"
    with open(fasta_file, "w") as outfile:
        outfile.write(handle.read())
    
    print(f"Downloaded {accession}: GenBank and FASTA files saved.")

# Read accession numbers from file
accessions_file = "../accessions.txt"
if not os.path.exists(accessions_file):
    print("Error: Accessions file not found.")
    exit()

with open(accessions_file, "r") as infile:
    accessions = [line.strip() for line in infile]

# Download files for each accession
for accession in accessions:
    download_accession(accession)

