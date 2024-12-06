import sys
import argparse
import os

def extract_sequence(fasta_file):
    sequence = ""
    with open(fasta_file) as file:
        for line in file:
            if not line.startswith('>'):
                sequence += line.strip()
    return sequence

def main():
    parser = argparse.ArgumentParser(description="Create a summary of Platon results for one sample.")
  ##
    # File path arguments (required within the script)
    parser.add_argument("-a", "--accession", help="Accesseion number to identify sample")
    parser.add_argument("-s", "--sequence", help="Path to the original Fasta file")
    parser.add_argument("-f", "--fasta_file", help="Path to the k-mer Fasta file")
    parser.add_argument("-p", "--platon_file", help="Path to the Platon Fasta file")
    parser.add_argument("-o", "--output_file", help="Path to the Output table file")
    parser.add_argument("-m", "--mode", help="Expect 'plasmid' or 'chromosome' sequences")
  ##
    # Optional arguments
    parser.add_argument("-k", type=int, default=5, help="Number of consecutive genes per k-mer")
  ##
    args = parser.parse_args()
    
    if not args.fasta_file or not args.platon_file or not args.output_file:
        print("Error: Missing required arguments.")
        print("Usage: parser.py -f <fasta_file> -p <platon_file> -o <output_file> [options]")
        exit(1)
    
    if not args.mode:
        print("Error: Please specify expected sequence type.")
        print("       --mode plasmid  or  --mode chromosome")
        exit(1)
    
    k = args.k; mode = args.mode
    if not os.path.exists(args.output_file):
        output_file = open(args.output_file, 'w')
        output_file.write(f"{mode}_ID\tlength_(bp)\t#{k}-mers\t#{mode}_{k}-mers\tlargest_non-{mode}_island_(#{k}-mers)\n")
        output_file.close()
    
    sequence = extract_sequence(args.fasta_file)
    length = len(sequence)
    
    platon_file = open(args.platon_file)
    fasta_file = open(args.fasta_file)
    
    num_platon_kmers = 0
    num_fasta_kmers = 0
    
    initial_island_length = 0
    current_island_length = 0
    largest_island_length = 0
    
    for platon_line in platon_file:
        if platon_line.startswith('>'):
            num_platon_kmers += 1
            
            for fasta_line in fasta_file:
                if fasta_line.startswith('>'):
                    num_fasta_kmers += 1
                    
                    if platon_line != fasta_line:
                        initial_island_length += 1
                    else:
                        break
            break
    
    for platon_line in platon_file:
        if platon_line.startswith('>'):
            num_platon_kmers += 1
            
            for fasta_line in fasta_file:
                if fasta_line.startswith('>'):
                    num_fasta_kmers += 1
                    
                    if platon_line != fasta_line:
                        current_island_length += 1
                    else:
                        largest_island_length = max(current_island_length, largest_island_length)
                        current_island_length = 0
                        break
    
    for fasta_line in fasta_file:
        if fasta_line.startswith('>'):
            num_fasta_kmers += 1
            
            initial_island_length += 1
    largest_island_length = max(initial_island_length, largest_island_length)
    
    output_file = open(args.output_file, 'a')
    output_file.write(f"{args.accession}\t{length}\t{num_fasta_kmers}\t{num_platon_kmers}\t{largest_island_length}\n")
    output_file.close()

if __name__ == "__main__":
    main()

