import sys
import argparse
import statistics

def extract_sequence(fasta_file):
    sequence = ""
    with open(fasta_file) as file:
        for line in file:
            if not line.startswith('>'):
                sequence += line.strip()
    return sequence

global_start = None; global_end = None
def extract_genes(genbank_file):
    global global_start; global global_end
    
    genes = []
    with open(genbank_file) as file:
        for line in file:
            if line.startswith("FEATURES"):
                break
        for line in file:
            if line.startswith("     source"):
                coords = line.strip().split()[1].split("..")
                global_start = int(coords[0])-1; global_end = int(coords[-1])
                break
        
        for line in file:
            if line.startswith("     gene"):
                coords = line.strip().split()[1]
                
                if coords.startswith("complement"):
                      coords = coords[11:-1]
                      orientation = "-"
                else: orientation = "+"
                
                if coords.startswith("join"):
                      coords = coords[5:-1]
                coords = coords.split("..")
                
                if coords[0].startswith("<") or coords[1].startswith(">"):
                    continue
                start = int(coords[0].split(",")[0])-1; end = int(coords[-1].split(",")[-1])
                
                if start < end:
                      length = end - start
                else: length = (global_end - start) + (end - global_start)
                
                if length < 150:
                    continue
                
                name = None
                for subline in file:
                    if subline.startswith("                     /locus_tag"):
                        name = subline.strip().split("=")[1][1:-1]
                        break
                
                genes.append({"name": name, "start": start, "end": end, "length": length, "orientation": orientation})
    return genes

def main():
    parser = argparse.ArgumentParser(description="Extract gene k-mers from a Fasta+GenBank file.")
  ##
    # File path arguments (required within the script)
    parser.add_argument("-f", "--fasta_file", help="Path to the Fasta file")
    parser.add_argument("-g", "--genbank_file", help="Path to the GenBank file")
    parser.add_argument("-o", "--output_file", help="Path to the Output file")
  ##
    # Optional arguments
    parser.add_argument("-k", type=int, default=5, help="Number of consecutive genes per k-mer")
    parser.add_argument("-i", action='store_true', default=False, help="Include gaps between genes in output")
  ##
    args = parser.parse_args()
    
    if not args.fasta_file or not args.genbank_file or not args.output_file:
        print("Error: Missing required arguments.")
        print("Usage: parser.py -f <fasta_file> -g <genbank_file> -o <output_file> [options]")
        exit(1)
    
    sequence = extract_sequence(args.fasta_file)
    gene_list = extract_genes(args.genbank_file)
    if global_start != 0 or global_end != len(sequence):
        print("Error: GenBank+FASTA file mismatch in length.")
        exit(1)
    
    #print()
    #for gene in gene_list:
    #    print(f'name: {gene["name"]}, start: {gene["start"]}, end: {gene["end"]}, length: {gene["length"]}, orientation: {gene["orientation"]}')
    #print()
    #print( "Statistics:")
    #print(f" + number of genes: {len(gene_list)}")
    
    gene_lengths = []
    if len(gene_list) == 0:
        gene_lengths.append(0)
    else:
        for gene in gene_list:
            gene_lengths.append(gene["length"])
    #print(f" + length of genes:")
    #print(f"    - min:   {int(min(gene_lengths))}")
    #print(f"    - max:   {int(max(gene_lengths))}")
    #print(f"    - mean:  {int(statistics.mean(gene_lengths))}")
    #print(f"    - stdev: {int(statistics.stdev(gene_lengths))}")
    
    gap_lengths = []
    if len(gene_list) == 0:
        gap_lengths.append(global_end - global_start)
    else:
        previous_end = gene_list[0]["end"]
        for next_gene in gene_list[1:]:
            gap_lengths.append(next_gene["start"] - previous_end)
            previous_end = next_gene["end"]
        if gene_list[0]["start"] < gene_list[0]["end"] and gene_list[-1]["start"] < gene_list[-1]["end"]:
            gap_lengths.append((global_end - gene_list[-1]["end"]) + (gene_list[0]["start"] - global_start))
        elif gene_list[0]["start"] > gene_list[0]["end"] and gene_list[-1]["start"] > gene_list[-1]["end"]:
            gap_lengths.append((gene_list[0]["start"] - global_end) + (global_start - gene_list[-1]["end"]))
        else:
            gap_lengths.append(gene_list[0]["start"] - gene_list[-1]["end"])
    #print(f" + length of gaps:")
    #print(f"    - min:   {int(min(gap_lengths))}")
    #print(f"    - max:   {int(max(gap_lengths))}")
    #print(f"    - mean:  {int(statistics.mean(gap_lengths))}")
    #print(f"    - stdev: {int(statistics.stdev(gap_lengths))}")
    
    #print()
    k = args.k; i = args.i
    if k > len(gene_list):
        print("Error: Not enough gene records for a k-mer.")
        exit(1)
    
    with open(args.output_file, "w") as file:
        if i:
            for pos in range(1-k, len(gene_list)-k+1):
                name = ""
                for q in range(k):
                    name += f'={gene_list[pos+q]["name"]}'
                
                start = gene_list[pos]["start"]
                end = gene_list[pos+k-1]["end"]
                if start <= end: kmer = sequence[start:end]
                else:            kmer = sequence[start:] + sequence[:end]
                
                file.write(f">{name[1:]}\n{kmer}\n")
        else:
            for pos in range(1-k, len(gene_list)-k+1):
                name = ""; kmer = ""
                for q in range(k):
                    name += f'|{gene_list[pos+q]["name"]}'
                    
                    start = gene_list[pos+q]["start"]
                    end = gene_list[pos+q]["end"]
                    if start <= end: kmer += sequence[start:end]
                    else:            kmer += sequence[start:] + sequence[:end]
                    
                file.write(f">{name[1:]}\n{kmer}\n")

if __name__ == "__main__":
    main()

