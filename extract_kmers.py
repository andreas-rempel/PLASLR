import sys, logging
import argparse
import json

level = logging.INFO

if not sys.stderr.isatty():
    logging.basicConfig( level=level, format='[%(asctime)s] %(message)s',
                                     datefmt='%Y-%m-%d %H:%M:%S' )
elif level == logging.DEBUG:
    logging.basicConfig( level=level, format='[%(asctime)s] %(message)s',
                                     datefmt='%H:%M:%S' )
else:
    logging.basicConfig( level=level, format='%(message)s' )


import networkx as nx
#from math import floor, ceil
#from math import sqrt, cbrt
#from matplotlib import pyplot as plt

G = nx.DiGraph()
node_IDs = {}
NODE_ID = 1


def build_gene_database(fastq_file):
    genes = {}
    with open(fastq_file) as file:
        
        while (line := file.readline()):
            if line.startswith("@"):
                 name = line.split()[0][1:]
            else:
                logging.error("Error: unexpected line in FASTQ file")
                exit(1)
            
            line = file.readline()
            sequence = line.strip()
            
            line = file.readline()
            if line.startswith("+"):
                 pass
            else:
                logging.error("Error: unexpected line in FASTQ file")
                exit(1)
            
            line = file.readline()
            genes[name] = sequence
    return genes


def reverse_complement(seq):
    translation_table = str.maketrans('ATCG', 'TAGC')
    return seq.translate(translation_table)[::-1]


def normalize_node(genes, sequence):
    fw = 0; rv = 0
    for gene in genes:
        if gene[0] == '+': fw += 1
        elif gene[0] == '-': rv += 1
    
    reverse = None
    if fw < rv:
        reverse = True
    elif fw > rv:
        reverse = False
    else:
        for i in range(len(genes)):
            if genes[i][0] == '-' and genes[-i-1][0] == '-':
                reverse = True; break
            elif genes[i][0] == '+' and genes[-i-1][0] == '+':
                reverse = False; break
    
    if reverse == None:
        fw_seq = sequence
        rv_seq = reverse_complement(sequence)
        if fw_seq > rv_seq:
            reverse = True
        elif fw_seq < rv_seq:
            reverse = False
        else:
            logging.warning("Warning: unexpected k-mer node detected")
    
    if reverse:
        canon = genes[::-1]
        for i in range(len(canon)):
            if canon[i][0] == '+': canon[i] = f"-{canon[i][1:]}"
            elif canon[i][0] == '-': canon[i] = f"+{canon[i][1:]}"
        sequence = reverse_complement(sequence)
    else:
        canon = genes[:]
    
    return tuple(canon), sequence


def normalize_edge(ID1, genes1, sequence1, ID2, genes2, sequence2):
    if genes1[1:] == genes2[:-1]:
        return (ID1, '+', ID2, '+')
    elif genes2[1:] == genes1[:-1]:
        return (ID2, '+', ID1, '+')
    else:
        for i in range(1, len(genes1)):
            if genes1[i][0] == genes2[-i][0] or genes1[i][1:] != genes2[-i][1:]:
                break
        else:
            fw_seq = sequence1 + reverse_complement(sequence2)
            rv_seq = sequence2 + reverse_complement(sequence1)
            if fw_seq < rv_seq:
                return (ID1, '+', ID2, '-')
            elif fw_seq > rv_seq:
                return (ID2, '+', ID1, '-')
            else:
                logging.error("Error: unexpected k-mer edge detected")
                exit(1)
        
        for i in range(len(genes2)-1):
            if genes1[-i-2][0] == genes2[i][0] or genes1[-i-2][1:] != genes2[i][1:]:
                break
        else:
            fw_seq = reverse_complement(sequence1) + sequence2
            rv_seq = reverse_complement(sequence2) + sequence1
            if fw_seq < rv_seq:
                return (ID1, '-', ID2, '+')
            elif fw_seq > rv_seq:
                return (ID2, '-', ID1, '+')
            else:
                logging.error("Error: unexpected k-mer edge detected")
                exit(1)


parser = argparse.ArgumentParser(description="Extract gene k-mers from a FASTQ+JSON file.")
##
# File path arguments (required within the script)
parser.add_argument("-f", "--fastq_file", help="Path to the FASTQ file")
parser.add_argument("-g", "--json_file", help="Path to the JSON file")
parser.add_argument("-o", "--output_file", help="Path to the Output file")
##
# Optional arguments
parser.add_argument("-k", type=int, default=5, help="Number of consecutive genes per k-mer")
##
args = parser.parse_args()

if not args.fastq_file or not args.json_file or not args.output_file:
    logging.error("Error: Missing required arguments.")
    logging.error("Usage: parser.py -f <fastq_file> -g <json_file> -o <output_file> [options]")
    exit(1)


gene_list = build_gene_database(args.fastq_file)

with open(args.json_file, 'r') as file:
    data = json.load(file)
    for read, genes in data.items():
        PREV_ID = 0; NEXT_ID = 0
        kmer_genes = []; kmer_sequence = []
        
        for gene in genes:
            gene_orientation = gene[0]
            gene_name = gene[1:]
            
            if gene_name in gene_list:
                kmer_genes.append(gene)
                
                gene_sequence = gene_list[gene_name]
                if gene_orientation == '+':
                    kmer_sequence.append(gene_sequence)
                else:
                    kmer_sequence.append(reverse_complement(gene_sequence))
            else:
                logging.warning(f"Warning: missing sequence for {gene_name}")
            
            if len(kmer_genes) == args.k:
                kmer, sequence = normalize_node(kmer_genes, kmer_sequence[args.k // 2])
                
                if kmer not in node_IDs:
                    G.add_node(NODE_ID)
                    G.nodes[NODE_ID]['SEQ'] = sequence
                    G.nodes[NODE_ID]['KMER'] = kmer
                    G.nodes[NODE_ID]['LEN'] = len(sequence)
                    G.nodes[NODE_ID]['LR'] = set()
                    node_IDs[kmer] = NODE_ID; NODE_ID += 1
                
                NEXT_ID = node_IDs[kmer]
                G.nodes[NEXT_ID]['LR'].add(read)
                
                if PREV_ID and not G.has_edge(PREV_ID, NEXT_ID)\
                           and not G.has_edge(NEXT_ID, PREV_ID):
                    
                    prev_kmer = G.nodes[PREV_ID]['KMER']
                    prev_sequence = G.nodes[PREV_ID]['SEQ']
                    ID1, ori1, ID2, ori2 = normalize_edge(PREV_ID, prev_kmer, prev_sequence, NEXT_ID, kmer, sequence)
                    G.add_edge(ID1, ID2)
                    G.edges[ID1, ID2]["From"] = ori1
                    G.edges[ID1, ID2]["To"] = ori2
                
                PREV_ID = NEXT_ID; NEXT_ID = 0
                del kmer_genes[0]; del kmer_sequence[0]


L = 0; Cov = 0
for node_id, attr in G.nodes(data=True):
    LEN = attr.get('LEN', 0)
    LR = attr.get('LR', set())
    XCov = len(LR) * LEN
    L += LEN; Cov += XCov
AvgCov = Cov / L

for node_id, attr in G.nodes(data=True):
    LEN = attr.get('LEN', 0)
    LR = attr.get('LR', set())
    XCov = len(LR) * LEN
    XAvgCov = XCov / LEN if LEN > 0 else 0
    XNormCov = XAvgCov / AvgCov
    G.nodes[node_id]['DP'] = XNormCov


with open(args.output_file, 'w') as fout:
    fout.write("H\n")
    
    # Write sequences
    for node_id, attr in G.nodes(data=True):
        SEQ = attr.get('SEQ', '')
        KMER = attr.get('KMER', tuple())
        LEN = attr.get('LEN', '')
        DP = attr.get('DP', '')
        LR = attr.get('LR', set())
        
        KMER = str(list(KMER)).replace("'", "").replace(' ', '')
        LR = str(list(LR)).replace("'", "").replace(' ', '')
        fout.write(f"S\t{node_id}\t{SEQ}\tGN:Z:{KMER}\tLN:i:{LEN}\tdp:f:{DP}\tLR:Z:{LR}\n")
    
    # Write links
    for u, v, attr in G.edges(data=True):
        FROM = attr.get('From', '')
        TO = attr.get('To', '')
        CIGAR = attr.get('CIGAR', '')
        fout.write(f"L\t{u}\t{FROM}\t{v}\t{TO}\t{CIGAR}\n")


# Position nodes
components = list(nx.weakly_connected_components(G))
components = sorted(components, key=lambda x: -len(x))
num_components = len(components)
logging.info(f"Number of connected components: {num_components}")

"""
num_pages = floor(cbrt(num_components))
num_plots = ceil(num_components/num_pages)
num_rows = floor(sqrt(num_plots))
num_columns = ceil(num_plots/num_rows)
logging.info('\n')

for i in range(num_pages):
    logging.info(f"Page: {i+1}/{num_pages}")
    fig, axes = plt.subplots(num_rows, num_columns)
    axes = axes.flatten()[:num_plots]

    for j in range(num_plots):
        if i*num_plots+j == num_components: break
        subgraph = nx.induced_subgraph(G, components[i*num_plots+j])
        layout = nx.kamada_kawai_layout(subgraph)
        nx.draw(subgraph, layout, ax=axes[j], with_labels=False)

    plt.show()
    logging.info('\r')
logging.info('\n')
"""
