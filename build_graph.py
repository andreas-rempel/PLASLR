import networkx as nx
#import layout as nx2
from matplotlib import pyplot as plt

from sys import argv, stderr
from time import time
from math import sqrt, floor, ceil
from threading import Thread as thread

# Create a directed graph
G = nx.DiGraph(); index = {}
node_IDs = {}; NODE_ID = 1
L = 0; R = 1; ___ = ','

# Add nodes and edges
with open(argv[1], 'r') as file:
    while line := file.readline():
        if line.startswith('>'):
            kmer, pos, length, cläss, read = line[1:-1].split('|')
            accession, coordinates = pos.split('@')
            start, end = coordinates.split('-')
            sequence = file.readline()[:-1]
            
            if kmer not in node_IDs:
                G.add_node(NODE_ID)
                G.nodes[NODE_ID]['SEQ'] = sequence
                G.nodes[NODE_ID]['KMER'] = kmer
                G.nodes[NODE_ID]['ACC'] = accession
                G.nodes[NODE_ID]['START'] = int(start)
                G.nodes[NODE_ID]['END'] = int(end)
                G.nodes[NODE_ID]['LEN'] = int(length)
                G.nodes[NODE_ID]['CLASS'] = cläss
                G.nodes[NODE_ID]['LR'] = set()
                
                L_link, ___, R_core = kmer.partition(___)
                L_core, ___, R_link = kmer.rpartition(___)
                
                if R_core not in index:
                    index[R_core] = (set(), set())
                if L_link not in index[R_core][L]:
                    index[R_core][L].add(L_link)
                
                if L_core not in index:
                    index[L_core] = (set(), set())
                if R_link not in index[L_core][R]:
                    index[L_core][R].add(R_link)
                
                for L_link in index[L_core][L]:
                    Lmer = f"{L_link}{___}{L_core}"
                    LEFT_ID = node_IDs[Lmer]
                    if not G.has_edge(LEFT_ID, NODE_ID):
                        
                        BSTART = G.nodes[NODE_ID]['START']
                        AEND = G.nodes[LEFT_ID]['END']
                        if BSTART <= AEND:
                            LENGTH = AEND - BSTART + 1
                        else:
                            ASTART = G.nodes[LEFT_ID]['START']
                            ALENGTH = G.nodes[LEFT_ID]['LEN']
                            LENGTH = ALENGTH - (BSTART-ASTART)
                        
                        G.add_edge(LEFT_ID, NODE_ID)
                        G.edges[LEFT_ID, NODE_ID]["From"] = '+'
                        G.edges[LEFT_ID, NODE_ID]["To"] = '+'
                        G.edges[LEFT_ID, NODE_ID]["CIGAR"] = f"{LENGTH}M"
                
                for R_link in index[R_core][R]:
                    Rmer = f"{R_core}{___}{R_link}"
                    RIGHT_ID = node_IDs[Rmer]
                    if not G.has_edge(NODE_ID, RIGHT_ID):
                        
                        BSTART = G.nodes[RIGHT_ID]['START']
                        AEND = G.nodes[NODE_ID]['END']
                        if BSTART <= AEND:
                            LENGTH = AEND - BSTART + 1
                        else:
                            ASTART = G.nodes[NODE_ID]['START']
                            ALENGTH = G.nodes[NODE_ID]['LEN']
                            LENGTH = ALENGTH - (BSTART-ASTART)
                        
                        G.add_edge(NODE_ID, RIGHT_ID)
                        G.edges[NODE_ID, RIGHT_ID]["From"] = '+'
                        G.edges[NODE_ID, RIGHT_ID]["To"] = '+'
                        G.edges[NODE_ID, RIGHT_ID]["CIGAR"] = f"{LENGTH}M"
                
                node_IDs[kmer] = NODE_ID; NODE_ID += 1
            G.nodes[node_IDs[kmer]]['LR'].add(read)


for node_id in G.nodes:
    LCUT = 0; RCUT = 0
    in_edges = [attr for _, _, attr in G.in_edges(node_id, data=True)]
    out_edges = [attr for _, _, attr in G.out_edges(node_id, data=True)]
    
    if len(in_edges) > 0:
        LCUT = floor(int(in_edges[0]["CIGAR"][:-1]) / 2)
    if len(out_edges) > 0:
        RCUT = ceil(int(out_edges[0]["CIGAR"][:-1]) / 2)
    
    RCUT = len(G.nodes[node_id]['SEQ']) - RCUT
    G.nodes[node_id]['SEQ'] = G.nodes[node_id]['SEQ'][LCUT:RCUT]
    G.nodes[node_id]['LEN'] = len(G.nodes[node_id]['SEQ'])

for edge_u, edge_v in G.edges():
    G.edges[edge_u, edge_v]["CIGAR"] = '0M'


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
    XAvgCov = XCov / LEN
    XNormCov = XAvgCov / AvgCov
    G.nodes[node_id]['DP'] = XNormCov


with open(argv[2], 'w') as fout:
    fout.write("H\n")

    # Write sequences
    for node_id, attr in G.nodes(data=True):
        SEQ = attr.get('SEQ', '')
        KMER = attr.get('KMER', '')
        ACC = attr.get('ACC', '')
        START = attr.get('START', '')
        END = attr.get('END', '')
        LEN = attr.get('LEN', '')
        DP = attr.get('DP', '')
        CLASS = attr.get('CLASS', '')
        LR = attr.get('LR', set())
        LR = str(list(LR)).replace("'", "").replace(' ', '')
        fout.write(f"S\t{node_id}\t{SEQ}\tGN:Z:{KMER}\tSEG:Z:{ACC}@{START}-{END}\tLN:i:{LEN}\tdp:f:{DP}\tclass:Z:{CLASS}\tLR:Z:{LR}\n")

    # Write links
    for u, v, attr in G.edges(data=True):
        FROM = attr.get('From', '')
        TO = attr.get('To', '')
        CIGAR = attr.get('CIGAR', '')
        fout.write(f"L\t{u}\t{FROM}\t{v}\t{TO}\t{CIGAR}\n")
exit()

DIST = None
def update_dist(value):
    global DIST
    DIST = value
    return True

POS = None; PRIORITY = 0
def update_pos(value, priority=0):
    global POS, PRIORITY
    if priority >= PRIORITY:
        POS = value
        PRIORITY = priority
        return True
    else:
        return False

# Position nodes
plt.ion()
dist = dict(nx.shortest_path_length(G, weight=None))
layout = nx2.draw_circular_layout(G, plt=plt)
layout = nx2.draw_spring_layout(G, pos=layout, plt=plt)
layout = nx2.draw_kamada_kawai_layout(G, dist=dist, pos=layout, plt=plt)
input("done"); exit()

#components = nx.strongly_connected_components(G)
#for component in components:
#    subgraph = nx.induced_subgraph(G, component)
#    layout = nx.spring_layout(subgraph, iterations=100)
#    nx.draw(subgraph, layout, with_labels=False)
#    plt.show()


# Position nodes
components = list(nx.weakly_connected_components(G))
components = sorted(components, key=lambda x: -len(x))
num_components = len(components)

stderr.write(f"Number of connected components: {num_components}"); stderr.flush()
num_pages = floor(cbrt(num_components))
num_plots = ceil(num_components/num_pages)
num_rows = floor(sqrt(num_plots))
num_columns = ceil(num_plots/num_rows)
stderr.write('\n'); stderr.flush()

for i in range(num_pages):
    stderr.write(f"Page: {i+1}/{num_pages}"); stderr.flush()
    fig, axes = plt.subplots(num_rows, num_columns)
    axes = axes.flatten()[:num_plots]

    for j in range(num_plots):
        if i*num_plots+j == num_components: break
        subgraph = nx.induced_subgraph(G, components[i*num_plots+j])
        layout = nx2.draw_kamada_kawai_layout(subgraph)
        nx.draw(subgraph, layout, ax=axes[j], with_labels=False)

    plt.show()
    stderr.write('\r'); stderr.flush()
stderr.write('\n'); stderr.flush()
