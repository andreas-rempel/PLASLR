# PLASLR: Plasmid Assembly from Long Reads

In this project, we would like to explore a novel method for assembling plasmids from long read sequencing data. The core idea is to identify genes within the long reads and represent each long read as a sequence of gene identifiers. Using these gene identifiers, we partition each long read into overlapping *k*-mers over the gene identifier alphabet, construct a de Bruijn graph from these *k*-mers, and try to find paths or cycles in the graph that might resemble plasmid structures.

1. Run [Amira](https://github.com/Danderson123/amira) on the long reads using the appropriate species model, e.g. [*Klebsiella pneumoniae*](https://github.com/Danderson123/amira?tab=readme-ov-file#pre-built-species-specific-panrgs):
   ```
   singularity exec amira.img amira
     --reads reads.fastq.gz
     --species Klebsiella_pneumoniae
     --panRG-path Klebsiella.pneumoniae.panidx.zip
     --no-trim
     --debug
     --cores 28
   ```
   This will generate as output (i) the gene calls, and (ii) for each gene the consensus sequence:
   ```
   gunzip amira_output/pandora_output/pandora.consensus.fq.gz
   ```

2. Extract the gene *k*-mers for the long reads from the Amira output and build a de Bruijn graph:
   ```
   python3 extract_kmers.py
     -f amira_output/pandora_output/pandora.consensus.fq
     -g amira_output/corrected_gene_calls_after_filtering.json
     -o graph.gfa
   ```

3. Run [HyPlAss](https://github.com/f0t1h/HyPlAss) on the long reads using the generated de Bruijn graph as a short-read assembly:
   ```
   python3 src/hyplass.py
     --sr-assembly graph.gfa
     --long-reads amira_output/reads.fastq.gz
     --platon-db platon/db
     -o HyPlAss_output
     -t 28
   ```
