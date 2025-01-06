# PLASLR: Plasmid Assembly from Long Reads

In this project, we would like to explore a novel method for assembling plasmids from long read sequencing data. The core idea is to identify genes within the long reads and represent each long read as a sequence of gene identifiers. Using these gene identifiers, we partition each long read into overlapping *k*-mers over the gene identifier alphabet, construct a de Bruijn graph from these *k*-mers, and try to find paths or cycles in the graph that might resemble plasmid structures.

To facilitate the assembly, especially given the uncertainty in real-world sequencing settings as to whether a long read originates from a plasmid or the bacterial chromosome, we introduce a preprocessing step where each *k*-mer is classified as either *chromosomal* or *plasmidic*. We want to assess the reliability of this preprocessing step over a series of experiments where we gradually increase the potential for errors and evaluate the effects based on some statistical measures.

## Project Overview

The core of the pipeline is written in **Bash** and **Python3**. Below is a high-level overview of the different steps:

1. **Gene Identification**: Identify genes on long reads and spell each long read as sequence of gene identifiers.
2. ***k*-mer Extraction**: Split each long read into a set of overlapping *k*-mers over the gene identifier alphabet.
3. **Preprocessing**: Classify each *k*-mer as either *chromosomal* or *plasmidic* using [Platon](https://github.com/oschwengers/platon "Platon") to aid in later assembly.
4. **Graph Construction**: Construct a de Bruijn graph from the *k*-mers and perform some basic error correction.
5. **Genome Assembly**: Identify potential plasmid-like structures by finding paths or cycles in the de Bruijn graph.

## Experimental Setup

We evaluate the preprocessing step across four different experimental setups with increasing error conditions:

1. **Experiment 1**: We download perfectly assembled plasmids from the NCBI database and use existing information about the genes provided in the associated GenBank files. The *k*-mers are extracted directly from the error-free plasmid sequence, and in the ideal case we would expect all *k*-mers to be classified as *plasmidic*.
2. **Experiment 2**: We simulate long reads from the assembled plasmid sequences using [Badread](https://github.com/rrwick/Badread "Badread") and look up the read coordinates in the original FASTA and GenBank files. The *k*-mers are extracted again from the original plasmid sequence, but compared to Exp. 1, some *k*-mers might be missing due to uncovered regions or erroneous in the case of chimeric reads.
3. **Experiment 3**: We simulate long reads from the assembled plasmid sequences and map the gene sequences from the original FASTA and GenBank files to the reads. The *k*-mers are extracted this time from the simulated long reads, thus compared to Exp. 2, some k-mers might be affected by sequencing errors and there could be multiple variants of a *k*-mer.
4. **Experiment 4**: We simulate long reads from the assembled plasmid sequences and perform the gene prediction on the long reads ourselves using [Prodigal](https://github.com/hyattpd/Prodigal "Prodigal") instead of relying on the information provided in the GenBank file. This is closer to real-world settings where the annotation is not available for newly sequenced reads, but further increases the potential for errors due to lower accuracy.

## Dependencies

To run the pipeline, make sure that you have the necessary dependencies installed and your environment is set up correctly.
The pipeline uses several tools to carry out its tasks. You can install the tools by following their respective installation guides:

* [**Platon**](https://github.com/oschwengers/platon "Platon"): Used for plasmid classification.
* [**Badread**](https://github.com/rrwick/Badread "Badread"): Used for simulating long reads.
* [**minimap2**](https://github.com/lh3/minimap2 "minimap2"): Used for gene-to-read mapping.
* [**Prodigal**](https://github.com/hyattpd/Prodigal "Prodigal"): Used for gene finding & annotation.

If you have any questions or suggestions for bug fixes or new features, feel free to contact the project maintainers.
