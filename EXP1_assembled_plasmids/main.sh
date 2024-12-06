#!/bin/bash
#SBATCH --time=05:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=32

module load StdEnv/2020 gcc/9.3.0 prodigal/2.6.3 diamond/2.0.15 blast+/2.12.0 mummer/4.0.0beta2 hmmer/3.3.2 infernal/1.1.4 python/3.9
PLATON="../../../tools/platon/bin/platon"

mkdir -p kmers
while IFS= read -r line; do
    echo -n "${line} "
    python3 extract_kmers.py -f "data/${line}.fasta" -g "data/${line}.gb" -o "kmers/${line}.fasta" && \
    "${PLATON}" "kmers/${line}.fasta" --output platon >/dev/null && \
    python3 evaluate_results.py -a "${line}" -s "data/${line}.fasta" -f "kmers/${line}.fasta" -p "platon/${line}.plasmid.fasta" -o "summary.tsv" -m "plasmid" && \
    echo
done < accessions.txt
