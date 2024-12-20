# EXP1_assembled_plasmids


### Augmenting summary files

Adding in summary files new fields recording BioProject, BioSampl, SRA accessions.

```
python update_summary.py accessions.txt EXP1a_genes-only/summary.tsv data/ EXP1a_genes-only/summary_2.tsv
python update_summary.py accessions.txt EXP1b_genes+gaps/summary.tsv data/ EXP1b_genes+gaps/summary_2.tsv
python update_summary.py accessions.txt EXP1c_genes-only__no-filters/summary.tsv data/ EXP1c_genes-only__no-filters/summary_2.tsv
python update_summary.py accessions.txt EXP1d_genes+gaps__no-filters/summary.tsv data/ EXP1d_genes+gaps__no-filters/summary_2.tsv
```