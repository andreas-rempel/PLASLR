#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Updating the summary files by adding BioSample, BioProject and Read accession ID

Usage: python update_summary.py <plasmids_accession_file> <summary_file> <data_dir> <updated_summary_file>
       plamids_accession_file: file containing the plasmids accession IDs
       summary_file: initial TSV summary file
       data_dir: directory wih GenBank files, one per plasmid
       updated_summary_file: generated summary file with additional columns recording BioSampe, BioProject, Reads

"""

"""
@Author: Cedric Chauve
@Email:  cedric.chauve@sfu.ca
@Date:   2024-12-20
"""

import sys
import os
import csv
from Bio import GenBank

def read_accession_file(in_accession_file):
    """
    Read an accession IDs file

    input:
    - in_accession_file: text file with one accession ID per line

    output:
    - list(st): list of accession IDs    
    """
    plasmids_list = []
    with open(in_accession_file) as in_file:
        for line in in_file:
            plasmids_list.append(line.rstrip())
    return plasmids_list

def read_dblink_from_genbank_file(in_gb_file):
    """
    Read a GenBank file and extracts the DBLINK entries

    input:
    - in_gb_file: text file with one accession ID per line

    output:
    - dict(str,str): database -> accession
    """
    dblink_dict = {}
    with open(in_gb_file) as handle:
        plasmid_record = GenBank.read(handle)
        for dblink in plasmid_record.dblinks:
            key = dblink.split(": ")[0]
            val = dblink.split(": ")[1]
            dblink_dict[key] = val
    return dblink_dict


def update_summary_file(in_summary_file, in_dblink_dict, out_summary_file):
    """
    Read a TSV summary file and augments it with sample information

    input:
    - in_summary_file: input summary TSV file
      ASSUMPTION: contains a field plasmid_ID with value corresponding to a GenBank ACCESSION key
    - in_dblink_dict: dict(str->(dict(str,str))
      key = plasmid ID, val= (key = DBLINK key (dataas), val = DBLINK accession)
    - out_summary_file: output TSV summary file

    output: None
    """
    # Adding all keys of in_dblink_dict entries
    summary_columns = [
        "plasmid_ID",
        "length_(bp)",
        "#5-mers",
        "#plasmid_5-mers",
        "largest_non-plasmid_island_(#5-mers)"
    ]
    additional_columns = []
    for dblink_dict in in_dblink_dict.values():
        for key in dblink_dict.keys():
            if key not in additional_columns:
                additional_columns.append(key)
    summary_columns += additional_columns
    # Adding DBLINK information to summary file
    with open(in_summary_file, newline="") as in_file, open(out_summary_file, "w") as out_file:
        reader = csv.DictReader(in_file, delimiter = "\t")
        writer = csv.DictWriter(out_file, fieldnames=summary_columns, delimiter="\t")
        writer.writeheader()
        for plasmid_row in reader:
            for key in additional_columns:
                plasmid_row[key] = "NA"
            plasmid_ID = plasmid_row["plasmid_ID"]
            if plasmid_ID in in_dblink_dict.keys():
                for key,val in in_dblink_dict[plasmid_ID].items():
                    plasmid_row[key] = val
            writer.writerow(plasmid_row)

def main(in_plasmids_accession_file, in_summary_file, in_data_dir, out_summary_file):
    plasmids_list = read_accession_file(in_plasmids_accession_file)
    plasmids_dblink_dict = {}
    for plasmid in plasmids_list:
        gb_file = os.path.join(in_data_dir,f"{plasmid}.gb")
        plasmids_dblink_dict[plasmid] = read_dblink_from_genbank_file(gb_file)
    update_summary_file(in_summary_file, plasmids_dblink_dict, out_summary_file)



if __name__ == '__main__':
    in_plasmids_accession_file = sys.argv[1]
    in_summary_file = sys.argv[2]
    in_data_dir = sys.argv[3]
    out_summary_file = sys.argv[4]

    main(in_plasmids_accession_file, in_summary_file, in_data_dir, out_summary_file)
    
    
