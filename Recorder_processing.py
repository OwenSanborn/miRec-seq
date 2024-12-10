import os
import matplotlib.pyplot as plt
import re
import pandas as pd
import pyfastx
import csv
from collections import Counter

# Function will extract the 3 barcodes and UMI from the reverse read. A vector is returned containing the read ID and each barcode.
def process_rec_BC(seq, read_id):
    BC1_pattern = re.compile(r'ATCCACGTGCTTGAGACTGTGG(.{8})')
    BC2_pattern = re.compile(r'(.{8})ATCCACGTGCTTGAGACTGTGG')
    BC3_pattern = re.compile(r'(.{8})GTGGCCGATGTTTCGCATCGGCGTACGACT')
    UMI = seq[1:10]

    match_BC1 = BC1_pattern.search(seq)
    BC1 = match_BC1.group(1) if match_BC1 else None

    match_BC2 = BC2_pattern.search(seq)
    BC2 = match_BC2.group(1) if match_BC2 else None
    
    match_BC3 = BC3_pattern.search(seq)
    BC3 = match_BC3.group(1) if match_BC3 else None

    return {
        'read_id':read_id,
        'BC1': BC1,
        'BC2': BC2,
        'BC3': BC3,
        'UMI': UMI
    }

# Function will extract recorder motif and hairpin from the forward read. The hairpin is verified as being 34 positions long and the number of edits is returned
def process_recorder(seq, read_id):
    motif_pattern = re.compile(r'CTATTCTGGCTG(.*?)TCCAACGCAAT')
    hp_pattern = re.compile(r'TTAAATT(.{34})')

    match_motif = motif_pattern.search(seq)
    motif = match_motif.group(1) if match_motif else None

    match_hp = hp_pattern.search(seq)
    hp = match_hp.group(1) if match_hp else None

    hp_length = len(hp) if hp else None
    edits_count = (hp.count('G')-7) if hp and hp_length == 34 else None

    return {
        'read_id':read_id,
        'motif': motif,
        'hairpin': hp,
        'hp_length': hp_length,
        'edits_count': edits_count
    }


# Function for iterating over a fastq file and extracting the motif, hairpin, and barcodes by calling process_recorder and process_recorder_BC.
def extract_paired_rec_reads(input_r1, input_r2):
    id = input_r1[:5]  
    rec = "CCAATCCAATCC"
    bc_adapter = "CGTGCTTGAG"

    # Define paths for the input FASTQ files (forward and reverse)
    input_path_r1 = os.path.join("/Users/osanborn/Documents/0.Obsidian Vault/Projects/01 - miRec/Data/sc-miRec-seq_03Oct24/raw_data", input_r1)
    input_path_r2 = os.path.join("/Users/osanborn/Documents/0.Obsidian Vault/Projects/01 - miRec/Data/sc-miRec-seq_03Oct24/raw_data", input_r2)
    
    # Output file paths for CSV and FASTQ (output both R1 and R2 reads)
    output_csv_1 = os.path.join("/Users/osanborn/Documents/0.Obsidian Vault/Projects/01 - miRec/Data/sc-miRec-seq_03Oct24/Rec_Analysis/rec_reads", f'{id}_rec_data_1.csv')
    output_csv_2 = os.path.join("/Users/osanborn/Documents/0.Obsidian Vault/Projects/01 - miRec/Data/sc-miRec-seq_03Oct24/Rec_Analysis/rec_reads", f'{id}_rec_data_2.csv')


    # Define CSV headers
    csv_headers_1 = ['read_id', 'motif', 'hairpin', 'hp_length', 'edits_count']

    # Check if the input files exist
    if not os.path.exists(input_path_r1) or not os.path.exists(input_path_r2):
        raise FileNotFoundError("Input FASTQ files do not exist")
    
    # Initialize total reads counter and list to store read information
    total_reads = 0

    # Open the two FASTQ files using pyfastx
    fq_r1 = pyfastx.Fastx(input_path_r1)
    fq_r2 = pyfastx.Fastx(input_path_r2)

    # Open the output files (CSV and FASTQ for both directions)
    with open(output_csv_1, 'w', newline='') as out_csv:
        # Initialize CSV writer
        csv_writer = csv.DictWriter(out_csv, fieldnames=csv_headers_1)
        csv_writer.writeheader()

        # Iterate over paired reads from both R1 and R2 simultaneously
        for name, seq, qual in fq_r1:
            total_reads += 1
            if re.search(rec, seq):
                rec_data_r1 = process_recorder(seq, name)
                combined_rec_data = {
                    'read_id': rec_data_r1['read_id'],
                    'motif': rec_data_r1['motif'],
                    'hairpin': rec_data_r1['hairpin'],
                    'hp_length': rec_data_r1['hp_length'],
                    'edits_count': rec_data_r1['edits_count'],
                }
                csv_writer.writerow(combined_rec_data)


    # Open the output files (CSV and FASTQ for both directions)
    csv_headers_2 = ['read_id', 'BC1', 'BC2', 'BC3', 'UMI']
    with open(output_csv_2, 'w', newline='') as out_csv:
        # Initialize CSV writer
        csv_writer_2 = csv.DictWriter(out_csv, fieldnames=csv_headers_2)
        csv_writer_2.writeheader()

        # Iterate over paired reads from both R1 and R2 simultaneously
        for name, seq, qual in fq_r2:

            if re.search(bc_adapter, seq):
                rec_data_r2 = process_rec_BC(seq, name)
                combined_bc_data = {
                    'read_id': rec_data_r2['read_id'],
                    'BC1': rec_data_r2['BC1'],
                    'BC2': rec_data_r2['BC2'],
                    'BC3': rec_data_r2['BC3'],
                    'UMI': rec_data_r2['UMI']
                }

                csv_writer_2.writerow(combined_bc_data)

    # Merge to one CSV
    df1 = pd.read_csv(output_csv_1)  # This is the CSV for forward reads
    df2 = pd.read_csv(output_csv_2)  # This is the CSV for reverse reads

    merged_df = pd.merge(df1, df2, on='read_id', how='inner')

    # Save the merged DataFrame to a new CSV file
    output_merged_csv = os.path.join("/Users/osanborn/Documents/0.Obsidian Vault/Projects/01 - miRec/Data/sc-miRec-seq_03Oct24/Rec_Analysis/rec_reads", f'{id}_rec_data_combined.csv')
    merged_df.to_csv(output_merged_csv, index=False)
