import os
import re
import pandas as pd
import pyfastx
import os
import csv
import pandas as pd

# Function will extract the 3 barcodes and UMI from the reverse read. A vector is returned containing the read ID and each barcode.
def process_rec_BC(seq, read_id):
    BC1_pattern = re.compile(r'AGACTGTGG(.{8})')
    BC2_pattern = re.compile(r'(.{8})ATCCACGT')
    BC3_pattern = re.compile(r'(.{8})GTGGCCGA')
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

# Handles each read to extract recorder or hairpin components
def process_recorder(seq, read_id):
    motif_pattern = re.compile(r'CTATTCTGGCTG(.*?)TCCAACGCAAT')
    #hp_pattern = re.compile(r'TTAAATT(.*?)AGGCCTG')
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


def extract_paired_rec_reads(raw_path, reads_out_path, input_r1, input_r2):
    id = input_r1[:5]  
    rec = "CCAATCCAATCC"
    bc_adapter = "CGTGCTTGAG"

    # Define paths for the input FASTQ files (forward and reverse)
    input_path_r1 = os.path.join(raw_path, input_r1)
    input_path_r2 = os.path.join(raw_path, input_r2)
    
    # Output file paths for CSV and FASTQ (output both R1 and R2 reads)
    output_csv_1 = os.path.join(reads_out_path, f'{id}_rec_data_1.csv')
    output_csv_2 = os.path.join(reads_out_path, f'{id}_rec_data_2.csv')


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
    output_merged_csv = os.path.join(reads_out_path, f'{id}_rec_data_combined.csv')
    merged_df.to_csv(output_merged_csv, index=False)

# Function will map the barcode sequences to wells to provide the cell ID
def map_barcodes_to_wells(bc1, bc2, bc3, bc1_dict, bc2_bc3_dict):

    bc1_well = bc1_dict.get(bc1, None)
    bc2_well = bc2_bc3_dict.get(bc2, None)
    bc3_well = bc2_bc3_dict.get(bc3, None)

    # Return the formatted barcode string in {BC1_well}_{BC2_well}_{BC3_well} format
    if bc1_well and bc2_well and bc3_well:
        return f"{bc1_well}_{bc2_well}_{bc3_well}"
    else:
        return None

# Processes reads with barcodes
def process_reads_with_barcodes(read_data, bc1_data_path, bc2_bc3_data_path):

    # Load the barcode data
    bc1_data = pd.read_csv(bc1_data_path)  # For BC1
    bc2_bc3_data = pd.read_csv(bc2_bc3_data_path)  # For BC2 and BC3

    # Create dictionaries for fast lookup
    bc1_dict = bc1_data.set_index('sequence')['well'].to_dict()
    bc2_bc3_dict = bc2_bc3_data.set_index('sequence')['well'].to_dict()

    # Ensure read_data is a DataFrame
    if not isinstance(read_data, pd.DataFrame):
        raise ValueError("read_data should be a pandas DataFrame.")

    # Apply the mapping function to each row using DataFrame.apply for vectorized operation
    read_data['well_barcode'] = read_data.apply(
        lambda row: map_barcodes_to_wells(
            row['BC1'], row['BC2'], row['BC3'], bc1_dict, bc2_bc3_dict
        ),
        axis=1
    )
    return read_data

# Collapses reads by UMI
def collapse_umi(read_data):
    # Use only required columns for grouping and aggregation
    grouping_columns = ['well_barcode', 'UMI']

    # Keep the first occurrence of all columns except 'read_id'
    agg_dict = {col: 'first' for col in read_data.columns if col not in ['read_id', 'UMI', 'well_barcode']}
    agg_dict['UMI'] = 'size'  # Count the occurrences of UMI

    # Perform groupby and aggregation
    collapsed_data = read_data.groupby(grouping_columns, as_index=False).agg(agg_dict)

    # Rename the 'UMI' column to 'read_count'
    collapsed_data.rename(columns={'UMI': 'read_count'}, inplace=True)

    return collapsed_data

# Cleans data to remove NA's and edit counts not between 0 and 7
def remove_rows_with_missing_values(df):
    cleaned_df = df.replace("", pd.NA).dropna()
    cleaned_df = cleaned_df[cleaned_df['edits_count'] >= 0]
    cleaned_df = cleaned_df[cleaned_df['edits_count'] < 8]
    return cleaned_df

# Map miRNA name to a motif
def find_miRNA(motif, motifs_df):
    # Check if any motif in the motifs dataframe appears as a substring within the given motif
    matching_miRNAs = motifs_df[motifs_df['Motif'].apply(lambda x: x in motif)]['miRNA'].unique()
    if len(matching_miRNAs) > 0:
        return matching_miRNAs
    else:
        return ["no match"]
