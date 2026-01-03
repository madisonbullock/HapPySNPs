import os
import sys
import pandas as pd
import argparse
import re

def parse_arguments():
    """Parses command-line arguments."""
    parser = argparse.ArgumentParser(description="Assign numeric allele IDs and generate dictionary file.")
    parser.add_argument("-i", "--input", required=True, help="Input CSV file with filtered microhaplotypes.")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file with assigned allele IDs.")
    parser.add_argument("-d", "--dictionary", required=True, help="Dictionary output CSV file mapping Gene, Block, Microhaplotype, Allele_ID.")
    
    return parser.parse_args()

def extract_gene_id(gene):
    """Extracts the last 4-digit number from the gene name."""
    match = re.findall(r'\d{4}', gene)  # Find all 4-digit numbers
    return match[-1] if match else "0000"  # Take the last one

def extract_block_id(block):
    """Extracts the numeric part from the block name."""
    match = re.search(r'\d+', block)  # Find all numbers
    return match.group() if match else "0000"

def assign_allele_ids(input_file, output_file, dictionary_file):
    """Assigns numeric allele IDs and generates the final output file and dictionary."""
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        sys.exit(1)

    print(f"DEBUG: Reading input file: {input_file}")
    df = pd.read_csv(input_file)

    required_columns = {'Gene', 'Sample', 'Haplotype_Block', 'Microhaplotype_1', 'Microhaplotype_2'}
    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        print(f"ERROR: Missing required columns {missing_columns} in {input_file}. Available columns: {list(df.columns)}")
        sys.exit(1)

    print("DEBUG: Extracting unique alleles for ID assignment...")
    
    # Dictionary to store unique allele IDs per gene and block
    allele_dict = {}
    dictionary_data = []

    # Generate unique numeric allele IDs per gene and block
    for _, row in df.iterrows():
        gene = row['Gene']
        block = row['Haplotype_Block']
        microhap_1 = row['Microhaplotype_1']
        microhap_2 = row['Microhaplotype_2']

        # Extract numeric gene ID and block ID
        gene_id = extract_gene_id(gene)
        block_id = extract_block_id(block)

        # Initialize allele counter for each gene if not already initialized
        if gene not in allele_dict:
            allele_dict[gene] = {}

        for microhap in [microhap_1, microhap_2]:
            if microhap not in allele_dict[gene]:
                # Assign a new numeric allele ID (1-based numbering)
                allele_dict[gene][microhap] = len(allele_dict[gene]) + 1

                # Add the dictionary entry with padded allele IDs
                dictionary_data.append([gene, block, microhap, str(allele_dict[gene][microhap]).zfill(3)])

    print(f"DEBUG: Assigned {sum(len(d) for d in allele_dict.values())} unique allele IDs.")

    # Add allele ID columns to the original dataframe based on the new numeric IDs with padding
    df['allele_ID_1'] = df.apply(lambda row: str(allele_dict[row['Gene']][row['Microhaplotype_1']]).zfill(3), axis=1)
    df['allele_ID_2'] = df.apply(lambda row: str(allele_dict[row['Gene']][row['Microhaplotype_2']]).zfill(3), axis=1)

    print(f"DEBUG: Saving final output to {output_file}")
    df.to_csv(output_file, index=False)
    
    # Write dictionary file
    print(f"DEBUG: Saving dictionary to {dictionary_file}")
    dictionary_df = pd.DataFrame(dictionary_data, columns=['Gene', 'Block', 'Microhaplotype', 'Allele_ID'])
    dictionary_df.to_csv(dictionary_file, index=False)

    print(f"Output saved to {output_file}")
    print(f"Dictionary saved to {dictionary_file}")

if __name__ == "__main__":
    args = parse_arguments()
    assign_allele_ids(args.input, args.output, args.dictionary)
