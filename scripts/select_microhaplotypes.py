import os
import sys
import pandas as pd
import argparse

def parse_arguments():
    """Parses command-line arguments for selecting blocks and input/output files."""
    parser = argparse.ArgumentParser(description="Extract gene and block ID based on selection, then retrieve matching microhaplotypes.")
    parser.add_argument("-I", "--Informativeness", action="store_true", help="Select Most_Informative_In column.")
    parser.add_argument("-A", "--Allele", action="store_true", help="Select Highest_Ae column.")
    parser.add_argument("-L", "--Longest_Block", action="store_true", help="Filter by Longest_Block.")
    parser.add_argument("-H", "--Highest_Het", action="store_true", help="Filter by Highest_Het.")
    parser.add_argument("-c", "--combined_stats", required=True, help="Input CSV file with combined gene-block statistics.")
    parser.add_argument("-m", "--microhaplotypes", required=True, help="Input CSV file with microhaplotypes per sample.")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file for filtered microhaplotypes.")
    
    args = parser.parse_args()
    
    # Ensure only one of the filter flags is chosen
    if sum([args.Informativeness, args.Allele, args.Longest_Block, args.Highest_Het]) != 1:
        parser.error("You must select exactly one filter flag: -I, -A, -L, or -H.")
    
    return args

def extract_selected_blocks(combined_stats_file, args):
    """Extracts Gene and the selected Block ID from the combined statistics file."""
    if not os.path.exists(combined_stats_file):
        print(f"ERROR: Combined stats file not found: {combined_stats_file}")
        sys.exit(1)
    
    print(f"DEBUG: Reading combined stats file: {combined_stats_file}")
    df = pd.read_csv(combined_stats_file)

    if args.Informativeness:
        column_to_use = "Most_Informative_In"
    elif args.Allele:
        column_to_use = "Highest_Ae"
    elif args.Longest_Block:
        column_to_use = "Longest_Block"
    elif args.Highest_Het:
        column_to_use = "Highest_Het"

    if column_to_use not in df.columns:
        print(f"ERROR: Column '{column_to_use}' not found in {combined_stats_file}. Available columns: {list(df.columns)}")
        sys.exit(1)

    print(f"DEBUG: Extracting column '{column_to_use}' for block selection")
    return df[['Gene', column_to_use]].rename(columns={column_to_use: 'Selected_Block'})

def filter_microhaplotypes(microhap_file, selected_blocks, output_file):
    """Filters microhaplotypes based on the selected blocks and writes to output CSV."""
    if not os.path.exists(microhap_file):
        print(f"ERROR: Microhaplotypes file not found: {microhap_file}")
        sys.exit(1)

    print(f"DEBUG: Reading microhaplotypes file: {microhap_file}")
    microhap_df = pd.read_csv(microhap_file)

    required_columns = {'Gene', 'Sample', 'Haplotype_Block', 'Microhaplotype_1', 'Microhaplotype_2'}
    missing_columns = required_columns - set(microhap_df.columns)
    
    if missing_columns:
        print(f"ERROR: Missing required columns {missing_columns} in {microhap_file}. Available columns: {list(microhap_df.columns)}")
        sys.exit(1)

    print("DEBUG: Merging selected blocks with microhaplotypes...")
    merged_df = pd.merge(microhap_df, selected_blocks, on="Gene", how="inner")

    print("DEBUG: Filtering matching Haplotype_Blocks...")
    filtered_df = merged_df[merged_df['Haplotype_Block'] == merged_df['Selected_Block']]

    if filtered_df.empty:
        print("WARNING: No matching microhaplotypes found after filtering. Output will be empty.")

    filtered_df = filtered_df[['Gene', 'Sample', 'Haplotype_Block', 'Microhaplotype_1', 'Microhaplotype_2']]
    
    print(f"DEBUG: Saving filtered microhaplotypes to {output_file}")
    filtered_df.to_csv(output_file, index=False)
    print(f"Filtered microhaplotypes saved to {output_file}")

if __name__ == "__main__":
    args = parse_arguments()
    print("DEBUG: Arguments parsed successfully.")

    selected_blocks = extract_selected_blocks(args.combined_stats, args)
    print(f"DEBUG: Selected blocks extracted: {selected_blocks.shape[0]} rows")

    filter_microhaplotypes(args.microhaplotypes, selected_blocks, args.output)
    print("DEBUG: Script execution completed.")
