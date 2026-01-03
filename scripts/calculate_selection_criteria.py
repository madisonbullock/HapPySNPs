#!/usr/bin/env python3
import os
import argparse
import pandas as pd
import numpy as np

def parse_arguments():
    """Parses command-line arguments for the combined microhaplotype input file and output summary file."""
    parser = argparse.ArgumentParser(description="Compute combined gene and block statistics summary from microhaplotype data.")
    parser.add_argument("-i", "--input", required=True, help="Input CSV file from the microhaplotype combination script.")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file for the gene summary.")
    return parser.parse_args()

def read_haplotype_blocks(file_path):
    """Reads haplotype block data from a file and returns a dataframe."""
    return pd.read_csv(file_path)

def calculate_effective_allele_number(freqs):
    """Calculates Effective Allele Number (Ae) given a list of frequencies."""
    s = sum(f**2 for f in freqs)
    return 1 / s if s > 0 else 0

def calculate_rosenberg_informativeness(freqs):
    """Calculates Rosenberg's Measure of Informativeness (In)."""
    n = len(freqs)  # Number of unique alleles
    return 1 - sum(f**2 for f in freqs) - (1 / n) if n > 1 else 0

def compute_block_statistics(df):
    """Computes block-level statistics and returns a dataframe with one row per block."""
    # Create a combined microhaplotype identifier
    df['combined_microhaplotypes'] = df['Microhaplotype_1'].astype(str) + "_" + df['Microhaplotype_2'].astype(str)
    
    # Calculate heterozygosity: proportion of rows where Microhaplotype_1 != Microhaplotype_2
    df['heterozygous'] = df['Microhaplotype_1'] != df['Microhaplotype_2']
    het = df.groupby(['Gene', 'Haplotype_Block'])['heterozygous'].mean().reset_index()
    het.columns = ['Gene', 'Haplotype_Block', 'Heterozygosity']
    
    # Calculate Effective Allele Number (Ae) and Informativeness (In) per block
    ae_in_list = []
    for (gene, block), group in df.groupby(['Gene', 'Haplotype_Block']):
        hap_counts = group['combined_microhaplotypes'].value_counts(normalize=True)
        Ae = calculate_effective_allele_number(hap_counts.values)
        In = calculate_rosenberg_informativeness(hap_counts.values)
        ae_in_list.append([gene, block, Ae, In])
    ae_in_df = pd.DataFrame(ae_in_list, columns=['Gene', 'Haplotype_Block', 'Effective_Allele_Number', 'Informativeness'])
    
    # Get average block length per block (assuming block length is consistent within each block)
    block_lengths = df.groupby(['Gene', 'Haplotype_Block'])['Haplotype_Block_Length'].mean().reset_index()
    
    # Merge block-level stats
    block_stats = pd.merge(het, ae_in_df, on=['Gene', 'Haplotype_Block'], how='outer')
    block_stats = pd.merge(block_stats, block_lengths, on=['Gene', 'Haplotype_Block'], how='outer')
    
    return block_stats

def compute_gene_summary(block_stats):
    """Creates a summary per gene:
       - Longest_Block: block with maximum Haplotype_Block_Length
       - Highest_Het: block with maximum Heterozygosity
       - Highest_Ae: block with maximum Effective_Allele_Number
       - Most_Informative_In: block with maximum Informativeness
       - Block_Comparison_Total: 'Same' if the above four blocks are identical, else 'Different'
       - Ae_vs_In: 'Same' if Highest_Ae block equals Most_Informative_In block, else 'Different'
    """
    summaries = []
    for gene, group in block_stats.groupby("Gene"):
        # Identify blocks by the metric
        longest = group.loc[group['Haplotype_Block_Length'].idxmax()]['Haplotype_Block']
        highest_het = group.loc[group['Heterozygosity'].idxmax()]['Haplotype_Block']
        highest_ae = group.loc[group['Effective_Allele_Number'].idxmax()]['Haplotype_Block']
        most_informative = group.loc[group['Informativeness'].idxmax()]['Haplotype_Block']
        
        # Compare the four block IDs
        blocks = [longest, highest_het, highest_ae, most_informative]
        block_comparison_total = "Same" if all(b == blocks[0] for b in blocks) else "Different"
        ae_vs_in = "Same" if highest_ae == most_informative else "Different"
        
        summaries.append({
            "Gene": gene,
            "Longest_Block": longest,
            "Highest_Het": highest_het,
            "Highest_Ae": highest_ae,
            "Most_Informative_In": most_informative,
            "Block_Comparison_Total": block_comparison_total,
            "Ae_vs_In": ae_vs_in
        })
    return pd.DataFrame(summaries)

if __name__ == "__main__":
    args = parse_arguments()
    
    # Read the combined microhaplotype data
    df = read_haplotype_blocks(args.input)
    
    # Compute block-level statistics
    block_stats = compute_block_statistics(df)
    
    # Compute per-gene summary from block stats
    gene_summary = compute_gene_summary(block_stats)
    
    # Write the summary CSV to the specified output file
    gene_summary.to_csv(args.output, index=False)
    print(f"Summary statistics saved to: {args.output}")

