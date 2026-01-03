#!/usr/bin/env python3

import csv
import os
import argparse
from collections import defaultdict

def combine_microhaplotypes_from_csv(input_csv, output_file):
    """Combines phase-contiguous microhaplotypes from a CSV file and outputs the result."""
    
    if not os.path.exists(input_csv):
        print(f"ERROR: Input CSV file not found: {input_csv}")
        return

    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    haplotypes = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: {
        "hap1": "",
        "hap2": "",
        "reads1": [],
        "depth1": [],
        "reads2": [],
        "depth2": [],
        "variant_positions": []
    })))

    with open(input_csv, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            gene_name = row["Gene"]
            sample_id = row["Sample"]
            block_id = row["Haplotype_Block"]

            # Skip variants not assigned to a block
            if block_id == "Unknown":
                continue
            
            # Skip variants flagged as unknown reason (e.g., "Variant not in GTF")
            if row.get("Reason_Unknown", "").strip():
                continue
            
            genotype = row["Genotype"].split("|")
            if len(genotype) != 2:
                continue  # Skip unphased genotypes
            
            allele1 = row["Allele1"]
            allele2 = row["Allele2"]

            try:
                total_reads1 = int(row["Total_Reads_Allele1"])
                read_depth1 = int(row["Read_Depth_Allele1"])
                total_reads2 = int(row["Total_Reads_Allele2"])
                read_depth2 = int(row["Read_Depth_Allele2"])
                variant_pos = int(row["Variant_Pos"])
            except ValueError:
                continue  # Skip rows with invalid numbers

            haplotype_data = haplotypes[gene_name][sample_id][block_id]

            haplotype_data["hap1"] += allele1
            haplotype_data["hap2"] += allele2
            haplotype_data["reads1"].append(total_reads1)
            haplotype_data["depth1"].append(read_depth1)
            haplotype_data["reads2"].append(total_reads2)
            haplotype_data["depth2"].append(read_depth2)
            haplotype_data["variant_positions"].append(variant_pos)

    with open(output_file, "w", newline="") as csvfile:
        fieldnames = [
            "Gene", "Sample", "Haplotype_Block", 
            "First_SNP_Position", "Last_SNP_Position", "Haplotype_Block_Length",
            "Microhaplotype_1", "Microhaplotype_2", "Microhaplotype_Length",
            "Avg_Reads_Per_SNP_Hap1", "Avg_Read_Depth_Per_SNP_Hap1",
            "Avg_Reads_Per_SNP_Hap2", "Avg_Read_Depth_Per_SNP_Hap2"
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for gene in haplotypes:
            for sample in haplotypes[gene]:
                for block in haplotypes[gene][sample]:
                    hap1 = haplotypes[gene][sample][block]["hap1"]
                    hap2 = haplotypes[gene][sample][block]["hap2"]
                    reads1 = haplotypes[gene][sample][block]["reads1"]
                    depth1 = haplotypes[gene][sample][block]["depth1"]
                    reads2 = haplotypes[gene][sample][block]["reads2"]
                    depth2 = haplotypes[gene][sample][block]["depth2"]
                    variant_positions = haplotypes[gene][sample][block]["variant_positions"]

                    microhap_length = (len(hap1) + len(hap2)) / 2 if hap1 and hap2 else 0

                    if variant_positions:
                        first_snp_pos = min(variant_positions)
                        last_snp_pos = max(variant_positions)
                        haplotype_block_length = last_snp_pos - first_snp_pos if last_snp_pos != first_snp_pos else 1
                    else:
                        first_snp_pos = last_snp_pos = haplotype_block_length = None

                    avg_reads1 = sum(reads1) / len(reads1) if reads1 else 0
                    avg_depth1 = sum(depth1) / len(depth1) if depth1 else 0
                    avg_reads2 = sum(reads2) / len(reads2) if reads2 else 0
                    avg_depth2 = sum(depth2) / len(depth2) if depth2 else 0

                    if round(microhap_length, 2) == 1:
                        continue
                    
                    writer.writerow({
                        "Gene": gene,
                        "Sample": sample,
                        "Haplotype_Block": block,
                        "First_SNP_Position": first_snp_pos,
                        "Last_SNP_Position": last_snp_pos,
                        "Haplotype_Block_Length": round(haplotype_block_length, 2),
                        "Microhaplotype_1": hap1,
                        "Microhaplotype_2": hap2,
                        "Microhaplotype_Length": round(microhap_length, 2),
                        "Avg_Reads_Per_SNP_Hap1": round(avg_reads1, 2),
                        "Avg_Read_Depth_Per_SNP_Hap1": round(avg_depth1, 2),
                        "Avg_Reads_Per_SNP_Hap2": round(avg_reads2, 2),
                        "Avg_Read_Depth_Per_SNP_Hap2": round(avg_depth2, 2)
                    })
    
    print(f"✅ Microhaplotype combination complete. Output saved to: {output_file}")
    print("ℹ️ Microhaplotypes assembled successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine phase-contiguous microhaplotypes from CSV.")
    parser.add_argument("-i", "--input_csv", required=True, help="Path to the input CSV file with phased genotypes.")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output CSV file for combined microhaplotypes.")

    args = parser.parse_args()
    combine_microhaplotypes_from_csv(args.input_csv, args.output_file)

