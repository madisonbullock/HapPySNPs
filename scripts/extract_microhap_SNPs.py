#!/usr/bin/env python3

import pysam
import argparse
import csv
import os
import subprocess
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Extract microhaplotype SNPs from VCF and BAMs with imputation awareness.")
    parser.add_argument("-v", "--vcf", required=True, help="Input phased VCF file (can be uncompressed or bgzipped)")
    parser.add_argument("-b", "--bam_list", required=True, help="BAM input file list: either CSV (Sample,BAM) or plain list of BAM paths")
    parser.add_argument("-g", "--gtf", required=True, help="GTF file containing block_id annotations")
    parser.add_argument("-i", "--imputation_log", required=True, help="CSV log of imputed genotypes")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
    parser.add_argument("-d", "--min_depth", type=int, default=6, help="Minimum depth to use for imputed SNPs (default: 6)")
    parser.add_argument("--bam_list_is_csv", action="store_true",
                        help="Specify if BAM list (-b) is a CSV with Sample,BAM columns. If not set, assumes plain BAM path list with sampleID extracted from filename.")
    return parser.parse_args()

def bgzip_and_index(vcf_path):
    """Compress and index a VCF if not already compressed."""
    if not vcf_path.endswith(".gz"):
        print(f"Compressing and indexing {vcf_path} ...")
        compressed_vcf = vcf_path + ".gz"
        with open(compressed_vcf, "wb") as f_out:
            subprocess.run(["bgzip", "-c", vcf_path], stdout=f_out, check=True)
        subprocess.run(["tabix", "-p", "vcf", compressed_vcf], check=True)
        return compressed_vcf
    else:
        index_path = vcf_path + ".tbi"
        if not os.path.exists(index_path):
            print(f"Index not found. Indexing {vcf_path} ...")
            subprocess.run(["tabix", "-p", "vcf", vcf_path], check=True)
        return vcf_path

def load_blocks_from_gtf(gtf_file):
    blocks = defaultdict(list)
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"): continue
            fields = line.strip().split("\t")
            chrom, _, feature, start, end, _, _, _, attributes = fields
            if feature != "block":
                continue
            start, end = int(start), int(end)
            block_id = attributes.split('block_id "')[1].split('"')[0]
            for pos in range(start, end + 1):
                blocks[chrom].append((pos, block_id))
    pos_to_block = defaultdict(str)
    for chrom in blocks:
        for pos, block_id in blocks[chrom]:
            pos_to_block[(chrom, pos)] = block_id
    return pos_to_block

def load_imputed_positions(log_file):
    imputed = set()
    with open(log_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample = row["Sample"]
            chrom_pos = row["SNP"]
            if ":" in chrom_pos:
                chrom, pos = chrom_pos.rsplit(":", 1)
                try:
                    pos = int(pos)
                except ValueError:
                    continue
                imputed.add((sample, chrom, pos))
    return imputed

def load_bam_paths_from_csv(csv_file):
    bam_paths = {}
    with open(csv_file) as f:
        reader = csv.reader(f)
        header = next(reader)
        sample_col = 0
        bam_col = 1
        if "Sample" in header and "BAM" in header:
            sample_col = header.index("Sample")
            bam_col = header.index("BAM")
        for row in reader:
            if len(row) < 2:
                continue
            sample, bam = row[sample_col], row[bam_col]
            bam_paths[sample] = bam
    return bam_paths

def load_bam_paths_from_list(bam_list_file):
    bam_paths = {}
    with open(bam_list_file) as f:
        for line in f:
            bam_path = line.strip()
            if not bam_path:
                continue
            filename = os.path.basename(bam_path)
            sample_id = filename.split('.')[0]
            bam_paths[sample_id] = bam_path
    return bam_paths

def fetch_alleles(vcf, record, sample):
    gt = record.samples[sample].get("GT")
    if gt is None or "." in str(gt):
        return None, None, "./."
    alleles = [record.alleles[i] if i < len(record.alleles) else "N" for i in gt]
    genotype_str = "|".join(map(str, gt))
    return alleles[0], alleles[1], genotype_str

def count_reads(bam_path, chrom, pos, alleles):
    samfile = pysam.AlignmentFile(bam_path, "rb")
    allele_counts = {alleles[0]: 0, alleles[1]: 0}
    for pileupcolumn in samfile.pileup(chrom, pos - 1, pos, truncate=True):
        if pileupcolumn.pos != pos - 1:
            continue
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del or pileupread.is_refskip:
                continue
            base = pileupread.alignment.query_sequence[pileupread.query_position]
            if base in allele_counts:
                allele_counts[base] += 1
    samfile.close()
    return allele_counts

def main():
    args = parse_args()
    vcf_path = bgzip_and_index(args.vcf)
    vcf = pysam.VariantFile(vcf_path)
    pos_to_block = load_blocks_from_gtf(args.gtf)
    imputed = load_imputed_positions(args.imputation_log)

    if args.bam_list_is_csv:
        bam_paths = load_bam_paths_from_csv(args.bam_list)
    else:
        bam_paths = load_bam_paths_from_list(args.bam_list)

    with open(args.output, "w", newline="") as out_f:
        writer = csv.writer(out_f)
        header = [
            "Gene", "Variant_Pos", "Sample", "Genotype", "Allele1", "Allele2",
            "Allele_Balance", "Total_Reads_Allele1", "Read_Depth_Allele1",
            "Total_Reads_Allele2", "Read_Depth_Allele2", "Haplotype_Block", "Imputed"
        ]
        writer.writerow(header)

        for record in vcf.fetch():
            chrom, pos = record.chrom, record.pos
            block_id = pos_to_block.get((chrom, pos), "NA")
            for sample in record.samples:
                if sample not in bam_paths:
                    continue
                allele1, allele2, genotype_str = fetch_alleles(vcf, record, sample)
                if allele1 is None:
                    continue
                imputed_flag = "Yes" if (sample, chrom, pos) in imputed else "No"

                if imputed_flag == "Yes":
                    ab = 0.0
                    dp = args.min_depth
                    writer.writerow([
                        chrom, pos, sample, genotype_str, allele1, allele2, ab,
                        dp, dp, dp, dp, block_id, imputed_flag
                    ])
                else:
                    allele_counts = count_reads(bam_paths[sample], chrom, pos, (allele1, allele2))
                    total1 = allele_counts.get(allele1, 0)
                    total2 = allele_counts.get(allele2, 0)
                    total_reads = total1 + total2
                    ab = round(total2 / total_reads, 3) if total_reads > 0 else 0
                    writer.writerow([
                        chrom, pos, sample, genotype_str, allele1, allele2, ab,
                        total1, total1, total2, total2, block_id, imputed_flag
                    ])

if __name__ == "__main__":
    main()
