#!/usr/bin/env python3

import pysam
import argparse
import os
import sys

def write_gtf(blocks, output_filename):
    """Writes detected phase blocks to a GTF file."""
    with open(output_filename, 'w') as gtf_file:
        for block in blocks:
            gtf_file.write(
                f"{block['chromosome']}\tWhatshap\tblock\t"
                f"{block['start']}\t{block['end']}\t.\t+\t.\t"
                f"block_id \"{block['block_id']}\";\n"
            )

def detect_phase_blocks(vcf_filename, output_filename):
    """Parses the VCF and detects phase blocks, then outputs a GTF file."""
    if not os.path.exists(vcf_filename):
        print(f"ERROR: VCF file not found: {vcf_filename}")
        sys.exit(1)

    vcf_reader = pysam.VariantFile(vcf_filename)

    blocks = []
    block_id = 1
    last_phase = None
    start_pos = None
    current_chromosome = None
    last_pos = None

    for record in vcf_reader:
        if not record.samples:
            continue

        sample = list(record.samples.keys())[0]
        phase = record.samples[sample].get("GT")

        if phase is None:
            continue

        new_chromosome = record.chrom

        if new_chromosome != current_chromosome:
            if start_pos is not None:
                blocks.append({
                    'chromosome': current_chromosome,
                    'start': start_pos,
                    'end': last_pos,
                    'block_id': f"block{block_id}"
                })
                block_id += 1
            current_chromosome = new_chromosome
            start_pos = record.pos
            last_pos = record.pos
            last_phase = phase
            continue

        if phase != last_phase:
            blocks.append({
                'chromosome': current_chromosome,
                'start': start_pos,
                'end': last_pos,
                'block_id': f"block{block_id}"
            })
            block_id += 1
            start_pos = record.pos

        last_phase = phase
        last_pos = record.pos

    if start_pos is not None:
        blocks.append({
            'chromosome': current_chromosome,
            'start': start_pos,
            'end': last_pos,
            'block_id': f"block{block_id}"
        })

    write_gtf(blocks, output_filename)
    print(f"✅ GTF file successfully written: {output_filename}")

def main():
    parser = argparse.ArgumentParser(
        description="Detect phase blocks from a phased VCF and output a GTF file."
    )
    parser.add_argument("-i", "--input_vcf", required=True, help="Path to the input phased VCF file.")
    parser.add_argument("-o", "--output_gtf", required=True, help="Path to the output GTF file.")

    args = parser.parse_args()
    detect_phase_blocks(args.input_vcf, args.output_gtf)

if __name__ == "__main__":
    main()

