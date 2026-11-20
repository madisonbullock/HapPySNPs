#!/usr/bin/env python3
"""
Summary:
This script generates consensus phase blocks from a multi-sample, WhatsHap-phased VCF.

Description:
Within each CHROM/gene, this script scans variants in position order. A new
consensus block begins whenever adding a SNP would introduce a different
nonmissing PS value for any sample already represented in the current block.

Only complete, phased, heterozygous diploid genotypes with nonmissing
PS values contribute phase information. Homozygous genotypes, missing
genotypes, unphased genotypes, and records without PS do not independently
create consensus block boundaries.

The resulting consensus blocks are written in GTF format and can
optionally be accompanied by a CSV summary containing block coordinates,
SNP positions, phase-supporting samples, and sample-specific PS values.

General Usage:
python3 01_determine_phase_blocks.py -i phased_variants.vcf.gz -o consensus_blocks.gtf --block-summary consensus_blocks_summary.csv
"""
## Import required packages
from __future__ import annotations

import argparse
import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple
import pandas as pd
import pysam

## Define the data structure used to store the variants and consensus phase blocks
@dataclass
class VariantInfo:
    chrom: str
    pos: int
    sample_ps: Dict[str, Optional[str]]


@dataclass
class ConsensusBlock:
    chrom: str
    positions: List[int]
    supporting_samples: Set[str]
    ps_values_by_sample: Dict[str, Set[str]]

    @property
    def start(self) -> int:
        return self.positions[0]

    @property
    def end(self) -> int:
        return self.positions[-1]

## Open output files
def open_output(path: Path, mode: str = "wt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return path.open(mode, encoding="utf-8")

## Parse arguments
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create consensus phase blocks from WhatsHap FORMAT/PS tags."
    )
    parser.add_argument(
        "-i", "--input-vcf", "--input_vcf",
        dest="input_vcf", required=True, type=Path,
        help="WhatsHap-phased VCF or VCF.GZ",
    )
    parser.add_argument(
        "-o", "--output-gtf", "--output_gtf",
        dest="output_gtf", required=True, type=Path,
        help="Output consensus-block GTF",
    )
    parser.add_argument(
        "--block-summary", type=Path, default=None,
        help="Optional per-block summary CSV",
    )
    return parser.parse_args()

## Normalize PS values and identify missing values
def normalize_ps(value) -> Optional[str]:
    if value is None:
        return None
    if isinstance(value, tuple):
        if not value:
            return None
        value = ",".join(str(item) for item in value)
    value = str(value)
    return None if value in {"", ".", "None"} else value

## Extract PS information only from complete, phased, heterozygous genotypes
def phased_heterozygote_ps(sample_data) -> Optional[str]:
    gt = sample_data.get("GT")
    if gt is None or len(gt) != 2 or any(a is None for a in gt):
        return None
    if not sample_data.phased or gt[0] == gt[1]:
        return None
    return normalize_ps(sample_data.get("PS"))

## Read input VCF and extract sample-specific PS information
def read_variants(
    vcf_path: Path,
) -> Tuple[List[str], Dict[str, List[VariantInfo]], int]:
    if not vcf_path.exists():
        raise FileNotFoundError(f"VCF not found: {vcf_path}")

    vcf = pysam.VariantFile(str(vcf_path))
    samples = list(vcf.header.samples)

    if not samples:
        raise ValueError("VCF contains no sample columns")
    if "PS" not in vcf.header.formats:
        raise ValueError(
            "VCF header does not define FORMAT/PS. Run WhatsHap with --tag=PS."
        )

    by_chrom: Dict[str, List[VariantInfo]] = {}
    total_records = 0

    for record in vcf:
        sample_ps = {
            sample: phased_heterozygote_ps(record.samples[sample])
            for sample in samples
        }
        by_chrom.setdefault(record.chrom, []).append(
            VariantInfo(record.chrom, record.pos, sample_ps)
        )
        total_records += 1

    vcf.close()

    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda x: x.pos)

    return samples, by_chrom, total_records

## Determine whether adding a variant to a block would introduce a conflicting PS
def would_conflict(
    current_ps: Dict[str, Set[str]],
    variant: VariantInfo,
) -> bool:
    for sample, new_ps in variant.sample_ps.items():
        if new_ps is None:
            continue
        existing = current_ps.get(sample, set())
        if existing and new_ps not in existing:
            return True
    return False

## Add a variant's phase information to current consensus block
def add_variant(
    variant: VariantInfo,
    current_ps: Dict[str, Set[str]],
    supporting_samples: Set[str],
) -> None:
    for sample, ps in variant.sample_ps.items():
        if ps is None:
            continue
        current_ps.setdefault(sample, set()).add(ps)
        supporting_samples.add(sample)

## Build consensus phase blocks independently for each chromosome/gene
def build_consensus_blocks(
    by_chrom: Dict[str, List[VariantInfo]],
) -> List[ConsensusBlock]:
    blocks: List[ConsensusBlock] = []

    for chrom, variants in by_chrom.items():
        positions: List[int] = []
        current_ps: Dict[str, Set[str]] = {}
        supporting_samples: Set[str] = set()

        for variant in variants:
            if positions and would_conflict(current_ps, variant):
                blocks.append(
                    ConsensusBlock(
                        chrom=chrom,
                        positions=list(positions),
                        supporting_samples=set(supporting_samples),
                        ps_values_by_sample={
                            sample: set(values)
                            for sample, values in current_ps.items()
                        },
                    )
                )
                positions = []
                current_ps = {}
                supporting_samples = set()

            positions.append(variant.pos)
            add_variant(variant, current_ps, supporting_samples)

        if positions:
            blocks.append(
                ConsensusBlock(
                    chrom=chrom,
                    positions=list(positions),
                    supporting_samples=set(supporting_samples),
                    ps_values_by_sample={
                        sample: set(values)
                        for sample, values in current_ps.items()
                    },
                )
            )

    return blocks

## Write consensus blocks to GTF
def write_gtf(
    blocks: Sequence[ConsensusBlock],
    output_path: Path,
) -> List[Tuple[str, ConsensusBlock]]:
    written: List[Tuple[str, ConsensusBlock]] = []
    block_number = 0

    with open_output(output_path, "wt") as handle:
        for block in blocks:

            block_number += 1
            block_id = f"block{block_number}"
            attributes = (
                f'block_id "{block_id}"; '
                f'snp_count "{len(block.positions)}"; '
                f'phase_support_samples "{len(block.supporting_samples)}";'
            )

            handle.write(
                f"{block.chrom}\tWhatsHapConsensus\tblock\t"
                f"{block.start}\t{block.end}\t.\t+\t.\t"
                f"{attributes}\n"
            )
            written.append((block_id, block))

    return written

## Generate an optional CSV summary containing information for all blocks
def write_summary(
    output_path: Path,
    written_blocks: Sequence[Tuple[str, ConsensusBlock]],
) -> None:
    """Write one row per consensus block using pandas."""
    rows = []

    for block_id, block in written_blocks:
        ps_text = ";".join(
            f"{sample}={','.join(sorted(values))}"
            for sample, values in sorted(block.ps_values_by_sample.items())
        )
        rows.append(
            {
                "Block_ID": block_id,
                "Gene": block.chrom,
                "Start": block.start,
                "End": block.end,
                "Span_bp": block.end - block.start + 1,
                "SNP_Count": len(block.positions),
                "SNP_Positions": ";".join(map(str, block.positions)),
                "Phase_Support_Sample_Count": len(block.supporting_samples),
                "Phase_Support_Samples": ";".join(
                    sorted(block.supporting_samples)
                ),
                "Sample_PS_Values": ps_text,
            }
        )

    columns = [
        "Block_ID", "Gene", "Start", "End", "Span_bp", "SNP_Count",
        "SNP_Positions", "Phase_Support_Sample_Count",
        "Phase_Support_Samples", "Sample_PS_Values",
    ]
    pd.DataFrame(rows, columns=columns).to_csv(
        output_path, index=False, encoding="utf-8"
    )

## Main execution
def main() -> None:
    args = parse_args()

    samples, by_chrom, total_records = read_variants(args.input_vcf)
    blocks = build_consensus_blocks(by_chrom)
    written = write_gtf(
        blocks,
        args.output_gtf,
    )

    if args.block_summary is not None:
        write_summary(args.block_summary, written)

    print("==============================================")
    print("Consensus phase-block detection complete")
    print("==============================================")
    print(f"Input VCF records:         {total_records}")
    print(f"Samples considered:        {len(samples)}")
    print(f"Genes/CHROM values:        {len(by_chrom)}")
    print(f"Consensus blocks detected: {len(blocks)}")
    print(f"Blocks written:            {len(written)}")
    print(f"Output GTF:                {args.output_gtf}")
    if args.block_summary is not None:
        print(f"Block summary:             {args.block_summary}")


if __name__ == "__main__":
    main()
