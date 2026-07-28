#!/usr/bin/env python3

import pysam
import argparse
import bisect
import csv
import os
import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed


class BlockIndex:
    """Per-chromosome interval lookup for block IDs, by bisect on the starts.

    Intervals are kept sorted and non-overlapping, so a position maps to at most
    one block. Where the GTF has overlapping blocks, the one later in the file
    wins over the overlap.
    """

    def __init__(self, resolved_by_chrom):
        """resolved_by_chrom: {chrom: sorted [(start, end_excl, block_id)]}"""
        self._starts = {}
        self._spans = {}
        for chrom, intervals in resolved_by_chrom.items():
            if intervals:
                self._starts[chrom] = [iv[0] for iv in intervals]
                self._spans[chrom] = [(iv[1], iv[2]) for iv in intervals]

    def lookup(self, chrom, pos):
        starts = self._starts.get(chrom)
        if not starts:
            return "NA"
        i = bisect.bisect_right(starts, pos) - 1
        if i < 0:
            return "NA"
        end_excl, block_id = self._spans[chrom][i]
        return block_id if pos < end_excl else "NA"


def _resolve_last_wins(intervals):
    """Return a sorted, non-overlapping list of (start, end_excl, block_id).

    Where inputs overlap, whichever interval comes later in the list keeps the
    overlap. Input that doesn't overlap (the usual case) just gets sorted.
    """
    if not intervals:
        return []
    result = []
    for s_new, e_new, bid_new in intervals:
        if s_new >= e_new:
            continue  # skip empty/malformed intervals
        # Clip existing intervals that overlap [s_new, e_new); the new block wins there
        clipped = []
        for s, e, bid in result:
            if s >= e_new or e <= s_new:
                clipped.append((s, e, bid))         # no overlap, keep intact
            else:
                if s < s_new:
                    clipped.append((s, s_new, bid)) # keep left portion
                if e > e_new:
                    clipped.append((e_new, e, bid)) # keep right portion
        clipped.append((s_new, e_new, bid_new))
        result = clipped
    result.sort(key=lambda x: x[0])
    return result


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract microhaplotype SNPs from VCF and BAMs with imputation awareness."
    )
    parser.add_argument("-v", "--vcf", required=True, help="Input phased VCF file (can be uncompressed or bgzipped)")
    parser.add_argument("-b", "--bam_list", required=True, help="BAM input file list: CSV (Sample,BAM) or plain list of BAM paths")
    parser.add_argument("-g", "--gtf", required=True, help="GTF file containing block_id annotations")
    parser.add_argument("-i", "--imputation_log", required=True, help="CSV log of imputed genotypes")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
    parser.add_argument("-d", "--min_depth", type=int, default=6, help="Minimum depth to use for imputed SNPs (default: 6)")
    parser.add_argument("--bam_list_is_csv", action="store_true",
                        help="Specify if BAM list (-b) is a CSV with Sample,BAM columns. If not set, assumes plain BAM path list with sampleID extracted from filename.")
    parser.add_argument("-t", "--threads", type=int, default=None,
                        help="Worker processes for parallel read counting (default: all CPUs).")
    return parser.parse_args()

def bgzip_and_index(vcf_path): # compress vcf
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

def load_blocks_from_gtf(gtf_file): # use GTF block file to identify blocks in sequences
    """Read the block features into a BlockIndex keyed by chromosome.

    If two blocks cover the same position, the one later in the file wins.
    """
    raw = defaultdict(list)
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            chrom, _, feature, start, end, _, _, _, attributes = fields
            if feature != "block":
                continue
            start, end = int(start), int(end)
            block_id = attributes.split('block_id "')[1].split('"')[0]
            raw[chrom].append((start, end + 1, block_id))  # GTF end is inclusive, store half-open

    resolved = {chrom: _resolve_last_wins(intervals) for chrom, intervals in raw.items()}
    return BlockIndex(resolved)

def load_imputed_positions(log_file): # load imputation log to determine which SNPs were imputed
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

def load_bam_paths_from_csv(csv_file): # load in BAM file paths
    bam_paths = {}
    with open(csv_file) as f:
        reader = csv.reader(f)
        header = next(reader)
        sample_col = header.index("Sample") if "Sample" in header else 0
        bam_col = header.index("BAM") if "BAM" in header else 1
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

def match_vcf_samples_to_bams(vcf_samples, bam_paths): # flexible matching of BAM files to VCF
    """Match VCF sample names to BAM paths robustly."""
    sample_to_bam = {}
    for sample in vcf_samples:
        # Direct match
        if sample in bam_paths:
            sample_to_bam[sample] = bam_paths[sample]
            continue

        # Match ignoring case
        sample_lower = sample.lower()
        for s, bam in bam_paths.items():
            if s.lower() == sample_lower:
                sample_to_bam[sample] = bam
                break
        else:
            # Match by partial name (sample appears in BAM filename)
            for s, bam in bam_paths.items():
                if sample in s or s in sample:
                    sample_to_bam[sample] = bam
                    break

        if sample not in sample_to_bam:
            print(f"⚠️ Warning: No BAM found for sample '{sample}'")
    return sample_to_bam

def fetch_alleles(vcf, record, sample): # extract alelles per individual
    gt = record.samples[sample].get("GT")
    if gt is None or "." in str(gt):
        return None, None, "./."
    alleles = [record.alleles[i] if i < len(record.alleles) else "N" for i in gt]
    genotype_str = "|".join(map(str, gt))
    return alleles[0], alleles[1], genotype_str

# nothing calls this; it's the plain one-position version, kept for reference.
# the live path is _count_reads_contig, which pileups a whole contig in one go
def count_reads(samfile, chrom, pos, alleles):
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
    return allele_counts


def _count_reads_contig(samfile, chrom, chrom_tasks):
    # counts every position we need on one contig in a single pileup pass and
    # returns {pos: {base: count}}. inner loop is the same as count_reads.
    needed = {}
    for t in chrom_tasks:
        needed.setdefault(t["pos"], {})
    if not needed:
        return needed
    lo, hi = min(needed), max(needed)
    for pileupcolumn in samfile.pileup(chrom, lo - 1, hi, truncate=True):
        hist = needed.get(pileupcolumn.pos + 1)  # pysam pos is 0-based
        if hist is None:
            continue
        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del or pileupread.is_refskip:
                continue
            base = pileupread.alignment.query_sequence[pileupread.query_position]
            hist[base] = hist.get(base, 0) + 1
    return needed


# open BAM handles per worker process, keyed by path
_BAM_HANDLES = {}


def _get_bam(bam_path):
    h = _BAM_HANDLES.get(bam_path)
    if h is None:
        h = pysam.AlignmentFile(bam_path, "rb")
        _BAM_HANDLES[bam_path] = h
    return h



def _process_sample(bam_path, sample, tasks):
    """Not called anywhere. Does a whole sample in one worker instead of one
    (sample, contig) like _process_unit, kept as a reference.
    Returns (variant_idx, sample_slot, row).
    """
    results = []
    samfile = pysam.AlignmentFile(bam_path, "rb")
    try:
        tasks_by_chrom = defaultdict(list)                               # this sample's tasks grouped by contig
        for t in tasks:
            tasks_by_chrom[t["chrom"]].append(t)
        hist_by_chrom = {c: _count_reads_contig(samfile, c, ts)          # one pileup per contig
                         for c, ts in tasks_by_chrom.items()}
        for t in tasks:
            allele1, allele2 = t["allele1"], t["allele2"]
            allele_counts = hist_by_chrom[t["chrom"]][t["pos"]]
            total1 = allele_counts.get(allele1, 0)
            total2 = allele_counts.get(allele2, 0)
            total_reads = total1 + total2
            ab = round(total2 / total_reads, 3) if total_reads > 0 else 0
            row = [
                t["chrom"], t["pos"], sample, t["genotype_str"],
                allele1, allele2, ab,
                total1, total1, total2, total2,
                t["block_id"], "No"
            ]
            results.append((t["variant_idx"], t["sample_slot"], row))
    finally:
        samfile.close()
    return results


def _process_unit(bam_path, sample, chrom, tasks):
    # worker for a single (sample, contig). rows come back in whatever order the
    # pool finishes; pass 3 puts them back in order.
    samfile = _get_bam(bam_path)              # reuse cached handle
    base_hist = _count_reads_contig(samfile, chrom, tasks)
    results = []
    for t in tasks:
        allele1, allele2 = t["allele1"], t["allele2"]
        allele_counts = base_hist[t["pos"]]
        total1 = allele_counts.get(allele1, 0)
        total2 = allele_counts.get(allele2, 0)
        total_reads = total1 + total2
        ab = round(total2 / total_reads, 3) if total_reads > 0 else 0
        row = [
            t["chrom"], t["pos"], sample, t["genotype_str"],
            allele1, allele2, ab,
            total1, total1, total2, total2,
            t["block_id"], "No"
        ]
        results.append((t["variant_idx"], t["sample_slot"], row))
    return results


def main():
    args = parse_args()
    vcf_path = bgzip_and_index(args.vcf)
    vcf = pysam.VariantFile(vcf_path)

    block_index = load_blocks_from_gtf(args.gtf)
    imputed = load_imputed_positions(args.imputation_log)

    if args.bam_list_is_csv:
        bam_paths_raw = load_bam_paths_from_csv(args.bam_list)
    else:
        bam_paths_raw = load_bam_paths_from_list(args.bam_list)

    bam_paths = match_vcf_samples_to_bams(vcf.header.samples, bam_paths_raw)
    sample_order = list(vcf.header.samples)

    print("VCF samples:", list(vcf.header.samples))
    print("Matched BAMs:", list(bam_paths.keys()))
    print(f"Total imputed positions loaded: {len(imputed)}")

    # --- Pass 1: walk the VCF, write imputed rows, queue the rest for the BAMs ---
    # samples with no matching BAM and missing GTs get dropped here. imputed rows
    # never touch a BAM, so they go straight into results.
    n_variants = 0
    work_by_unit = defaultdict(list)    # (sample, contig) -> task dicts
    results = {}                        # (variant_idx, sample_slot) -> final CSV row

    for variant_idx, record in enumerate(vcf.fetch()):
        chrom, pos = record.chrom, record.pos
        block_id = block_index.lookup(chrom, pos)
        n_variants = variant_idx + 1

        for sample_slot, sample in enumerate(sample_order):
            if sample not in bam_paths:
                continue
            allele1, allele2, genotype_str = fetch_alleles(vcf, record, sample)
            if allele1 is None:
                continue
            if (sample, chrom, pos) in imputed:
                dp = args.min_depth
                results[(variant_idx, sample_slot)] = [
                    chrom, pos, sample, genotype_str, allele1, allele2,
                    0.0, dp, dp, dp, dp, block_id, "Yes"
                ]
            else:
                work_by_unit[(sample, chrom)].append({
                    "variant_idx": variant_idx,
                    "sample_slot": sample_slot,
                    "chrom": chrom,
                    "pos": pos,
                    "allele1": allele1,
                    "allele2": allele2,
                    "genotype_str": genotype_str,
                    "block_id": block_id,
                })

    # --- Pass 2: read counting in parallel, one task per (sample, contig) ---
    # splitting per contig rather than per sample stops one huge sample from
    # holding up the pool at the end. results arrive out of order.
    n_workers = max(1, min(
        args.threads if args.threads is not None else (os.cpu_count() or 1),
        len(work_by_unit) if work_by_unit else 1
    ))

    if work_by_unit:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            futures = {
                pool.submit(_process_unit, bam_paths[s], s, c, tasks): (s, c)
                for (s, c), tasks in work_by_unit.items()
            }
            for fut in as_completed(futures):
                for vi, ss, row in fut.result():
                    results[(vi, ss)] = row

    # --- Pass 3: write rows in VCF variant order, samples in header order ---
    # walking (variant_idx, sample_slot) upward gets that order back no matter
    # how the workers finished
    header = [
        "Gene", "Variant_Pos", "Sample", "Genotype", "Allele1", "Allele2",
        "Allele_Balance", "Total_Reads_Allele1", "Read_Depth_Allele1",
        "Total_Reads_Allele2", "Read_Depth_Allele2", "Haplotype_Block", "Imputed"
    ]
    with open(args.output, "w", newline="") as out_f:
        writer = csv.writer(out_f)
        writer.writerow(header)
        for vi in range(n_variants):
            for ss in range(len(sample_order)):
                row = results.get((vi, ss))
                if row is not None:
                    writer.writerow(row)


if __name__ == "__main__":
    main()
