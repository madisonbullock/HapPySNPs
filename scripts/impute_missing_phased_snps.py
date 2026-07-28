import pandas as pd
import numpy as np
from sklearn.impute import SimpleImputer
import argparse
from io import StringIO
from collections import Counter
#from scipy.stats import mode

def read_vcf(file_path, all_lines=None):
    """Parse a VCF into a DataFrame, optionally reusing already-read lines."""
    if all_lines is None:
        with open(file_path, 'r') as f:
            all_lines = f.readlines()
    lines = [line for line in all_lines if not line.startswith('##')]
    header = lines[0].strip().lstrip('#').split('\t')
    df = pd.read_csv(StringIO(''.join(lines[1:])), sep='\t', names=header)
    return df

def is_missing(geno):
    return geno.startswith('./.') or geno.startswith('.|.')

def _missing_mask(df):
    """Boolean DataFrame, True where a genotype cell is missing."""
    return df.apply(lambda c: c.astype(str).str.startswith(('./.', '.|.')))

def extract_genotype_only(geno):
    return geno.split(':')[0] if ':' in geno else geno

def encode_genotypes(df):
    mapping = {
        '0/0': 0, '0|0': 0,
        '0/1': 1, '1/0': 1, '0|1': 1, '1|0': 1,
        '1/1': 2, '1|1': 2
    }
    numeric = df.map(lambda g: mapping.get(extract_genotype_only(g), np.nan))
    return numeric

def most_common_phased_het(genotypes):
    counter = Counter()
    for gt in genotypes:
        if gt in {'0|1', '1|0'}:
            counter[gt] += 1
    if not counter:
        return '0|1'
    most_common = counter.most_common()
    if len(most_common) == 1 or most_common[0][1] != most_common[1][1]:
        return most_common[0][0]
    return '0|1'

def decode_genotypes(imputed, original, ref_df, method, detailed_log_path):
    inv_map = {0: '0|0', 1: 'HET', 2: '1|1'}
    decoded = original.copy()
    detailed_entries = []

    for col in imputed.columns:
        phased_het = most_common_phased_het(original[col].map(extract_genotype_only))
        for idx in imputed.index:
            orig_val = original.at[idx, col]
            imp_val = imputed.at[idx, col]

            if is_missing(orig_val):
                if not np.isnan(imp_val):
                    snp_id = f"{ref_df.at[idx, 'CHROM']}:{ref_df.at[idx, 'POS']}"
                    confidence = 1.0  # mode imputation gives no per-genotype score; column kept so the log format doesn't change
                    gt_code = int(round(imp_val))
                    if inv_map[gt_code] == 'HET':
                        new_gt = phased_het
                    else:
                        new_gt = inv_map[gt_code]
                    new_gt += (orig_val[3:] if ':' in orig_val else '')
                    decoded.at[idx, col] = new_gt
                    detailed_entries.append(f"{col},{snp_id},{method},{confidence:.2f}")
                else:
                    decoded.at[idx, col] = './.' + (orig_val[3:] if ':' in orig_val else '')

    with open(detailed_log_path, 'w') as log:
        log.write("Sample,SNP,Method,Confidence\n")
        log.write('\n'.join(detailed_entries) + '\n')

    return decoded

def count_missing_genotypes(df):
    return int(_missing_mask(df).sum().sum())

#def calculate_mode_confidence(encoded_col):
#    vals = encoded_col[~np.isnan(encoded_col)]
#    if len(vals) == 0:
#        return 0.0
#    m, count = mode(vals, nan_policy='omit')
#    return count[0] / len(vals)

def impute_vcf_phased(input_vcf, output_vcf, log_file, detailed_log):
    # read the file once, the same lines feed the header scan and read_vcf
    with open(input_vcf, 'r') as f:
        all_lines = f.readlines()

    header_lines = [l for l in all_lines if l.startswith('##')]
    column_header_line = next(l for l in all_lines if l.startswith('#CHROM'))

    df = read_vcf(input_vcf, all_lines)
    ref_df = df[['CHROM', 'POS']].copy()
    original_geno = df.iloc[:, 9:].copy()

    encoded = encode_genotypes(original_geno)

    # A sample with no usable biallelic genotype anywhere encodes to an all-NaN
    # column. SimpleImputer silently drops those, which leaves the frames a column
    # short and breaks the reconstruction below, so drop them here instead and keep
    # every frame and the #CHROM header lined up. Logged as well, so a dead sample
    # is visible rather than quietly filled in as 0|0.
    dropped_samples = [c for c in encoded.columns if bool(encoded[c].isna().all())]
    if dropped_samples:
        print(f"WARNING: dropping {len(dropped_samples)} sample(s) with no usable "
              f"biallelic genotypes: {', '.join(map(str, dropped_samples))}")
        df = df.drop(columns=dropped_samples)
        original_geno = original_geno.drop(columns=dropped_samples)
        encoded = encoded.drop(columns=dropped_samples)

    missing_before = count_missing_genotypes(original_geno)

    method = "Mode"
    imputer = SimpleImputer(strategy="most_frequent")
    imputed_array = imputer.fit_transform(encoded)
    imputed_df = pd.DataFrame(imputed_array, index=encoded.index, columns=encoded.columns)

    decoded_df = decode_genotypes(imputed_df, original_geno, ref_df, method, detailed_log)

    df.iloc[:, 9:] = decoded_df

    missing_after = count_missing_genotypes(decoded_df)
    snps_before = df.shape[0]
    df_clean = df[~_missing_mask(df.iloc[:, 9:]).any(axis=1)]
    snps_removed = snps_before - df_clean.shape[0]

    with open(log_file, 'w') as log:
        if dropped_samples:
            log.write(f"Samples dropped (no usable biallelic genotypes): {len(dropped_samples)}\n")
            log.write(f"Dropped sample names: {', '.join(map(str, dropped_samples))}\n")
        log.write(f"Missing genotypes before imputation: {missing_before}\n")
        log.write(f"Missing genotypes after imputation: {missing_after}\n")
        log.write(f"SNPs removed due to remaining missing data: {snps_removed}\n")

    with open(output_vcf, 'w') as out_f:
        for line in header_lines:
            out_f.write(line)
        # rebuild #CHROM from the current columns so the header width always matches
        # the body. same as the input line unless samples were dropped.
        out_f.write('#' + '\t'.join(map(str, df_clean.columns)) + '\n')
        df_clean.to_csv(out_f, sep='\t', index=False, header=False)

def main():
    parser = argparse.ArgumentParser(description="Impute missing phased genotypes in a VCF.")
    parser.add_argument("-i", "--input", required=True, help="Input VCF file")
    parser.add_argument("-o", "--output", required=True, help="Output VCF file")
    parser.add_argument("-l", "--log", required=True, help="Summary log file")
    parser.add_argument("-d", "--detailed-log", required=True, help="Detailed per-genotype log")
    args = parser.parse_args()

    impute_vcf_phased(args.input, args.output, args.log, args.detailed_log)

if __name__ == "__main__":
    main()
