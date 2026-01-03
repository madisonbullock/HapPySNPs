import pandas as pd
import numpy as np
from sklearn.impute import KNNImputer, SimpleImputer
import argparse
from io import StringIO
from collections import Counter
from scipy.stats import mode

def read_vcf(file_path):
    with open(file_path, 'r') as f:
        lines = [line for line in f if not line.startswith('##')]
    header = lines[0].strip().lstrip('#').split('\t')
    df = pd.read_csv(StringIO(''.join(lines[1:])), sep='\t', names=header)
    return df

def is_missing(geno):
    return geno.startswith('./.') or geno.startswith('.|.')

def extract_genotype_only(geno):
    return geno.split(':')[0] if ':' in geno else geno

def encode_genotypes(df):
    mapping = {
        '0/0': 0, '0|0': 0,
        '0/1': 1, '1/0': 1, '0|1': 1, '1|0': 1,
        '1/1': 2, '1|1': 2
    }
    numeric = df.applymap(lambda g: mapping.get(extract_genotype_only(g), np.nan))
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

def decode_genotypes(imputed, original, ref_df, method, detailed_log_path, confidence_threshold):
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
                    confidence = 1.0  # Currently placeholder (confidence filtering done earlier)
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
            else:
                decoded.at[idx, col] = orig_val

    with open(detailed_log_path, 'w') as log:
        log.write("Sample,SNP,Method,Confidence\n")
        log.write('\n'.join(detailed_entries) + '\n')

    return decoded

def count_missing_genotypes(df):
    return sum(is_missing(val) for col in df.columns for val in df[col])

def calculate_mode_confidence(encoded_col):
    vals = encoded_col[~np.isnan(encoded_col)]
    if len(vals) == 0:
        return 0.0
    m, count = mode(vals, nan_policy='omit')
    return count[0] / len(vals)

def calculate_knn_confidence(imputer, encoded_df, col_idx):
    # This is an approximation: use variance of imputed values as inverse confidence
    imputed_values = imputer.transform(encoded_df)[:, col_idx]
    var = np.var(imputed_values)
    confidence = 1 - min(var, 1)  # Clamp variance max to 1
    return confidence

def impute_vcf_phased(input_vcf, output_vcf, log_file, detailed_log, knn=True, confidence_threshold=0.8):
    with open(input_vcf, 'r') as f:
        all_lines = f.readlines()

    header_lines = [l for l in all_lines if l.startswith('##')]
    column_header_line = next(l for l in all_lines if l.startswith('#CHROM'))

    df = read_vcf(input_vcf)
    ref_df = df[['CHROM', 'POS']].copy()
    original_geno = df.iloc[:, 9:].copy()
    missing_before = count_missing_genotypes(original_geno)

    encoded = encode_genotypes(original_geno)

    if knn:
        imputer = KNNImputer(n_neighbors=5)
        method = "KNN"
        imputed_array = imputer.fit_transform(encoded)
        imputed_df = pd.DataFrame(imputed_array, index=encoded.index, columns=encoded.columns)

        confidences = []
        for col_idx, col in enumerate(encoded.columns):
            conf = calculate_knn_confidence(imputer, encoded, col_idx)
            confidences.append(conf)

        # Filter variants by confidence threshold
        high_conf_variants = [col for col, conf in zip(encoded.columns, confidences) if conf >= confidence_threshold]
        low_conf_variants = [col for col in encoded.columns if col not in high_conf_variants]

        imputed_df = imputed_df[high_conf_variants]
        encoded = encoded[high_conf_variants]
        original_geno = original_geno[high_conf_variants]
        ref_df = ref_df.loc[encoded.index]

    else:
        method = "Mode"
        imputer = SimpleImputer(strategy="most_frequent")
        imputed_array = imputer.fit_transform(encoded)
        imputed_df = pd.DataFrame(imputed_array, index=encoded.index, columns=encoded.columns)

        # For mode imputation, keep all variants (no confidence filtering)
        high_conf_variants = list(encoded.columns)
        low_conf_variants = []

    decoded_df = decode_genotypes(imputed_df, original_geno, ref_df, method, detailed_log, confidence_threshold)

    df.iloc[:, 9:] = decoded_df

    missing_after = count_missing_genotypes(decoded_df)
    snps_before = df.shape[0]
    df_clean = df[~df.iloc[:, 9:].applymap(is_missing).any(axis=1)]
    snps_removed = snps_before - df_clean.shape[0] + len(low_conf_variants)

    with open(log_file, 'w') as log:
        log.write(f"Missing genotypes before imputation: {missing_before}\n")
        log.write(f"Missing genotypes after imputation: {missing_after}\n")
        log.write(f"SNPs removed due to remaining missing data: {snps_removed}\n")
        if knn:
            log.write(f"SNPs removed due to low confidence (<{confidence_threshold}): {len(low_conf_variants)}\n")
        else:
            log.write("No confidence filtering applied for mode imputation.\n")
        log.write(f"Confidence threshold used: {confidence_threshold}\n")
        log.write(f"KNN used: {knn}\n")

    with open(output_vcf, 'w') as out_f:
        for line in header_lines:
            out_f.write(line)
        out_f.write(column_header_line)
        df_clean.to_csv(out_f, sep='\t', index=False, header=False)

def main():
    parser = argparse.ArgumentParser(description="Impute missing phased genotypes in a VCF.")
    parser.add_argument("-i", "--input", required=True, help="Input VCF file")
    parser.add_argument("-o", "--output", required=True, help="Output VCF file")
    parser.add_argument("-l", "--log", required=True, help="Summary log file")
    parser.add_argument("-d", "--detailed-log", required=True, help="Detailed per-genotype log")
    parser.add_argument("--knn", action="store_true", help="Use KNN imputation (default is mode imputation)")
    parser.add_argument("-c", "--confidence-threshold", type=float, default=0.8, help="Confidence threshold for filtering KNN-imputed variants")
    args = parser.parse_args()

    impute_vcf_phased(args.input, args.output, args.log, args.detailed_log, args.knn, args.confidence_threshold)

if __name__ == "__main__":
    main()
