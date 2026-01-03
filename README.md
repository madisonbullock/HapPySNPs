# HapPy SNPs

⚠️ Under development ⚠️

## Purpose
HapPy SNPs is a batch of Python scripts in a pipeline designed to process plant target capture sequencing data to extract phased microhaplotypes for downstream population genetic analyses. It integrates preparation steps such as assembly, variant calling, and phasing with custom imputation and microhaplotype summarization to generate high‑quality microhaplotype datasets from targeted nuclear gene regions, specifically [Angiosperms353](https://academic.oup.com/sysbio/article/68/4/594/5237557).

## Initial File Preparation
### Software
- HybPiper: https://github.com/mossmatters/HybPiper
- R: https://www.r-project.org/
- HybSep-SNP-Extraction: https://github.com/lindsawi/HybSeq-SNP-Extraction
  - BWA: https://bio-bwa.sourceforge.net/
  - GATK4: https://gatk.broadinstitute.org/hc/en-us
  - Samtools: https://www.htslib.org/
- FreeBayes: https://github.com/freebayes/freebayes
- VCFtools: https://vcftools.github.io/index.html
- WhatsHap: https://whatshap.readthedocs.io/en/latest/#

### Taxon-Specific Pseudo-Reference File Generation
First, run HybPiper on all samples to assemble target genes from paired-end trimmed FASTQ files using the mega353 target file. This can be done for each sample or in a loop, such as:
```bash
while read sample
do
hybpiper assemble -r ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz -t_dna mega353.fasta --bwa --prefix ${sample}
done < samples.txt
```
Then, generate statistics on the HybPiper output to determine per-sample, per-gene recovery:
```bash
hybpiper stats -t_dna mega353.fasta gene --hybpiper_dir /path/to/hybpiper/output
samples.txt
```
The sequence lengths file is used as an input to the R code `persample_pergene_recoverystats.R`, which helps determine which sample(s) will be used to create a taxon-specific pseudo-reference file. In general, it is best to limit the number of individuals who contribute to the pseudo-reference file. So, a good rule of thumb is to select the sample with the highest total number of genes recovered, identify the genes not recovered for that sample, and determine the smallest group of additional individuals to fill in the genes not recovered for the first sample.

Once the samples that will be used in the pseudo-reference file have been determined, retrieve supercontig sequences for each individual via:
```bash
hybpiper retrieve_sequences -t_dna mega353.fasta --single_sample_name SampleID --hybpiper_dir /path/to/hybpiper/output/ --fasta_dir /path/to/fasta/files supercontig
```
Then, manually copy the missing supercontig sequences from the additional individuals into the supercontig file of the sample with the highest total number of genes recovered. There should be only one sequence per gene in the final taxon-specific pseudo-reference file.

Finally, index the pseudo-reference file with:
```bash
bwa index pseudo_reference.fasta
samtools index pseudo_reference.fasta
```

### Variant Calling
Next, this protocol utilizes a BAM file generation script adapted from the [HybSep-SNP-Extraction](https://github.com/lindsawi/HybSeq-SNP-Extraction) script `variantcall.sh`. The adapted script is called `BAM_generation.sh`. This script will map paired-end reads to supercontigs, replace read groups for mapped BAM files, mark duplicate reads, and remove intermediate BAM files.
- Usage for single samples: `bash BAM_generation.sh pseudo_reference.fasta SampleID`
- Usage for multiple samples: 
```bash
while read sample
do
bash BAM_generation.sh pseudo_reference.fasta $sample
done < samples.txt
```
To efficiently generate a variant file, run FreeBayes in parallel to call simple SNPs and indels from a list of BAM file paths (`bam_list.txt`) against the taxon-specific pseudo-reference file. This process splits the genome into 250 bp chunks, ignores complex variants and priors, saves results to a VCF file, and logs any errors:
```bash
nohup freebayes-parallel <(fasta_generate_regions.py pseudo_reference.fasta.fai 250) 10 -f pseudo_reference.fasta -L bam_list.txt --haplotype-length 0 -kwVaXi > variants.markedDups.noMNP.noComplex.noPriors.vcf 2> freebayes_errors.log
```
Filter the VCF output from FreeBayes using VCFtools:
```bash
vcftools --gzvcf variants.markedDups.noMNP.noComplex.noPriors.vcf --max-missing 0.8 --remove-indels --minQ 36 --minDP 6 --recode --out variants.markedDups.noMNP.noComplex.noPriors.0.8Missing.minQ36.minDP6
```
The code above, in particular, keeps SNPs that are called in ≥80% of samples, have a minimum quality greater than or equal to 36, have a coverage depth of greater than or equal to 6 reads, and are not indels. These filters can be changed depending on the researcher's preference. However, in preliminary analysis using the HapPySNPs optional imputation script, we found that decreasing `--max-missing` below 0.8 does not substantially increase the total number of SNPs retained. If the imputation script is not used, `--max-missing` must be set to `1.0`.

### Variant Phasing
WhatsHap is used to perform read-backed phasing of variants in the filtered VCF file using aligned reads from the BAM files and the pseudo-reference file, subsequently outputting a phased VCF that will be the input to the HapPySNPs workflow. WhatsHap was specifically chosen as it can leverage paired-end reads spanning multiple variants to infer haplotypes, and read-backed phasing is often more accurate than statistical phasing for small sample sizes or diverse populations.
```bash
whatshap phase -o variants.markedDups.noMNP.noComplex.noPriors.0.8Missing.minQ36.minDP6.phased.vcf --reference pseudo_reference.fasta variants.markedDups.noMNP.noComplex.noPriors.0.8Missing.minQ36.minDP6.vcf /path/to/BAM/files/*.markedDups.bam
```

## Workflow
### Software
- Python
  - pandas: https://pandas.pydata.org/
  - numpy: https://numpy.org/
  - scikit-learn: https://scikit-learn.org/stable/
  - pysam: https://pysam.readthedocs.io/en/stable/

### 0. Variant Imputation (Optional) - impute_missing_phased_snps.py
This script will:
- Read a phased VCF file with missing genotypes
- Fill in missing genotypes using KNN or mode imputation
- Keep the phasing information intact
- Remove low-confidence variants when using KNN imputation
- Remove variants that still have missing data after imputation
- Save a new VCF file with imputed genotypes
- Create logs summarizing the imputation results

Usage:
```bash
python impute_missing_phased_snps.py -i variants.markedDups.noMNP.noComplex.noPriors.0.8Missing.minQ36.minDP6.phased.vcf -o variants.markedDups.noMNP.noComplex.noPriors.0.8Missing.minQ36.minDP6.phased.imputed.vcf -l imputation_summary.log -d imputation_detailed.log --knn -c 0.9
```

Output:
- Imputed phased VCF file
- Summary log with missing data and filtering info
- Detailed log of each imputed genotype with confidence scores

### 1. Determine Phase Blocks - determine_phase_blocks.py
This script will:
- Read a phased VCF file
- Find blocks where the genotype phase is consistent
- Create a GTF file listing these phase blocks with their positions

Usage:
```bash
python determine_phase_blocks.py -i variants.markedDups.noMNP.noComplex.noPriors.0.8Missing.minQ36.minDP6.phased.imputed.vcf -o output_phase_blocks.gtf
```

Output:
- GTF file listing detected phase blocks

### 2. Extract Microhaplotype SNPs - extract_microhap_SNPs.py
This script will:
- Read a phased VCF, BAM file list, phase block GTF, and imputation log
- Extract genotypes and count allele-supporting reads per sample
- Mark which genotypes were imputed and assign default values for them
- Output a CSV summarizing genotype, read counts, phasing blocks, and imputation status

Usage:
```bash
python extract_microhap_SNPs.py -v variants.markedDups.noMNP.noComplex.noPriors.0.8Missing.minQ36.minDP6.phased.imputed.vcf -b bam_list.txt -g output_phase_blocks.gtf -i imputation_detailed.log -o microhap_snps.csv
```

Output:
- CSV file containing microhaplotype SNP data per sample, including genotype, allele counts, allele balance, phase block IDs, and imputation status

### 3. Assemble Microhaplotypes - assemble_microhaplotypes.py
This script will:
- Read a CSV file with phased microhaplotype SNP data per sample and variant
= Group variants by gene, sample, and haplotype block
- Combine phased SNP alleles into continuous microhaplotype sequences
- Calculate average read counts and depths per SNP for each haplotype
- Filter out single-SNP haplotypes to keep only multi-SNP microhaplotypes
- Save the combined microhaplotype information with summary stats to a new CSV file

Usage:
```bash
python assemble_microhaplotypes.py -i microhap_snps.csv -o combined_microhaps.csv
```

Output:
- CSV file summarizing combined microhaplotypes per gene, sample, and block, including haplotype sequences, SNP positions, block length, and read depth averages

### 4. Calculate Selection Criteria - calculate_selection_criteria.py
This script will:
- Read a CSV file with combined microhaplotype data per gene and block
- Calculate and summarize block-level statistics, including:
  - Heterozygosity
  - Effective allele number (Ae)
  - Rosenberg’s Informativeness Index (In)
  - Average haplotype block length
- Compare blocks across metrics to assess concordance
- Output a CSV file with gene-level summary statistics

Usage:
```bash
python calculate_selection_criteria.py -i combined_microhaps.csv -o selection_criteria_summary.csv
```

Output:
- CSV file summarizing diversity and informativeness statistics per gene and haplotype block

### 5. Select Microhaplotypes - select_microhaplotypes.py
This script will:
- Read combined gene-block statistics from a CSV file
- Select blocks based on one chosen metric:
  - Highest heterozygosity
  - Highest effective allele number (Ae)
  - Highest Rosenberg’s Informativeness index (In)
  - Highest average haplotype block length
- Read microhaplotype genotype data per sample from a CSV file and filter microhaplotypes to keep only those matching the selected blocks per gene
- Output a filtered CSV file with the selected microhaplotypes

Usage:
```bash
python select_microhaplotypes.py -I -c selection_criteria_summary.csv -m combined_microhaps.csv -o selected_microhaps.csv
```
*Replace -I with one of -A (Highest effective allele number), -L (Highest average haplotype block length), or -H (Highest heterozygosity) to select a different filter*

Output:
- Filtered microhaplotype CSV file containing only microhaplotypes from the selected blocks based on the specified metric

### 6. Assign Allele IDs - assign_allele_IDs.py
This script will:
- Read a CSV file with filtered microhaplotypes per gene, sample, and block
- Assign unique numeric allele IDs to each distinct microhaplotype sequence per gene and block, similar to what is used with microsatellites
- Add these numeric allele IDs as new columns for each microhaplotype in the CSV
- Create a dictionary CSV mapping genes, blocks, microhaplotype sequences, and their assigned allele IDs
- Save the updated microhaplotype data with allele IDs and the dictionary file

Usage:
```bash
python assign_allele_IDs.py -i selected_microhaps.csv -o selected_microhaps_with_alleleIDs.csv -d allele_dictionary.csv
```

Output:
- CSV file with microhaplotypes and their assigned numeric allele IDs per haplotype
- CSV dictionary mapping genes, blocks, microhaplotype sequences, and allele IDs
