# HapPy SNPs

⚠️ Under development ⚠️

## Purpose

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
Next, this protocol utilizes a BAM file generation script adapted from the [HybSep-SNP-Extraction](https://github.com/lindsawi/HybSeq-SNP-Extraction) script `variantcall.sh`. The adapted script is called `variantcall_HapPySNP.sh`. This script will map paired-end reads to supercontigs, replace read groups for mapped BAM files, mark duplicate reads, and remove intermediate BAM files.

Usage for single samples: `bash variantcall_HapPySNP.sh pseudo_reference.fasta SampleID`
Usage for multiple samples: 
```bash
while read sample
do
bash variantcall_HapPySNP.sh pseudo_reference.fasta $sample
done < samples.txt
```

## Workflow
### Software
- Python
  - pandas: https://pandas.pydata.org/
  - numpy: https://numpy.org/
  - scikit-learn: https://scikit-learn.org/stable/
  - pysam: https://pysam.readthedocs.io/en/stable/


