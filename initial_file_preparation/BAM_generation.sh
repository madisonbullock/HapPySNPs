#!/bin/bash

## Stop script if an error happens:
set -euo pipefail

## General usage:
### bash variantcall_microhaps.sh reference.fasta sampleID /path/to/fastq 4
### Where 4 is the number of threads to run on and can be changed to user preferences

## Set variables:
reference=$1
prefix=$2
fastq_dir=$3
threads=${4:-4}
read1fn="${fastq_dir}/${prefix}_trimmed_R1.fastq.gz" # Change to match naming format of the file
read2fn="${fastq_dir}/${prefix}_trimmed_R2.fastq.gz" # Change to match naming format of the file

## Set outputs:
sorted_bam="${prefix}.sorted.bam"
rg_bam="${prefix}.sorted-RG.bam"
final_bam="${prefix}.marked.noDups.bam"
final_bai="${prefix}.marked.noDups.bai"
metrics_file="${prefix}.marked_dup_metrics.txt"

## Align reads, filter alignments, and sort BAM:
bwa mem -t "$threads" "$reference" "$read1fn" "$read2fn" |
samtools view -b -q 20 -F 2304 - |
samtools sort -@ "$threads" -o "$sorted_bam" -

## Validate sorted BAM:
samtools quickcheck "$sorted_bam"

## Add read-group information and validate:
gatk AddOrReplaceReadGroups -I "$sorted_bam" -O "$rg_bam" -RGID "$prefix" -RGLB "$prefix" -RGPL ILLUMINA -RGPU "$prefix" -RGSM "$prefix"
samtools quickcheck "$rg_bam"

## Remove duplicate reads, validate, and index final BAM:
gatk MarkDuplicates -I "$rg_bam" -O "$final_bam" -M "$metrics_file" --REMOVE_DUPLICATES true
samtools quickcheck "$final_bam"
samtools index -@ "$threads" -o "$final_bai" "$final_bam"

## Remove intermediate files:
rm -f "$sorted_bam" "$rg_bam"
