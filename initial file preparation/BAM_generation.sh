#!/bin/bash
set -eo pipefail

#### Usage: bash variantcall.sh ref.fasta SampleID
#### Where .fasta is reference sequence + ' ' + samplename 

## Set variables
reference=$1
prefix=$2
read1fn="${prefix}_trimmed_R1.fastq.gz"
read2fn="${prefix}_trimmed_R2.fastq.gz"

#Align read files to reference sequence and map
bwa mem $reference $read1fn $read2fn | samtools view -bS -q 20 -F 2048 -F 256 - | samtools sort - -o "$prefix.sorted.bam"
gatk FastqToSam -F1 $read1fn -F2 $read2fn -O $prefix.unmapped.bam -SM $prefix.sorted.bam

#Replace read groups to mapped and unmapped bam files using library prep and sequencing information
gatk AddOrReplaceReadGroups -I  $prefix.sorted.bam -O $prefix.sorted-RG.bam -RGID $prefix -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM $prefix

#Mark duplicate sequences
gatk MarkDuplicates -I $prefix.sorted-RG.bam -O $prefix.markedDups.bam -M $prefix.metrics.txt
samtools index $prefix.markedDups.bam

#Remove intermediate files, more intermediate files can be added if prefered
rm $prefix.unmapped.bam $prefix.unmapped-RG.bam
