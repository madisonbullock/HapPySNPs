#!/bin/bash
set -eo pipefail

#### Usage: bash variantcall.sh ref.fasta SampleID
#### Where .fasta is reference sequence + ' ' + samplename 

## Set variables
reference=$1
prefix=$2
THREADS=${3:-$(nproc)}
read1fn="${prefix}_trimmed_R1.fastq.gz"
read2fn="${prefix}_trimmed_R2.fastq.gz"

#Read group is set at alignment time, so no separate AddOrReplaceReadGroups pass is needed
RG="@RG\tID:${prefix}\tSM:${prefix}\tPL:illumina\tLB:lib1\tPU:unit1"

#Align read files to reference and sort; -F 2304 drops secondary (256) and supplementary (2048) alignments
#-K pins the batch size so the insert size estimate, and the proper-pair flags that
#follow from it, don't change with -t
#-u leaves the pipe uncompressed for sort; -l 1 because the sorted BAM is an
#intermediate that gets deleted below
bwa mem -K 100000000 -t $THREADS -R "$RG" $reference $read1fn $read2fn \
  | samtools view -@ $THREADS -u -q 20 -F 2304 - \
  | samtools sort -@ $THREADS -l 1 - -o "${prefix}.sorted.bam"

#Mark duplicate sequences on the RG-tagged BAM
#--CREATE_INDEX writes the .bai as the BAM is written, so no separate indexing pass
#--TMP_DIR keeps spill files off the shared filesystem
gatk MarkDuplicates \
  -I "${prefix}.sorted.bam" \
  -O "${prefix}.marked.noDups.bam" \
  -M "${prefix}.metrics.txt" \
  --REMOVE_DUPLICATES true \
  --CREATE_INDEX true \
  --TMP_DIR "${TMPDIR:-/tmp}"

#picard names the index <base>.bai, htslib wants <base>.bam.bai
mv -f "${prefix}.marked.noDups.bai" "${prefix}.marked.noDups.bam.bai"

#Remove intermediate files, more intermediate files can be added if prefered
rm -f "${prefix}.sorted.bam"
