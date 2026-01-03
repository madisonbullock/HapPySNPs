# Target Capture Gene Recovery per Sample
# The intention of this code is to determine which sample(s) recover the most genes, and to select samples for pseudo-reference file creation.

# Load required libraries
library(tibble)   # For data frame manipulation (e.g., relocate)
library(dplyr)    # For data manipulation using pipes and verbs

# Read sequence length table
# Rows = samples (plus reference), columns = genes
sample.data = as.matrix(
  read.table(
    "seq_lengths.tsv",
    header = TRUE,
    row.names = 1,
    sep = "\t"
  )
)

# Remove the first row (reference) to keep only sample sequence lengths
sample.len = sample.data[2:nrow(sample.data), ]

# Extract reference sequence lengths (first row)
reference.len = as.numeric(sample.data[1, ])

# Calculate percent recovery by dividing sample lengths by reference lengths
percent.len = sweep(sample.len, 2, as.numeric(reference.len), "/")

# Alternative calculation of percent recovery (overwrites previous line)
percent.len <- sample.len / as.numeric(reference.len)

# Cap percent recovery at 1 (100%)
percent.len = ifelse(percent.len > 1, 1, percent.len)

# Transpose matrix so rows = genes, columns = samples
pl.transposed <- t(percent.len)

# Identify the sample with the best recovery for each gene
bestrecovery <- apply(
  pl.transposed,
  1,
  function(x) colnames(pl.transposed)[which(x == max(x))]
)

# Count number of genes unrecovered (0 length) per sample
genesUnrecovered <- rowSums(percent.len == 0)

# Count number of genes recovered (non-zero length) per sample
genesRecovered <- rowSums(percent.len != 0)

# Add recovery statistics as new columns
percent.len.rec <- cbind(percent.len, genesUnrecovered, genesRecovered)

# Convert matrix to data frame
percent.len.rec.df <- as.data.frame(percent.len.rec)

# Move recovery summary columns to the front
seqlen <- percent.len.rec.df %>%
  relocate(genesRecovered, genesUnrecovered, .before = X4471)

# Subset samples with at least 50 genes recovered
subseqlen <- subset(seqlen, genesRecovered >= 50)
subseqlen <- as.data.frame(subseqlen)

# Subset samples with fewer than 50 genes recovered
dropped <- subset(seqlen, genesRecovered < 50)
dropped <- as.data.frame(dropped)

# Count unrecovered genes per gene across retained samples
samplesUnrecovered <- colSums(subseqlen == 0)

# Count recovered genes per gene across retained samples
samplesRecovered <- colSums(subseqlen != 0)

# Append sample-level recovery statistics as rows
subseqlen.srec <- rbind(
  subseqlen,
  samplesUnrecovered,
  samplesRecovered
)

# Rename the appended rows
# NOTE: Indices must match the size of 'subseqlen.srec'
row.names(subseqlen.srec)[49] <- "samplesUnrecovered"
row.names(subseqlen.srec)[50] <- "samplesRecovered"

# Reorder rows to place summary rows at the top
# NOTE: Adjust indices to match your data frame
subseqlen.srec <- subseqlen.srec[c(49, 50, 1:48), ]

# Write final gene recovery table to CSV
write.csv(
  subseqlen.srec,
  file = "GeneRecovery_withstats.csv"
)