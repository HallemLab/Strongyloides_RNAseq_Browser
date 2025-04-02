## mRNA Extraction Script
# This script takes as input a GFF3 file and a genome fasta file and generates a fasta file containing the mRNA sequences

# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GenomicFeatures", "Biostrings", "txdbmaker"))

library(GenomicFeatures)
library(Biostrings)
library(txdbmaker)

# Define input and output file paths
gff_file <- "sster_liftoff_WBPS18.gff3"        
genome_fasta <- "strongyloides_stercoralis.PRJNA930454.WBPS19.genomic.fa"  
output_fasta <- "strongyloides_stercoralis.PRJNA930454.liftoff.fasta"

# Load genome sequences
genome <- readDNAStringSet(genome_fasta, format = "fasta")

# Clean chromosome names by removing "length=" and any numbers following it
names(genome) <- sub(" length=.*", "", names(genome))

# Create a TxDb object from the GFF3 file
txdb <- makeTxDbFromGFF(gff_file, format = "gff3")

# Extract exons grouped by transcripts
exons_by_transcripts <- exonsBy(txdb, by = "tx", use.names = TRUE)

names(exons_by_transcripts)<- paste0(names(exons_by_transcripts), ' gene=', names(exons_by_transcripts))
names(exons_by_transcripts)<- sub("\\.[0-9]$", "", test)

# Extract mRNA sequences using the exons
mRNA_sequences <- extractTranscriptSeqs(genome, exons_by_transcripts)

# Write the mRNA sequences to a FASTA file
writeXStringSet(mRNA_sequences, filepath = output_fasta)