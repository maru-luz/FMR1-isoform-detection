#!/bin/bash
# ==============================================================================
# SCRIPT 1: DEMULTIPLEXING, MAPPING AND FORMAT CONVERSION
# ==============================================================================
# Description: 
# This script performs basecalling demultiplexing, aligns ONT long-reads to the 
# human reference genome (GRCh38), and converts output files (SAM -> BAM -> BED12).
# 
# Biological Note on Reference Genome:
# We explicitly use the PRIMARY ASSEMBLY (unmasked) rather than 'toplevel' to 
# exclude alternative haplotypes and patches. This avoids multi-mapping ambiguities 
# and secondary assignments, which is crucial for characterizing FMR1 (ENSG00000102081) 
# located on chrX (Xq27.3). This also ensures strict compatibility with the standard 
# Ensembl GTF annotation.
#
# Requirements: 
# - dorado (v0.9.1)
# - conda environment containing minimap2 and bedtools (e.g., 'minimap2' env)
# ==============================================================================

# ------------------------------------------------------------------------------
# VARIABLES
# ------------------------------------------------------------------------------
# Change this variable to process different samples (e.g., "01", "11", etc.)
BARCODE="01"

# Paths to reference genome and annotation (Ensembl Release 115)
REF_FASTA="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
REF_MMI="Homo_sapiens.GRCh38.dna.primary_assembly.mmi"
REF_GTF="Homo_sapiens.GRCh38.115.gtf"

# Paths to data
RAW_FASTQ="barcode${BARCODE}-trimmed-ok.fastq"
MAPPED_SAM="bc${BARCODE}-mapped.sam"
SORTED_BAM="bc${BARCODE}-mapped.sorted.bam"
OUTPUT_BED="bc${BARCODE}.bed12"

echo "Starting Pipeline for Barcode: ${BARCODE}"

# ------------------------------------------------------------------------------
# 1. DEMULTIPLEXING (Dorado)
# ------------------------------------------------------------------------------
# Note: Ensure dorado is in your system path or specify the full path.
echo "[1/4] Demultiplexing..."
~/dorado-0.9.1-linux-x64/bin/dorado demux \
    -o barcode${BARCODE} \
    --kit-name SQK-NBD114-24 \
    --emit-fastq ${RAW_FASTQ}

# ------------------------------------------------------------------------------
# 2. GENOME INDEXING
# ------------------------------------------------------------------------------
# This step only needs to be run once. Uncomment if the .mmi index does not exist.
# echo "[*] Indexing reference genome..."
# minimap2 -d ${REF_MMI} ${REF_FASTA}

# ------------------------------------------------------------------------------
# 3. READ MAPPING (Minimap2)
# ------------------------------------------------------------------------------
echo "[2/4] Mapping reads to the reference genome..."
# --splice-flank yes: assumes biological conservation of splice sites
minimap2 -ax splice ${REF_MMI} ${RAW_FASTQ} \
    --splice-flank yes \
    -o ${MAPPED_SAM}

# ------------------------------------------------------------------------------
# 4. FORMAT CONVERSIONS (Samtools & Bedtools)
# ------------------------------------------------------------------------------
echo "[3/4] Converting SAM to sorted BAM..."
samtools view -bS ${MAPPED_SAM} | samtools sort -o ${SORTED_BAM}

echo "[*] Indexing BAM file..."
samtools index ${SORTED_BAM}

echo "[4/4] Converting BAM to BED12 for FLAIR..."
# The -bed12 flag is mandatory to retain splice junction structures for FLAIR correct
bamToBed -bed12 -i ${SORTED_BAM} > ${OUTPUT_BED}

echo "✔ Script 1 completed successfully for Barcode ${BARCODE}."
