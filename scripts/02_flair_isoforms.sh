#!/bin/bash
# ==============================================================================
# SCRIPT 2: ISOFORM CORRECTION, COLLAPSE, AND COMBINE (FLAIR)
# ==============================================================================
# Description: 
# This script uses FLAIR to correct splice junctions based on genome annotation,
# collapses reads into high-confidence isoforms, and merges sample-specific 
# transcriptomes. 
# ==============================================================================

# ------------------------------------------------------------------------------
# VARIABLES (Adjust paths accordingly)
# ------------------------------------------------------------------------------
BARCODE="01" # Change this per sample for steps 1 and 2

# Reference files (MUST MATCH SCRIPT 1 - Ensembl Release 115)
REF_FASTA="Homo_sapiens.GRCh38.dna.primary_assembly.fa" # Unzipped for FLAIR
REF_GTF="Homo_sapiens.GRCh38.115.gtf"

# Inputs per sample
INPUT_BED12="bc${BARCODE}.bed12"
RAW_FASTQ="barcode${BARCODE}.fastq"

# ------------------------------------------------------------------------------
# PART A: SAMPLE-SPECIFIC PROCESSING (Run for each barcode)
# ------------------------------------------------------------------------------

echo "[1/3] Correcting splice junctions for bc${BARCODE}..."
# Output will be bcXX_all_corrected.bed
flair correct \
    -q ${INPUT_BED12} \
    -g ${REF_FASTA} \
    -f ${REF_GTF} \
    --output bc${BARCODE}

echo "[2/3] Collapsing reads into high-confidence isoforms for bc${BARCODE}..."
flair collapse \
    --gtf ${REF_GTF} \
    -g ${REF_FASTA} \
    -q bc${BARCODE}_all_corrected.bed \
    -r ${RAW_FASTQ} \
    --stringent \
    --check_splice \
    --generate_map \
    --output bc${BARCODE}

# ==============================================================================
# PART B: GLOBAL PROCESSING (Run ONLY ONCE after all samples are collapsed)
# ==============================================================================
# IMPORTANT: Ensure 'manifest.tsv' is properly formatted:
# sample_id    condition    batch    path/to/reads.fastq
# ==============================================================================

# Uncomment the following lines to run global step:

# echo "[3/3] Combining all sample-specific transcriptomes..."
# flair combine \
#     -m manifest.tsv \
#     -o combined \
#     -p 0 \
#     --convert_gtf

echo "✔ Script 2 sample-specific steps completed for Barcode ${BARCODE}."
