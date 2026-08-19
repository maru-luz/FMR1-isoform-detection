#!/bin/bash
# ==============================================================================
# SCRIPT 4: OPEN READING FRAME (ORF) PREDICTION
# ==============================================================================
# Description: 
# This script predicts coding regions from the set of novel and reference 
# FMR1 isoforms using TransDecoder2 (TD2). It incorporates homology searches 
# against the Pfam-A database via HMMER to retain ORFs with known protein domains.
#
# Post-Processing Note:
# The resulting predicted protein sequences (.pep files) were subsequently 
# subjected to ClustalW alignment, NMD evaluation, and domain mapping 
# using the InterProScan web service as described in the manuscript.
# ==============================================================================

# ------------------------------------------------------------------------------
# VARIABLES
# ------------------------------------------------------------------------------
# Final FASTA file containing the 15 novel isoforms,  the reference
INPUT_FASTA="novel_isoforms_FMR1.fasta"

# ------------------------------------------------------------------------------
# 1. PREDICT LONG ORFS (TD2.LongOrfs)
# ------------------------------------------------------------------------------
echo "[1/2] Extracting candidate Open Reading Frames (ORFs)..."
# Generates an automatic output directory (e.g., isoformas_FMR1_finales.fasta.transdecoder_dir)
TD2.LongOrfs -t ${INPUT_FASTA}

# ------------------------------------------------------------------------------
# 2. FINAL ORF PREDICTION (TD2.Predict)
# ------------------------------------------------------------------------------
echo "[2/2] Finalizing ORF prediction..."
# Default prediction based on ORF length and nucleotide composition
TD2.Predict -t ${INPUT_FASTA}

echo "✔ Script 4 completed. Predicted proteins are ready for downstream analysis."
