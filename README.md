# FMR1 Isoform Detection and Differential Expression Pipeline

## Overview
This repository contains the bioinformatics pipeline and custom scripts used for the characterization of alternative splicing and isoform diversity of the *FMR1* gene using Oxford Nanopore Technologies (ONT) long-read sequencing. 

## Data Availability
Raw long-read sequencing data generated for this study has been deposited in the NCBI Sequence Read Archive (SRA) under BioProject accession **PRJNA1509465** *(Note: Data access will be made public upon peer-reviewed publication)*.

## Experimental Metadata
The metadata required to reproduce the analysis is provided in the `metadata/` directory:
* [Execution Manifest](./metadata/manifest.tsv): Contains sample IDs and FASTQ file paths required for the FLAIR pipeline.
* [Experimental Design](./metadata/deseq2_design.csv): Maps individual barcodes to their corresponding biological groups (B1, B2, and GC) for differential expression analysis.

## Bioinformatics Workflow
The analysis is divided into four main modular scripts. They must be executed sequentially:

1. **Read Mapping & Conversion:** 
   [`01_mapping_and_conversion.sh`](./scripts/01_mapping_and_conversion.sh)
   Performs demultiplexing, aligns ONT long-reads to the Ensembl human reference genome (GRCh38 primary assembly) using Minimap2, and generates BED12 files.

2. **Isoform Assembly & Quantification:** 
   [`02_flair_isoforms.sh`](./scripts/02_flair_isoforms.sh)
   Uses FLAIR to correct splice junctions, collapse reads into high-confidence transcript models, and generate a unified count matrix for the entire dataset.

3. **Differential Transcript Expression (DTE):** 
   [`03_differential_expression.R`](./scripts/03_differential_expression.R)
   R script utilizing DESeq2 to identify differentially expressed *FMR1* isoforms across B1, B2, and GC groups.

4. **ORF Prediction:** 
   [`04_orf_prediction.sh`](./scripts/04_orf_prediction.sh)
   Predicts coding regions from the novel and reference *FMR1* isoforms using TransDecoder2, anchoring predictions to the canonical start codon (ATG) prior to downstream functional annotation.

## Citation
*(Citation details will be updated upon publication of the manuscript).*
