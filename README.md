# *FMR1* Isoform Detection and Differential Expression Pipeline

## Overview
This repository contains the bioinformatics pipeline and custom scripts used for the characterization of alternative splicing and isoform diversity of the *FMR1* gene using Oxford Nanopore Technologies (ONT) long-read sequencing.

## Data Availability
Raw long-read sequencing data generated for this study has been deposited in the NCBI Sequence Read Archive (SRA) under BioProject accession **PRJNA1509465** *(Note: Data access will be made public upon peer-reviewed publication)*.

## Experimental Metadata
The metadata required to reproduce the analysis is provided in the `metadata/` directory:
* [Execution Manifest](./metadata/manifest.tsv): Contains sample IDs and FASTQ file paths required for the FLAIR pipeline.
* [Experimental Design](./metadata/deseq2_design.csv): Maps individual barcodes to their corresponding biological groups for differential expression analysis.

## Bioinformatics Workflow

**Script 1 — Read Mapping & Conversion**
[`01_mapping_and_conversion.sh`](./scripts/01_mapping_and_conversion.sh)
Performs demultiplexing, aligns ONT long-reads to a reference genome using Minimap2, and generates BED12 files. Run once per barcode by setting the `BARCODE` variable at the top of the script.

**Script 2 — Isoform Assembly & Quantification**
[`02_flair_isoforms.sh`](./scripts/02_flair_isoforms.sh)
Uses FLAIR to correct splice junctions against the reference annotation, collapse reads into high-confidence transcript models, and merge all sample-specific transcriptomes into a unified isoform set with count matrix. Run the sample-specific steps (Part A) once per barcode; then uncomment and run the global `flair combine` step (Part B) once after all barcodes have been processed.

**Script 3 — Differential Transcript Expression (DTE)**
[`03_differential_expression.R`](./scripts/03_differential_expression.R)
R script utilizing DESeq2 to identify differentially expressed *FMR1* isoforms across biological groups. Requires the filtered count matrix (`filtered_isoforms.tsv`) and the experimental design file (`metadata/deseq2_design.csv`).

**Script 4 — ORF Prediction**
[`04_orf_prediction.sh`](./scripts/04_orf_prediction.sh)
Predicts coding regions from *FMR1* isoforms using TransDecoder2, anchoring predictions to the canonical start codon (ATG) prior to downstream functional annotation.

## Study-Specific Parameters
The following parameters are specific to the study described in the associated manuscript:

* **Reference genome:** Ensembl GRCh38 primary assembly, release 115
* **Biological groups:** Blood pre-stimulation (B1), Blood post-stimulation (B2), and Granulosa Cells (GC)
* **Outlier samples excluded:** barcodes 02, 03, and 09 (identified by PCA-based unsupervised clustering)

## Software Requirements
| Tool | Version |
|------|---------|
| Dorado | v0.9.1 |
| Minimap2 | v2.24-r1122 |
| SAMtools | v1.18 |
| BEDtools | v2.31.1 |
| FLAIR correct/collapse | v2.0.0 |
| FLAIR combine | v2.2.0 |
| R | 4.6.0 |
| DESeq2 | 1.52.0 |
| TransDecoder | v2 |

## Citation
*(Citation details will be updated upon publication of the manuscript).*
