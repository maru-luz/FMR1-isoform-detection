# Long-Read RNA-seq Analysis Pipeline

This repository contains the bioinformatics workflow used for the analysis of Oxford Nanopore Technologies (ONT) long-read RNA-seq data. The pipeline covers read demultiplexing, mapping, isoform identification, differential expression analysis and coding potential prediction.

## Repository Structure

*   `scripts/`: Contains the `.sh` and `.R` scripts used for data processing and statistical analysis.
*   `metadata/`: Contains metadata files required for the analysis (e.g., sample manifest).
*   `data/`: Information regarding the availability of raw and processed datasets.

## Pipeline Overview

> **Note:** The bash commands provided below are representative examples for a single sample (`bc01`). During the study, these steps were executed for all biological replicates.

### 1. Demultiplexing
Raw reads were demultiplexed using [Dorado](https://github.com/nanoporetech/dorado) with the Native Barcoding Kit.

```bash
# Example for barcode 01
dorado demux -o barcode01 --kit-name SQK-NBD114-24 --emit-fastq barcode01.fastq
```

### 2. Genome Indexing, Mapping and Format Conversion
Reads were aligned to the human reference genome (GRCh38, primary assembly) using minimap2 in splice-aware mode. Alignments were sorted and indexed using samtools, and subsequently converted to BED12 format using bedtools.

```bash
# Map reads
minimap2 -ax splice Homo_sapiens.GRCh38.dna.primary_assembly.mmi barcode01.fastq \
    --splice-flank yes -o bc01-mapped.sam

# Convert to BAM, sort, and index
samtools view -bS bc01-mapped.sam | samtools sort -o bc01-mapped.sorted.bam
samtools index bc01-mapped.sorted.bam

# Convert BAM to BED12
bamToBed -bed12 -i bc01-mapped.sorted.bam > bc01.bed12
````

### 3. Read Correction
Splice junctions from the mapped reads were corrected against the reference genome and annotation using FLAIR correct.

````bash
flair correct \
    -q bc01.bed12 \
    -g Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -f Homo_sapiens.GRCh38.115.gtf \
    --output bc01
````

### 4. Isoform Identification
Novel and known isoforms were identified using FLAIR. Reads were corrected prior to the collapse step. All sample-specific flairomes were subsequently merged into a single comprehensive transcriptome.

````bash
# Collapse identical reads into isoforms
flair collapse \
    --gtf Homo_sapiens.GRCh38.115.gtf \
    -g Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -q bc01_all_corrected.bed \
    -r barcode01.fastq \
    --stringent --check_splice --generate_map --output bc01

# Combine all transcriptomes using a manifest file
flair_combine -m metadata/manifest.txt -o combined -p 0 --convert_gtf
````

### 5. Differential Expression Analysis
Statistical analysis was performed in R using the DESeq2 package. The provided R script includes low-count filtering, basic quality control (PCA), and differential expression testing.

Script: scripts/differential_expression.R

Input: Count matrices generated from the combined FLAIR transcriptome.

### 6. Open Reading Frame (ORF) Prediction
To identify candidate coding regions within the identified transcripts, we utilized TransDecoder2. The prediction pipeline incorporates homology searches against the Pfam database to improve accuracy. The workflow consists of:

Extraction of candidate long ORFs.

Domain identification using hmmscan against the Pfam-A database.

Final ORF prediction retaining candidates with biological evidence (Pfam hits).

Script: scripts/orf_prediction.sh

Input: Fasta file of the combined novel and reference isoforms.


