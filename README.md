# Somatic TEI Detection Pipeline

## Overview

This pipeline is designed to identify and refine **somatic transposable element insertions (TEIs)** with high confidence. It integrates a series of tools and custom scripts to detect, correct, and annotate somatic structural variants (SVs), particularly large-scale insertions (INSs) and duplications (DUPs), focusing on transposable elements (TEs). The pipeline leverages **Sniffles2**, **Iris**, **Racon**, **RepeatMasker**, **Tandem Repeat Finder**, and **Sdust** to ensure high accuracy and reliability in somatic TEI identification.
## Pipeline Workflow

The pipeline consists of the following key steps:
### Step 1: Merging Alignments and Detect Somatic Large-Scale Insertions

- Merge the alignments of **tumor samples** with their paired **blood samples**.
- Use **Sniffles2** to detect large-scale insertions (INSs) (>50 bp) from the merged alignments.
- A somatic INS is considered valid if:
    - It has **>3 supporting reads** from the tumor sample.
    - It has **no supporting reads** from the paired blood sample.
### Step 2: Additional Screening for High-confidence SVs

- SV positions were restricted to a standard deviation of less than 150 bp
-  Compare alignment-supported lengths with Sniffles-reported lengths for INSs and DUPs:
    - Retain SVs within:
        ```
        A.mean(Lengtha) - 2 * STRI(Lengtha) <= Lengths <= mean(Lengtha) + 2 * STRI(Lengtha)
        B.Q1 - 1.5 * IQR <= Length <= Q3 + 1.5 * IQR
        ```
- Assess SVs in repeat regions for significant length differences between tumor and blood alignments, retaining only those exhibiting clear distinctions.
### Step 3: SV Correction And Annotation

- Use **Iris** (version 1.0.4) to correct the positions and sequences of somatic SVs (INSs and DUPs)
- Use **RepeatMasker** (version 4.1.2) to annotate corrected insertion sequences
- Re-annotate sequences by RepeatMasker using **Tandem Repeat Finder** (TRF, version 4.10.0) and **Sdust** (version 0.1). This part of the script is in a folder called script.
### Step 4: Polish Reference Genome and Realign Reads

- To reduce potential false positives caused by differences between the reference and individual genomes:
    1. Polish the reference genome at each insertion locus using **Racon**, based on alignments from the paired blood samples.
    2. Realign tumor reads to the polished reference genome.
    3. Retain somatic INSs after realignment as candidate somatic INSs.
### Step 5: Classify Somatic TEIs
-  Base on Step 3 and Step 4 classify somatic TEIs based on the following criteria:
    - If the 50 bp flanking sequences and the insertion sequence are not annotated as the same type of repeat.
    - The insertion carries a transposable element (TE) sequence.
## Installation
### Prerequisites
Ensure the following tools are installed:
- **Samtools** (version 2.8.0)
- **minimap2** (version 2.2.17)
- **Sniffles2** (version 2.0.6)
- **Iris** (version 1.0.4)
- **Racon** (version 1.4.16)
- **RepeatMasker** (version 4.1.2)
- **Tandem Repeat Finder** (version 4.10.0)
- **Sdust** (version 0.1)
- **Snakemake** (version 3.13.3)

## Usage

### Input Data
- Tumor sample fastq filte (fastq files).
- Paired blood sample fastq filte (fastq files).
- Reference genome (FASTA format).
### Output

- **Sample_merge_minimap2_sniffles_v2_tumor_somatic.bed**: This file contains the **initial set of somatic large-scale insertions (INSs)** detected in Step 1 of the pipeline.
- **Sample_iris_out.vcf**: This file contains high-confidence insertion (INS) calls refined by **Iris** in Step 3 of the pipeline.
- **all.INS.sdust.trf.replaced.cor.type.TE_TD_de_novo_type.tsv**: Annotation results,Includes RepeatMasker, TRF, and Sdust annotations in Step 3 of the pipeline.
- **Sample_tumor_somatic_TE.bed**: A list of validated somatic TEIs in BED format in Step 5 of the pipeline.


