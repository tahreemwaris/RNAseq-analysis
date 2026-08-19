# Arabidopsis Heat Stress RNA-seq Analysis

## Overview

This project implements a reproducible RNA-seq analysis pipeline for studying gene-expression changes in *Arabidopsis thaliana* under heat-stress conditions.

The pipeline starts from publicly available sequencing data and processes the data through quality control, adapter/quality trimming, genome alignment, gene-level read counting, and differential expression analysis using DESeq2.

The overall workflow is:

```text
Raw sequencing data
        ↓
FASTQ conversion
        ↓
Quality Control (FastQC)
        ↓
Read trimming (fastp)
        ↓
Quality Control after trimming
        ↓
MultiQC
        ↓
HISAT2 alignment
        ↓
SAM/BAM processing
        ↓
featureCounts
        ↓
Gene × Sample count matrix
        ↓
DESeq2
        ↓
Differentially Expressed Genes (DEGs)
```

---

## Biological Objective

The main objective is to identify genes whose expression changes in *Arabidopsis thaliana* in response to heat stress.

The differential-expression analysis compares:

* **Control / normal condition**
* **Heat-stress condition**

The resulting differentially expressed genes can subsequently be used for downstream biological interpretation, including:

* Gene Ontology (GO) enrichment
* KEGG pathway analysis
* Functional annotation
* Identification of heat-stress-responsive genes
* Comparison with other stress-response datasets

---

## Dataset

The project uses publicly available *Arabidopsis thaliana* RNA-seq data obtained from the NCBI Sequence Read Archive (SRA).

The current analysis uses biological replicates from control and heat-stress conditions.

The sequencing data are intentionally **not stored in this GitHub repository** because raw FASTQ/SRA files are large and are excluded through `.gitignore`.

---

## Pipeline

### 1. SRA Data Acquisition

SRA accession data are downloaded using the NCBI SRA Toolkit.

Tools:

* `prefetch`
* `fasterq-dump`

Script:

```text
scripts/00_prefetch_heatstress.sh
scripts/01_sra_to_fastq.sh
```

The SRA files are converted into paired-end FASTQ files for downstream analysis.

---

### 2. Raw Read Quality Control

FastQC is used to evaluate the quality of the raw sequencing reads.

Script:

```text
scripts/02_fastqc_raw.sh
```

Important quality metrics include:

* Per-base sequence quality
* Sequence quality scores
* Adapter contamination
* Sequence duplication
* GC content
* Overrepresented sequences

---

### 3. Read Trimming

Low-quality bases and adapter sequences are removed using `fastp`.

Script:

```text
scripts/03_trim_fastp.sh
```

The purpose of trimming is to improve the quality of reads before genome alignment.

---

### 4. Post-trimming Quality Control

FastQC is performed again after trimming.

Script:

```text
scripts/04_fastqc_trimmed.sh
```

The results are compared with the raw-read QC to determine whether trimming improved read quality.

---

### 5. MultiQC

MultiQC is used to combine quality-control reports from multiple samples into a single report.

Script:

```text
scripts/05_multiqc_trimmed.sh
```

This makes it easier to compare sequencing quality across all samples.

---

### 6. Genome Alignment

The cleaned paired-end reads are aligned to the *Arabidopsis thaliana* reference genome using HISAT2.

Script:

```text
scripts/06_hisat2_alignment.sh
```

The alignment produces SAM/BAM files containing the genomic locations of sequencing reads.

HISAT2 indexes are generated locally but are not stored in this repository.

---

### 7. Gene-level Read Counting

After alignment, reads are assigned to genes using `featureCounts`.

Script:

```text
scripts/06_featurecounts.sh
```

The output is a gene-by-sample count matrix.

Conceptually:

```text
              Sample1   Sample2   Sample3   Heat1   Heat2   Heat3
Gene1            120       135       110      450     490     470
Gene2             80        75        90       30      25      35
Gene3            500       520       480      700     680     720
...
```

Rows represent genes and columns represent biological samples.

This count matrix is the primary input for DESeq2.

---

### 8. Differential Expression Analysis

DESeq2 is used to statistically identify genes whose expression differs between control and heat-stress conditions.

Script:

```text
scripts/06_deseq2.R
```

The analysis includes:

* Count normalization
* Estimation of dispersion
* Statistical testing
* Log2 fold-change estimation
* Multiple-testing correction
* Identification of differentially expressed genes

Typical outputs include:

* Differential-expression table
* Log2 fold changes
* Adjusted p-values
* MA plot
* PCA plot
* Volcano plot
* Heatmaps of selected genes

---

## Statistical Interpretation

For a gene, the expression change between heat stress and control is represented using the log2 fold change.

```text
log2FC > 0
    ↓
Higher expression under heat stress
```

```text
log2FC < 0
    ↓
Lower expression under heat stress
```

Statistical significance is assessed using the DESeq2-adjusted p-value.

A typical DEG-selection criterion may be:

```text
|log2FC| ≥ 1
adjusted p-value < 0.05
```

The exact thresholds should be defined according to the objectives of the final analysis.

---

## Repository Structure

```text
RNAseq-analysis/
│
├── README.md
│
├── .gitignore
│
└── scripts/
    │
    ├── 00_prefetch_heatstress.sh
    ├── 01_sra_to_fastq.sh
    ├── 02_fastqc_raw.sh
    ├── 03_trim_fastp.sh
    ├── 04_fastqc_trimmed.sh
    ├── 05_multiqc_trimmed.sh
    ├── 06_hisat2_alignment.sh
    ├── 06_featurecounts.sh
    └── 06_deseq2.R
```

Large sequencing and generated files are excluded from version control.

Excluded file types include:

```text
FASTQ / FASTQ.GZ
SRA
SAM
BAM
BAI
FASTA
GTF
HISAT2 indexes
large analysis outputs
```

---

## Main Software

The pipeline uses the following major tools:

| Tool             | Purpose                          |
| ---------------- | -------------------------------- |
| NCBI SRA Toolkit | Download and convert SRA data    |
| FastQC           | Read-quality assessment          |
| fastp            | Adapter and quality trimming     |
| MultiQC          | Aggregate QC reports             |
| HISAT2           | Genome alignment                 |
| SAMtools         | BAM processing                   |
| featureCounts    | Gene-level read counting         |
| R                | Statistical analysis             |
| DESeq2           | Differential expression analysis |

---

## Reproducibility

The purpose of this repository is to maintain a reproducible record of the computational workflow.

Large sequencing files and generated intermediate files are excluded from GitHub.

The analysis scripts are version-controlled so that changes to the computational pipeline can be tracked over time.

---

## Future Analysis

After differential-expression analysis, the project can be extended to downstream biological interpretation.

Planned analyses include:

```text
Differentially Expressed Genes
          ↓
GO enrichment
          ↓
KEGG pathway analysis
          ↓
Heat-stress response pathways
          ↓
Candidate heat-stress-responsive genes
```

The resulting genes can also be compared with other *Arabidopsis* stress-response datasets to identify conserved stress-responsive genes.

---

## Project Status

Current stage:

```text
✓ SRA data acquisition
✓ FASTQ generation
✓ Raw read QC
✓ Read trimming
✓ Trimmed-read QC
✓ MultiQC
✓ Reference preparation
✓ HISAT2 alignment
✓ Gene-level counting
✓ DESeq2 analysis
```

The pipeline is under active development as downstream analyses are added.

---

## Author

**Tahreem Waris**

Bioinformatics / RNA-seq Analysis

GitHub:

`tahreemwaris`
