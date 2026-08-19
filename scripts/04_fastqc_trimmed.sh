#!/bin/bash

set -euo pipefail

PROJECT_DIR="$HOME/Arabidopsis_HeatStress_RNAseq"
TRIMMED_DIR="$PROJECT_DIR/results/trimmed"
FASTQC_DIR="$PROJECT_DIR/results/fastqc_trimmed"
LOG_DIR="$PROJECT_DIR/logs"

mkdir -p "$FASTQC_DIR"
mkdir -p "$LOG_DIR"

echo "============================================================"
echo "        FASTQC QUALITY CONTROL - TRIMMED READS"
echo "============================================================"
echo "Start time: $(date)"
echo

echo "Input directory:"
echo "$TRIMMED_DIR"
echo

echo "Output directory:"
echo "$FASTQC_DIR"
echo

echo "Running FastQC..."

fastqc \
    "$TRIMMED_DIR"/*.trimmed.fastq.gz \
    --outdir "$FASTQC_DIR" \
    --threads 4

echo
echo "FastQC completed successfully."
echo "End time: $(date)"
echo
echo "Reports saved in:"
echo "$FASTQC_DIR"
echo "============================================================"
