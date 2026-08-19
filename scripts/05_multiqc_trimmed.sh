#!/bin/bash

set -e

PROJECT_DIR="/home/tahreem/Arabidopsis_HeatStress_RNAseq"

INPUT_DIR="$PROJECT_DIR/results/fastqc_trimmed"
OUTPUT_DIR="$PROJECT_DIR/results/multiqc_trimmed"

MULTIQC="$HOME/multiqc_env/bin/multiqc"

echo "============================================================"
echo "        MULTIQC - TRIMMED READS"
echo "============================================================"
echo "Start time: $(date)"
echo
echo "Input:  $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo

# Check MultiQC
if [ ! -x "$MULTIQC" ]; then
    echo "ERROR: MultiQC 1.35 not found at:"
    echo "$MULTIQC"
    exit 1
fi

echo "Using:"
"$MULTIQC" --version
echo

# Check input directory
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory does not exist:"
    echo "$INPUT_DIR"
    exit 1
fi

# Skip if MultiQC report already exists
if [ -f "$OUTPUT_DIR/multiqc_report.html" ]; then
    echo "MultiQC report already exists."
    echo "Skipping MultiQC analysis."
    echo
    echo "Existing report:"
    echo "$OUTPUT_DIR/multiqc_report.html"
    echo
    echo "End time: $(date)"
    exit 0
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "Running MultiQC..."
echo

"$MULTIQC" \
    "$INPUT_DIR" \
    --outdir "$OUTPUT_DIR" \
    --force

echo
echo "============================================================"
echo "MultiQC completed successfully."
echo "Report:"
echo "$OUTPUT_DIR/multiqc_report.html"
echo "End time: $(date)"
echo "============================================================"
