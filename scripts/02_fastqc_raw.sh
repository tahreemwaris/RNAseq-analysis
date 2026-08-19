#!/bin/bash

set -euo pipefail

# ============================================================
# Arabidopsis thaliana Heat-Stress RNA-seq
# Step 02: Quality Control of Raw FASTQ Files
# Tool: FastQC
#
# Processes raw paired-end FASTQ files
# ============================================================

PROJECT_DIR="$HOME/Arabidopsis_HeatStress_RNAseq"
FASTQ_DIR="$PROJECT_DIR/data/fastq"
OUTPUT_DIR="$PROJECT_DIR/results/fastqc_raw"
LOG_DIR="$PROJECT_DIR/logs"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

LOG_FILE="$LOG_DIR/02_fastqc_raw.log"

# Save output to terminal and log file
exec > >(tee -a "$LOG_FILE") 2>&1

echo "============================================================"
echo "Arabidopsis Heat-Stress RNA-seq"
echo "Step 02: Raw FASTQ Quality Control"
echo "Started: $(date)"
echo "============================================================"

# ------------------------------------------------------------
# Check FastQC installation
# ------------------------------------------------------------

if ! command -v fastqc >/dev/null 2>&1; then
    echo "ERROR: FastQC is not installed."
    exit 1
fi

echo
echo "FastQC location:"
which fastqc

echo
echo "FastQC version:"
fastqc --version

# ------------------------------------------------------------
# Check FASTQ files
# ------------------------------------------------------------

echo
echo "Checking input FASTQ files..."

FASTQ_FILES=(
    "$FASTQ_DIR/Normal_1_R1.fastq.gz"
    "$FASTQ_DIR/Normal_1_R2.fastq.gz"
    "$FASTQ_DIR/Normal_2_R1.fastq.gz"
    "$FASTQ_DIR/Normal_2_R2.fastq.gz"
    "$FASTQ_DIR/Normal_3_R1.fastq.gz"
    "$FASTQ_DIR/Normal_3_R2.fastq.gz"
    "$FASTQ_DIR/HeatStress_1_R1.fastq.gz"
    "$FASTQ_DIR/HeatStress_1_R2.fastq.gz"
    "$FASTQ_DIR/HeatStress_2_R1.fastq.gz"
    "$FASTQ_DIR/HeatStress_2_R2.fastq.gz"
    "$FASTQ_DIR/HeatStress_3_R1.fastq.gz"
    "$FASTQ_DIR/HeatStress_3_R2.fastq.gz"
)

for FILE in "${FASTQ_FILES[@]}"; do
    if [[ ! -f "$FILE" ]]; then
        echo "ERROR: FASTQ file not found:"
        echo "$FILE"
        exit 1
    fi
done

echo "All 12 FASTQ files found."

echo
echo "Input FASTQ files:"
ls -lh "${FASTQ_FILES[@]}"

# ------------------------------------------------------------
# Run FastQC
# ------------------------------------------------------------

echo
echo "============================================================"
echo "Running FastQC..."
echo "============================================================"

fastqc \
    --threads 2 \
    --outdir "$OUTPUT_DIR" \
    "${FASTQ_FILES[@]}"

# ------------------------------------------------------------
# Verify output
# ------------------------------------------------------------

echo
echo "============================================================"
echo "FastQC completed"
echo "Finished: $(date)"
echo "============================================================"

echo
echo "FastQC output files:"
ls -lh "$OUTPUT_DIR"

echo
echo "Output directory:"
echo "$OUTPUT_DIR"

echo
echo "Log file:"
echo "$LOG_FILE"

echo
echo "Remaining disk space:"
df -h "$PROJECT_DIR"
