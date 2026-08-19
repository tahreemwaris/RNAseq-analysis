#!/bin/bash

set -euo pipefail

# ============================================================
# Arabidopsis thaliana Heat-Stress RNA-seq
# Step 03: Adapter and Quality Trimming
# Tool: fastp
#
# Input:
#   data/fastq/*_R1.fastq.gz
#   data/fastq/*_R2.fastq.gz
#
# Output:
#   results/trimmed/*_R1.trimmed.fastq.gz
#   results/trimmed/*_R2.trimmed.fastq.gz
#   results/trimmed/*.html
#   results/trimmed/*.json
# ============================================================

PROJECT_DIR="$HOME/Arabidopsis_HeatStress_RNAseq"
INPUT_DIR="$PROJECT_DIR/data/fastq"
OUTPUT_DIR="$PROJECT_DIR/results/trimmed"
LOG_DIR="$PROJECT_DIR/logs"

mkdir -p "$OUTPUT_DIR" "$LOG_DIR"

LOG_FILE="$LOG_DIR/03_trim_fastp.log"

exec > >(tee -a "$LOG_FILE") 2>&1

echo "============================================================"
echo "Arabidopsis Heat-Stress RNA-seq"
echo "Step 03: fastp trimming"
echo "Started: $(date)"
echo "============================================================"

# ------------------------------------------------------------
# Check fastp
# ------------------------------------------------------------

if ! command -v fastp >/dev/null 2>&1; then
    echo "ERROR: fastp is not installed."
    exit 1
fi

echo
echo "fastp:"
which fastp

echo
fastp --version

# ------------------------------------------------------------
# Samples in fixed order
# ------------------------------------------------------------

SAMPLES=(
    Normal_1
    Normal_2
    Normal_3
    HeatStress_1
    HeatStress_2
    HeatStress_3
)

# ------------------------------------------------------------
# Process one sample at a time
# ------------------------------------------------------------

for SAMPLE in "${SAMPLES[@]}"; do

    R1="$INPUT_DIR/${SAMPLE}_R1.fastq.gz"
    R2="$INPUT_DIR/${SAMPLE}_R2.fastq.gz"

    OUT_R1="$OUTPUT_DIR/${SAMPLE}_R1.trimmed.fastq.gz"
    OUT_R2="$OUTPUT_DIR/${SAMPLE}_R2.trimmed.fastq.gz"

    HTML="$OUTPUT_DIR/${SAMPLE}_fastp.html"
    JSON="$OUTPUT_DIR/${SAMPLE}_fastp.json"

    echo
    echo "============================================================"
    echo "Processing: $SAMPLE"
    echo "Started: $(date)"
    echo "============================================================"

    # --------------------------------------------------------
    # Check input files
    # --------------------------------------------------------

    if [[ ! -f "$R1" ]]; then
        echo "ERROR: Missing input file:"
        echo "$R1"
        exit 1
    fi

    if [[ ! -f "$R2" ]]; then
        echo "ERROR: Missing input file:"
        echo "$R2"
        exit 1
    fi

    # --------------------------------------------------------
    # Skip completed sample
    # --------------------------------------------------------

    if [[ -f "$OUT_R1" &&
          -f "$OUT_R2" &&
          -f "$HTML" &&
          -f "$JSON" ]]; then

        echo "$SAMPLE already trimmed. Skipping."
        continue
    fi

    # --------------------------------------------------------
    # Run fastp
    #
    # --detect_adapter_for_pe
    #     Automatically detects adapters in paired-end data
    #
    # --qualified_quality_phred 20
    #     Bases below Q20 are considered low quality
    #
    # --length_required 50
    #     Discard reads shorter than 50 bp after trimming
    #
    # --thread 2
    #     Conservative CPU usage
    # --------------------------------------------------------

    fastp \
        --in1 "$R1" \
        --in2 "$R2" \
        --out1 "$OUT_R1" \
        --out2 "$OUT_R2" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --thread 2 \
        --html "$HTML" \
        --json "$JSON" \
        --report_title "fastp - $SAMPLE"

    # --------------------------------------------------------
    # Verify outputs
    # --------------------------------------------------------

    if [[ ! -s "$OUT_R1" || ! -s "$OUT_R2" ]]; then
        echo "ERROR: Trimmed FASTQ output missing for $SAMPLE."
        exit 1
    fi

    echo
    echo "Output:"
    ls -lh "$OUT_R1" "$OUT_R2"

    echo
    echo "fastp reports:"
    ls -lh "$HTML" "$JSON"

    echo
    echo "Remaining storage:"
    df -h "$PROJECT_DIR"

    echo
    echo "Finished: $SAMPLE"
    echo "Time: $(date)"

done

echo
echo "============================================================"
echo "Trimming completed for all six samples"
echo "Finished: $(date)"
echo "============================================================"

echo
echo "Trimmed FASTQ files:"
ls -lh "$OUTPUT_DIR"/*.trimmed.fastq.gz

echo
echo "Total trimmed data:"
du -sh "$OUTPUT_DIR"

echo
echo "Final storage:"
df -h "$PROJECT_DIR"

echo
echo "Log:"
echo "$LOG_FILE"

echo
echo "============================================================"
