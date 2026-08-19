#!/bin/bash

set -euo pipefail

# ============================================================
# Arabidopsis thaliana Heat-Stress RNA-seq
# Step 00: Download SRA files using prefetch
#
# Memory-efficient:
#   - Downloads ONE accession at a time
#   - Does not load FASTQ data into RAM
#   - Stops immediately if a download fails
# ============================================================

PROJECT_DIR="$HOME/Arabidopsis_HeatStress_RNAseq"
SRA_DIR="$PROJECT_DIR/data/sra"
LOG_DIR="$PROJECT_DIR/logs"

mkdir -p "$SRA_DIR" "$LOG_DIR"

LOG_FILE="$LOG_DIR/00_prefetch_heatstress.log"

# ------------------------------------------------------------
# SRA accessions
# ------------------------------------------------------------

SAMPLES=(
    SRR25485412
    SRR25485411
    SRR25485410
    SRR25485409
    SRR25485408
    SRR25485407
)

# ------------------------------------------------------------
# Start log
# ------------------------------------------------------------

exec > >(tee -a "$LOG_FILE") 2>&1

echo "============================================================"
echo "Arabidopsis Heat-Stress RNA-seq"
echo "SRA download started"
echo "Started: $(date)"
echo "============================================================"

echo
echo "Project directory:"
echo "$PROJECT_DIR"

echo
echo "SRA directory:"
echo "$SRA_DIR"

echo
echo "Available disk space:"
df -h "$PROJECT_DIR"

echo
echo "Number of samples: ${#SAMPLES[@]}"

# ------------------------------------------------------------
# Check prefetch
# ------------------------------------------------------------

if ! command -v prefetch >/dev/null 2>&1; then
    echo "ERROR: prefetch is not installed."
    exit 1
fi

echo
echo "prefetch:"
which prefetch

# ------------------------------------------------------------
# Download one sample at a time
# ------------------------------------------------------------

for SRA in "${SAMPLES[@]}"; do

    echo
    echo "------------------------------------------------------------"
    echo "Downloading: $SRA"
    echo "Time: $(date)"
    echo "------------------------------------------------------------"

    # Skip if the SRA file already exists
    if find "$SRA_DIR/$SRA" -type f -name "$SRA.sra" -print -quit 2>/dev/null | grep -q .; then
        echo "$SRA already exists. Skipping download."
        continue
    fi

    # Download one accession
    prefetch \
        --output-directory "$SRA_DIR" \
        --max-size 100G \
        "$SRA"

    echo
    echo "Finished: $SRA"
    echo "Time: $(date)"

    echo
    echo "Remaining disk space:"
    df -h "$PROJECT_DIR"

done

# ------------------------------------------------------------
# Final summary
# ------------------------------------------------------------

echo
echo "============================================================"
echo "All SRA downloads completed"
echo "Finished: $(date)"
echo "============================================================"

echo
echo "Downloaded SRA files:"
find "$SRA_DIR" -type f -name "*.sra" -exec ls -lh {} \;

echo
echo "Final disk space:"
df -h "$PROJECT_DIR"

echo
echo "Log:"
echo "$LOG_FILE"

echo
echo "============================================================"
