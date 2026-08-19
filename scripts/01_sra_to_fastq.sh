#!/bin/bash

set -euo pipefail

# ============================================================
# Arabidopsis thaliana Heat-Stress RNA-seq
# Step 01: SRA → paired-end FASTQ
#
# Processes ONE sample at a time to reduce resource usage.
# ============================================================

PROJECT_DIR="$HOME/Arabidopsis_HeatStress_RNAseq"
SRA_DIR="$PROJECT_DIR/data/sra"
FASTQ_DIR="$PROJECT_DIR/data/fastq"
LOG_DIR="$PROJECT_DIR/logs"

mkdir -p "$FASTQ_DIR" "$LOG_DIR"

LOG_FILE="$LOG_DIR/01_sra_to_fastq.log"

exec > >(tee -a "$LOG_FILE") 2>&1

echo "============================================================"
echo "SRA → FASTQ conversion"
echo "Started: $(date)"
echo "============================================================"

# ------------------------------------------------------------
# Check software
# ------------------------------------------------------------

if ! command -v fasterq-dump >/dev/null 2>&1; then
    echo "ERROR: fasterq-dump is not installed."
    exit 1
fi

if ! command -v gzip >/dev/null 2>&1; then
    echo "ERROR: gzip is not installed."
    exit 1
fi

echo "fasterq-dump: $(which fasterq-dump)"
fasterq-dump --version

# ------------------------------------------------------------
# Sample mapping
# ------------------------------------------------------------

declare -A SAMPLES

SAMPLES["SRR25485412"]="Normal_1"
SAMPLES["SRR25485411"]="Normal_2"
SAMPLES["SRR25485410"]="Normal_3"

SAMPLES["SRR25485409"]="HeatStress_1"
SAMPLES["SRR25485408"]="HeatStress_2"
SAMPLES["SRR25485407"]="HeatStress_3"

# ------------------------------------------------------------
# Check storage
# ------------------------------------------------------------

echo
echo "Available storage:"
df -h "$PROJECT_DIR"

# ------------------------------------------------------------
# Process one sample at a time
# ------------------------------------------------------------

for SRA in "${!SAMPLES[@]}"; do

    SAMPLE="${SAMPLES[$SRA]}"

    echo
    echo "============================================================"
    echo "Processing: $SAMPLE"
    echo "SRA: $SRA"
    echo "Started: $(date)"
    echo "============================================================"

    SRA_PATH="$SRA_DIR/$SRA/$SRA.sra"

    if [[ ! -f "$SRA_PATH" ]]; then
        echo "ERROR: SRA file not found:"
        echo "$SRA_PATH"
        exit 1
    fi

    # --------------------------------------------------------
    # Skip already-compressed samples
    # --------------------------------------------------------

    if [[ -f "$FASTQ_DIR/${SAMPLE}_R1.fastq.gz" &&
          -f "$FASTQ_DIR/${SAMPLE}_R2.fastq.gz" ]]; then

        echo "$SAMPLE already converted. Skipping."
        continue
    fi

    # --------------------------------------------------------
    # Convert SRA → FASTQ
    # --------------------------------------------------------

    echo "Converting $SRA to FASTQ..."

    fasterq-dump \
        --split-files \
        --threads 2 \
        --outdir "$FASTQ_DIR" \
        "$SRA_PATH"

    # --------------------------------------------------------
    # Rename FASTQ files
    # --------------------------------------------------------

    if [[ ! -f "$FASTQ_DIR/${SRA}_1.fastq" ||
          ! -f "$FASTQ_DIR/${SRA}_2.fastq" ]]; then

        echo "ERROR: Expected paired FASTQ files were not produced."
        exit 1
    fi

    mv "$FASTQ_DIR/${SRA}_1.fastq" \
       "$FASTQ_DIR/${SAMPLE}_R1.fastq"

    mv "$FASTQ_DIR/${SRA}_2.fastq" \
       "$FASTQ_DIR/${SAMPLE}_R2.fastq"

    # --------------------------------------------------------
    # Compress FASTQ files
    # --------------------------------------------------------

    echo "Compressing R1..."

    gzip -1 "$FASTQ_DIR/${SAMPLE}_R1.fastq"

    echo "Compressing R2..."

    gzip -1 "$FASTQ_DIR/${SAMPLE}_R2.fastq"

    # --------------------------------------------------------
    # Verify output
    # --------------------------------------------------------

    if [[ ! -f "$FASTQ_DIR/${SAMPLE}_R1.fastq.gz" ||
          ! -f "$FASTQ_DIR/${SAMPLE}_R2.fastq.gz" ]]; then

        echo "ERROR: FASTQ compression failed."
        exit 1
    fi

    echo
    echo "FASTQ files created:"
    ls -lh \
        "$FASTQ_DIR/${SAMPLE}_R1.fastq.gz" \
        "$FASTQ_DIR/${SAMPLE}_R2.fastq.gz"

    echo
    echo "Remaining storage:"
    df -h "$PROJECT_DIR"

    echo
    echo "Finished: $SAMPLE"
    echo "Time: $(date)"

done

echo
echo "============================================================"
echo "All six samples converted successfully."
echo "Finished: $(date)"
echo "============================================================"

echo
echo "FASTQ files:"
ls -lh "$FASTQ_DIR"/*.fastq.gz

echo
echo "Total FASTQ storage:"
du -sh "$FASTQ_DIR"

echo
echo "Final storage:"
df -h "$PROJECT_DIR"

echo
echo "Log:"
echo "$LOG_FILE"
