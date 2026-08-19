#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# FeatureCounts - Arabidopsis thaliana RNA-seq
# ============================================================

PROJECT_DIR="$HOME/Arabidopsis_HeatStress_RNAseq"

BAM_DIR="$PROJECT_DIR/results/alignment"
GTF="$PROJECT_DIR/data/reference/Arabidopsis_thaliana.TAIR10.63.gtf"
OUT_DIR="$PROJECT_DIR/results/counts"
LOG_DIR="$PROJECT_DIR/logs"

# Exactly 2 threads
THREADS=2

# Output files
COUNTS="$OUT_DIR/gene_counts.txt"
SUMMARY="$OUT_DIR/gene_counts.txt.summary"
LOG="$LOG_DIR/06_featurecounts.log"

# ------------------------------------------------------------
# Create required directories
# ------------------------------------------------------------

mkdir -p "$OUT_DIR"
mkdir -p "$LOG_DIR"

# ------------------------------------------------------------
# Logging
# ------------------------------------------------------------

exec > >(tee -a "$LOG") 2>&1

echo "============================================================"
echo "FeatureCounts"
echo "Started: $(date)"
echo "============================================================"

# ------------------------------------------------------------
# Check featureCounts
# ------------------------------------------------------------

if ! command -v featureCounts >/dev/null 2>&1; then
    echo "ERROR: featureCounts is not installed."
    exit 1
fi

echo "featureCounts version:"
featureCounts -v

# ------------------------------------------------------------
# Check GTF
# ------------------------------------------------------------

if [[ ! -f "$GTF" ]]; then
    echo "ERROR: GTF file not found:"
    echo "$GTF"
    exit 1
fi

echo "GTF:"
echo "$GTF"

# ------------------------------------------------------------
# Skip if already completed
# ------------------------------------------------------------

if [[ -s "$COUNTS" && -s "$SUMMARY" ]]; then
    echo ""
    echo "FeatureCounts output already exists."
    echo "Skipping featureCounts."
    echo ""
    echo "Output:"
    echo "$COUNTS"
    echo "$SUMMARY"
    echo ""
    echo "Finished: $(date)"
    exit 0
fi

# ------------------------------------------------------------
# Check BAM files
# ------------------------------------------------------------

SAMPLES=(
    "Normal_1"
    "Normal_2"
    "Normal_3"
    "HeatStress_1"
    "HeatStress_2"
    "HeatStress_3"
)

BAMS=()

for SAMPLE in "${SAMPLES[@]}"; do

    BAM="$BAM_DIR/${SAMPLE}.sorted.bam"
    BAI="$BAM.bai"

    if [[ ! -s "$BAM" ]]; then
        echo "ERROR: BAM file missing:"
        echo "$BAM"
        exit 1
    fi

    if [[ ! -s "$BAI" ]]; then
        echo "ERROR: BAM index missing:"
        echo "$BAI"
        exit 1
    fi

    BAMS+=("$BAM")

done

echo ""
echo "BAM files to process:"
printf '  %s\n' "${BAMS[@]}"

# ------------------------------------------------------------
# Run featureCounts
# ------------------------------------------------------------

echo ""
echo "Running featureCounts..."
echo "Threads: $THREADS"
echo "Paired-end mode: YES"
echo ""

featureCounts \
    -T "$THREADS" \
    -p \
    -B \
    -C \
    -s 0 \
    -t exon \
    -g gene_id \
    -a "$GTF" \
    -o "$COUNTS" \
    "${BAMS[@]}"

# ------------------------------------------------------------
# Verify output
# ------------------------------------------------------------

if [[ ! -s "$COUNTS" ]]; then
    echo "ERROR: featureCounts did not produce the expected output."
    exit 1
fi

if [[ ! -s "$SUMMARY" ]]; then
    echo "ERROR: featureCounts summary file was not created."
    exit 1
fi

# ------------------------------------------------------------
# Completion
# ------------------------------------------------------------

echo ""
echo "============================================================"
echo "FeatureCounts completed successfully."
echo "Finished: $(date)"
echo ""
echo "Count matrix:"
echo "$COUNTS"
echo ""
echo "Summary:"
echo "$SUMMARY"
echo "============================================================"
