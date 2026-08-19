#!/bin/bash

set -euo pipefail

# ============================================================
# HISAT2 Alignment + Sorted BAM + BAM Index
# Arabidopsis Heat Stress RNA-seq
# ============================================================

PROJECT="$HOME/Arabidopsis_HeatStress_RNAseq"

TRIMMED="$PROJECT/results/trimmed"
INDEX="$PROJECT/data/reference/hisat2_index/TAIR10"
OUTDIR="$PROJECT/results/alignment"
LOGDIR="$PROJECT/logs"

# Use exactly 2 threads
THREADS=2

mkdir -p "$OUTDIR"
mkdir -p "$LOGDIR"

echo "============================================================"
echo "HISAT2 Alignment Pipeline"
echo "Started: $(date)"
echo "Project: $PROJECT"
echo "Threads: $THREADS"
echo "============================================================"

for R1 in "$TRIMMED"/*_R1.trimmed.fastq.gz
do

    SAMPLE=$(basename "$R1" _R1.trimmed.fastq.gz)

    R2="$TRIMMED/${SAMPLE}_R2.trimmed.fastq.gz"

    BAM="$OUTDIR/${SAMPLE}.sorted.bam"
    BAI="$OUTDIR/${SAMPLE}.sorted.bam.bai"
    LOG="$LOGDIR/${SAMPLE}_hisat2.log"

    echo
    echo "------------------------------------------------------------"
    echo "Sample: $SAMPLE"
    echo "------------------------------------------------------------"

    # --------------------------------------------------------
    # Check R2 exists
    # --------------------------------------------------------

    if [[ ! -f "$R2" ]]; then
        echo "ERROR: R2 file missing for $SAMPLE"
        echo "Expected: $R2"
        continue
    fi

    # --------------------------------------------------------
    # Skip completed samples
    # --------------------------------------------------------

    if [[ -f "$BAM" && -f "$BAI" ]]; then
        echo "SKIPPING $SAMPLE"
        echo "Sorted BAM already exists:"
        echo "$BAM"
        echo "BAM index already exists:"
        echo "$BAI"
        continue
    fi

    # --------------------------------------------------------
    # Remove incomplete BAM if present
    # --------------------------------------------------------

    if [[ -f "$BAM" && ! -f "$BAI" ]]; then
        echo "WARNING: BAM exists but index is missing."
        echo "Creating BAM index..."
        samtools index -@ "$THREADS" "$BAM"
        echo "Index created."
        continue
    fi

    echo "R1: $R1"
    echo "R2: $R2"
    echo "Output BAM: $BAM"
    echo "Log: $LOG"
    echo "Started: $(date)"

    # --------------------------------------------------------
    # HISAT2 → SAMtools sort → sorted BAM
    # --------------------------------------------------------

    hisat2 \
        -p "$THREADS" \
        -x "$INDEX" \
        -1 "$R1" \
        -2 "$R2" \
        2> "$LOG" \
    | samtools sort \
        -@ "$THREADS" \
        -o "$BAM" -

    # --------------------------------------------------------
    # Create BAM index
    # --------------------------------------------------------

    samtools index \
        -@ "$THREADS" \
        "$BAM"

    # --------------------------------------------------------
    # Verify output
    # --------------------------------------------------------

    if [[ -f "$BAM" && -f "$BAI" ]]; then

        echo "SUCCESS: $SAMPLE"
        echo "BAM: $BAM"
        echo "BAI: $BAI"

        echo
        echo "BAM statistics:"
        samtools quickcheck "$BAM"

        if [[ $? -eq 0 ]]; then
            echo "BAM quickcheck: PASSED"
        else
            echo "BAM quickcheck: FAILED"
            exit 1
        fi

    else

        echo "ERROR: Alignment failed for $SAMPLE"
        exit 1

    fi

    echo "Finished: $(date)"

done

echo
echo "============================================================"
echo "ALL ALIGNMENT JOBS COMPLETED"
echo "Finished: $(date)"
echo "============================================================"
