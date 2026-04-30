#!/bin/bash
# run_saw_batch.sh - Run SAW pipeline across multiple samples
# ============================================================
# Designed for autonomous overnight execution of STOmics SAW
# pipeline to generate cellbin data for chips that need it.
#
# Usage:
#   1. Edit CONFIGURATION section below
#   2. Populate SAMPLES array (after running saw_data_audit.py)
#   3. Run in tmux/screen:
#      tmux new -s saw_batch
#      bash run_saw_batch.sh 2>&1 | tee saw_batch_full.log
#
# ============================================================

set -euo pipefail

# ============================================================
# CONFIGURATION - Edit these paths before running
# ============================================================
SAW_SIF="/opt/stomics/SAW_7.1.sif"
REF_INDEX="/data/references/STAR_SJ100"
ANNOTATION="/data/references/Homo_sapiens.GRCh38.109.gtf"
SPECIES="Homo_sapiens"
THREADS=16
SPLIT_COUNT=1    # 1 for Q40 data, 16 for Q4 data
OUTPUT_BASE="/data/saw_output"
LOG_DIR="${OUTPUT_BASE}/logs"

# Singularity bind mounts (comma-separated source:dest pairs)
# Adjust if your data is on network drives or non-standard locations
BIND_MOUNTS="/mnt/t:/mnt/t,/data:/data"

mkdir -p "${LOG_DIR}"

# ============================================================
# SAMPLE DEFINITIONS
# ============================================================
# Populate after running saw_data_audit.py
# Format: "SAMPLE_ID|CHIP_ID|FQ1|FQ2|MASK|TISSUE|IPR|TAR"
# - FQ1/FQ2: comma-separated for multi-lane data
# - IPR and TAR: use "NONE" if not available
# - TISSUE: endometrium, ovary, etc.
#
# Example:
# SAMPLES=(
#   "SO1|D03453A6|/mnt/t/Sammy Data/FASTQ/D03453A6_R1.fq.gz|/mnt/t/Sammy Data/FASTQ/D03453A6_R2.fq.gz|/mnt/t/Sammy Data/masks/D03453A6.h5|ovary|/mnt/t/Sammy Data/D03453A6/D03453A6.ipr|/mnt/t/Sammy Data/D03453A6/D03453A6_image.tar.gz"
# )

SAMPLES=(
    # ---- UNCOMMENT AND FILL IN AFTER DATA AUDIT ----
    # "SO1|D03453A6|...|...|...|ovary|...|..."
)

# ============================================================
# FUNCTIONS
# ============================================================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

run_saw_sample() {
    local SAMPLE="$1" CHIP="$2" FQ1="$3" FQ2="$4" MASK="$5" TISSUE="$6" IPR="$7" TAR="$8"

    local OUTDIR="${OUTPUT_BASE}/${CHIP}"
    local LOGFILE="${LOG_DIR}/${CHIP}_$(date +%Y%m%d_%H%M%S).log"

    log "Processing ${SAMPLE} (Chip: ${CHIP})"

    # Skip if cellbin already exists in output
    if [ -f "${OUTDIR}/03.ssDNA_analysis/${CHIP}.cellbin.gef" ]; then
        log "  SKIP: cellbin.gef already exists at ${OUTDIR}"
        return 2  # return code 2 = skipped
    fi

    # Validate inputs
    local MISSING=""
    [ ! -f "$FQ1" ] && MISSING="${MISSING} FQ1(${FQ1})"
    [ ! -f "$FQ2" ] && MISSING="${MISSING} FQ2(${FQ2})"
    [ ! -f "$MASK" ] && MISSING="${MISSING} MASK(${MASK})"

    if [ -n "$MISSING" ]; then
        log "  ERROR: Missing inputs:${MISSING}"
        echo "$(date) | FAILED | ${SAMPLE} | ${CHIP} | Missing:${MISSING}" >> "${LOG_DIR}/failures.log"
        return 1
    fi

    mkdir -p "${OUTDIR}"

    # Build SAW command
    local SAW_CMD="singularity exec"

    # Add bind mounts if configured
    if [ -n "${BIND_MOUNTS}" ]; then
        SAW_CMD="${SAW_CMD} --bind ${BIND_MOUNTS}"
    fi

    SAW_CMD="${SAW_CMD} ${SAW_SIF} stereoPipeline.sh \
        -splitCount ${SPLIT_COUNT} \
        -maskFile ${MASK} \
        -fq1 ${FQ1} \
        -fq2 ${FQ2} \
        -refIndex ${REF_INDEX} \
        -annotationFile ${ANNOTATION} \
        -speciesName ${SPECIES} \
        -tissueType ${TISSUE} \
        -outDir ${OUTDIR} \
        -threads ${THREADS} \
        -doCellBin Y \
        -rRNAremove Y"

    # Add image files if available
    if [ "$IPR" != "NONE" ] && [ -f "$IPR" ]; then
        SAW_CMD="${SAW_CMD} -imageRecordFile ${IPR}"
    fi
    if [ "$TAR" != "NONE" ] && [ -f "$TAR" ]; then
        SAW_CMD="${SAW_CMD} -imageCompressedFile ${TAR}"
    fi

    log "  Output: ${OUTDIR}"
    log "  Log: ${LOGFILE}"

    local START_TIME=$(date +%s)

    if eval "${SAW_CMD}" > "${LOGFILE}" 2>&1; then
        local END_TIME=$(date +%s)
        local DURATION=$(( (END_TIME - START_TIME) / 60 ))
        log "  SUCCESS: Completed in ${DURATION} minutes"
        echo "$(date) | SUCCESS | ${SAMPLE} | ${CHIP} | ${DURATION} min" >> "${LOG_DIR}/summary.log"
        return 0
    else
        local END_TIME=$(date +%s)
        local DURATION=$(( (END_TIME - START_TIME) / 60 ))
        log "  FAILED: Check ${LOGFILE} (ran for ${DURATION} min)"
        echo "$(date) | FAILED | ${SAMPLE} | ${CHIP} | ${DURATION} min | See ${LOGFILE}" >> "${LOG_DIR}/failures.log"
        return 1
    fi
}

# ============================================================
# MAIN EXECUTION
# ============================================================

TOTAL=${#SAMPLES[@]}

if [ "$TOTAL" -eq 0 ]; then
    echo "ERROR: No samples defined. Edit the SAMPLES array in this script."
    echo "Run saw_data_audit.py first to identify chips and their input file paths."
    exit 1
fi

COMPLETED=0
FAILED=0
SKIPPED=0

log "============================================================"
log "SAW Batch Pipeline - Starting"
log "Total samples: ${TOTAL}"
log "SAW container: ${SAW_SIF}"
log "Reference: ${REF_INDEX}"
log "Output: ${OUTPUT_BASE}"
log "Threads: ${THREADS}"
log "============================================================"

# Initialize summary log
echo "$(date) | SAW Batch Pipeline Started | ${TOTAL} samples" > "${LOG_DIR}/summary.log"

for entry in "${SAMPLES[@]}"; do
    IFS='|' read -r SAMPLE CHIP FQ1 FQ2 MASK TISSUE IPR TAR <<< "$entry"

    echo ""
    log "────────────────────────────────────────────────────────────"
    log "Sample ${SAMPLE} (${CHIP}) - $((COMPLETED + FAILED + SKIPPED + 1))/${TOTAL}"
    log "────────────────────────────────────────────────────────────"

    run_saw_sample "$SAMPLE" "$CHIP" "$FQ1" "$FQ2" "$MASK" "$TISSUE" "$IPR" "$TAR"
    RC=$?

    case $RC in
        0) COMPLETED=$((COMPLETED + 1)) ;;
        1) FAILED=$((FAILED + 1)) ;;
        2) SKIPPED=$((SKIPPED + 1)) ;;
    esac
done

# ============================================================
# FINAL SUMMARY
# ============================================================
echo ""
log "============================================================"
log "SAW Batch Pipeline - COMPLETE"
log "============================================================"
log "  Total:     ${TOTAL}"
log "  Completed: ${COMPLETED}"
log "  Failed:    ${FAILED}"
log "  Skipped:   ${SKIPPED}"
log ""
log "  Outputs: ${OUTPUT_BASE}/"
log "  Logs:    ${LOG_DIR}/"

if [ -f "${LOG_DIR}/failures.log" ]; then
    echo ""
    log "FAILURES:"
    cat "${LOG_DIR}/failures.log"
fi

echo "$(date) | SAW Batch Pipeline Complete | Done=${COMPLETED} Failed=${FAILED} Skipped=${SKIPPED}" >> "${LOG_DIR}/summary.log"
