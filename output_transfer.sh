#!/bin/bash

# ============================================================
# Copy selected UFEMISM output files from Arrhenius HPC
#
# Usage:
#   bash copy_results.sh <folder_name>
#
# Example:
#   bash copy_results.sh results_R-LIS_init_BCAP_H2_O1_SMB1_GHF1_2km
# ============================================================

FOLDER="$1"

REMOTE_BASE="/nobackup/proj/disk/bolinc/personal/frare/outputs"
REMOTE="frare@login.hpc.arrhenius.naiss.se:${REMOTE_BASE}/${FOLDER}"

# ------------------------------------------------------------
# Check input
# ------------------------------------------------------------

if [ -z "$FOLDER" ]; then
    echo "Error: no folder specified."
    echo
    echo "Usage:"
    echo "  bash copy_results.sh <folder_name>"
    exit 1
fi

# ------------------------------------------------------------
# Files to copy
# ------------------------------------------------------------

INCLUDE_FILES=(
    "scalar_output_ANT_00001.nc"
    "main_output_ANT_00001.nc"
    # Add more files here
)

# ------------------------------------------------------------
# Create local directory
# ------------------------------------------------------------

mkdir -p "$FOLDER"

# ------------------------------------------------------------
# Build rsync include arguments
# ------------------------------------------------------------

RSYNC_ARGS=()

for FILE in "${INCLUDE_FILES[@]}"; do
    RSYNC_ARGS+=(--include="$FILE")
done

RSYNC_ARGS+=(--exclude='*')

# ------------------------------------------------------------
# Transfer
# ------------------------------------------------------------

echo "Copying selected files from:"
echo "  $REMOTE"
echo
echo "Files:"

for FILE in "${INCLUDE_FILES[@]}"; do
    echo "  $FILE"
done

echo

rsync -av \
    "${RSYNC_ARGS[@]}" \
    "$REMOTE/" \
    "$FOLDER/"

echo
echo "Transfer complete."
