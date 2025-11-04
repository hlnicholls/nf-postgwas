#!/bin/bash
set -euo pipefail

# Determine paths
SCRIPT_DIR=$(dirname "$0")
BASE_DIR=$(dirname "$SCRIPT_DIR")
CONFIG_FILE="$BASE_DIR/config_shell.sh"

# Load configuration variables (e.g., TRAITS, MTAG_TRAITS, MAGMA_PATH, etc.)
source "$CONFIG_FILE"

# Ensure Results directory exists
mkdir -p "$MAGMA_PATH/Results"

# Define portable MAGMA binary inside the repo
PORTABLE_MAGMA="$BASE_DIR/external/magma/magma"

# Check if MAGMA binary exists and is executable
if [[ ! -x "$PORTABLE_MAGMA" ]]; then
    echo "ERROR: MAGMA binary not found or not executable at: $PORTABLE_MAGMA" >&2
    exit 1
fi

# Run MAGMA for each primary phenotype
for phenotype in "${TRAITS[@]}"; do
    echo "Running MAGMA for trait: $phenotype"
    "$PORTABLE_MAGMA" \
        --bfile "$MAGMA_1000G" \
        --gene-annot "$MAGMA_annot" \
        --pval "$MAGMA_PATH/${phenotype}_38_rsid_magma.txt" ncol=N \
        --gene-model snp-wise=mean \
        --out "$MAGMA_PATH/Results/magma_${phenotype}"
done

echo "✓ MAGMA analysis completed successfully
