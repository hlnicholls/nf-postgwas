#!/bin/bash

SCRIPT_DIR=$(dirname "$0")
BASE_DIR=$(dirname "$SCRIPT_DIR")
CONFIG_FILE="$BASE_DIR/config_shell.sh"
source "$CONFIG_FILE"

for phenotype in "${TRAITS[@]}"
    do
    echo "Processing phenotype: $phenotype"
    python $POPS_SCRIPT \
    --gene_annot_path $POPS_ANNOT \
    --feature_mat_prefix $POPS_FEATURES_ALL \
    --num_feature_chunks 2 \
    --magma_prefix $MAGMA_PATH/Results/magma_${phenotype} \
    --control_features_path $POPS_CONTROL \
    --out_prefix $POPS_PATH/${phenotype}_pops_all_features
 done
