#!/bin/bash

SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

# Navigate three levels up from the current script directory
PROJECT_ROOT=$(dirname "$(dirname "$SCRIPT_DIR")")

# Define the path to the config file in the project root
CONFIG_FILE="$PROJECT_ROOT/config_shell.sh"
source "$CONFIG_FILE"

mkdir -p $LDSC_OUT

for trait in "${TRAITS[@]}"
do
    echo "Running LDSC for ${trait} with original allele order..."
    ORIGINAL_FILE="${OUTPUT_PATH}/GWAS_Preprocessing/${trait}_GWAS_37_corr.txt"
    OUTPUT_FILE="${LDSC_OUT}/${trait}_ldsc_res_ukb"

    python2.7 $LDSC_PATH \
    --h2 $ORIGINAL_FILE \
    --ref-ld ${LDSC_UKBB_DATA} \
    --w-ld ${LDSC_UKBB_DATA} \
    --out $OUTPUT_FILE

done
