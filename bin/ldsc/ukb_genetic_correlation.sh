#!/bin/bash

SCRIPT_DIR=$(dirname "$0")
BASE_DIR=$(dirname "$SCRIPT_DIR")/..
CONFIG_FILE=$(realpath "$BASE_DIR/config_shell.sh")
echo "SCRIPT_DIR: $SCRIPT_DIR"
echo "BASE_DIR: $BASE_DIR"
echo "CONFIG_FILE: $CONFIG_FILE"
source "$CONFIG_FILE"

mkdir -p $LDSC_OUT

for phenotype in "${TRAITS[@]}"
do
    other_phenotypes=()
    for item in "${TRAITS[@]}"; do
        if [[ "$item" != "$phenotype" ]]; then
            other_phenotypes+=("$item")
        fi
    done

    other_phenotypes_string=$(printf ",${OUTPUT_PATH}/GWAS_Preprocessing/%s_GWAS_37_corr.txt" "${other_phenotypes[@]}")
    other_phenotypes_string=${other_phenotypes_string:1}

    python2.7 $LDSC_PATH \
    --rg ${OUTPUT_PATH}/GWAS_Preprocessing/${phenotype}_GWAS_37_corr.txt,${other_phenotypes_string} \
    --ref-ld $LDSC_UKBB_DATA \
    --w-ld $LDSC_UKBB_DATA \
        if [ "${#TRAITS[@]}" -le 1 ]; then
            echo "Genetic correlation requires at least 2 traits. Skipping."
            exit 0
        fi
    --out $LDSC_OUT/genetic_correlation_ldsc_res_ukb_${phenotype}
done
