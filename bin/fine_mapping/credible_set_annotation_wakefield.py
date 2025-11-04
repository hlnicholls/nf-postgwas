import pandas as pd
import sys
import os

# Import config from current working directory (runroot)
sys.path.insert(0, os.getcwd())
from config_python import var_path, wakefield_cs, susie_cs, carma_cs, annotated_credible_set_path

var_annots = pd.read_csv(os.path.join(var_path, 'all_variant_annotations.csv'))
var_annots = var_annots.drop_duplicates(subset=['SNP', 'Phenotype'])

cs = pd.read_csv(wakefield_cs, sep='\t')
cs['credible_set_99'] = cs['credible_set']
cs_to_merge = cs[['SNP', 'Phenotype', 'credible_set_99', 'credible_set_95']]
cs_to_merge = cs_to_merge.drop_duplicates(subset=['SNP', 'Phenotype'])


merged_df = pd.merge(var_annots, cs_to_merge, on=['SNP', 'Phenotype'], how='left')
merged_df['R2'] = pd.to_numeric(merged_df['R2'], errors='coerce')


filtered_df = merged_df[
    (merged_df['credible_set_99'] == 'Yes') |
    ((merged_df['credible_set_99'] == 'No') & (merged_df['R2'] > 0.8))
]

filtered_df['RegulomeDB_numeric'] = filtered_df['RegulomeDB_ranking'].str.extract(r'(\d+\.?\d*)').astype(float)
filtered_df['Damaging_in_CADD_or_RegulomeDB'] = filtered_df.apply(
    lambda row: 'Yes' if row['RegulomeDB_numeric'] == 1 or row['CADD_RawScore'] > 20 else 'No',
    axis=1
)

column_order = [
    "SNP", "Phenotype", "Method", "Gene_Symbol", "Locus_name", "lead_snp", "R2", 
    'credible_set_99',	'credible_set_95',
    "CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", 
    "A1FREQ", "N", "BETA", "SE", "P", "MAF", "GENPOS_hg19", "rsid_1kg", 
    "CADD_PHRED", "CADD_RawScore", "HiC_gene", "HiC_tissue", 
    "RegulomeDB_ranking", "RegulomeDB_probability_score",
    "Damaging_in_CADD_or_RegulomeDB"
]

filtered_df = filtered_df[column_order]
filtered_df.to_csv(os.path.join(annotated_credible_set_path, 'Annotated_CredibleSets_and_LD08.csv'), index=False)

total_rows = filtered_df.shape[0]
unique_lead_snps = filtered_df['lead_snp'].nunique()
unique_locus_names = filtered_df['Locus_name'].nunique()
non_na_r2_count = filtered_df['R2'].notna().sum()

lead_snp_groups = filtered_df.groupby('lead_snp')

damaging_snp_groups_cadd_regul = lead_snp_groups['Damaging_in_CADD_or_RegulomeDB'].apply(lambda x: (x == "Yes").any()).sum()
total_damaging_snps_cadd_regul = (filtered_df['Damaging_in_CADD_or_RegulomeDB'] == "Yes").sum()


report_path = os.path.join(annotated_credible_set_path, 'CredibleSet_report.txt')
with open(report_path, 'w') as report_file:
    report_file.write(f"Total number of rows: {total_rows}\n")
    report_file.write(f"Number of unique lead SNPs: {unique_lead_snps}\n")
    report_file.write(f"Number of unique Locus names: {unique_locus_names}\n")
    report_file.write(f"Number of rows with non-NA R2 values: {non_na_r2_count}\n")
    report_file.write(f"Number of lead SNP groups with damaging CADD > 20 or RegulomeDB < 1f variants: {damaging_snp_groups_cadd_regul}\n")
    report_file.write(f"Total of SNPs with damaging CADD > 20 or RegulomeDB < 1f across all rows: {total_damaging_snps_cadd_regul}\n")
print(f"Report generated at: {report_path}")