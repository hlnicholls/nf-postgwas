import pandas as pd
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
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
    "RegulomeDB_ranking", "RegulomeDB_probability_score", "IMPACT", "REVEL", 
    "Consequence", "SYMBOL", "SIFT", "PolyPhen", "CLIN_SIG", "ClinVar_CLNDN", 
    "Nonsynonymous_consequence_by_impact", "Nonsynonymous_Damaging_sift_or_polyphen",
    "Damaging_in_CADD_or_RegulomeDB"
]

filtered_df = filtered_df[column_order]
filtered_df.to_csv(os.path.join(annotated_credible_set_path, 'Annotated_CredibleSets_and_LD08.csv'), index=False)

total_rows = filtered_df.shape[0]
unique_lead_snps = filtered_df['lead_snp'].nunique()
unique_locus_names = filtered_df['Locus_name'].nunique()
non_na_r2_count = filtered_df['R2'].notna().sum()

lead_snp_groups = filtered_df.groupby('lead_snp')

damaging_snp_groups = lead_snp_groups['Nonsynonymous_Damaging_sift_or_polyphen'].apply(lambda x: (x == "Damaging").any()).sum()
total_damaging_snps = (filtered_df['Nonsynonymous_Damaging_sift_or_polyphen'] == "Damaging").sum()

damaging_snp_groups_cadd_regul = lead_snp_groups['Damaging_in_CADD_or_RegulomeDB'].apply(lambda x: (x == "Yes").any()).sum()
total_damaging_snps_cadd_regul = (filtered_df['Damaging_in_CADD_or_RegulomeDB'] == "Yes").sum()


exonic_conseq = {"synonymous_variant", "missense_variant", "inframe_insertion", "inframe_deletion", "stop_gained", 
                 "frameshift_variant", "coding_sequence_variant", "stop_lost", "stop_retained_variant", 
                 "incomplete_terminal_codon_variant", "start_retained_variant", "start_lost"}

intronic_conseq = {"splice_acceptor_variant", "splice_donor_variant", "splice_donor_5th_base_variant", 
                   "splice_region_variant", "splice_donor_region_variant", "splice_polypyrimidine_tract_variant", 
                   "intron_variant"}

other_conseq = {"transcript_ablation", "transcript_amplification", "feature_elongation", "feature_truncation", 
                "protein_altering_variant", "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", 
                "non_coding_transcript_exon_variant", "NMD_transcript_variant", "non_coding_transcript_variant", 
                "coding_transcript_variant", "upstream_gene_variant", "downstream_gene_variant", "TFBS_ablation", 
                "TFBS_amplification", "TF_binding_site_variant", "regulatory_region_ablation", 
                "regulatory_region_amplification", "regulatory_region_variant", "sequence_variant"}


damaging_exonic_snps = filtered_df[
    (filtered_df['Nonsynonymous_Damaging_sift_or_polyphen'] == "Damaging") &
    (filtered_df['Consequence'].apply(lambda x: any(c in exonic_conseq for c in str(x).split(',')) if isinstance(x, str) else False))
].shape[0]

damaging_exonic_snps_cadd_regul = filtered_df[
    (filtered_df['Damaging_in_CADD_or_RegulomeDB'] == "Yes") &
    (filtered_df['Consequence'].apply(lambda x: any(c in exonic_conseq for c in str(x).split(',')) if isinstance(x, str) else False))
].shape[0]


filtered_df['Consequence'] = filtered_df['Consequence'].astype(str).str.split(',')
all_consequences = [item for sublist in filtered_df['Consequence'].dropna() for item in sublist if item.lower() != 'nan']
consequence_counts = pd.Series(all_consequences).value_counts()

total_count = consequence_counts.sum()

exonic_count = sum(consequence_counts.get(c, 0) for c in exonic_conseq)
intronic_count = sum(consequence_counts.get(c, 0) for c in intronic_conseq)
other_count = sum(consequence_counts.get(c, 0) for c in other_conseq)
intergenic_count = consequence_counts.get("intergenic_variant", 0)

exonic_percentage = (exonic_count / total_count) * 100
intronic_percentage = (intronic_count / total_count) * 100
other_percentage = (other_count / total_count) * 100
intergenic_percentage = (intergenic_count / total_count) * 100

total_percentage = exonic_percentage + intronic_percentage + other_percentage + intergenic_percentage

report_path = os.path.join(annotated_credible_set_path, 'CredibleSet_report.txt')
with open(report_path, 'w') as report_file:
    report_file.write(f"Total number of rows: {total_rows}\n")
    report_file.write(f"Number of unique lead SNPs: {unique_lead_snps}\n")
    report_file.write(f"Number of unique Locus names: {unique_locus_names}\n")
    report_file.write(f"Number of rows with non-NA R2 values: {non_na_r2_count}\n")
    report_file.write(f"Number of lead SNP groups with damaging non-synonymous variants: {damaging_snp_groups}\n")
    report_file.write(f"Total 'Damaging' in Nonsynonymous_Damaging_sift_or_polyphen across all rows: {total_damaging_snps}\n")
    report_file.write(f"Number of 'Damaging' SNPs that are also exonic consequences: {damaging_exonic_snps}\n")
    report_file.write(f"Number of lead SNP groups with damaging CADD > 20 or RegulomeDB < 1f variants: {damaging_snp_groups_cadd_regul}\n")
    report_file.write(f"Total of SNPs with damaging CADD > 20 or RegulomeDB < 1f across all rows: {total_damaging_snps_cadd_regul}\n")
    report_file.write(f"Number of CADD > 20 or RegulomeDB < 1f SNPs that are also exonic consequences: {damaging_exonic_snps_cadd_regul}\n")
    report_file.write(f"Consequences (VEP):\n")
    report_file.write(f"Exonic: {exonic_count} ({exonic_percentage:.2f}%)\n")
    report_file.write(f"Intronic: {intronic_count} ({intronic_percentage:.2f}%)\n")
    report_file.write(f"Intergenic: {intergenic_count} ({intergenic_percentage:.2f}%)\n")
    report_file.write(f"Other: {other_count} ({other_percentage:.2f}%)\n")
    report_file.write(f"Other group consists of non-coding RNA, untranslated, and upstream/downstream variants\n")
    report_file.write(f"Total percentages sum to: {total_percentage:.2f}%\n")
    if total_percentage < 99.9999999999:
        report_file.write("Warning: Consequence categories do not sum to 100%! Manual checking needed.\n")
        report_file.write("Unique Consequences found:\n")
        report_file.write("\n".join(consequence_counts.index) + "\n")
        missing_consequences = set(consequence_counts.index) - (exonic_conseq | intronic_conseq | other_conseq | {"intergenic_variant"})
        report_file.write("Missing categorized consequences:\n")
        report_file.write("\n".join(missing_consequences) + "\n")

print(f"Report generated at: {report_path}")