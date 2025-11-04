#!/usr/bin/env python3
# collate_LD.py
"""
Collate PLINK LD outputs across files, merge with loci and per-trait GWAS rsID maps,
and produce per-trait and combined locus+LD tables.

- Reads config from config_python.py in the current working directory (Nextflow runroot).
- Collates all *.ld files from ld_path and writes:
    <output_path>/Loci_Preprocessing/collated_ld_all_pheno_all_loci.txt
- Loads loci_path (lead SNPs / phenotype table), drops plotting/QC extras, and merges with LD.
- Discovers traits from Phenotype column if present; otherwise falls back to config `traits`.
- For each trait, merges with <output_path>/GWAS_Preprocessing/<trait>_38_37_rsids.txt
    and writes <output_path>/Loci_Preprocessing/<trait>_38_loci_with_ld.txt
- Finally writes the combined all-traits file:
    <output_path>/Loci_Preprocessing/all_traits_loci_38_with_ld.txt
"""

from __future__ import annotations

import os
import sys
import math
from typing import List, Sequence

import pandas as pd


class LDCollator:
    """
    Orchestrates:
      1) Collation of PLINK *.ld files
      2) Merge with loci table
      3) Trait discovery
      4) Per-trait merge with GWAS rsID maps
      5) Final combined export
    """

    LD_COLS = ['CHR_A', 'BP_A', 'SNP_A', 'CHR_B', 'BP_B', 'SNP_B', 'R2']
    DROP_LOCI_COLS = [
        "CHROM", "GENPOS", "ALLELE1", "ALLELE0", "A1FREQ",
        "MAF", "INFO", "BETA", "SE", "P", "rsid_1kg"
    ]
    DROP_AFTER_MERGE_COLS = ["CHR_A", "BP_A", "CHR_B", "BP_B", "GENPOS_hg19"]

    def __init__(self, traits: Sequence[str], output_path: str, ld_path: str, loci_path: str):
        self.cfg_traits = list(traits or [])
        self.output_path = output_path
        self.ld_path = ld_path
        self.loci_path = loci_path

        # Output locations
        self.out_dir_pre = os.path.join(self.output_path, "Loci_Preprocessing")
        self.out_dir_gwas = os.path.join(self.output_path, "GWAS_Preprocessing")
        os.makedirs(self.out_dir_pre, exist_ok=True)

        # Dataframes populated during run
        self.dat_ld: pd.DataFrame | None = None
        self.ld_pheno: pd.DataFrame | None = None
        self.actual_traits: List[str] = []

    # -------------------- Stage 1: LD collation --------------------

    def _list_ld_files(self) -> List[str]:
        """Return list of *.ld files in ld_path; raise if none."""
        files = [f for f in os.listdir(self.ld_path) if f.endswith(".ld")]
        if not files:
            raise FileNotFoundError(f"No .ld files found in {self.ld_path}")
        return files

    def collate_ld(self) -> pd.DataFrame:
        """
        Read and vertically concatenate all PLINK *.ld files (space-separated),
        enforcing the expected column schema.
        """
        files = self._list_ld_files()
        first_path = os.path.join(self.ld_path, files[0])
        dat_ld = pd.read_csv(first_path, sep=r"\s+", header=0, names=self.LD_COLS)

        for fname in files[1:]:
            tmp = pd.read_csv(os.path.join(self.ld_path, fname), sep=r"\s+", header=0, names=self.LD_COLS)
            dat_ld = pd.concat([dat_ld, tmp], ignore_index=True)

        # Write collated LD for reference (same as original behavior)
        collated_path = os.path.join(self.out_dir_pre, "collated_ld_all_pheno_all_loci.txt")
        dat_ld.to_csv(collated_path, index=False, sep="\t")
        self.dat_ld = dat_ld
        return dat_ld

    # -------------------- Stage 2: Loci loading & merge --------------------

    def _load_loci_table(self) -> pd.DataFrame:
        """
        Load loci/lead SNPs table from loci_path (CSV with header),
        rename SNP->ID, drop plotting/QC extras (if present).
        """
        pheno_file = pd.read_csv(self.loci_path)
        pheno_file = pheno_file.rename(columns={"SNP": "ID"})
        # Drop optional columns if present (errors="ignore" equivalent via selective list comp)
        to_drop = [c for c in self.DROP_LOCI_COLS if c in pheno_file.columns]
        if to_drop:
            pheno_file = pheno_file.drop(columns=to_drop)
        return pheno_file

    def merge_ld_with_loci(self) -> pd.DataFrame:
        """
        Merge concatenated LD (SNP_A as ID) with loci table on lead SNP ID.
        Produces ld_pheno with columns: lead_snp, SNP (partner in LD), R2, Phenotype, etc.
        """
        if self.dat_ld is None:
            raise RuntimeError("LD data is not loaded. Call collate_ld() first.")

        pheno_file = self._load_loci_table()

        dat_ld = self.dat_ld.rename(columns={'SNP_A': 'ID', 'SNP_B': 'SNP'})
        ld_pheno = pd.merge(dat_ld, pheno_file, on='ID', how='left').rename(columns={'ID': 'lead_snp'})

        # Clean up unused columns if present
        drop_cols = [c for c in self.DROP_AFTER_MERGE_COLS if c in ld_pheno.columns]
        if drop_cols:
            ld_pheno = ld_pheno.drop(columns=drop_cols)

        self.ld_pheno = ld_pheno
        return ld_pheno

    # -------------------- Stage 3: Trait discovery --------------------

    @staticmethod
    def _is_valid_trait(t) -> bool:
        """Return True if trait value is a non-empty, non-NA string."""
        if t is None:
            return False
        s = str(t).strip()
        if not s:
            return False
        low = s.lower()
        if low in {"na", "nan"}:
            return False
        # guard for float NaN
        if isinstance(t, float) and math.isnan(t):
            return False
        return True

    def discover_traits(self) -> List[str]:
        """
        Discover unique traits from ld_pheno['Phenotype'] if present; otherwise use config traits.
        Apply the same filtering as the original script.
        """
        if self.ld_pheno is not None and "Phenotype" in self.ld_pheno.columns:
            raw = self.ld_pheno["Phenotype"].unique().tolist()
        else:
            raw = self.cfg_traits

        actual = [t for t in raw if self._is_valid_trait(t)]
        # Preserve ordering: keep first occurrence order
        seen = set()
        ordered: List[str] = []
        for t in actual:
            key = str(t)
            if key not in seen:
                seen.add(key)
                ordered.append(key)

        self.actual_traits = ordered
        print(f"Filtered traits: {self.actual_traits}")
        return self.actual_traits

    # -------------------- Stage 4: Per-trait merge with GWAS rsIDs --------------------

    def _gwas_rsids_path(self, trait: str) -> str:
        """Return path to the per-trait rsIDs file produced by GET_RSIDS."""
        return os.path.join(self.out_dir_gwas, f"{trait}_38_37_rsids.txt")

    def _per_trait_out_path(self, trait: str) -> str:
        """Return path to the per-trait merged LD+loci+GWAS output."""
        return os.path.join(self.out_dir_pre, f"{trait}_38_loci_with_ld.txt")

    def merge_per_trait(self, trait: str, ld_pheno: pd.DataFrame) -> pd.DataFrame:
        """
        For a given trait:
          - Read GWAS rsID map: <output>/GWAS_Preprocessing/<trait>_38_37_rsids.txt
          - Merge with rows from ld_pheno where Phenotype == trait on 'SNP'
          - Write <output>/Loci_Preprocessing/<trait>_38_loci_with_ld.txt
          - Return the merged DataFrame (can be empty)
        """
        gwas_file = self._gwas_rsids_path(trait)
        gwas = pd.read_csv(gwas_file, sep="\t")

        trait_rows = ld_pheno[ld_pheno["Phenotype"] == trait] if "Phenotype" in ld_pheno.columns else ld_pheno
        merged = pd.merge(trait_rows, gwas, on="SNP", how="left")

        out_path = self._per_trait_out_path(trait)
        if merged.empty:
            # Preserve header behavior like the original
            header = list(merged.columns) if len(merged.columns) > 0 else list(ld_pheno.columns) + list(gwas.columns)
            pd.DataFrame(columns=header).to_csv(out_path, index=False, sep="\t")
        else:
            merged.to_csv(out_path, index=False, sep="\t")
        return merged

    # -------------------- Stage 5: Final combined file --------------------

    def write_final_all_traits(self, frames: List[pd.DataFrame]) -> str:
        """
        Concatenate all per-trait frames, filter rows with non-null CHROM,
        sort by CHROM, and write the final combined output.
        """
        if not frames:
            # Create an empty file with no rows if absolutely nothing merged
            final_out = os.path.join(self.out_dir_pre, "all_traits_loci_38_with_ld.txt")
            pd.DataFrame().to_csv(final_out, index=False, sep="\t")
            return final_out

        all_traits_df = pd.concat(frames, ignore_index=True, sort=False)
        if "CHROM" in all_traits_df.columns:
            all_traits_df = all_traits_df[all_traits_df["CHROM"].notna()].sort_values(by="CHROM")

        final_out = os.path.join(self.out_dir_pre, "all_traits_loci_38_with_ld.txt")
        all_traits_df.to_csv(final_out, index=False, sep="\t")
        return final_out

    # -------------------- Orchestration --------------------

    def run(self) -> None:
        """Run the complete pipeline in the same order as the original script."""
        # 1) Collate LD across *.ld files
        self.collate_ld()

        # 2) Merge with loci table
        ld_pheno = self.merge_ld_with_loci()

        # 3) Discover traits (print for parity with original)
        traits = self.discover_traits()

        # 4) Per-trait merge with GWAS rsIDs and write outputs
        per_trait_frames: List[pd.DataFrame] = []
        for trait in traits:
            merged = self.merge_per_trait(trait, ld_pheno)
            per_trait_frames.append(merged)

        # 5) Write final combined file
        self.write_final_all_traits(per_trait_frames)


def main() -> None:
    """
    Entry-point: import config from runroot, instantiate LDCollator, and run.
    Keeps the original behavior of preferring config_python.py from CWD.
    """
    # Ensure current working directory is first on sys.path (Nextflow runroot)
    sys.path.insert(0, os.getcwd())

    # Import config from runroot (raises if missing)
    from config_python import traits, output_path, ld_path, loci_path  # type: ignore

    collator = LDCollator(
        traits=traits,
        output_path=output_path,
        ld_path=ld_path,
        loci_path=loci_path,
    )
    collator.run()


if __name__ == "__main__":
    main()
