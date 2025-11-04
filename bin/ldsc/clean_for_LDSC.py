#!/usr/bin/env python3
from __future__ import annotations
import os
import sys
import logging
from pathlib import Path
from typing import Iterable

import pandas as pd
import numpy as np

sys.path.insert(0, os.getcwd())
from config_python import output_path, traits  # type: ignore


class GWASPreprocessor:
    """
    Convert *_38_37_rsids.txt into LDSC-friendly GWAS files:
      <trait>_GWAS_37_corr.txt
    """

    def __init__(self, output_dir: Path, trait_names: Iterable[str]) -> None:
        self.output_dir = Path(output_dir)
        self.trait_names = list(trait_names)
        self.prep_dir = self.output_dir / "GWAS_Preprocessing"

        # I/O column expectations
        self.input_cols = [
            "CHROM", "GENPOS", "ALLELE0", "ALLELE1",
            "A1FREQ", "N", "BETA", "SE", "P", "GENPOS_hg19"
        ]

    def setup(self) -> None:
        self.prep_dir.mkdir(parents=True, exist_ok=True)
        logging.debug("Ensured output directory exists: %s", self.prep_dir)

    def _read_input(self, pheno: str) -> pd.DataFrame:
        in_path = self.prep_dir / f"{pheno}_38_37_rsids.txt"
        logging.info("Reading: %s", in_path)
        df = pd.read_csv(in_path, sep="\t", dtype=str)
        # Convert known numeric columns safely after read
        for col in ["GENPOS", "GENPOS_hg19", "A1FREQ", "N", "BETA", "SE", "P"]:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")
        return df

    @staticmethod
    def _drop_inf_nan(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
        for c in cols:
            if c in df.columns:
                df = df.dropna(subset=[c])
                df = df[~df[c].isin([np.inf, -np.inf])]
        return df

    def _clean(self, df: pd.DataFrame) -> pd.DataFrame:
        df = df[self.input_cols].copy()

        # Drop bad rows for both coordinates
        df = self._drop_inf_nan(df, ["GENPOS_hg19", "GENPOS"])

        # Cast coordinates to int (after cleaning)
        df["GENPOS_hg19"] = df["GENPOS_hg19"].astype(int)
        df["GENPOS"] = df["GENPOS"].astype(int)

        # Use lifted (hg19) positions
        df["GENPOS"] = df["GENPOS_hg19"]

        # Z-score
        df["Z"] = df["BETA"] / df["SE"]
        return df

    @staticmethod
    def _with_snp_id(df: pd.DataFrame) -> pd.DataFrame:
        out = df.rename(columns={"ALLELE0": "A1", "ALLELE1": "A2"}).copy()
        out["SNP"] = (
            out["CHROM"].astype(str)
            + ":" + out["GENPOS"].astype(str)
            + ":" + out["A1"].astype(str)
            + ":" + out["A2"].astype(str)
        )
        # Put SNP first
        cols = ["SNP"] + [c for c in out.columns if c != "SNP"]
        return out[cols]

    def _write(self, pheno: str, df: pd.DataFrame) -> Path:
        out_path = self.prep_dir / f"{pheno}_GWAS_37_corr.txt"
        logging.info("Writing: %s", out_path)
        df.to_csv(out_path, sep="\t", index=False)
        return out_path

    def process_trait(self, pheno: str) -> Path:
        df = self._read_input(pheno)
        df = self._clean(df)
        df = self._with_snp_id(df)
        return self._write(pheno, df)

    def run(self) -> None:
        self.setup()
        for pheno in self.trait_names:
            try:
                self.process_trait(pheno)
            except Exception as e:
                logging.exception("Failed processing trait %s: %s", pheno, e)
                raise  # preserve original fail-fast behaviour


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    proc = GWASPreprocessor(output_dir=Path(output_path), trait_names=traits)
    proc.run()


if __name__ == "__main__":
    main()
