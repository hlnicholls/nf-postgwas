#!/usr/bin/env python3
# step2_initial_QC.py 
import argparse
import os
import logging
from typing import Optional

import pandas as pd
import numpy as np


class Step1InitialQC:
    """
    Single-file, trait-agnostic QC:
      - Reads a TSV with all columns as strings
      - Uppercases column names
      - If LOG10P exists, computes P = 10 ** (-LOG10P)
      - Casts A1FREQ to float and computes MAF = min(A1FREQ, 1 - A1FREQ); keeps MAF >= 0.01
      - Drops EXTRA if present
      - If INFO exists, coerces to numeric and keeps INFO >= 0.3
      - Constructs SNP = CHROM:GENPOS:ALLELE0:ALLELE1 (as strings)
      - Writes TSV to --out
    """

    def __init__(self, inp: str, out: str, logger: Optional[logging.Logger] = None) -> None:
        self.inp = inp
        self.out = out
        self.logger = logger or logging.getLogger(__name__)
        self.df: Optional[pd.DataFrame] = None

    # ---------- high-level orchestration ----------

    def run(self) -> None:
        """Execute the full QC pipeline with identical semantics to the original script."""
        self._prepare_output_dir()
        self._load_input()
        self._uppercase_columns()
        self._maybe_compute_p_from_log10p()
        self._compute_maf_and_filter()
        self._maybe_drop_extra()
        self._maybe_filter_info()
        self._build_snp_column()
        self._write_output()

    # ---------- steps (match original behavior) ----------

    def _prepare_output_dir(self) -> None:
        os.makedirs(os.path.dirname(self.out), exist_ok=True)
        self.logger.debug("Ensured output directory exists: %s", os.path.dirname(self.out))

    def _load_input(self) -> None:
        # EXACT: read as strings, tab-separated
        self.df = pd.read_csv(self.inp, sep="\t", dtype=str)
        self.logger.debug("Loaded input with shape %s from %s", self.df.shape, self.inp)

    def _uppercase_columns(self) -> None:
        assert self.df is not None
        self.df.columns = self.df.columns.str.upper()
        self.logger.debug("Uppercased columns: %s", list(self.df.columns))

    def _maybe_compute_p_from_log10p(self) -> None:
        assert self.df is not None
        if 'LOG10P' in self.df.columns:
            # EXACT: P = 10 ** (-LOG10P.astype(float))
            self.df['P'] = 10 ** (-self.df['LOG10P'].astype(float))
            self.logger.debug("Computed P from LOG10P")

    def _compute_maf_and_filter(self) -> None:
        assert self.df is not None
        # EXACT: A1FREQ -> float; MAF = where(A1FREQ > 0.5, 1-A1FREQ, A1FREQ); keep MAF >= 0.01
        self.df['A1FREQ'] = self.df['A1FREQ'].astype(float)
        self.df['MAF'] = np.where(self.df['A1FREQ'] > 0.5, 1 - self.df['A1FREQ'], self.df['A1FREQ'])
        before = len(self.df)
        self.df = self.df[self.df['MAF'] >= 0.01]
        self.logger.debug("Filtered on MAF>=0.01: %d -> %d rows", before, len(self.df))

    def _maybe_drop_extra(self) -> None:
        assert self.df is not None
        if 'EXTRA' in self.df.columns:
            self.df = self.df.drop(columns=['EXTRA'])
            self.logger.debug("Dropped EXTRA column")

    def _maybe_filter_info(self) -> None:
        assert self.df is not None
        if 'INFO' in self.df.columns:
            # EXACT: coerce to numeric (NaN if bad), keep INFO >= 0.3
            self.df['INFO'] = pd.to_numeric(self.df['INFO'], errors='coerce')
            before = len(self.df)
            self.df = self.df[self.df['INFO'] >= 0.3]
            self.logger.debug("Filtered on INFO>=0.3: %d -> %d rows", before, len(self.df))

    def _build_snp_column(self) -> None:
        assert self.df is not None
        # EXACT: concatenate as strings with ':'
        self.df['SNP'] = (
            self.df['CHROM'] + ":" +
            self.df['GENPOS'] + ":" +
            self.df['ALLELE0'] + ":" +
            self.df['ALLELE1']
        )
        self.logger.debug("Constructed SNP column")

    def _write_output(self) -> None:
        assert self.df is not None
        self.df.to_csv(self.out, sep="\t", index=False)
        self.logger.debug("Wrote output to %s with shape %s", self.out, self.df.shape)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument('--in', dest='inp', required=True, help="Input TSV file")
    p.add_argument('--out', dest='out', required=True, help="Output TSV file")
    return p.parse_args()


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s"
    )


def main() -> None:
    configure_logging()
    args = parse_args()
    qc = Step1InitialQC(inp=args.inp, out=args.out)
    qc.run()


if __name__ == "__main__":
    main()
