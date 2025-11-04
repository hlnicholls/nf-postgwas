#!/usr/bin/env python3
# step3_LiftOver_38to37.py

import argparse
import os
import logging
from typing import Optional, Dict, Tuple

import pandas as pd
from pyliftover import LiftOver


class Step3LiftOver38to37:
    """
    Performs liftover from hg38 (GENPOS) to hg19 (GENPOS_hg19) with identical behavior to the original script:

    - Reads input TSV (all columns as string), uppercases column names.
    - Requires columns: CHROM, GENPOS, ALLELE0, ALLELE1.
    - Normalizes chromosome labels to 'chrN' (accepts '1' or 'chr1').
    - Coerces GENPOS to numeric (nullable Int64).
    - Deduplicates unique (CHROM, GENPOS) pairs for liftover to avoid repeated calls.
    - Uses LiftOver("hg38", "hg19") to map positions and applies back to full DataFrame.
    - Writes TSV with GENPOS_hg19 as nullable Int64, same as original script.
    """

    REQUIRED = ["CHROM", "GENPOS", "ALLELE0", "ALLELE1"]

    def __init__(self, inp: str, out: str, logger: Optional[logging.Logger] = None) -> None:
        self.inp = inp
        self.out = out
        self.logger = logger or logging.getLogger(__name__)
        self.df: Optional[pd.DataFrame] = None
        self.lo = LiftOver("hg38", "hg19")

    # -------- Orchestration --------

    def run(self) -> None:
        """Run the full liftover pipeline with identical semantics to the original script."""
        self._prepare_output_dir()
        self._load_input_uppercase()
        self._verify_required_columns()
        chrom_series, pos_series = self._normalize_chrom_and_pos()
        uniq_keys = self._build_unique_keys(chrom_series, pos_series)
        mapping = self._liftover_unique_keys(uniq_keys)
        self._apply_mapping_to_dataframe(mapping)
        self._write_output()

    # -------- Steps (mirror original behavior) --------

    def _prepare_output_dir(self) -> None:
        os.makedirs(os.path.dirname(self.out), exist_ok=True)
        self.logger.debug("Ensured output directory: %s", os.path.dirname(self.out))

    def _load_input_uppercase(self) -> None:
        # EXACT: read as strings, tab-separated
        self.df = pd.read_csv(self.inp, sep="\t", dtype=str)
        self.df.columns = self.df.columns.str.upper()
        self.logger.debug("Loaded %s with shape %s", self.inp, self.df.shape)

    def _verify_required_columns(self) -> None:
        assert self.df is not None
        missing = [c for c in self.REQUIRED if c not in self.df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")
        self.logger.debug("Required columns present")

    def _normalize_chrom_and_pos(self) -> Tuple[pd.Series, pd.Series]:
        """
        - CHROM: accept both '1' and 'chr1'; convert to 'chrN'
        - GENPOS: to numeric (nullable), as in original
        """
        assert self.df is not None
        chrom_series = self.df["CHROM"].astype(str)
        chrom_series = chrom_series.map(lambda c: c if c.startswith("chr") else f"chr{c}")
        pos_series = pd.to_numeric(self.df["GENPOS"], errors="coerce").astype("Int64")
        self.logger.debug("Normalized CHROM and GENPOS")
        return chrom_series, pos_series

    def _build_unique_keys(self, chrom_series: pd.Series, pos_series: pd.Series) -> pd.DataFrame:
        """
        Build unique (CHROM, GENPOS) combinations after dropping rows with NA positions.
        This replicates:
            uniq_keys = DataFrame({"CHROM": chrom_series, "GENPOS": pos_series}).dropna().drop_duplicates()
        """
        uniq_keys = pd.DataFrame({"CHROM": chrom_series, "GENPOS": pos_series}).dropna().drop_duplicates()
        self.logger.debug("Unique key count for liftover: %d", len(uniq_keys))
        return uniq_keys

    def _liftover_one(self, chrom: str, pos: int) -> Optional[int]:
        """
        Apply pyliftover to a single (chrom, pos). Returns new position or None if unmapped.
        Behavior identical to original: take first hit if multiple mappings.
        """
        result = self.lo.convert_coordinate(chrom, int(pos))
        return result[0][1] if result else None

    def _liftover_unique_keys(self, uniq_keys: pd.DataFrame) -> Dict[Tuple[str, int], Optional[int]]:
        """
        For each unique (CHROM, GENPOS) pair, compute hg19 position and store in a dict.
        """
        mapping: Dict[Tuple[str, int], Optional[int]] = {}
        for row in uniq_keys.itertuples(index=False):
            mapping[(row.CHROM, int(row.GENPOS))] = self._liftover_one(row.CHROM, int(row.GENPOS))
        self.logger.debug("Completed liftover for %d unique keys", len(mapping))
        return mapping

    def _apply_mapping_to_dataframe(self, mapping: Dict[Tuple[str, int], Optional[int]]) -> None:
        """
        Apply the liftover mapping back to the full DataFrame, exactly as the original list comprehension:
            df["GENPOS_hg19"] = [
                mapping.get((normalized_chrom, int(p)) if notna(p) else None, None)
                for c, p in zip(df["CHROM"].astype(str), to_numeric(df["GENPOS"]))
            ]
        Then coerce to nullable Int64.
        """
        assert self.df is not None
        # Rebuild normalized chrom and numeric pos exactly as original apply step
        norm_chrom = self.df["CHROM"].astype(str).map(lambda c: c if c.startswith("chr") else f"chr{c}")
        numeric_pos = pd.to_numeric(self.df["GENPOS"], errors="coerce")

        self.df["GENPOS_hg19"] = [
            mapping.get((c, int(p)) if pd.notna(p) else None, None)
            for c, p in zip(norm_chrom, numeric_pos)
        ]
        self.df["GENPOS_hg19"] = pd.to_numeric(self.df["GENPOS_hg19"], errors="coerce").astype("Int64")
        self.logger.debug("Applied mapping to full DataFrame")

        # (Optional duplicate handling from original remains commented)
        # dups_hg19 = self.df.duplicated(subset=["CHROM", "GENPOS_hg19", "ALLELE0", "ALLELE1"], keep=False)

    def _write_output(self) -> None:
        assert self.df is not None
        self.df.to_csv(self.out, sep="\t", index=False)
        self.logger.debug("Wrote output to %s with shape %s", self.out, self.df.shape)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--in", dest="inp", required=True, help="Input QC file: <trait>_assoc_regenie_allchr.txt")
    parser.add_argument("--out", dest="out", required=True, help="Output file: <trait>_38_37.txt")
    return parser.parse_args()


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )


def main() -> None:
    configure_logging()
    args = parse_args()
    lifter = Step3LiftOver38to37(inp=args.inp, out=args.out)
    lifter.run()


if __name__ == "__main__":
    main()
