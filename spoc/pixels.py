"""This part of spoc is responsible for binned higher order contacts in the form of 'genomic pixels'"""
from typing import Union
import pandas as pd
import dask.dataframe as dd
from cooler.util import binnify
import pyranges as pr


class GenomicBinner:
    """Bins higher order contacts into genomic bins of fixed size.
    Is capable of sorting genomic bins along columns based on sister chromatid
    identity and adding"""

    def __init__(
        self, bin_size: int, chrom_sizes: pd.DataFrame, use_dask: bool = True
    ) -> None:
        self._use_dask = use_dask
        self.bins = self._create_bins(chrom_sizes, bin_size)

    def _create_bins(self, chrom_sizes: pd.DataFrame, bin_size: int) -> pr.PyRanges:
        bins_df = binnify(chrom_sizes, bin_size)
        bins_df.index.name = "bin_id"
        return pr.PyRanges(
            bins_df.reset_index().rename(
                columns={"start": "Start", "end": "End", "chrom": "Chromosome"}
            )
        )

    def bin_contacts(self, contacts: Union[pd.DataFrame, dd.DataFrame]) -> pd.DataFrame:
        """Bins genomic contacts"""
