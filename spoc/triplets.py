"""Managing multi-way contacts in particular generating binned triplets."""

import pandas as pd
import dask.dataframe as dd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sbn
from cooler.util import binnify
import pyranges as pr
import uuid


class TripletBinner:
    """Bins contacts over sequencing reads"""

    def __init__(self, use_dask: bool = False) -> None:
        if use_dask:
            self._reader_func = dd.read_parquet
        else:
            self._reader_func = pd.read_parquet

    def load_triplets(self, path: str):
        """Loads triplets as provided from the the sister-pore-c-snakemake pipeline"""
        triplets = self._reader_func(path)
        return triplets

    # def sort_triplets():
    #     """Sorts triplets such that the trans """
