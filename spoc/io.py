"""Persisting functionality of spoc that manages writing to and reading from the filesystem."""

import pickle
from typing import Dict, Union
import pandas as pd
import dask.dataframe as dd
from .contacts import Contacts
from .dataframe_models import FragmentSchema, ContactSchema
from .fragments import Fragments


class FileManager:
    """Is responsible for loading and writing files"""

    def __init__(
        self, verify_schemas_on_load: bool = True, use_dask: bool = False
    ) -> None:
        self._verify_schemas = verify_schemas_on_load
        if use_dask:
            self._parquet_reader_func = dd.read_parquet
        else:
            self._parquet_reader_func = pd.read_parquet

    @staticmethod
    def write_label_library(path: str, data: Dict[str, bool]) -> None:
        """Writes label library to file"""
        with open(path, "wb") as handle:
            pickle.dump(data, handle)

    @staticmethod
    def load_label_library(path: str):
        """Load label library"""
        with open(path, "rb") as handle:
            label_library = pickle.load(handle)
        return label_library

    def load_fragments(self, path: str):
        """Load annotated fragments"""
        data = self._parquet_reader_func(path)
        return Fragments(data)

    @staticmethod
    def write_fragments(path: str, fragments: Fragments) -> None:
        """Write annotated fragments"""
        fragments.data.to_parquet(path)

    @staticmethod
    def write_multiway_contacts(path: str, contacts: Contacts) -> None:
        """Write multiway contacts"""
        contacts.data.to_parquet(path)


    # TODO; find a solution to hold information about symmetry flipping, label_sorting and binary_labeels equating
    def load_multiway_contacts(self, path: str, number_fragments: int) -> Contacts:
        """Load multiway contacts"""
        return Contacts(
            self._parquet_reader_func(path), number_fragments=number_fragments
        )

    @staticmethod
    def load_chromosome_sizes(path: str):
        """Load chromosome sizes"""
        # TODO: validate schema for this
        return pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["chrom", "size"],
            index_col=["chrom"],
            squeeze=True,
        )

    @staticmethod
    def load_pixels(path: str):
        """Loads pixels"""
        raise NotImplementedError

    @staticmethod
    def write_pixels(path: str, data: Union[pd.DataFrame, dd.DataFrame]) -> None:
        """Write pixels"""
        data.to_parquet(path)
