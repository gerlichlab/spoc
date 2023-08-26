"""Persisting functionality of spoc that manages writing to and reading from the filesystem."""

import pickle
from typing import Dict, Union
import pandas as pd
import dask.dataframe as dd
import pyarrow as pa
import pyarrow.parquet as pq
import json
from pydantic import BaseModel
from .contacts import Contacts
from .pixels import Pixels
from .dataframe_models import FragmentSchema, ContactSchema, DataFrame
from .fragments import Fragments


class FileManager:
    """Is responsible for loading and writing files"""

    def __init__(
        self, use_dask: bool = False
    ) -> None:
        if use_dask:
            self._parquet_reader_func = dd.read_parquet
        else:
            self._parquet_reader_func = pd.read_parquet

    @staticmethod
    def _update_parquet_metadata(path:str, update_metadata: BaseModel):
        """Update parquet metadata"""
        md = pa.parquet.read_metadata(path)
        metadata_dict = md.to_dict()
        metadata_dict.update(update_metadata.dict())
        md.set_metadata(metadata_dict)
        pq.write_metadata(md, path)

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
        # Write fragments
        fragments.data.to_parquet(path, row_group_size=1024*1024)

    def write_multiway_contacts(self, path: str, contacts: Contacts) -> None:
        """Write multiway contacts"""
        # Write contacts
        contacts.data.to_parquet(path, row_group_size=1024*1024)
        # Write metadata
        self._update_parquet_metadata(path, contacts.get_global_parameters())
        


    # TODO; find a solution to hold information about symmetry flipping, label_sorting and binary_labeels equating
    # -> do it with metadata of parquet file
    def load_contacts(self, path: str, *args, **kwargs) -> Contacts:
        """Load multiway contacts"""
        return Contacts(
            self._parquet_reader_func(path), *args, **kwargs
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
        # TODO:make this conform with new file format
        raise NotImplementedError

    @staticmethod
    def write_pixels(path: str, pixels: Pixels) -> None:
        """Write pixels"""
        pixels.data.to_parquet(path, row_group_size=1024*1024)
