"""Persisting functionality of spoc that manages writing to and reading from the filesystem."""

import pickle
from typing import Dict, Union
import pandas as pd
import dask.dataframe as dd
from .dataframe_models import FragmentSchema, AnnotatedFragmentSchema, HigherOrderContactSchema


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

    def load_porec_fragments(self, path: str):
        """Load porec fragments"""
        data = self._parquet_reader_func(path)
        if self._verify_schemas:
            return FragmentSchema.validate(data)
        return data

    def load_annotated_fragments(self, path: str):
        """Load annotated fragments"""
        data = self._parquet_reader_func(path)
        if self._verify_schemas:
            return AnnotatedFragmentSchema.validate(data)
        return data

    @staticmethod
    def write_annotated_fragments(
        path: str, data: Union[pd.DataFrame, dd.DataFrame]
    ) -> None:
        """Write annotated fragments"""
        data.to_parquet(path)

    @staticmethod
    def write_multiway_contacts(
        path: str, data: Union[pd.DataFrame, dd.DataFrame]
    ) -> None:
        """Write multiway contacts"""
        data.to_parquet(path)

    def load_multiway_contacts(
        self, path: str, number_fragments: Union[int, None] = None
    ) -> Union[pd.DataFrame, dd.DataFrame]:
        """Load multiway contacts"""
        data = self._parquet_reader_func(path)
        if self._verify_schemas:
            if number_fragments is None:
                raise ValueError(
                    "Number of fragments needs to be specified if schema is validated!"
                )
            schema = HigherOrderContactSchema(number_fragments)
            schema.validate_header(
                data
            )  # this is needed to catch problems that cause task graph construction failures
            return schema.validate(data)
        return data

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
