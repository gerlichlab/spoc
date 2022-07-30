"""Persisting functionality of spoc that manages writing to and reading from the filesystem."""

import pickle
from typing import Dict, Union
import pandas as pd
import dask.dataframe as dd
from .models import fragment_schema, annotated_fragment_schema, HigherOrderContactSchema


class FileManager:
    """Is responsible for loading and writing files"""

    def __init__(self, verify_schemas_on_load: bool = True, use_dask: bool = False) -> None:
        self._verify_schemas = verify_schemas_on_load
        if use_dask:
            self._reader_func = dd.read_parquet
        else:
            self._reader_func = pd.read_parquet

    def write_label_library(self, path: str, data: Dict[str, bool]) -> None:
        with open(path, "wb") as handle:
            pickle.dump(data, handle)

    def load_label_library(self, path: str):
        with open(path, "rb") as handle:
            label_library = pickle.load(handle)
        return label_library

    def load_porec_fragments(self, path: str):
        data = self._reader_func(path)
        if self._verify_schemas:
            return fragment_schema.validate(data)
        return data

    def load_annotated_fragments(self, path: str):
        data = self._reader_func(path)
        if self._verify_schemas:
            return annotated_fragment_schema.validate(data)
        return data

    def write_annotated_fragments(self, path: str, data: Union[pd.DataFrame, dd.DataFrame]) -> None:
        data.to_parquet(path)

    def write_multiway_contacts(
        self, path: str, data: Union[pd.DataFrame, dd.DataFrame]
    ) -> None:
        data.to_parquet(path)

    def load_multiway_contacts(
        self, path: str, number_fragments: Union[int, None] = None
    ):
        data = self._reader_func(path)
        if self._verify_schemas:
            if number_fragments is None:
                raise ValueError(
                    "Number of fragments need to be specified is schema is validated!"
                )
            return HigherOrderContactSchema(number_fragments).validate(data)
        return data
