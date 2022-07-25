"""Persisting functionality of spoc that manages writing to and reading from the filesystem."""

import pickle
from typing import Dict, Union
import pandas as pd
import dask.dataframe as dd
from .models import fragment_schema, HigherOrderContactSchema


class FileManager:
    """Is responsible for loading and writing files"""

    def __init__(self, verify_schemas_on_load: bool = True) -> None:
        self._verify_schemas = verify_schemas_on_load

    def write_label_library(self, path: str, data: Dict[str, bool]) -> None:
        with open(path, "wb") as handle:
            pickle.dump(data, handle)

    def load_label_library(self, path: str):
        with open(path, "rb") as handle:
            label_library = pickle.load(handle)
        return label_library

    def load_porec_fragments(self, path: str):
        data = pd.read_parquet(path)
        if self._verify_schemas:
            fragment_schema.validate(data)
        return data

    def write_multiway_contacts(
        self, path: str, data: Union[pd.DataFrame, dd.DataFrame]
    ) -> Union[pd.DataFrame, dd.DataFrame]:
        data.to_parquet(path)

    def load_multiway_contacts(
        self, path: str, number_fragments: Union[int, None] = None
    ):
        data = pd.read_parquet(path)
        if self._verify_schemas:
            if number_fragments is None:
                raise ValueError(
                    "Number of fragments need to be specified is schema is validated!"
                )
            HigherOrderContactSchema(number_fragments).validate(data)
        return data
