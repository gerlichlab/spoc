"""Persisting functionality of spoc that manages writing to and reading from the filesystem."""

import pickle
from typing import Dict, Optional
import os
import json
import pandas as pd
import dask.dataframe as dd
import pyarrow as pa
import pyarrow.parquet as pq
from pydantic import BaseModel
from .contacts import Contacts
from .pixels import Pixels
from .dataframe_models import FragmentSchema, ContactSchema, DataFrame
from .file_parameter_models import ContactsParameters
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

    def _write_parquet_dask(self, path: str, df: dd.DataFrame, global_parameters: BaseModel) -> None:
        """Write parquet file using dask"""
        custom_meta_data = {
            'spoc'.encode(): json.dumps(global_parameters.dict()).encode()
        }
        dd.to_parquet(
            df,
            path,
            custom_metadata=custom_meta_data
        )
    
    def _write_parquet_pandas(self, path: str, df: pd.DataFrame ,global_parameters: BaseModel) -> None:
        """Write parquet file using pandas. Pyarrow is needed because
        the pandas .to_parquet method does not support writing custom metadata."""
        table = pa.Table.from_pandas(df)
        # Add metadata
        custom_meta_key = "spoc"
        existing_meta = table.schema.metadata
        combined_meta = {
            custom_meta_key.encode() : json.dumps(global_parameters.dict()).encode(),
            **existing_meta
        }
        table = table.replace_schema_metadata(combined_meta)
        # write table
        pq.write_table(table, path)


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
        if contacts.is_dask:
            self._write_parquet_dask(path, contacts.data, contacts.get_global_parameters())
        else:
            self._write_parquet_pandas(path, contacts.data, contacts.get_global_parameters())
        

    def load_contacts(self, path: str, global_parameters: Optional[ContactsParameters] = None) -> Contacts:
        """Load multiway contacts"""
        if global_parameters is None:
            # check if path is a directory, if so, we need to read the schema from one of the partitioned files
            if os.path.isdir(path):
                path = path + '/' + os.listdir(path)[0]
            global_parameters = pa.parquet.read_schema(path).metadata.get("spoc".encode())
            if global_parameters is not None:
                global_parameters = json.loads(global_parameters.decode())
            else:
                global_parameters = ContactsParameters()
        else:
            global_parameters = global_parameters.dict()
        return Contacts(
            self._parquet_reader_func(path), **global_parameters
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
