"""Persisting functionality of spoc that manages writing to and reading from the filesystem."""

import pickle
from typing import Dict, Optional
from hashlib import md5
import os
import json
from pathlib import Path
import pandas as pd
import dask.dataframe as dd
import pyarrow as pa
import pyarrow.parquet as pq
from pydantic import BaseModel
from spoc.contacts import Contacts
from spoc.pixels import Pixels
from spoc.file_parameter_models import ContactsParameters, PixelParameters
from spoc.fragments import Fragments


class FileManager:
    """Is responsible for loading and writing files"""

    def __init__(self, use_dask: bool = False) -> None:
        if use_dask:
            self._parquet_reader_func = dd.read_parquet
        else:
            self._parquet_reader_func = pd.read_parquet

    def _write_parquet_dask(
        self, path: str, df: dd.DataFrame, global_parameters: BaseModel
    ) -> None:
        """Write parquet file using dask"""
        custom_meta_data = {
            "spoc".encode(): json.dumps(global_parameters.dict()).encode()
        }
        dd.to_parquet(df, path, custom_metadata=custom_meta_data)

    def _write_parquet_pandas(
        self, path: str, df: pd.DataFrame, global_parameters: BaseModel
    ) -> None:
        """Write parquet file using pandas. Pyarrow is needed because
        the pandas .to_parquet method does not support writing custom metadata."""
        table = pa.Table.from_pandas(df)
        # Add metadata
        custom_meta_key = "spoc"
        existing_meta = table.schema.metadata
        combined_meta = {
            custom_meta_key.encode(): json.dumps(global_parameters.dict()).encode(),
            **existing_meta,
        }
        table = table.replace_schema_metadata(combined_meta)
        # write table
        pq.write_table(table, path)

    @staticmethod
    def _load_parquet_global_parameters(path: str) -> BaseModel:
        """Load global parameters from parquet file"""
        # check if path is a directory, if so, we need to read the schema from one of the partitioned files
        if os.path.isdir(path):
            path = path + "/" + os.listdir(path)[0]
        global_parameters = pa.parquet.read_schema(path).metadata.get("spoc".encode())
        if global_parameters is not None:
            global_parameters = json.loads(global_parameters.decode())
        else:
            # use default parameters
            global_parameters = ContactsParameters().dict()
        return global_parameters

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
        fragments.data.to_parquet(path, row_group_size=1024 * 1024)

    def write_multiway_contacts(self, path: str, contacts: Contacts) -> None:
        """Write multiway contacts"""
        if contacts.is_dask:
            self._write_parquet_dask(
                path, contacts.data, contacts.get_global_parameters()
            )
        else:
            self._write_parquet_pandas(
                path, contacts.data, contacts.get_global_parameters()
            )

    def load_contacts(
        self, path: str, global_parameters: Optional[ContactsParameters] = None
    ) -> Contacts:
        """Load multiway contacts"""
        if global_parameters is None:
            global_parameters = self._load_parquet_global_parameters(path)
        else:
            global_parameters = global_parameters.dict()
        return Contacts(self._parquet_reader_func(path), **global_parameters)

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
    def _load_pixel_metadata(path: str):
        """Load metadata"""
        metadata_path = Path(path) / "metadata.json"
        if metadata_path.exists():
            with open(metadata_path, "r", encoding="UTF-8") as f:
                metadata = json.load(f)
        else:
            raise ValueError(f"Metadata file not found at {metadata_path}")
        return metadata

    @staticmethod
    def list_pixels(path: str):
        """List available pixels"""
        # read metadata.json
        metadata = FileManager._load_pixel_metadata(path)
        # instantiate pixel parameters
        pixels = [PixelParameters(**params) for params in metadata.values()]
        return pixels

    def load_pixels(
        self, path: str, global_parameters: PixelParameters, load_dataframe: bool = True
    ) -> Pixels:
        """Loads specific pixels instance based on global parameters.
        load_dataframe specifies whether the dataframe should be loaded, or whether pixels
         should be instantiated based on the path alone."""
        metadata = self._load_pixel_metadata(path)
        # find matching pixels
        for pixel_path, value in metadata.items():
            param_value = PixelParameters(**value)
            if param_value == global_parameters:
                selected_value = param_value
                break
        else:
            raise ValueError(f"No matching pixels found for {global_parameters}")
        # rewrite path to contain parent folder
        pixel_path = Path(path) / pixel_path
        if load_dataframe:
            df = self._parquet_reader_func(pixel_path)
        else:
            df = pixel_path
        return Pixels(df, **selected_value.dict())

    @staticmethod
    def _get_pixel_hash_path(path: str, pixels: Pixels) -> str:
        hash_string = path + json.dumps(pixels.get_global_parameters().dict())
        return md5(hash_string.encode()).hexdigest() + ".parquet"

    def write_pixels(self, path: str, pixels: Pixels) -> None:
        """Write pixels"""
        # check whether path exists
        metadata_path = Path(path) / "metadata.json"
        if not Path(path).exists():
            # we need to create the directory first
            os.mkdir(path)
            current_metadata = {}
        else:
            current_metadata = FileManager._load_pixel_metadata(path)
        # create new file path -> hash of directory path and parameters
        write_path = Path(path) / self._get_pixel_hash_path(path, pixels)
        # write pixels
        if pixels.data is None:
            raise ValueError(
                "Writing pixels only suppported for pixels hodling dataframes!"
            )
        pixels.data.to_parquet(write_path, row_group_size=1024 * 1024)
        # write metadata
        current_metadata[write_path.name] = pixels.get_global_parameters().dict()
        with open(metadata_path, "w", encoding="UTF-8") as f:
            json.dump(current_metadata, f)
