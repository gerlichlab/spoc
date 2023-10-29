"""Persisting functionality of spoc that manages writing to and reading from the filesystem."""

import pickle
from typing import Dict, Optional, Union, List
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
from spoc.dataframe_models import FragmentSchema, ContactSchema, DataFrame
from spoc.file_parameter_models import ContactsParameters, PixelParameters
from spoc.fragments import Fragments


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
    def _load_metadata(path: str):
        """Load metadata"""
        metadata_path = Path(path) / "metadata.json"
        if metadata_path.exists():
            with open(metadata_path, "r") as f:
                metadata = json.load(f)
        else:
            raise ValueError(f"Metadata file not found at {metadata_path}")
        return metadata
    
    @staticmethod
    def list_pixels(path: str) -> List[PixelParameters]:
        """List available pixels"""
        # read metadata.json
        metadata = FileManager._load_metadata(path)
        # instantiate pixel parameters
        pixels = [
            PixelParameters(**params) for params in metadata.values()
        ]
        return pixels

    @staticmethod
    def list_contacts(path: str) -> List[ContactsParameters]:
        """List available contacts"""
        # read metadata.json
        metadata = FileManager._load_metadata(path)
        # instantiate pixel parameters
        contacts = [
            ContactsParameters(**params) for params in metadata.values()
        ]
        return contacts

    def load_pixels(self, path: str, global_parameters: PixelParameters, load_dataframe:bool = True) -> Pixels:
        """Loads specific pixels instance based on global parameters.
        load_dataframe specifies whether the dataframe should be loaded, or whether pixels
         should be instantiated based on the path alone. """
        metadata = self._load_metadata(path)
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

    def load_contacts(self, path: str, global_parameters: ContactsParameters) -> Contacts:
        """Loads specific contacts instance based on global parameters.
        load_dataframe specifies whether the dataframe should be loaded"""
        metadata = self._load_metadata(path)
        # find matching contacts
        for contacts_path, value in metadata.items():
            param_value = ContactsParameters(**value)
            if param_value == global_parameters:
                selected_value = param_value
                break
        else:
            raise ValueError(f"No matching contacts found for {global_parameters}")
        # rewrite path to contain parent folder
        contacts_path = Path(path) / contacts_path
        df = self._parquet_reader_func(contacts_path)
        return Contacts(df, **selected_value.dict())

    @staticmethod
    def _get_object_hash_path(path: str, data_object: Union[Pixels, Contacts]) -> str:
        hash_string = path + json.dumps(data_object.get_global_parameters().dict())
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
            current_metadata = FileManager._load_metadata(path)
        # create new file path -> hash of directory path and parameters
        write_path = Path(path) / self._get_object_hash_path(path, pixels)
        # write pixels
        if pixels.data is None:
            raise ValueError("Writing pixels only suppported for pixels hodling dataframes!")
        pixels.data.to_parquet(write_path, row_group_size=1024*1024)
        # write metadata
        current_metadata[write_path.name] = pixels.get_global_parameters().dict()
        with open(metadata_path, "w") as f:
            json.dump(current_metadata, f)


    def write_contacts(self, path: str, contacts: Contacts) -> None:
        """Write contacts"""
        # check whether path exists
        metadata_path = Path(path) / "metadata.json"
        if not Path(path).exists():
            # we need to create the directory first
            os.mkdir(path)
            current_metadata = {}
        else:
            current_metadata = self._load_metadata(path)
        # create new file path -> hash of directory path and parameters
        write_path = Path(path) / self._get_object_hash_path(path, contacts)
        # write contacts
        if contacts.data is None:
            raise ValueError("Writing contacts only suppported for contacts hodling dataframes!")
        contacts.data.to_parquet(write_path, row_group_size=1024*1024)
        # write metadata
        current_metadata[write_path.name] = contacts.get_global_parameters().dict()
        with open(metadata_path, "w") as f:
            json.dump(current_metadata, f)