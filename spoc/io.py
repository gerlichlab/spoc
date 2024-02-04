"""Persisting functionality of spoc that manages writing to and reading from the filesystem."""

import pickle
from typing import Dict, Union, List, Optional, Tuple
from hashlib import md5
import os
import json
from pathlib import Path
import pandas as pd
import dask.dataframe as dd
from functools import partial
import duckdb
from spoc.contacts import Contacts
from spoc.models.dataframe_models import DataMode
from spoc.pixels import Pixels
from spoc.models.file_parameter_models import (
    ContactsParameters,
    PixelParameters,
    GlobalParameters,
)
from spoc.fragments import Fragments

# Instantiate one duckdb connection to be used for all duckdb relations
DUCKDB_CONNECTION = duckdb.connect(database=":memory:")


class FileManager:
    """Is responsible for loading and writing files

    Args:
        data_mode (DataMode, optional): Data mode. Defaults to DataMode.PANDAS.
    """

    def __init__(self, data_mode: DataMode = DataMode.PANDAS) -> None:
        if data_mode == DataMode.DUCKDB:
            self._parquet_reader_func = partial(duckdb.read_parquet, connection=DUCKDB_CONNECTION)
        elif data_mode == DataMode.DASK:
            self._parquet_reader_func = dd.read_parquet
        elif data_mode == DataMode.PANDAS:
            self._parquet_reader_func = pd.read_parquet
        else:
            raise ValueError(f"Data mode {data_mode} not supported!")

    @staticmethod
    def write_label_library(path: str, data: Dict[str, bool]) -> None:
        """Writes label library to file

        Args:
            path (str): Path to write the file to.
            data (Dict[str, bool]): Label library data.

        Returns:
            None
        """
        with open(path, "wb") as handle:
            pickle.dump(data, handle)

    @staticmethod
    def load_label_library(path: str) -> Dict:
        """Load label library

        Args:
            path (str): Path to the label library file.

        Returns:
            Dict: Label library data.
        """
        with open(path, "rb") as handle:
            label_library = pickle.load(handle)
        return label_library

    def load_fragments(self, path: str) -> Fragments:
        """Load annotated fragments

        Args:
            path (str): Path to the fragments file.

        Returns:
            Fragments: Fragments object containing the fragment data.

        """
        data = self._parquet_reader_func(path)
        return Fragments(data)

    @staticmethod
    def write_fragments(path: str, fragments: Fragments) -> None:
        """Write annotated fragments

        Args:
            path (str): Path to write the file to.
            fragments (Fragments): Fragments object containing the fragment data.

        Returns:
            None

        """
        # Write fragments
        fragments.data.to_parquet(path, row_group_size=1024 * 1024)

    @staticmethod
    def load_chromosome_sizes(path: str) -> pd.DataFrame:
        """Load chromosome sizes

        Args:
            path (str): Path to the chromosome sizes file.

        Returns:
            pd.DataFrame: DataFrame containing the chromosome sizes.
        """
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
            with open(metadata_path, "r", encoding="UTF-8") as f:
                metadata = json.load(f)
        else:
            raise ValueError(f"Metadata file not found at {metadata_path}")
        return metadata

    @staticmethod
    def list_pixels(path: str) -> List[PixelParameters]:
        """List available pixels

        Args:
            path (str): Path to the pixel data.

        Returns:
            List[PixelParameters]: List of PixelParameters objects.
        """
        # read metadata.json
        metadata = FileManager._load_metadata(path)
        # instantiate pixel parameters
        pixels = [PixelParameters(**params) for params in metadata.values()]
        return pixels

    @staticmethod
    def list_contacts(path: str) -> List[ContactsParameters]:
        """List available contacts

        Args:
            path (str): Path to the contacts data.

        Returns:
            List[ContactsParameters]: List of ContactsParameters objects.
        """
        # read metadata.json
        metadata = FileManager._load_metadata(path)
        # instantiate pixel parameters
        contacts = [ContactsParameters(**params) for params in metadata.values()]
        return contacts

    def _parse_uri(
        self, uri: str, uri_parameters: List[str], min_fields: int = 2
    ) -> Tuple[str, Dict[str, str]]:
        """Parse URI

        Args:
            uri (str): URI to parse.
            uri_parameters (List[str]): List of URI parameters.
            min_fields (int, optional): Minimum number of fields that the URI should contain. Defaults to 2.
        Returns:
            Tuple(str, Dict[str, str]): Tuple containing the path and a dictionary of parameters.
        """
        # parse uri
        uri = uri.split("::")
        # validate uri
        if len(uri) < min_fields:
            raise ValueError(
                f"Uri: {uri} is not valid. Must contain at least Path, number_fragments and binsize"
            )
        params = dict(zip(uri_parameters, uri[1:]))
        # rewrite metadata_combi parameter
        if "metadata_combi" in params.keys() and params["metadata_combi"] != "None":
            params["metadata_combi"] = str(tuple(params["metadata_combi"]))
        return uri[0], params

    def _fuzzy_match_parameters(
        self,
        target_parameters: Dict[str, str],
        candidate_parameters: Dict[str, GlobalParameters],
    ) -> Tuple[str, GlobalParameters]:
        """Fuzzy match parameters

        Args:
            target_parameters (Dict[str,str]): Target parameters.
            candidate_parameters (Dict[str,GlobalParameters]): Candidate parameters.
        Returns:
            Tuple[str,GlobalParameters]: Tuple containing the path and a dictionary of parameters.
        """
        # get fuzzy matched parameters
        matched_parameters = [
            (path, param)
            for path, param in candidate_parameters.items()
            if all(
                str(value) == str(param.dict()[key])
                for key, value in target_parameters.items()
            )
        ]
        # check whether there was a unique match
        if len(matched_parameters) == 0:
            raise ValueError(f"No matches found for parameters: {target_parameters}!")
        if len(matched_parameters) > 1:
            raise ValueError(
                f"Multiple matches found for parameters: {target_parameters}!"
            )
        return matched_parameters[0]

    def load_pixels(
        self, path: str, global_parameters: Optional[PixelParameters] = None
    ) -> Pixels:
        """Loads specific pixels instance based on global parameters.
        load_dataframe specifies whether the dataframe should be loaded, or whether pixels
         should be instantiated based on the path alone.

        Args:
            path (str): Path to the pixel data.
            global_parameters (PixelParameters): Global parameters.

        Returns:
            Pixels: Pixels object containing the pixel data.

        """
        # if global parameters is None, path is assumed to be a uri
        if global_parameters is None:
            # parse uri
            path, parsed_parameters = self._parse_uri(
                path, PixelParameters.get_uri_fields(), min_fields=3
            )
        else:
            parsed_parameters = global_parameters.dict()
        # get fuzzy matched parameters
        metadata = {
            path: PixelParameters(**values)
            for path, values in self._load_metadata(path).items()
        }
        pixel_path, matched_parameters = self._fuzzy_match_parameters(
            parsed_parameters, metadata
        )
        # rewrite path to contain parent folder
        pixel_path = Path(path) / pixel_path
        df = self._parquet_reader_func(pixel_path)
        return Pixels(df, **matched_parameters.dict())

    def load_contacts(
        self, path: str, global_parameters: Optional[ContactsParameters] = None
    ) -> Contacts:
        """Loads specific contacts instance based on global parameters.
        load_dataframe specifies whether the dataframe should be loaded

        Args:
            path (str): Path to the contacts data.
            global_parameters (ContactsParameters): Global parameters.
        Returns:
            Contacts: Contacts object containing the contacts data.
        """
        # if global parameters is None, path is assumed to be a uri
        if global_parameters is None:
            # parse uri
            path, parsed_parameters = self._parse_uri(
                path, ContactsParameters.get_uri_fields(), min_fields=2
            )
        else:
            parsed_parameters = global_parameters.dict()
        # get fuzzy matched parameters
        metadata = {
            path: ContactsParameters(**values)
            for path, values in self._load_metadata(path).items()
        }
        contacts_path, matched_parameters = self._fuzzy_match_parameters(
            parsed_parameters, metadata
        )
        # rewrite path to contain parent folder
        contacts_path = Path(path) / contacts_path
        df = self._parquet_reader_func(contacts_path)
        return Contacts(df, **matched_parameters.dict())

    @staticmethod
    def _get_object_hash_path(path: str, data_object: Union[Pixels, Contacts]) -> str:
        hash_string = path + json.dumps(
            data_object.get_global_parameters().dict(),
            skipkeys=False,
            ensure_ascii=True,
            check_circular=True,
            allow_nan=True,
            indent=None,
            sort_keys=False,
            separators=None,
            default=None,
        )
        return md5(hash_string.encode(encoding="utf-8")).hexdigest() + ".parquet"

    def write_pixels(self, path: str, pixels: Pixels) -> None:
        """Write pixels

        Args:
            path (str): Path to write the pixel data to.
            pixels (Pixels): Pixels object containing the pixel data.

        Returns:
            None

        """
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
            raise ValueError(
                "Writing pixels only suppported for pixels hodling dataframes!"
            )
        pixels.data.to_parquet(write_path, row_group_size=1024 * 1024)
        # write metadata
        current_metadata[write_path.name] = pixels.get_global_parameters().dict()
        with open(metadata_path, "w", encoding="UTF-8") as f:
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
            raise ValueError(
                "Writing contacts only suppported for contacts hodling dataframes!"
            )
        contacts.data.to_parquet(write_path, row_group_size=1024 * 1024)
        # write metadata
        current_metadata[write_path.name] = contacts.get_global_parameters().dict()
        with open(metadata_path, "w") as f:
            json.dump(current_metadata, f)
