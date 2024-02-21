"""Dataframe models"""
import copy
from enum import auto
from enum import Enum
from typing import Dict
from typing import Iterable
from typing import List
from typing import Optional
from typing import Protocol
from typing import Union

import dask.dataframe as dd
import duckdb
import pandas as pd
import pandera as pa

# Define dataframe type

DataFrame = Union[pd.DataFrame, dd.DataFrame, duckdb.DuckDBPyRelation]

FragmentSchema = pa.DataFrameSchema(
    {
        "chrom": pa.Column(str),
        "start": pa.Column(int),
        "end": pa.Column(int),
        "strand": pa.Column(bool),
        "read_name": pa.Column(str),
        "read_start": pa.Column(int),
        "read_end": pa.Column(int),
        "read_length": pa.Column(int),
        "mapping_quality": pa.Column(int),
        "align_score": pa.Column(int),
        "align_base_qscore": pa.Column(int),
        "pass_filter": pa.Column(bool, required=False),
        "metadata": pa.Column(str, required=False),
    },
    coerce=True,
)

RegionSchema = pa.DataFrameSchema(
    {
        "region_id": pa.Column(),
        "region_chrom": pa.Column(str),
        "region_start": pa.Column(int),
        "region_end": pa.Column(int),
    },
    coerce=True,
    unique=["region_id"],
)

# Protocol for genomic data


class GenomicDataSchema(Protocol):
    """Protocol for genomic data schema
    to be used in the query engine"""

    def get_position_fields(self) -> Dict[int, List[str]]:
        """Returns the position fields as a dictionary
        of framgent index to the respective fields"""

    def get_contact_order(self) -> int:
        """Returns the order of the genomic data"""

    def get_schema(self) -> pa.DataFrameSchema:
        """Return the schema of the underlying data"""

    def get_binsize(self) -> Optional[int]:
        """Returns the binsize of the genomic data"""

    def get_region_number(self) -> Optional[int]:
        """Returns the number of regions in the genomic data
        if present."""

    def get_half_window_size(self) -> Optional[int]:
        """Returns the window size of the genomic data
        if present."""


class QueryStepDataSchema:
    """Implements GenomicDataSchema for query steps
    with generic columns"""

    # pylint: disable=too-many-arguments
    # arguments needed to define the schema
    def __init__(
        self,
        columns: List[str],
        position_fields: Dict[int, List[str]],
        contact_order: int,
        binsize: Optional[int] = None,
        region_number: Optional[int] = None,
        half_window_size: Optional[int] = None,
    ) -> None:
        self._columns = columns
        self._contact_order = contact_order
        self._position_fields = position_fields
        self._binsize = binsize
        self._region_number = region_number
        self._half_window_size = half_window_size
        self._schema = pa.DataFrameSchema(
            {column: pa.Column() for column in columns},
            coerce=True,
        )

    def get_position_fields(self) -> Dict[int, List[str]]:
        """
        Returns the position fields as a dictionary.

        Returns:
            A dictionary where the keys are integers representing positions
            and the values are lists of strings representing the fields.
        """
        return self._position_fields

    def get_contact_order(self) -> int:
        """
        Returns the contact order of the object.
        """
        return self._contact_order

    def get_schema(self) -> pa.DataFrameSchema:
        """Return the schema of the underlying data"""
        return self._schema

    def get_binsize(self) -> Optional[int]:
        """Returns the binsize of the genomic data"""
        return self._binsize

    def get_region_number(self) -> Optional[int]:
        """Returns the number of regions in the genomic data
        if present."""
        return self._region_number

    def get_half_window_size(self) -> Optional[int]:
        """Returns the half window size of the genomic data
        if present."""
        return self._half_window_size


# schemas for higher order contacts


class ContactSchema:
    """Dynamic schema for N-way contacts

    Args:
        number_fragments (int, optional): Number of fragments. Defaults to 3.
        contains_metadata (bool, optional): Whether the contact data contains metadata. Defaults to True.
    """

    # field groups

    common_fields = {
        "read_name": pa.Column(str),
        "read_length": pa.Column(int),
    }

    contact_fields = {
        "chrom": pa.Column(str),
        "start": pa.Column(int),
        "end": pa.Column(int),
        "mapping_quality": pa.Column(int),
        "align_score": pa.Column(int),
        "align_base_qscore": pa.Column(int),
        "metadata": pa.Column(str, required=False),
    }

    def __init__(
        self, number_fragments: int = 3, contains_metadata: bool = True
    ) -> None:
        self._number_fragments = number_fragments
        self._schema = pa.DataFrameSchema(
            dict(
                self.common_fields,
                **self._expand_contact_fields(
                    range(1, number_fragments + 1), contains_metadata
                ),
            ),
            coerce=True,
        )

    @classmethod
    def get_contact_fields(cls, contains_metadata: bool) -> Dict:
        """returns contact fields

        Args:
            contains_metadata (bool): Whether the contact data contains metadata.

        Returns:
            Dict: Dictionary containing the contact fields.
        """
        if contains_metadata:
            return copy.deepcopy(cls.contact_fields)
        return {
            key: value
            for key, value in copy.deepcopy(cls.contact_fields).items()
            if key not in ["metadata"]
        }

    def _expand_contact_fields(
        self, expansions: Iterable = (1, 2, 3), contains_metadata: bool = True
    ) -> dict:
        """adds suffixes to fields"""
        output = {}
        for i in expansions:
            for key, value in self.get_contact_fields(contains_metadata).items():
                output[key + f"_{i}"] = value
        return output

    def validate_header(self, data_frame: DataFrame) -> None:
        """Validates only header, needed to validate that dask taskgraph can be built before
        evaluation.

        Args:
            data_frame (DataFrame): The DataFrame to validate.
        """
        for column in data_frame.columns:
            if column not in self._schema.columns:
                raise pa.errors.SchemaError(
                    self._schema, data_frame, "Header is invalid!"
                )

    def get_schema(self) -> pa.DataFrameSchema:
        """
        Get the schema of the DataFrame.

        Returns:
            pa.DataFrameSchema: The schema of the DataFrame.
        """
        return self._schema

    def get_position_fields(self) -> Dict[int, List[str]]:
        """Returns the position fields as a dictionary
        of framgent index to the respective fields"""
        return {
            i: [f"chrom_{i}", f"start_{i}", f"end_{i}"]
            for i in range(1, self._number_fragments + 1)
        }

    def get_contact_order(self) -> int:
        """Returns the order of the genomic data"""
        return self._number_fragments

    def validate(self, data_frame: DataFrame) -> DataFrame:
        """Validate multiway contact dataframe

        Args:
            data_frame (DataFrame): The DataFrame to validate.
        """
        self.validate_header(data_frame)
        if isinstance(data_frame, duckdb.DuckDBPyRelation):
            # duckdb does not support schema validation
            return data_frame
        return self._schema.validate(data_frame)

    def get_binsize(self) -> Optional[int]:
        """Returns the binsize of the genomic data"""
        return None

    def get_region_number(self) -> Optional[int]:
        """Returns the number of regions in the genomic data
        if present."""
        return None

    def get_half_window_size(self) -> Optional[int]:
        """Returns the window size of the genomic data
        if present."""
        return None


class PixelSchema:
    """Dynamic schema for N-way pixels

    Args:
        number_fragments (int, optional): Number of fragments. Defaults to 3.
        same_chromosome (bool, optional): Whether the fragments are on the same chromosome. Defaults to True.
    """

    def __init__(
        self,
        number_fragments: int = 3,
        same_chromosome: bool = True,
        binsize: Optional[int] = None,
    ) -> None:
        self._number_fragments = number_fragments
        self._same_chromosome = same_chromosome
        self._binsize = binsize
        self._schema = pa.DataFrameSchema(
            dict(
                self._get_constant_fields(),
                **self._expand_contact_fields(range(1, number_fragments + 1)),
            ),
            coerce=True,
        )

    def _get_contact_fields(self):
        if self._same_chromosome:
            return {
                "start": pa.Column(int),
            }
        return {"chrom": pa.Column(str), "start": pa.Column(int)}

    def _get_constant_fields(self):
        if self._same_chromosome:
            return {
                "chrom": pa.Column(str),
                "count": pa.Column(int),
                "corrected_count": pa.Column(float, required=False),
            }
        return {
            "count": pa.Column(int),
            "corrected_count": pa.Column(float, required=False),
        }

    def _expand_contact_fields(self, expansions: Iterable = (1, 2, 3)) -> dict:
        """adds suffixes to fields"""
        output = {}
        for i in expansions:
            for key, value in self._get_contact_fields().items():
                output[key + f"_{i}"] = value
        return output

    def validate_header(self, data_frame: DataFrame) -> None:
        """Validates only header, needed to validate that dask taskgraph can be built before
        evaluation

        Args:
            data_frame (DataFrame): The DataFrame to validate.
        """
        for column in data_frame.columns:
            if column not in self._schema.columns:
                raise pa.errors.SchemaError(
                    self._schema, data_frame, "Header is invalid!"
                )

    def get_schema(self) -> pa.DataFrameSchema:
        """
        Get the schema of the DataFrame.

        Returns:
            pa.DataFrameSchema: The schema of the DataFrame.
        """
        return self._schema

    def get_position_fields(self) -> Dict[int, List[str]]:
        """Returns the position fields as a dictionary
        of framgent index to the respective fields"""
        if self._same_chromosome:
            return {
                i: ["chrom", f"start_{i}"] for i in range(1, self._number_fragments + 1)
            }
        else:
            return {
                i: [f"chrom_{i}", f"start_{i}"]
                for i in range(1, self._number_fragments + 1)
            }

    def get_binsize(self) -> Optional[int]:
        """Returns the binsize of the genomic data"""
        return self._binsize

    def get_region_number(self) -> Optional[int]:
        """Returns the number of regions in the genomic data
        if present."""
        return None

    def get_contact_order(self) -> int:
        """Returns the order of the genomic data"""
        return self._number_fragments

    def get_half_window_size(self) -> Optional[int]:
        """Returns the window size of the genomic data
        if present."""
        return None

    def validate(self, data_frame: DataFrame) -> DataFrame:
        """Validate multiway contact dataframe

        Args:
            data_frame (DataFrame): The DataFrame to validate.

        """
        self.validate_header(data_frame)
        if isinstance(data_frame, duckdb.DuckDBPyRelation):
            # duckdb does not support schema validation
            return data_frame
        return self._schema.validate(data_frame)


class DataMode(Enum):
    """Enum for data mode"""

    PANDAS = auto()
    DASK = auto()
    DUCKDB = auto()
