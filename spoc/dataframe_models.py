"""Dataframe models"""

from typing import Iterable, Union, Dict
import copy
import pandera as pa
import pandas as pd
import dask.dataframe as dd

# Define dataframe type

DataFrame = Union[pd.DataFrame, dd.DataFrame]

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

    def validate(self, data_frame: DataFrame) -> DataFrame:
        """Validate multiway contact dataframe
        
        Args:
            data_frame (DataFrame): The DataFrame to validate.
        """
        self.validate_header(data_frame)
        return self._schema.validate(data_frame)


class PixelSchema:
    """Dynamic schema for N-way pixels
    
    Args:
        number_fragments (int, optional): Number of fragments. Defaults to 3.
        same_chromosome (bool, optional): Whether the fragments are on the same chromosome. Defaults to True.
    """

    def __init__(self, number_fragments: int = 3, same_chromosome: bool = True) -> None:
        self._number_fragments = number_fragments
        self._same_chromosome = same_chromosome
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

    def validate(self, data_frame: DataFrame) -> DataFrame:
        """Validate multiway contact dataframe
        
        Args:
            data_frame (DataFrame): The DataFrame to validate.
        
        """
        return self._schema.validate(data_frame)
