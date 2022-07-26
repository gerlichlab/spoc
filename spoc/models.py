"""Dataframe models"""

from typing import Iterable, Union
import pandera as pa
import pandas as pd
import dask.dataframe as dd

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
        "pass_filter": pa.Column(bool),
    },
    coerce=True,
)

AnnotatedFragmentSchema = pa.DataFrameSchema(
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
        "is_labelled": pa.Column(bool),
        "pass_filter": pa.Column(bool),
        "sister_identity": pa.Column(
            str, checks=[pa.Check(lambda x: x.isin(["SisterA", "SisterB"]))]
        ),
    },
    coerce=True,
)

# schemas for higher order contacts


class HigherOrderContactSchema:
    """Dynamic schema for n-way contacts"""

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
        "is_labelled": pa.Column(bool),
        "sister_identity": pa.Column(
            str, checks=[pa.Check(lambda x: x.isin(["SisterA", "SisterB"]))]
        ),
    }

    def __init__(self, number_fragments: int = 3) -> None:
        self._number_fragments = number_fragments
        self._schema = pa.DataFrameSchema(
            dict(
                self.common_fields,
                **self._expand_contact_fields(range(1, number_fragments + 1)),
            ),
            coerce=True,
        )

    @staticmethod
    def _get_contact_fields():
        return {
            "chrom": pa.Column(str),
            "start": pa.Column(int),
            "end": pa.Column(int),
            "mapping_quality": pa.Column(int),
            "align_score": pa.Column(int),
            "align_base_qscore": pa.Column(int),
            "is_labelled": pa.Column(bool),
            "sister_identity": pa.Column(
                str, checks=[pa.Check(lambda x: x.isin(["SisterA", "SisterB"]))]
            ),
        }

    def _expand_contact_fields(self, expansions: Iterable = (1, 2, 3)) -> dict:
        """adds suffixes to fields"""
        output = {}
        for i in expansions:
            for key, value in self._get_contact_fields().items():
                output[key + f"_{i}"] = value
        return output

    def validate_header(self, data_frame: Union[pd.DataFrame, dd.DataFrame]) -> None:
        """Validates only header, needed to validate that dask taskgraph can be built before
        evaluation"""
        for column in data_frame.columns:
            if column not in self._schema.columns:
                raise pa.errors.SchemaError(
                    self._schema, data_frame, "Header is invalid!"
                )

    def validate(
        self, data_frame: Union[pd.DataFrame, dd.DataFrame]
    ) -> Union[pd.DataFrame, dd.DataFrame]:
        """Validate multiway contact dataframe"""
        return self._schema.validate(data_frame)


# schemas for higher order pixels TODO: add higher order pixels or generic order

TripletPixelSchema = pa.DataFrameSchema(
    {
        "chrom": pa.Column(str),
        "start_1": pa.Column(int),
        "start_2": pa.Column(int),
        "start_3": pa.Column(int),
        "count": pa.Column(int),
        "corrected_count": pa.Column(float, required=False),
    },
    coerce=True,
)
