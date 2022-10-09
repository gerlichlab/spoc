"""Dataframe models"""

import pandera as pa
import pandas as pd
from typing import Iterable

fragment_schema = pa.DataFrameSchema(
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
    coerce=True
)

annotated_fragment_schema = pa.DataFrameSchema(
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
    coerce=True
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
            coerce=True
        )

    def _get_contact_fields(self):
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

    def _expand_contact_fields(self, expansions: Iterable = [1, 2, 3]) -> dict:
        """adds suffixes to fields"""
        output = {}
        for i in expansions:
            for key, value in self._get_contact_fields().items():
                output[key + f"_{i}"] = value
        return output

    def validate(self, df: pd.DataFrame) -> None:
        """Validate multiway contact dataframe"""
        return self._schema.validate(df)

# schemas for higher order pixels

triplet_pixel_schema = pa.DataFrameSchema(
    {
        "chrom": pa.Column(str),
        "start_1": pa.Column(int),
        "start_2": pa.Column(int),
        "start_3": pa.Column(int),
        "count": pa.Column(int),
        "corrected_count": pa.Column(float, required=False)
    },
    coerce=True
)