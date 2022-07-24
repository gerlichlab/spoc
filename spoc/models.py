"""Dataframe models"""

import pandera as pa
import pandas as pd
from typing import Iterable

fragment_schema = pa.DataFrameSchema(
    {
        "chrom": pa.Column(str, checks=[pa.Check.str_startswith("chr")]),
        "start": pa.Column(int, checks=pa.Check.gt(0)),
        "end": pa.Column(int, checks=pa.Check.gt(0)),
        "strand": pa.Column(bool),
        "read_name": pa.Column(str),
        "mapping_quality": pa.Column(int),
        "align_score": pa.Column(int),
        "align_base_qscore": pa.Column(int),
    }
)

annotated_fragment_schema = pa.DataFrameSchema(
    {
        "chrom": pa.Column(str, checks=[pa.Check.str_startswith("chr")]),
        "start": pa.Column(int, checks=pa.Check.gt(0)),
        "end": pa.Column(int, checks=pa.Check.gt(0)),
        "strand": pa.Column(bool),
        "read_name": pa.Column(str),
        "mapping_quality": pa.Column(int),
        "align_score": pa.Column(int),
        "align_base_qscore": pa.Column(int),
        "is_labelled": pa.Column(bool),
        "sister_identity": pa.Column(
            str, checks=[pa.Check(lambda x: x.isin(["SisterA", "SisterB"]))]
        ),
    }
)

# schemas for higher order contacts


class HigherOrderContactSchema:
    """Dynamic schema for n-way contacts"""

    # field groups

    _common_fields = {
        "read_name": pa.Column(str),
        "read_length": pa.Column(int, checks=pa.Check.gt(0)),
    }

    _contact_fields = {
        "chrom": pa.Column(str, checks=[pa.Check.str_startswith("chr")]),
        "start": pa.Column(int, checks=pa.Check.gt(0)),
        "end": pa.Column(int, checks=pa.Check.gt(0)),
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
                self._common_fields,
                **self._expand_fields(
                    self._contact_fields, range(1, number_fragments + 1)
                ),
            )
        )

    def _expand_fields(self, fields: dict, expansions: Iterable = [1, 2, 3]) -> dict:
        """adds suffixes to fields"""
        output = {}
        for i in expansions:
            for key, value in fields.items():
                output[key + f"_{i}"] = value
        return output

    def validate(self, df: pd.DataFrame) -> None:
        """Validate multiway contact dataframe"""
        self._schema.validate(df)
