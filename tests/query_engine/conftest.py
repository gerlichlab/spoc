"""Shared fixtures for query engine tests."""
from __future__ import annotations

import pandas as pd
import pytest


@pytest.fixture(name="example_2d_df")
def example_2d_df_fixture():
    """Example 2d contacts"""
    return pd.DataFrame(
        {
            "chrom_1": ["chr1", "chr1", "chr1", "chr1"],
            "start_1": [100, 100, 750, 400],
            "end_1": [200, 200, 780, 500],
            "mapping_quality_1": [30, 30, 30, 30],
            "align_score_1": [100, 100, 100, 100],
            "align_base_qscore_1": [100, 100, 100, 100],
            "chrom_2": ["chr1", "chr1", "chr1", "chr1"],
            "start_2": [300, 100, 300, 750],
            "end_2": [400, 200, 400, 780],
            "mapping_quality_2": [30, 30, 30, 30],
            "align_score_2": [100, 100, 100, 100],
            "align_base_qscore_2": [100, 100, 100, 100],
            # read name serves as id
            "read_name": ["read1", "read2", "read3", "read4"],
            "read_length": [100, 100, 100, 100],
        },
    )


@pytest.fixture(name="pixel_dataframe")
def pixel_dataframe_fixture():
    """A dataframe containing pixels"""
    return pd.DataFrame(
        {
            "chrom": ["chr1"] * 4,
            "start_1": [180, 180, 750, 400],
            "start_2": [300, 180, 300, 750],
            "count": [1, 2, 3, 4],  # contact count is id
        },
    )


@pytest.fixture(name="single_region")
def single_region_fixture():
    """Single region"""
    return pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [150],
            "end": [200],
        },
    )


@pytest.fixture(name="single_region_2")
def single_region_2_fixture():
    """Single region"""
    return pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [700],
            "end": [800],
        },
    )


@pytest.fixture(name="multi_region")
def multi_region_fixture():
    """Multi region"""
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [150, 700],
            "end": [200, 800],
        },
    )


@pytest.fixture(name="multi_region_2")
def multi_region_2_fixture():
    """Multi region"""
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [150, 180],
            "end": [200, 220],
        },
    )
