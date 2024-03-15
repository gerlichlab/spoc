"""Shared fixtures for query engine tests."""
from __future__ import annotations

import pandas as pd
import pytest

from spoc.contacts import Contacts
from spoc.pixels import Pixels
from spoc.query_engine import Anchor
from spoc.query_engine import Overlap


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


@pytest.fixture(name="single_region_3")
def single_region_3_fixture():
    """Single region"""
    return pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [750],
            "end": [850],
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


@pytest.fixture(name="contacts_without_regions")
def contacts_without_regions_fixture(example_2d_df):
    """Example 2d contacts"""
    return Contacts(example_2d_df)


@pytest.fixture(name="pixels_without_regions")
def pixels_wihtout_regions_fixture(pixel_dataframe):
    """Pixels without regions"""
    return Pixels(pixel_dataframe, number_fragments=2, binsize=10)


@pytest.fixture(name="contacts_with_single_region")
def contacts_with_single_region_fixture(contacts_without_regions, single_region):
    """Contacts with single region"""
    return Overlap(single_region, anchor_mode=Anchor(fragment_mode="ANY"))(
        contacts_without_regions,
    )


@pytest.fixture(name="contacts_with_multiple_regions")
def contacts_with_multiple_regions_fixture(contacts_without_regions, multi_region):
    """Contacts with multiple regions"""
    return Overlap(multi_region, anchor_mode=Anchor(fragment_mode="ANY"))(
        contacts_without_regions,
    )


@pytest.fixture(name="pixels_with_single_region")
def pixels_with_single_region_fixture(pixels_without_regions, single_region):
    """Pixels with single region"""
    return Overlap(single_region, anchor_mode=Anchor(fragment_mode="ANY"))(
        pixels_without_regions,
    )


@pytest.fixture(name="pixels_with_multiple_regions")
def pixels_with_multiple_regions_fixture(pixels_without_regions, multi_region):
    """Pixels with multiple regions"""
    return Overlap(multi_region, anchor_mode=Anchor(fragment_mode="ANY"))(
        pixels_without_regions
    )
