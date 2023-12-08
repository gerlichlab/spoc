"""Tests for the pixel selection"""
import pytest
import pandas as pd
import dask.dataframe as dd
import duckdb
from spoc.query_engine import BasicQuery, Anchor, Snipper
from spoc.pixels import Pixels
from spoc.io import DUCKDB_CONNECTION

@pytest.fixture(name="pixel_dataframe")
def pixel_dataframe_fixture():
    """A dataframe containing pixels"""
    return pd.DataFrame(
        {
            "chrom": ["chr1"] * 4,
            "start_1": [100_000, 10_000_000, 100_000, 10_000_000],
            "start_2": [500_000, 25_000_000, 500_000, 25_000_000],
            "start_3": [600_000, 6_000_000, 600_000, 6_000_000],
            "contact_count": [1, 2, 3, 4],
        }
    )

@pytest.fixture(name="pixels_dask")
def pixels_dask_fixture(pixel_dataframe):
    """A dask dataframe containing pixels"""
    return Pixels(dd.from_pandas(pixel_dataframe, npartitions=2))

@pytest.fixture(name="pixels_pandas")
def pixels_pandas_fixture(pixel_dataframe):
    """A pandas dataframe containing pixels"""
    return Pixels(pixel_dataframe)

@pytest.fixture(name="pixels_duckdb")
def pixels_duckdb_fixture(pixel_dataframe):
    """A duckdb dataframe containing pixels"""
    return Pixels(duckdb.from_df(pixel_dataframe, DUCKDB_CONNECTION))

# happy path

def test_no_filter_returns_all_pixels(contact_fixture, request):
    """Test that no filter returns all pixels"""

def test_any_anchor_region_returns_correct_pixels(
    contact_fixture, single_region, request
):
    """Test that any anchor region returns correct pixels"""


def test_all_anchor_regions_returns_correct_pixels(
    contact_fixture, single_region, request
):
    """Test that all anchor regions returns correct pixels"""

def test_specific_anchor_regions_returns_correct_pixels(
    contact_fixture, anchors, expected_reads, single_region_2, request
):
    """Test that specific anchor regions returns correct pixels"""

def test_any_anchor_region_returns_correct_pixels_multi_region(
    contact_fixture, multi_region, request
):
    """Test that any anchor region returns correct pixels"""


def test_all_anchor_regions_returns_correct_pixels_multi_region(
    contact_fixture, multi_region, request
):
    """Test that all anchor regions returns correct pixels"""

def test_pixels_duplicated_for_multiple_overlapping_regions(
    contact_fixture, multi_region_2, request
):
    """
    This test verifies that when multiple overlapping regions are specified as anchor regions,
    the query returns duplicated pixels for each overlapping region.
    """

def test_specific_anchor_regions_returns_correct_pixels_multi_region(
    contact_fixture, anchors, expected_reads, multi_region, request
):
    """Test that specific anchor regions returns correct pixels"""

# validation problems

def test_specific_anchor_region_not_in_pixels_raises_error(
    contact_fixture, single_region, request
):
    """Test that specific anchor region not in pixels raises error"""

