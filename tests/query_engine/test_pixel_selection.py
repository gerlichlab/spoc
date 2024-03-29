"""Tests for the pixel selection"""
import dask.dataframe as dd
import duckdb
import numpy as np
import pytest

from spoc.io import DUCKDB_CONNECTION
from spoc.pixels import Pixels
from spoc.query_engine import Anchor
from spoc.query_engine import Overlap
from spoc.query_engine import Query


@pytest.fixture(name="pixels_dask")
def pixels_dask_fixture(pixel_dataframe):
    """A dask dataframe containing pixels"""
    return Pixels(
        dd.from_pandas(pixel_dataframe, npartitions=2),
        number_fragments=2,
        binsize=10,
    )


@pytest.fixture(name="pixels_pandas")
def pixels_pandas_fixture(pixel_dataframe):
    """A pandas dataframe containing pixels"""
    return Pixels(pixel_dataframe, number_fragments=2, binsize=10)


@pytest.fixture(name="pixels_duckdb")
def pixels_duckdb_fixture(pixel_dataframe):
    """A duckdb dataframe containing pixels"""
    return Pixels(
        duckdb.from_df(pixel_dataframe, DUCKDB_CONNECTION),
        number_fragments=2,
        binsize=10,
    )


# happy path


@pytest.mark.parametrize(
    "pixels_fixture",
    [
        "pixels_pandas",
        "pixels_duckdb",
    ],
)
def test_no_filter_returns_all_pixels(pixels_fixture, request):
    """Test that no filter returns all pixels"""
    pixels = request.getfixturevalue(pixels_fixture)
    query = Query(query_steps=[])
    result = query.build(pixels)
    assert result.compute().shape[0] == 4


@pytest.mark.parametrize(
    "pixels_fixture",
    [
        "pixels_pandas",
        "pixels_dask",
        "pixels_duckdb",
    ],
)
def test_any_anchor_region_returns_correct_pixels(
    pixels_fixture, single_region, request
):
    """Test that any anchor region returns correct pixels"""
    # setup
    pixels = request.getfixturevalue(pixels_fixture)
    query_plan = [Overlap(regions=single_region, anchor_mode=Anchor(mode="ANY"))]
    # execution
    query = Query(query_steps=query_plan)
    result = query.build(pixels)
    # test correct selection
    assert result.compute().shape[0] == 2
    assert sorted(result.compute()["count"].tolist()) == sorted([1, 2])
    # test addition of end columns
    assert np.allclose(result.compute()["start_1"] + 10, result.compute()["end_1"])
    assert np.allclose(result.compute()["start_2"] + 10, result.compute()["end_2"])


@pytest.mark.parametrize(
    "pixels_fixture",
    [
        "pixels_pandas",
        "pixels_dask",
        "pixels_duckdb",
    ],
)
def test_all_anchor_regions_returns_correct_pixels(
    pixels_fixture, single_region, request
):
    """Test that all anchor regions returns correct pixels"""
    # setup
    pixels = request.getfixturevalue(pixels_fixture)
    query_plan = [Overlap(regions=single_region, anchor_mode=Anchor(mode="ALL"))]
    # execution
    query = Query(query_steps=query_plan)
    result = query.build(pixels)
    # test
    assert result.compute().shape[0] == 1
    assert sorted(result.compute()["count"].tolist()) == sorted([2])


@pytest.mark.parametrize(
    "pixel_fixture,anchors,expected_reads",
    [
        (source_data, anchors, expected_reads)
        for source_data, anchors, expected_reads in zip(
            [
                "pixels_pandas",
                "pixels_dask",
                "pixels_duckdb",
            ]
            * 3,
            [[1]] * 3 + [[2]] * 3 + [[1, 2]] * 3,
            [[3]] * 3 + [[4]] * 3 + [[]] * 3,
        )
    ],
)
def test_specific_anchor_regions_returns_correct_pixels(
    pixel_fixture, anchors, expected_reads, single_region_2, request
):
    """Test that specific anchor regions returns correct pixels"""
    # setup
    pixels = request.getfixturevalue(pixel_fixture)
    query_plan = [
        Overlap(
            regions=single_region_2, anchor_mode=Anchor(mode="ALL", anchors=anchors)
        )
    ]
    # execution
    query = Query(query_steps=query_plan)
    result = query.build(pixels)
    # test
    assert result.compute().shape[0] == len(expected_reads)
    assert sorted(result.compute()["count"].tolist()) == sorted(expected_reads)


@pytest.mark.parametrize(
    "pixels_fixture",
    [
        "pixels_pandas",
        "pixels_dask",
        "pixels_duckdb",
    ],
)
def test_any_anchor_region_returns_correct_pixels_multi_region(
    pixels_fixture, multi_region, request
):
    """Test that any anchor region returns correct pixels"""
    # setup
    pixels = request.getfixturevalue(pixels_fixture)
    query_plan = [Overlap(regions=multi_region, anchor_mode=Anchor(mode="ANY"))]
    # execution
    query = Query(query_steps=query_plan)
    result = query.build(pixels)
    # test
    assert result.compute().shape[0] == 4
    assert sorted(result.compute()["count"].tolist()) == sorted([1, 2, 3, 4])


@pytest.mark.parametrize(
    "pixels_fixture",
    [
        "pixels_pandas",
        "pixels_dask",
        "pixels_duckdb",
    ],
)
def test_all_anchor_regions_returns_correct_pixels_multi_region(
    pixels_fixture, multi_region, request
):
    """Test that all anchor regions returns correct pixels"""
    # setup
    pixels = request.getfixturevalue(pixels_fixture)
    query_plan = [Overlap(regions=multi_region, anchor_mode=Anchor(mode="ALL"))]
    # execution
    query = Query(query_steps=query_plan)
    result = query.build(pixels)
    # test
    assert result.compute().shape[0] == 1
    assert sorted(result.compute()["count"].tolist()) == sorted([2])


@pytest.mark.parametrize(
    "pixels_fixture",
    [
        "pixels_pandas",
        "pixels_dask",
        "pixels_duckdb",
    ],
)
def test_pixels_duplicated_for_multiple_overlapping_regions(
    pixels_fixture, multi_region_2, request
):
    """
    This test verifies that when multiple overlapping regions are specified as anchor regions,
    the query returns duplicated pixels for each overlapping region.
    """
    # setup
    pixels = request.getfixturevalue(pixels_fixture)
    query_plan = [Overlap(regions=multi_region_2, anchor_mode=Anchor(mode="ALL"))]
    # execution
    query = Query(query_steps=query_plan)
    result = query.build(pixels)
    # test
    assert result.compute().shape[0] == 2
    assert sorted(result.compute()["count"].tolist()) == sorted([2, 2])


@pytest.mark.parametrize(
    "pixels_fixture,anchors,expected_reads",
    [
        (source_data, anchors, expected_reads)
        for source_data, anchors, expected_reads in zip(
            [
                "pixels_pandas",
                "pixels_dask",
                "pixels_duckdb",
            ]
            * 3,
            [[1]] * 3 + [[2]] * 3 + [[1, 2]] * 3,
            [[1, 2, 3]] * 3 + [[2, 4]] * 3 + [[2]] * 3,
        )
    ],
)
def test_specific_anchor_regions_returns_correct_pixels_multi_region(
    pixels_fixture, anchors, expected_reads, multi_region, request
):
    """Test that specific anchor regions returns correct pixels"""
    # setup
    pixels = request.getfixturevalue(pixels_fixture)
    query_plan = [
        Overlap(regions=multi_region, anchor_mode=Anchor(mode="ALL", anchors=anchors))
    ]
    # execution
    query = Query(query_steps=query_plan)
    result = query.build(pixels)
    # test
    assert result.compute().shape[0] == len(expected_reads)
    assert sorted(result.compute()["count"].tolist()) == sorted(expected_reads)


# validation problems


@pytest.mark.parametrize(
    "pixels_fixture",
    [
        "pixels_pandas",
        "pixels_dask",
        "pixels_duckdb",
    ],
)
def test_specific_anchor_region_not_in_pixels_raises_error(
    pixels_fixture, single_region, request
):
    """Test that specific anchor region not in pixels raises error"""
    # setup
    pixels = request.getfixturevalue(pixels_fixture)
    query_plan = [
        Overlap(regions=single_region, anchor_mode=Anchor(mode="ALL", anchors=[3]))
    ]
    with pytest.raises(ValueError):
        query = Query(query_steps=query_plan)
        query.build(pixels)
