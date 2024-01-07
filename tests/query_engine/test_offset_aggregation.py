"""Tests for the aggregation functions in the query engine."""
from itertools import product

import numpy as np
import pandas as pd
import pytest

from spoc.pixels import Pixels
from spoc.query_engine import AggregationFunction
from spoc.query_engine import Anchor
from spoc.query_engine import BasicQuery
from spoc.query_engine import OffsetAggregation
from spoc.query_engine import OffsetMode
from spoc.query_engine import RegionOffsetTransformation
from spoc.query_engine import Snipper


@pytest.fixture(name="pixels_with_offset")
def pixels_with_offset_fixture(pixels_with_single_region):
    """Pixels with single region"""
    return RegionOffsetTransformation()(pixels_with_single_region)


@pytest.fixture(name="complete_synthetic_pixels_df")
def complete_synthetic_pixels_df_fixture():
    """Pixels that span two regions densely"""
    np.random.seed(42)
    # genomic region_1
    pixels_1 = [
        {
            "chrom": tup[0],
            "start_1": tup[1],
            "start_2": tup[2],
            "start_3": tup[3],
            "count": np.random.randint(0, 10),
        }
        for tup in product(
            ["chr1"],
            np.arange(900_000, 1_150_000, 50_000),
            np.arange(900_000, 1_150_000, 50_000),
            np.arange(900_000, 1_150_000, 50_000),
        )
    ]
    # genomic region_2
    pixels_2 = [
        {
            "chrom": tup[0],
            "start_1": tup[1],
            "start_2": tup[2],
            "start_3": tup[3],
            "count": np.random.randint(0, 10),
        }
        for tup in product(
            ["chr2"],
            np.arange(900_000, 1_150_000, 50_000),
            np.arange(900_000, 1_150_000, 50_000),
            np.arange(900_000, 1_150_000, 50_000),
        )
    ]
    return pd.concat((pd.DataFrame(pixels_1), pd.DataFrame(pixels_2)))


@pytest.fixture(name="incomplete_synthetic_pixels_df")
def incomplete_synthetic_pixels_df_fixture():
    """Pixels that span two regions sparsely"""
    np.random.seed(42)
    # genomic region 1
    pixels_1 = [
        {
            "chrom": tup[0],
            "start_1": tup[1],
            "start_2": tup[2],
            "start_3": tup[3],
            "count": np.random.randint(0, 10),
        }
        for tup in product(
            ["chr1"],
            np.arange(900_000, 1_000_000, 50_000),
            np.arange(900_000, 1_000_000, 50_000),
            np.arange(900_000, 1_000_000, 50_000),
        )
    ]
    # genomic region_2
    pixels_2 = [
        {
            "chrom": tup[0],
            "start_1": tup[1],
            "start_2": tup[2],
            "start_3": tup[3],
            "count": np.random.randint(0, 10),
        }
        for tup in product(
            ["chr2"],
            np.arange(1_000_000, 1_150_000, 50_000),
            np.arange(1_000_000, 1_150_000, 50_000),
            np.arange(1_000_000, 1_150_000, 50_000),
        )
    ]
    return pd.concat((pd.DataFrame(pixels_1), pd.DataFrame(pixels_2)))


@pytest.fixture
def single_region():
    """Single region"""
    return pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [900_000],
            "end": [1_100_000],
        },
        index=[0],
    )


@pytest.fixture
def single_region_not_binaligned():
    """Single region that is not aligned to the bin size"""
    return pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [920_000],
            "end": [1_120_000],
        },
        index=[0],
    )


@pytest.fixture
def two_regions():
    """Two regions"""
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "start": [900_000, 900_000],
            "end": [1_100_000, 1_100_000],
        }
    )


@pytest.mark.parametrize(
    "genomic_data_fixture",
    [
        "contacts_without_regions",
        "pixels_without_regions",
        "contacts_with_single_region",
        "contacts_with_multiple_regions",
        "pixels_with_single_region",
        "pixels_with_multiple_regions",
    ],
)
def test_input_wo_offset_rejected(genomic_data_fixture, request):
    """Test that the validation fails for incorrect inputs."""
    genomic_data = request.getfixturevalue(genomic_data_fixture)
    with pytest.raises(ValueError):
        query = BasicQuery(
            query_plan=[
                OffsetAggregation(
                    value_column="value", function=AggregationFunction.COUNT
                ),
            ],
        )
        query.query(genomic_data)


def test_input_wo_data_column_rejected(pixels_with_offset):
    """Test that the validation fails for incorrect inputs."""
    with pytest.raises(ValueError):
        query = BasicQuery(
            query_plan=[
                OffsetAggregation(
                    value_column="test_column_does_no_exist",
                    function=AggregationFunction.AVG,
                ),
            ],
        )
        query.query(pixels_with_offset)


@pytest.mark.parametrize(
    "aggregation_spoc, aggregation_pandas, region_fixture",
    [
        (AggregationFunction.COUNT, "count", "single_region"),
        (AggregationFunction.COUNT, "count", "two_regions"),
        (AggregationFunction.COUNT, "count", "single_region_not_binaligned"),
        (AggregationFunction.SUM, "sum", "single_region"),
        (AggregationFunction.SUM, "sum", "two_regions"),
        (AggregationFunction.SUM, "sum", "single_region_not_binaligned"),
        (AggregationFunction.AVG, "mean", "single_region"),
        (AggregationFunction.AVG, "mean", "two_regions"),
        (AggregationFunction.AVG, "mean", "single_region_not_binaligned"),
    ],
)
def test_aggregations_on_dense_input(
    complete_synthetic_pixels_df,
    aggregation_spoc,
    aggregation_pandas,
    region_fixture,
    request,
):
    """Test sum aggregation on dense input."""
    # setup (pixels here are points to make the test easier)
    pixels = Pixels(complete_synthetic_pixels_df, binsize=50_000, number_fragments=3)
    region = request.getfixturevalue(region_fixture)
    mapped_pixels = BasicQuery(
        query_plan=[
            Snipper(region, anchor_mode=Anchor(mode="ANY")),
            RegionOffsetTransformation(offset_mode=OffsetMode.LEFT),
        ],
    ).query(pixels)
    mapped_pixels_df = mapped_pixels.load_result()
    cat_dtype = pd.CategoricalDtype(range(-100_000, 150_000, 50_000))
    mapped_pixels_df["offset_1"] = mapped_pixels_df["offset_1"].astype(cat_dtype)
    mapped_pixels_df["offset_2"] = mapped_pixels_df["offset_2"].astype(cat_dtype)
    mapped_pixels_df["offset_3"] = mapped_pixels_df["offset_3"].astype(cat_dtype)
    expected_aggregation = (
        mapped_pixels_df.groupby(["offset_1", "offset_2", "offset_3"])
        .agg(count=("count", aggregation_pandas))
        .astype(float)
        .reset_index()
        .rename(columns={"count": f"count_{aggregation_pandas}"})
        .sort_values(["offset_1", "offset_2", "offset_3"])
    )
    # execute aggregation
    query = BasicQuery(
        query_plan=[
            OffsetAggregation(
                value_column="count", function=aggregation_spoc, densify_output=False
            ),
        ],
    )
    actual_aggregation = query.query(mapped_pixels).load_result()
    # test
    np.testing.assert_array_almost_equal(
        expected_aggregation.values,
        actual_aggregation.values,
    )


@pytest.mark.parametrize(
    "aggregation_spoc, aggregation_pandas, region_fixture",
    [
        (AggregationFunction.COUNT, "count", "single_region"),
        (AggregationFunction.COUNT, "count", "two_regions"),
        (AggregationFunction.COUNT, "count", "single_region_not_binaligned"),
        (AggregationFunction.SUM, "sum", "single_region"),
        (AggregationFunction.SUM, "sum", "two_regions"),
        (AggregationFunction.SUM, "sum", "single_region_not_binaligned"),
        (AggregationFunction.AVG, "mean", "single_region"),
        (AggregationFunction.AVG, "mean", "two_regions"),
        (AggregationFunction.AVG, "mean", "single_region_not_binaligned"),
    ],
)
def test_aggregations_on_sparse_input(
    incomplete_synthetic_pixels_df,
    aggregation_spoc,
    aggregation_pandas,
    region_fixture,
    request,
):
    """Test sum aggregation on dense input."""
    # setup
    pixels = Pixels(incomplete_synthetic_pixels_df, binsize=50_000, number_fragments=3)
    region = request.getfixturevalue(region_fixture)
    mapped_pixels = BasicQuery(
        query_plan=[
            Snipper(region, anchor_mode=Anchor(mode="ANY")),
            RegionOffsetTransformation(offset_mode=OffsetMode.LEFT),
        ],
    ).query(pixels)
    mapped_pixels_df = mapped_pixels.load_result()
    cat_dtype = pd.CategoricalDtype(range(-100_000, 150_000, 50_000))
    mapped_pixels_df["offset_1"] = mapped_pixels_df["offset_1"].astype(cat_dtype)
    mapped_pixels_df["offset_2"] = mapped_pixels_df["offset_2"].astype(cat_dtype)
    mapped_pixels_df["offset_3"] = mapped_pixels_df["offset_3"].astype(cat_dtype)
    expected_aggregation = (
        mapped_pixels_df.groupby(["offset_1", "offset_2", "offset_3"])
        .agg(count=("count", aggregation_pandas))
        .astype(float)
        .reset_index()
        .rename(columns={"count": f"count_{aggregation_pandas}"})
        .sort_values(["offset_1", "offset_2", "offset_3"])
    )
    # execute aggregation
    query = BasicQuery(
        query_plan=[
            OffsetAggregation(
                value_column="count", function=aggregation_spoc, densify_output=True
            ),
        ],
    )
    actual_aggregation = query.query(mapped_pixels).load_result()
    # test
    np.testing.assert_array_almost_equal(
        expected_aggregation.values,
        actual_aggregation.values,
    )
