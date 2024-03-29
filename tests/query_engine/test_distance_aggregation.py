"""Tests for the aggregation functions in the query engine."""
from itertools import product

import numpy as np
import pandas as pd
import pytest

from spoc.pixels import Pixels
from spoc.query_engine import AggregationFunction
from spoc.query_engine import Anchor
from spoc.query_engine import DistanceAggregation
from spoc.query_engine import DistanceMode
from spoc.query_engine import DistanceTransformation
from spoc.query_engine import Overlap
from spoc.query_engine import Query


@pytest.fixture(name="pixels_with_distance")
def pixels_with_distance_fixture(pixels_with_single_region):
    """Pixels with single region"""
    return DistanceTransformation()(pixels_with_single_region)


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


@pytest.fixture(name="incomplete_synthetic_pixels_dense_df")
def incomplete_synthetic_pixels_dense_df_fixture(
    complete_synthetic_pixels_df, incomplete_synthetic_pixels_df
):
    """Pixels that span two regions sparsely
    with missing pixels filled with 0."""
    return incomplete_synthetic_pixels_df.merge(
        complete_synthetic_pixels_df[["chrom", "start_1", "start_2", "start_3"]],
        on=["chrom", "start_1", "start_2", "start_3"],
        how="outer",
    ).fillna(0)


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
def test_input_wo_distance_rejected(genomic_data_fixture, request):
    """Test that the validation fails for incorrect inputs."""
    genomic_data = request.getfixturevalue(genomic_data_fixture)
    with pytest.raises(ValueError):
        query = Query(
            query_steps=[
                DistanceAggregation(
                    value_column="value", function=AggregationFunction.COUNT
                ),
            ],
        )
        query.build(genomic_data)


def test_input_wo_data_column_rejected(pixels_with_distance):
    """Test that the validation fails for incorrect inputs."""
    with pytest.raises(ValueError):
        query = Query(
            query_steps=[
                DistanceAggregation(
                    value_column="test_column_does_no_exist",
                    function=AggregationFunction.AVG,
                ),
            ],
        )
        query.build(pixels_with_distance)


def test_contacts_rejected(contacts_with_single_region):
    """Test that the validation fails for incorrect inputs."""
    with pytest.raises(ValueError):
        query = Query(
            query_steps=[
                DistanceAggregation(
                    value_column="value", function=AggregationFunction.COUNT
                ),
            ],
        )
        query.build(contacts_with_single_region)


@pytest.mark.parametrize(
    "aggregation_spoc, aggregation_pandas, region_fixture",
    [
        (AggregationFunction.COUNT, "count", "single_region"),
        (AggregationFunction.COUNT, "count", "two_regions"),
        (AggregationFunction.COUNT, "count", "single_region_not_binaligned"),
        ("COUNT", "count", "single_region_not_binaligned"),
        (AggregationFunction.SUM, "sum", "single_region"),
        (AggregationFunction.SUM, "sum", "two_regions"),
        (AggregationFunction.SUM, "sum", "single_region_not_binaligned"),
        ("SUM", "sum", "single_region_not_binaligned"),
        (AggregationFunction.AVG, "mean", "single_region"),
        (AggregationFunction.AVG, "mean", "two_regions"),
        (AggregationFunction.AVG, "mean", "single_region_not_binaligned"),
        ("AVG", "mean", "single_region_not_binaligned"),
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
    mapped_pixels = Query(
        query_steps=[
            Overlap(region, anchor_mode=Anchor(mode="ANY")),
            DistanceTransformation(distance_mode=DistanceMode.LEFT),
        ],
    ).build(pixels)
    mapped_pixels_df = mapped_pixels.compute()
    cat_dtype = pd.CategoricalDtype(range(-100_000, 150_000, 50_000))
    mapped_pixels_df["distance_1"] = mapped_pixels_df["distance_1"].astype(cat_dtype)
    mapped_pixels_df["distance_2"] = mapped_pixels_df["distance_2"].astype(cat_dtype)
    mapped_pixels_df["distance_3"] = mapped_pixels_df["distance_3"].astype(cat_dtype)
    expected_aggregation = (
        mapped_pixels_df.groupby(["distance_1", "distance_2", "distance_3"])
        .agg(count=("count", aggregation_pandas))
        .astype(float)
        .reset_index()
        .rename(columns={"count": f"count_{aggregation_pandas}"})
        .sort_values(["distance_1", "distance_2", "distance_3"])
    )
    # execute aggregation
    query = Query(
        query_steps=[
            DistanceAggregation(
                value_column="count", function=aggregation_spoc, densify_output=False
            ),
        ],
    )
    actual_aggregation = query.build(mapped_pixels).compute()
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
        ("COUNT", "count", "single_region_not_binaligned"),
        (AggregationFunction.SUM, "sum", "single_region"),
        (AggregationFunction.SUM, "sum", "two_regions"),
        (AggregationFunction.SUM, "sum", "single_region_not_binaligned"),
        ("SUM", "sum", "single_region_not_binaligned"),
        (AggregationFunction.AVG, "mean", "single_region"),
        (AggregationFunction.AVG, "mean", "two_regions"),
        (AggregationFunction.AVG, "mean", "single_region_not_binaligned"),
        ("AVG", "mean", "single_region_not_binaligned"),
    ],
)
def test_aggregations_on_dense_input_with_reduced_dimensionality(
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
    mapped_pixels = Query(
        query_steps=[
            Overlap(region, anchor_mode=Anchor(mode="ANY")),
            DistanceTransformation(distance_mode=DistanceMode.LEFT),
        ],
    ).build(pixels)
    mapped_pixels_df = mapped_pixels.compute()
    cat_dtype = pd.CategoricalDtype(range(-100_000, 150_000, 50_000))
    mapped_pixels_df["distance_1"] = mapped_pixels_df["distance_1"].astype(cat_dtype)
    mapped_pixels_df["distance_2"] = mapped_pixels_df["distance_2"].astype(cat_dtype)
    expected_aggregation = (
        mapped_pixels_df.groupby(["distance_1", "distance_2"])
        .agg(count=("count", aggregation_pandas))
        .astype(float)
        .reset_index()
        .rename(columns={"count": f"count_{aggregation_pandas}"})
        .sort_values(["distance_1", "distance_2"])
    )
    # execute aggregation
    query = Query(
        query_steps=[
            DistanceAggregation(
                value_column="count",
                function=aggregation_spoc,
                densify_output=False,
                position_list=[1, 2],
            ),
        ],
    )
    actual_aggregation = query.build(mapped_pixels).compute()
    # test
    np.testing.assert_array_almost_equal(
        expected_aggregation.values,
        actual_aggregation.values,
    )


# pylint: disable=too-many-arguments
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
        (AggregationFunction.AVG_WITH_EMPTY, "mean", "single_region"),
        (AggregationFunction.AVG_WITH_EMPTY, "mean", "two_regions"),
        (AggregationFunction.AVG_WITH_EMPTY, "mean", "single_region_not_binaligned"),
    ],
)
def test_aggregations_on_sparse_input(
    incomplete_synthetic_pixels_df,
    incomplete_synthetic_pixels_dense_df,
    aggregation_spoc,
    aggregation_pandas,
    region_fixture,
    request,
):
    """Test sum aggregation on dense input."""
    # setup
    incomplete_pixels = Pixels(
        incomplete_synthetic_pixels_df, binsize=50_000, number_fragments=3
    )
    incomplete_dense_pixels = Pixels(
        incomplete_synthetic_pixels_dense_df, binsize=50_000, number_fragments=3
    )
    query_plan = Query(
        query_steps=[
            Overlap(
                request.getfixturevalue(region_fixture), anchor_mode=Anchor(mode="ANY")
            ),
            DistanceTransformation(distance_mode=DistanceMode.LEFT),
        ],
    )
    mapped_pixels = query_plan.build(incomplete_pixels)
    mapped_incomplete_dense_pixels_df = query_plan.build(
        incomplete_dense_pixels
    ).compute()
    if aggregation_spoc == AggregationFunction.AVG_WITH_EMPTY:
        # when we test the AVG_WITH_EMPTY function, we need to use the dense pixels
        # where missing values with 0 count are filled in
        pixel_frame_for_expected = mapped_incomplete_dense_pixels_df
    else:
        pixel_frame_for_expected = mapped_pixels.compute()
    pixel_frame_for_expected["distance_1"] = pixel_frame_for_expected[
        "distance_1"
    ].astype(pd.CategoricalDtype(range(-100_000, 150_000, 50_000)))
    pixel_frame_for_expected["distance_2"] = pixel_frame_for_expected[
        "distance_2"
    ].astype(pd.CategoricalDtype(range(-100_000, 150_000, 50_000)))
    pixel_frame_for_expected["distance_3"] = pixel_frame_for_expected[
        "distance_3"
    ].astype(pd.CategoricalDtype(range(-100_000, 150_000, 50_000)))
    expected_aggregation = (
        pixel_frame_for_expected.groupby(["distance_1", "distance_2", "distance_3"])
        .agg(count=("count", aggregation_pandas))
        .astype(float)
        .reset_index()
        .rename(columns={"count": f"count_{aggregation_pandas}"})
        .sort_values(["distance_1", "distance_2", "distance_3"])
    )
    # execute aggregation
    query = Query(
        query_steps=[
            DistanceAggregation(
                value_column="count", function=aggregation_spoc, densify_output=True
            ),
        ],
    )
    # test
    np.testing.assert_array_almost_equal(
        expected_aggregation.values, query.build(mapped_pixels).compute()
    )


# pylint: disable=too-many-arguments
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
        # for aggregation function with empty, we need to sum up and divide by region number
        # just taking the mean with take the mean with respect to every triplet pixel
        (AggregationFunction.AVG_WITH_EMPTY, lambda s: s.sum(), "single_region"),
        (AggregationFunction.AVG_WITH_EMPTY, lambda s: s.sum() / 2, "two_regions"),
        (
            AggregationFunction.AVG_WITH_EMPTY,
            lambda s: s.sum(),
            "single_region_not_binaligned",
        ),
    ],
)
def test_aggregations_on_sparse_input_with_reduced_dimensionality(
    incomplete_synthetic_pixels_df,
    incomplete_synthetic_pixels_dense_df,
    aggregation_spoc,
    aggregation_pandas,
    region_fixture,
    request,
):
    """Test aggregation on sparse input with reduced dimensionality."""
    # setup
    incomplete_pixels = Pixels(
        incomplete_synthetic_pixels_df, binsize=50_000, number_fragments=3
    )
    incomplete_dense_pixels = Pixels(
        incomplete_synthetic_pixels_dense_df, binsize=50_000, number_fragments=3
    )
    query_plan = Query(
        query_steps=[
            Overlap(
                request.getfixturevalue(region_fixture), anchor_mode=Anchor(mode="ANY")
            ),
            DistanceTransformation(distance_mode=DistanceMode.LEFT),
        ],
    )
    mapped_pixels = query_plan.build(incomplete_pixels)
    mapped_incomplete_dense_pixels_df = query_plan.build(
        incomplete_dense_pixels
    ).compute()
    if aggregation_spoc == AggregationFunction.AVG_WITH_EMPTY:
        # when we test the AVG_WITH_EMPTY function, we need to use the dense pixels
        # where missing values with 0 count are filled in
        pixel_frame_for_expected = mapped_incomplete_dense_pixels_df
    else:
        pixel_frame_for_expected = mapped_pixels.compute()
    pixel_frame_for_expected["distance_1"] = pixel_frame_for_expected[
        "distance_1"
    ].astype(pd.CategoricalDtype(range(-100_000, 150_000, 50_000)))
    pixel_frame_for_expected["distance_2"] = pixel_frame_for_expected[
        "distance_2"
    ].astype(pd.CategoricalDtype(range(-100_000, 150_000, 50_000)))
    expected_aggregation = (
        pixel_frame_for_expected.groupby(["distance_1", "distance_2"])
        .agg(count=("count", aggregation_pandas))
        .astype(float)
        .reset_index()
        .rename(columns={"count": f"count_{aggregation_pandas}"})
        .sort_values(["distance_1", "distance_2"])
    )
    # execute aggregation
    query = Query(
        query_steps=[
            DistanceAggregation(
                value_column="count",
                function=aggregation_spoc,
                densify_output=True,
                position_list=[1, 2],
            ),
        ],
    )
    # test
    np.testing.assert_array_almost_equal(
        expected_aggregation.values, query.build(mapped_pixels).compute()
    )
