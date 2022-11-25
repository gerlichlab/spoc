import os
import shutil
import pandas as pd
import numpy as np
from itertools import product
import pytest
from spoc.snipping.snipping_strategies import TripletCCT1DSnippingStrategy, SnippingValues
from spoc.pixels import PersistedPixels

# Fixtures


@pytest.fixture
def complete_synthetic_pixels():
    np.random.seed(42)
    # region_1
    pixels_1 = [
        {
            "chrom": tup[0],
            "start_1": tup[1],
            "start_2": tup[2],
            "start_3": tup[3],
            "contact_count": np.random.randint(0, 10),
        }
        for tup in product(
            ["chr1"],
            np.arange(900_000, 1_150_000, 50_000),
            np.arange(900_000, 1_150_000, 50_000),
            np.arange(900_000, 1_150_000, 50_000),
        )
    ]
    # region_2
    pixels_2 = [
        {
            "chrom": tup[0],
            "start_1": tup[1],
            "start_2": tup[2],
            "start_3": tup[3],
            "contact_count": np.random.randint(0, 10),
        }
        for tup in product(
            ["chr2"],
            np.arange(900_000, 1_150_000, 50_000),
            np.arange(900_000, 1_150_000, 50_000),
            np.arange(900_000, 1_150_000, 50_000),
        )
    ]
    return pd.concat((pd.DataFrame(pixels_1), pd.DataFrame(pixels_2)))


@pytest.fixture
def incomplete_synthetic_pixels():
    np.random.seed(42)
    # region 1
    pixels_1 = [
        {
            "chrom": tup[0],
            "start_1": tup[1],
            "start_2": tup[2],
            "start_3": tup[3],
            "contact_count": np.random.randint(0, 10),
        }
        for tup in product(
            ["chr1"],
            np.arange(900_000, 1_000_000, 50_000),
            np.arange(900_000, 1_000_000, 50_000),
            np.arange(900_000, 1_000_000, 50_000),
        )
    ]
    # region_2
    pixels_2 = [
        {
            "chrom": tup[0],
            "start_1": tup[1],
            "start_2": tup[2],
            "start_3": tup[3],
            "contact_count": np.random.randint(0, 10),
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
def complete_persisted_pixels(complete_synthetic_pixels):
    # setup
    os.mkdir("tmp")
    complete_synthetic_pixels.to_parquet("tmp/csp.parquet")
    yield PersistedPixels("tmp/csp.parquet")
    # teardown
    shutil.rmtree("tmp")


@pytest.fixture
def expected_entire_single_region_complete(complete_synthetic_pixels):
    return (
        complete_synthetic_pixels.assign(
            offset_1=lambda df_: (df_.start_1 - 1_000_000) // 1000,
            offset_2=lambda df_: (df_.start_2 - 1_000_000) // 1000,
        )
        .query("chrom == 'chr1'")
        .pivot_table(
            index="offset_1", columns="offset_2", values="contact_count", aggfunc=np.sum
        )
        .astype(float)
    )


@pytest.fixture
def expected_entire_single_region_incomplete(incomplete_synthetic_pixels):
    return (
        incomplete_synthetic_pixels.assign(
            offset_1=lambda df_: (df_.start_1 - 1_000_000) // 1000,
            offset_2=lambda df_: (df_.start_2 - 1_000_000) // 1000,
        )
        .astype(
            {
                "offset_1": pd.CategoricalDtype(
                    np.arange(-100, 150, 50)
                ),  # needed to set possible values
                "offset_2": pd.CategoricalDtype(
                    np.arange(-100, 150, 50)
                ),  # needed to set possible values
            }
        )
        .query("chrom == 'chr1'")
        .pivot_table(
            index="offset_1", columns="offset_2", values="contact_count", aggfunc=np.sum
        )
        .astype(float)
    )


@pytest.fixture
def expected_two_entire_regions_complete(complete_synthetic_pixels):
    # region1
    first = (
        complete_synthetic_pixels.assign(
            offset_1=lambda df_: (df_.start_1 - 1_000_000) // 1000,
            offset_2=lambda df_: (df_.start_2 - 1_000_000) // 1000,
        )
        .query("chrom == 'chr1'")
        .groupby(["offset_1", "offset_2"])
        .contact_count.sum()
    )
    # region2
    second = (
        complete_synthetic_pixels.assign(
            offset_1=lambda df_: (df_.start_1 - 1_000_000) // 1000,
            offset_2=lambda df_: (df_.start_2 - 1_000_000) // 1000,
        )
        .query("chrom == 'chr2'")
        .groupby(["offset_1", "offset_2"])
        .contact_count.sum()
    )
    # average and return
    return (
        pd.concat((first, second))
        .reset_index()
        .pivot_table(
            index="offset_1",
            columns="offset_2",
            values="contact_count",
            aggfunc=np.mean,
        )
        .astype(float)
    )


@pytest.fixture
def expected_two_entire_regions_incomplete(incomplete_synthetic_pixels):
    # region1
    first = (
        incomplete_synthetic_pixels.assign(
            offset_1=lambda df_: (df_.start_1 - 1_000_000) // 1000,
            offset_2=lambda df_: (df_.start_2 - 1_000_000) // 1000,
        )
        .astype(
            {
                "offset_1": pd.CategoricalDtype(
                    np.arange(-100, 150, 50)
                ),  # needed to set possible values
                "offset_2": pd.CategoricalDtype(
                    np.arange(-100, 150, 50)
                ),  # needed to set possible values
            }
        )
        .query("chrom == 'chr1'")
        .groupby(["offset_1", "offset_2"])
        .contact_count.sum()
    )
    # region2
    second = (
        incomplete_synthetic_pixels.assign(
            offset_1=lambda df_: (df_.start_1 - 1_000_000) // 1000,
            offset_2=lambda df_: (df_.start_2 - 1_000_000) // 1000,
        )
        .astype(
            {
                "offset_1": pd.CategoricalDtype(
                    np.arange(-100, 150, 50)
                ),  # needed to set possible values
                "offset_2": pd.CategoricalDtype(
                    np.arange(-100, 150, 50)
                ),  # needed to set possible values
            }
        )
        .query("chrom == 'chr2'")
        .groupby(["offset_1", "offset_2"])
        .contact_count.sum()
    )
    # average and return
    return (
        pd.concat((first, second))
        .reset_index()
        .pivot_table(
            index="offset_1",
            columns="offset_2",
            values="contact_count",
            aggfunc=np.mean,
        )
        .astype(float)
    )


@pytest.fixture
def expected_one_region_center_selection(complete_synthetic_pixels):
    return (
        complete_synthetic_pixels.assign(
            offset_1=lambda df_: (df_.start_1 - 1_000_000) // 1000,
            offset_2=lambda df_: (df_.start_2 - 1_000_000) // 1000,
            offset_3=lambda df_: (df_.start_3 - 1_000_000) // 1000,
        )
        .query("chrom == 'chr1' and offset_3 == 0")
        .pivot_table(
            index="offset_1", columns="offset_2", values="contact_count", aggfunc=np.sum
        )
        .astype(float)
    )


@pytest.fixture
def expected_one_region_50k_offset_selection(complete_synthetic_pixels):
    return (
        complete_synthetic_pixels.assign(
            offset_1=lambda df_: (df_.start_1 - 1_000_000) // 1000,
            offset_2=lambda df_: (df_.start_2 - 1_000_000) // 1000,
            offset_3=lambda df_: (df_.start_3 - 1_000_000) // 1000,
        )
        .query("chrom == 'chr1' and offset_3 == 50")
        .pivot_table(
            index="offset_1", columns="offset_2", values="contact_count", aggfunc=np.sum
        )
        .astype(float)
    )


@pytest.fixture
def expected_two_regions_center_selection(complete_synthetic_pixels):
    # region1
    first = (
        complete_synthetic_pixels.assign(
            offset_1=lambda df_: (df_.start_1 - 1_000_000) // 1000,
            offset_2=lambda df_: (df_.start_2 - 1_000_000) // 1000,
            offset_3=lambda df_: (df_.start_3 - 1_000_000) // 1000,
        )
        .query("chrom == 'chr1' and offset_3 == 0")
        .groupby(["offset_1", "offset_2"])
        .contact_count.sum()
    )
    # region2
    second = (
        complete_synthetic_pixels.assign(
            offset_1=lambda df_: (df_.start_1 - 1_000_000) // 1000,
            offset_2=lambda df_: (df_.start_2 - 1_000_000) // 1000,
            offset_3=lambda df_: (df_.start_3 - 1_000_000) // 1000,
        )
        .query("chrom == 'chr2' and offset_3 == 0")
        .groupby(["offset_1", "offset_2"])
        .contact_count.sum()
    )
    # average and return
    return (
        pd.concat((first, second))
        .reset_index()
        .pivot_table(
            index="offset_1",
            columns="offset_2",
            values="contact_count",
            aggfunc=np.mean,
        )
        .astype(float)
    )


@pytest.fixture
def single_region():
    return pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [1_000_000],
            "end": [1_000_000],
        },
        index=[0],
    )


@pytest.fixture
def two_regions():
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr2"],
            "start": [1_000_000, 1_000_000],
            "end": [1_000_000, 1_000_000],
        }
    )


@pytest.fixture
def standard_snipping_strategy():
    return TripletCCT1DSnippingStrategy(
        bin_size=50_000,
        half_window_size=100_000,
        snipping_value=SnippingValues.ICCF,
        position_slack=1_000_000,
    )

@pytest.fixture
def standard_snipping_strategy_from_string():
    return TripletCCT1DSnippingStrategy(
        bin_size=50_000,
        half_window_size=100_000,
        snipping_value="iccf",
        position_slack=1_000_000,
    )


@pytest.fixture
def center_snipping_strategy():
    return TripletCCT1DSnippingStrategy(
        bin_size=50_000,
        half_window_size=100_000,
        snipping_value=SnippingValues.ICCF,
        position_slack=0,
    )


@pytest.fixture
def center_snipping_strategy_w_offset():
    return TripletCCT1DSnippingStrategy(
        bin_size=50_000,
        half_window_size=100_000,
        snipping_value=SnippingValues.ICCF,
        position_slack=0,
        relative_offset=50_000,
    )


@pytest.fixture
def standard_snipping_strategy_obs_exp(single_region):
    return TripletCCT1DSnippingStrategy(
        bin_size=50_000,
        half_window_size=100_000,
        snipping_value=SnippingValues.OBSEXP,
        position_slack=1_000_000,
        n_random_regions=100,
        genome="hg19",
    )


# Tests

## ICCF


def test_entire_region_with_complete_pixels(
    complete_synthetic_pixels,
    single_region,
    expected_entire_single_region_complete,
    standard_snipping_strategy,
):
    """Test whether snipping of an entire region produces correct results when pixels are supplied as dataframe."""
    # do the snipping
    result = standard_snipping_strategy.snip(complete_synthetic_pixels, single_region)
    # check result
    np.testing.assert_array_almost_equal(
        result.values, expected_entire_single_region_complete.values
    )


def test_entire_region_with_complete_pixels_strategy_from_string(
    complete_synthetic_pixels,
    single_region,
    expected_entire_single_region_complete,
    standard_snipping_strategy_from_string,
):
    """Test whether snipping of an entire region produces correct results when pixels are supplied as dataframe
    and strategy is instantiated with a string"""
    # do the snipping
    result = standard_snipping_strategy_from_string.snip(complete_synthetic_pixels, single_region)
    # check result
    np.testing.assert_array_almost_equal(
        result.values, expected_entire_single_region_complete.values
    )

def test_entire_region_with_complete_pixels_from_file(
    complete_persisted_pixels,
    single_region,
    expected_entire_single_region_complete,
    standard_snipping_strategy,
):
    """Test whether snipping of an entire region produces correct results when pixels are supplied as file."""
    # do the snipping
    result = standard_snipping_strategy.snip(complete_persisted_pixels, single_region)
    # check result
    np.testing.assert_array_almost_equal(
        result.values, expected_entire_single_region_complete.values
    )


def test_two_regions_with_complete_pixels(
    complete_synthetic_pixels,
    two_regions,
    expected_two_entire_regions_complete,
    standard_snipping_strategy,
):
    """Test whether snipping of two regions produces correct results when pixels are supplied as dataframe."""
    # do the snipping
    result = standard_snipping_strategy.snip(complete_synthetic_pixels, two_regions)
    # check result
    np.testing.assert_array_almost_equal(
        result.values, expected_two_entire_regions_complete.values
    )


def test_entire_region_with_incomplete_pixels(
    incomplete_synthetic_pixels,
    single_region,
    expected_entire_single_region_incomplete,
    standard_snipping_strategy,
):
    """Test whether snipping of an entire region with incomplete pixels (not dense) produces correct results."""
    # do the snipping
    result = standard_snipping_strategy.snip(incomplete_synthetic_pixels, single_region)

    # check result
    np.testing.assert_array_almost_equal(
        result.values, expected_entire_single_region_incomplete.values
    )


def test_two_regions_with_incomplete_pixels(
    incomplete_synthetic_pixels,
    two_regions,
    expected_two_entire_regions_incomplete,
    standard_snipping_strategy,
):
    """Test whether snipping of two entire regions produces correct result."""
    # do the snipping
    result = standard_snipping_strategy.snip(incomplete_synthetic_pixels, two_regions)
    # check result
    np.testing.assert_array_almost_equal(
        result.values, expected_two_entire_regions_incomplete.values
    )


def test_one_region_center_selection(
    complete_synthetic_pixels,
    single_region,
    expected_one_region_center_selection,
    center_snipping_strategy,
):
    """Tests whether snipping of complete region with center selection produces correct results"""
    # do the snipping
    result = center_snipping_strategy.snip(complete_synthetic_pixels, single_region)
    # check result
    np.testing.assert_array_almost_equal(
        result.values, expected_one_region_center_selection.values
    )


def test_two_regions_center_selection(
    complete_synthetic_pixels,
    two_regions,
    expected_two_regions_center_selection,
    center_snipping_strategy,
):
    """Tests whether snipping of two complete regions with center selection produces correct results"""
    # do the snipping
    result = center_snipping_strategy.snip(complete_synthetic_pixels, two_regions)
    # check result
    np.testing.assert_array_almost_equal(
        result.values, expected_two_regions_center_selection.values
    )


def test_relative_offset_works(
    complete_synthetic_pixels,
    single_region,
    expected_one_region_50k_offset_selection,
    center_snipping_strategy_w_offset,
):
    # do the snipping
    result = center_snipping_strategy_w_offset.snip(
        complete_synthetic_pixels, single_region
    )
    # check result
    np.testing.assert_array_almost_equal(
        result.values, expected_one_region_50k_offset_selection.values
    )


## Obs/Exp


def test_one_region_center_selection_obs_exp(
    complete_synthetic_pixels, single_region, standard_snipping_strategy_obs_exp, mocker
):
    """Tests whether snipping of complete region with center selection produces correct results"""
    # patch random coordinates to return same coordinates
    mocker.patch.object(
        TripletCCT1DSnippingStrategy, "_get_random_coordinates", return_value=single_region
    )
    # do the snipping
    result = standard_snipping_strategy_obs_exp.snip(
        complete_synthetic_pixels, single_region
    )
    # check result -> since observed and expected is the same, should be ones everywhere
    np.testing.assert_array_almost_equal(result.values, np.ones(result.values.shape))
