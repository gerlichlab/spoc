import os
import shutil
import pandas as pd
import numpy as np
from itertools import product
import pytest
from spoc.snipping.snipping_strategies import Triplet1DSnippingStrategy, SnippingValues
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
            "contact_count": np.random.randint(0, 10)
        }
        for tup in product(
                            ['chr1'],
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
            "contact_count": np.random.randint(0, 10)
        }
        for tup in product(
                            ['chr2'],
                            np.arange(900_000, 1_150_000, 50_000),
                            np.arange(900_000, 1_150_000, 50_000),
                            np.arange(900_000, 1_150_000, 50_000),
                        )
    ]
    return pd.concat(
        (
            pd.DataFrame(pixels_1),
            pd.DataFrame(pixels_2)
        )
    )

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
            "contact_count": np.random.randint(0, 10)
        }
        for tup in product(
                            ['chr1'],
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
            "contact_count": np.random.randint(0, 10)
        }
        for tup in product(
                            ['chr2'],
                            np.arange(900_000, 1_000_000, 50_000),
                            np.arange(900_000, 1_000_000, 50_000),
                            np.arange(900_000, 1_000_000, 50_000),
                        )
    ]
    return pd.concat(
        (
            pd.DataFrame(pixels_1),
            pd.DataFrame(pixels_2)
        )
    )

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
    return complete_synthetic_pixels.assign(
        offset_1=lambda df_: (df_.start_1 - 1_000_000)//1000,
        offset_2=lambda df_: (df_.start_2 - 1_000_000)//1000
    ).query("chrom == 'chr1'").pivot_table(index="offset_1", columns="offset_2", values="contact_count", aggfunc=np.sum).astype(float)


@pytest.fixture
def expected_entire_single_region_incomplete(incomplete_synthetic_pixels):
    return incomplete_synthetic_pixels.assign(
        offset_1=lambda df_: (df_.start_1 - 1_000_000)//1000,
        offset_2=lambda df_: (df_.start_2 - 1_000_000)//1000
    ).astype(
        {
            "offset_1": pd.CategoricalDtype(np.arange(-100, 150, 50)), # needed to set possible values
            "offset_2": pd.CategoricalDtype(np.arange(-100, 150, 50)), # needed to set possible values
        }
    ).query("chrom == 'chr1'").pivot_table(index="offset_1", columns="offset_2", values="contact_count", aggfunc=np.sum).astype(float)



@pytest.fixture
def expected_two_entire_regions_complete(complete_synthetic_pixels):
    # region1
    first = complete_synthetic_pixels.assign(
        offset_1=lambda df_: (df_.start_1 - 1_000_000)//1000,
        offset_2=lambda df_: (df_.start_2 - 1_000_000)//1000
    ).query("chrom == 'chr1'").groupby(["offset_1", "offset_2"]).contact_count.sum()
    # region2
    second = complete_synthetic_pixels.assign(
        offset_1=lambda df_: (df_.start_1 - 1_000_000)//1000,
        offset_2=lambda df_: (df_.start_2 - 1_000_000)//1000
    ).query("chrom == 'chr2'").groupby(["offset_1", "offset_2"]).contact_count.sum()
    # average and return
    return pd.concat((
        first, second
    )).reset_index().pivot_table(index="offset_1", columns="offset_2", values="contact_count", aggfunc=np.mean).astype(float)

@pytest.fixture
def single_region():
    return pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [1_000_000],
            "end": [1_000_000],
        },
        index=[0]
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


# Tests

def test_triplet_1d_entire_region_with_complete_pixels(complete_synthetic_pixels, single_region, expected_entire_single_region_complete):
    """Test whether snipping of an entire region produces correct results when pixels are supplied as dataframe."""
    # instantiate
    snip_strat = Triplet1DSnippingStrategy(
        bin_size=50_000,
        snipping_value=SnippingValues.ICCF
    )
    # do the snipping
    result = snip_strat.snip(
        complete_synthetic_pixels,
        single_region,
        100_000,
        position_slack=1_000_000 # this selects the entire region
    )
    # check result
    np.testing.assert_array_almost_equal(result.values, expected_entire_single_region_complete.values)

def test_triplet_1d_entire_region_with_complete_pixels_from_file(complete_persisted_pixels, single_region, expected_entire_single_region_complete):
    """Test whether snipping of an entire region produces correct results when pixels are supplied as file."""
    # instantiate
    snip_strat = Triplet1DSnippingStrategy(
        bin_size=50_000,
        snipping_value=SnippingValues.ICCF
    )
    # do the snipping
    result = snip_strat.snip(
        complete_persisted_pixels,
        single_region,
        100_000,
        position_slack=1_000_000 # this selects the entire region
    )
    # check result
    np.testing.assert_array_almost_equal(result.values, expected_entire_single_region_complete.values)

def test_triplet_1d_two_regions_with_complete_pixels(complete_synthetic_pixels, two_regions, expected_two_entire_regions_complete):
    """Test whether snipping of an entire region produces correct results when pixels are supplied as dataframe."""
    # instantiate
    snip_strat = Triplet1DSnippingStrategy(
        bin_size=50_000,
        snipping_value=SnippingValues.ICCF
    )
    # do the snipping
    result = snip_strat.snip(
        complete_synthetic_pixels,
        two_regions,
        100_000,
        position_slack=1_000_000 # this selects the entire region
    )
    # check result
    np.testing.assert_array_almost_equal(result.values, expected_two_entire_regions_complete.values)


def test_triplet_1d_entire_region_with_incomplete_pixels(incomplete_synthetic_pixels, single_region, expected_entire_single_region_incomplete):
    """Test whether snipping of an entire region produces correct results."""
    # instantiate
    snip_strat = Triplet1DSnippingStrategy(
        bin_size=50_000,
        snipping_value=SnippingValues.ICCF
    )
    # do the snipping
    result = snip_strat.snip(
        incomplete_synthetic_pixels,
        single_region,
        100_000,
        position_slack=1_000_000 # this selects the entire region
    )

    # check result
    np.testing.assert_array_almost_equal(result.values, expected_entire_single_region_incomplete.values)