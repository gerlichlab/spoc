from duckdb import aggregate
import pandas as pd
import numpy as np
from itertools import product
import pytest
from spoc.snipping.snipping_strategies import Triplet1DSnippingStrategy, SnippingValues

@pytest.fixture
def complete_synthetic_pixels():
    np.random.seed(42)
    pixels = [
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
    return pd.DataFrame(pixels)

@pytest.fixture
def incomplete_synthetic_pixels():
    np.random.seed(42)
    pixels = [
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
    return pd.DataFrame(pixels)

@pytest.fixture
def expected_incomplete_synthetic_region(incomplete_synthetic_pixels):
    return incomplete_synthetic_pixels.assign(
        offset_1=lambda df_: (df_.start_1 - 1_000_000)//1000,
        offset_2=lambda df_: (df_.start_2 - 1_000_000)//1000
    ).astype(
        {
            "offset_1": pd.CategoricalDtype(np.arange(-100, 150, 50)), # needed to set possible values
            "offset_2": pd.CategoricalDtype(np.arange(-100, 150, 50)),
        }
    )\
    .pivot_table(index="offset_1", columns="offset_2", values="contact_count", aggfunc=np.sum).astype(float)


@pytest.fixture
def expected_entire_synthetic_region(complete_synthetic_pixels):
    return complete_synthetic_pixels.assign(
        offset_1=lambda df_: (df_.start_1 - 1_000_000)//1000,
        offset_2=lambda df_: (df_.start_2 - 1_000_000)//1000
    ).pivot_table(index="offset_1", columns="offset_2", values="contact_count", aggfunc=np.sum).astype(float)


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


def test_triplet_1d_entire_region_with_complete_pixels(complete_synthetic_pixels, single_region, expected_entire_synthetic_region):
    """Test whether snipping of an entire region produces correct results."""
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
    np.array_equal(result.values, expected_entire_synthetic_region.values)


def test_triplet_1d_entire_region_with_incomplete_pixels(incomplete_synthetic_pixels, single_region, expected_incomplete_synthetic_region):
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
    np.array_equal(result.values, expected_incomplete_synthetic_region.values)