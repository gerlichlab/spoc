"""Tests for transformations"""
from __future__ import annotations

import numpy as np
import pytest

from spoc.contacts import Contacts
from spoc.pixels import Pixels
from spoc.query_engine import Anchor
from spoc.query_engine import BasicQuery
from spoc.query_engine import RegionOffsetTransformation
from spoc.query_engine import Snipper


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
    return Snipper(single_region, anchor_mode=Anchor(mode="ANY"))(
        contacts_without_regions,
    )


@pytest.fixture(name="contacts_with_multiple_regions")
def contacts_with_multiple_regions_fixture(contacts_without_regions, multi_region):
    """Contacts with multiple regions"""
    return Snipper(multi_region, anchor_mode=Anchor(mode="ANY"))(
        contacts_without_regions,
    )


@pytest.fixture(name="pixels_with_single_region")
def pixels_with_single_region_fixture(pixels_without_regions, single_region):
    """Pixels with single region"""
    return Snipper(single_region, anchor_mode=Anchor(mode="ANY"))(
        pixels_without_regions,
    )


@pytest.fixture(name="pixels_with_multiple_regions")
def pixels_with_multiple_regions_fixture(pixels_without_regions, multi_region):
    """Pixels with multiple regions"""
    return Snipper(multi_region, anchor_mode=Anchor(mode="ANY"))(pixels_without_regions)


@pytest.mark.parametrize(
    "genomic_data_fixture",
    [
        "contacts_without_regions",
        "pixels_without_regions",
    ],
)
def test_incompatible_input_rejected(genomic_data_fixture, request):
    """Tests that incompatible input raises a ValueError"""
    genomic_data = request.getfixturevalue(genomic_data_fixture)
    with pytest.raises(ValueError):
        query = BasicQuery(
            query_plan=[
                RegionOffsetTransformation(),
            ],
        )
        query.query(genomic_data)


@pytest.mark.parametrize(
    "genomic_data_fixture",
    [
        "contacts_with_single_region",
        "contacts_with_multiple_regions",
        "pixels_with_single_region",
        "pixels_with_multiple_regions",
    ],
)
def test_offset_calculated_correctly(genomic_data_fixture, request):
    """Tests that the offset is calculated correctly for contacts"""
    genomic_data = request.getfixturevalue(genomic_data_fixture)
    query = BasicQuery(
        query_plan=[
            RegionOffsetTransformation(),
        ],
    )
    result = query.query(genomic_data).load_result()
    # check that the offset is correct
    region_midpoints = (result["region_start"] + result["region_end"]) // 2
    position_1_midpoint = (result["start_1"] + result["end_1"]) // 2
    position_2_midpoint = (result["start_2"] + result["end_2"]) // 2
    assert np.allclose(
        position_1_midpoint - region_midpoints,
        result["offset_1"],
    )
    assert np.allclose(
        position_2_midpoint - region_midpoints,
        result["offset_2"],
    )
