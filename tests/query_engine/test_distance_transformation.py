"""Tests for transformations"""
from __future__ import annotations

import numpy as np
import pytest

from spoc.query_engine import DistanceMode
from spoc.query_engine import DistanceTransformation
from spoc.query_engine import Query


@pytest.mark.parametrize(
    "genomic_data_fixture",
    [
        "contacts_without_regions",
        "pixels_without_regions",
        "contacts_with_multiple_regions_overlapped",
    ],
)
def test_incompatible_input_rejected(genomic_data_fixture, request):
    """Tests that incompatible input raises a ValueError"""
    genomic_data = request.getfixturevalue(genomic_data_fixture)
    with pytest.raises(ValueError):
        query = Query(
            query_steps=[
                DistanceTransformation(DistanceMode.MIDPOINT),
            ],
        )
        query.build(genomic_data)


@pytest.mark.parametrize(
    "genomic_data_fixture",
    ["contacts_with_single_region", "contacts_with_multiple_regions"],
)
def test_distance_calculated_correctly_contacts(genomic_data_fixture, request):
    """Tests that the distance is calculated correctly for contacts"""
    genomic_data = request.getfixturevalue(genomic_data_fixture)
    query = Query(
        query_steps=[
            DistanceTransformation(DistanceMode.MIDPOINT),
        ],
    )
    result = query.build(genomic_data).compute()
    # check that the distance is correct
    region_midpoints = (result["region_start"] + result["region_end"]) // 2
    position_1_midpoint = (result["start_1"] + result["end_1"]) // 2
    position_2_midpoint = (result["start_2"] + result["end_2"]) // 2
    assert np.allclose(
        position_1_midpoint - region_midpoints,
        result["distance_1"],
    )
    assert np.allclose(
        position_2_midpoint - region_midpoints,
        result["distance_2"],
    )


@pytest.mark.parametrize(
    "genomic_data_fixture",
    ["pixels_with_single_region", "pixels_with_multiple_regions"],
)
def test_distance_midpoint_rejected_pixels(genomic_data_fixture, request):
    """Tests that the distance calculation is rejected with midpoint distance mode for pixels."""
    genomic_data = request.getfixturevalue(genomic_data_fixture)
    query = Query(
        query_steps=[
            DistanceTransformation(DistanceMode.MIDPOINT),
        ],
    )
    with pytest.raises(ValueError):
        query.build(genomic_data)


@pytest.mark.parametrize(
    "genomic_data_fixture",
    ["pixels_with_single_region", "pixels_with_multiple_regions"],
)
def test_distance_pixels(genomic_data_fixture, request):
    """Tests distance calculation succeeds for pixels with distancemode left."""
    genomic_data = request.getfixturevalue(genomic_data_fixture)
    query = Query(
        query_steps=[
            DistanceTransformation(DistanceMode.LEFT),
        ],
    )
    result = query.build(genomic_data).compute()
    # check that the distance is correct
    region_midpoints = (result["region_start"] + result["region_end"]) // 2
    position_1_midpoint = result["start_1"]
    position_2_midpoint = result["start_2"]
    assert np.allclose(
        position_1_midpoint - region_midpoints,
        result["distance_1"],
    )
    assert np.allclose(
        position_2_midpoint - region_midpoints,
        result["distance_2"],
    )


@pytest.mark.parametrize(
    "genomic_data_fixture",
    ["pixels_with_single_region", "pixels_with_multiple_regions"],
)
def test_distance_pixels_left_string_arg(genomic_data_fixture, request):
    """Tests distance calculation succeeds for pixels with distancemode left, specified as string."""
    genomic_data = request.getfixturevalue(genomic_data_fixture)
    query = Query(
        query_steps=[
            DistanceTransformation("left"),
        ],
    )
    result = query.build(genomic_data).compute()
    # check that the distance is correct
    region_midpoints = (result["region_start"] + result["region_end"]) // 2
    position_1_midpoint = result["start_1"]
    position_2_midpoint = result["start_2"]
    assert np.allclose(
        position_1_midpoint - region_midpoints,
        result["distance_1"],
    )
    assert np.allclose(
        position_2_midpoint - region_midpoints,
        result["distance_2"],
    )
