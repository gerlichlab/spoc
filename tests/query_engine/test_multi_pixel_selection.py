"""These set of tests test selection of pixels overlapping
with a list of mulitple regions"""
import pytest

from spoc.pixels import Pixels
from spoc.query_engine import Anchor
from spoc.query_engine import Overlap
from spoc.query_engine import Query


@pytest.fixture(name="pixels_pandas")
def pixels_pandas_fixture(pixel_dataframe):
    """A pandas dataframe containing pixels"""
    return Pixels(pixel_dataframe, number_fragments=2, binsize=10)


@pytest.mark.parametrize(
    "fragment_mode,region_mode,number_pixels",
    [
        ("ALL", "ALL", 0),
        ("ANY", "ALL", 0),
        ("ANY", "ANY", 4),
        ("ALL", "ANY", 1),
    ],
)
def test_overlap_without_position_subset(
    fragment_mode,
    region_mode,
    number_pixels,
    single_region,
    single_region_2,
    pixels_pandas,
):
    """Test that overlap without position subset"""
    # setup
    query_plan = [
        Overlap(
            regions=[single_region, single_region_2],
            anchor_mode=Anchor(fragment_mode=fragment_mode, region_mode=region_mode),
            half_window_size=100,
        )
    ]
    query = Query(query_steps=query_plan)

    # run
    result = query.build(pixels_pandas)

    # assert
    assert result.compute().shape[0] == number_pixels


@pytest.mark.parametrize(
    "fragment_mode,region_mode,number_pixels",
    [
        ("ALL", "ALL", 0),
        ("ANY", "ALL", 0),
        ("ANY", "ANY", 3),
        ("ALL", "ANY", 3),
    ],
)
def test_overlap_with_position_subset(
    fragment_mode,
    region_mode,
    number_pixels,
    single_region,
    single_region_2,
    pixels_pandas,
):
    """Test that overlap with position subset"""
    # setup
    query_plan = [
        Overlap(
            regions=[single_region, single_region_2],
            anchor_mode=Anchor(
                fragment_mode=fragment_mode, region_mode=region_mode, positions=[1]
            ),
            half_window_size=100,
        )
    ]
    query = Query(query_steps=query_plan)

    # run
    result = query.build(pixels_pandas)

    # assert
    assert result.compute().shape[0] == number_pixels


@pytest.mark.parametrize(
    "add_overlap_columns,number_pixels",
    [
        (True, 5),
        (False, 3),
    ],
)
def test_duplicates_after_overlap_handled_correctly(
    add_overlap_columns,
    number_pixels,
    multi_region_2,
    single_region,
    pixels_pandas,
):
    """Test that duplicates after overlap are handled correctly"""
    # setup
    query_plan = [
        Overlap(
            regions=[multi_region_2, single_region],
            anchor_mode=Anchor(fragment_mode="ANY", region_mode="ANY"),
            half_window_size=100,
            add_overlap_columns=add_overlap_columns,
        )
    ]
    query = Query(query_steps=query_plan)

    # run
    result = query.build(pixels_pandas)

    # assert
    assert result.compute().shape[0] == number_pixels
    if not add_overlap_columns:
        assert len(result.compute().filter(regex="region").columns) == 0
