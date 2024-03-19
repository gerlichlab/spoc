"""These set of tests test selection of contacts overlapping
with a list of mulitple regions"""
import dask.dataframe as dd
import duckdb
import pytest

from spoc.contacts import Contacts
from spoc.io import DUCKDB_CONNECTION
from spoc.query_engine import Anchor
from spoc.query_engine import Overlap
from spoc.query_engine import Query


@pytest.fixture(name="example_2d_contacts_pandas")
def example_2d_contacts_pandas_fixture(example_2d_df):
    """Example 2d contacts"""
    return Contacts(example_2d_df)


@pytest.fixture(name="example_2d_contacts_dask")
def example_2d_contacts_dask_fixture(example_2d_df):
    """Example 2d contacts"""
    return Contacts(dd.from_pandas(example_2d_df, npartitions=2))


@pytest.fixture(name="example_2d_contacts_duckdb")
def example_2d_contacts_duckdb_fixture(example_2d_df):
    """Example 2d contacts"""
    return Contacts(duckdb.from_df(example_2d_df, connection=DUCKDB_CONNECTION))


def test_different_half_window_size_throws_error(single_region, single_region_2):
    """Test that different half window size throws error"""
    with pytest.raises(ValueError):
        Overlap(
            regions=[single_region, single_region_2],
            anchor_mode=Anchor(fragment_mode="ANY"),
        )


def test_negative_half_window_size_throws_error(single_region, single_region_2):
    """Test that negative half window size throws error"""
    with pytest.raises(ValueError):
        Overlap(
            regions=[single_region, single_region_2],
            anchor_mode=Anchor(fragment_mode="ANY"),
            half_window_size=-1,
        )


@pytest.mark.parametrize(
    "fragment_mode,region_mode,number_contacts",
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
    number_contacts,
    single_region,
    single_region_2,
    example_2d_contacts_pandas,
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
    result = query.build(example_2d_contacts_pandas)

    # assert
    assert result.compute().shape[0] == number_contacts


@pytest.mark.parametrize(
    "fragment_mode,region_mode,number_contacts",
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
    number_contacts,
    single_region,
    single_region_2,
    example_2d_contacts_pandas,
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
    result = query.build(example_2d_contacts_pandas)

    # assert
    assert result.compute().shape[0] == number_contacts


@pytest.mark.parametrize(
    "add_overlap_columns,number_contacts",
    [
        (True, 5),
        (False, 3),
    ],
)
def test_duplicates_after_overlap_handled_correctly(
    add_overlap_columns,
    number_contacts,
    multi_region_2,
    single_region,
    example_2d_contacts_pandas,
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
    result = query.build(example_2d_contacts_pandas)

    # assert
    assert result.compute().shape[0] == number_contacts
    if not add_overlap_columns:
        assert len(result.compute().filter(regex="region").columns) == 0
