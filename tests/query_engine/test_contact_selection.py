"""Tests for contact selection"""
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


# happy path


@pytest.mark.parametrize(
    "contact_fixture",
    [
        "example_2d_contacts_pandas",
        "example_2d_contacts_duckdb",
    ],
)
def test_no_filter_returns_all_contacts(contact_fixture, request):
    """Test that no filter returns all contacts"""
    contacts = request.getfixturevalue(contact_fixture)
    query = Query(query_steps=[])
    result = query.build(contacts)
    assert result.compute().shape[0] == 4


@pytest.mark.parametrize(
    "contact_fixture",
    [
        "example_2d_contacts_pandas",
        "example_2d_contacts_dask",
        "example_2d_contacts_duckdb",
    ],
)
def test_any_anchor_region_returns_correct_contacts(
    contact_fixture, single_region, request
):
    """Test that any anchor region returns correct contacts"""
    # setup
    contacts = request.getfixturevalue(contact_fixture)
    query_plan = [Overlap(regions=single_region, anchor_mode=Anchor(mode="ANY"))]
    # execution
    query = Query(query_steps=query_plan)
    result = query.build(contacts)
    # test
    assert result.compute().shape[0] == 2
    assert sorted(result.compute().read_name.tolist()) == sorted(["read1", "read2"])


@pytest.mark.parametrize(
    "contact_fixture",
    [
        "example_2d_contacts_pandas",
        "example_2d_contacts_dask",
        "example_2d_contacts_duckdb",
    ],
)
def test_all_anchor_regions_returns_correct_contacts(
    contact_fixture, single_region, request
):
    """Test that all anchor regions returns correct contacts"""
    # setup
    contacts = request.getfixturevalue(contact_fixture)
    query_plan = [Overlap(regions=single_region, anchor_mode=Anchor(mode="ALL"))]
    # execution
    query = Query(query_steps=query_plan)
    result = query.build(contacts)
    # test
    assert result.compute().shape[0] == 1
    assert sorted(result.compute().read_name.tolist()) == sorted(["read2"])


@pytest.mark.parametrize(
    "contact_fixture,anchors,expected_reads",
    [
        (source_data, anchors, expected_reads)
        for source_data, anchors, expected_reads in zip(
            [
                "example_2d_contacts_pandas",
                "example_2d_contacts_dask",
                "example_2d_contacts_duckdb",
            ]
            * 3,
            [[1]] * 3 + [[2]] * 3 + [[1, 2]] * 3,
            [["read3"]] * 3 + [["read4"]] * 3 + [[]] * 3,
        )
    ],
)
def test_specific_anchor_regions_returns_correct_contacts(
    contact_fixture, anchors, expected_reads, single_region_2, request
):
    """Test that specific anchor regions returns correct contacts"""
    # setup
    contacts = request.getfixturevalue(contact_fixture)
    query_plan = [
        Overlap(
            regions=single_region_2, anchor_mode=Anchor(mode="ALL", anchors=anchors)
        )
    ]
    # execution
    query = Query(query_steps=query_plan)
    result = query.build(contacts)
    # test
    assert result.compute().shape[0] == len(expected_reads)
    assert sorted(result.compute().read_name.tolist()) == sorted(expected_reads)


@pytest.mark.parametrize(
    "contact_fixture,anchors,expected_reads",
    [
        (source_data, anchors, expected_reads)
        for source_data, anchors, expected_reads in zip(
            [
                "example_2d_contacts_pandas",
                "example_2d_contacts_dask",
                "example_2d_contacts_duckdb",
            ]
            * 3,
            [[1]] * 3 + [[2]] * 3 + [[1, 2]] * 3,
            [["read3"]] * 3 + [["read4"]] * 3 + [[]] * 3,
        )
    ],
)
def test_specific_anchor_regions_returns_correct_contacts_point_region(
    contact_fixture, anchors, expected_reads, single_region_3, request
):
    """Test that specific anchor regions returns correct contacts, when
    the region is a point region and windowsize is passed"""
    # setup
    contacts = request.getfixturevalue(contact_fixture)
    query_plan = [
        Overlap(
            regions=single_region_3,
            anchor_mode=Anchor(mode="ALL", anchors=anchors),
            window_size=50,
        )
    ]
    # execution
    query = Query(query_steps=query_plan)
    result = query.build(contacts)
    # test
    assert result.compute().shape[0] == len(expected_reads)
    assert sorted(result.compute().read_name.tolist()) == sorted(expected_reads)
    assert result.get_schema().get_window_size() == 50


@pytest.mark.parametrize(
    "contact_fixture",
    [
        "example_2d_contacts_pandas",
        "example_2d_contacts_dask",
        "example_2d_contacts_duckdb",
    ],
)
def test_any_anchor_region_returns_correct_contacts_multi_region(
    contact_fixture, multi_region, request
):
    """Test that any anchor region returns correct contacts"""
    # setup
    contacts = request.getfixturevalue(contact_fixture)
    query_plan = [Overlap(regions=multi_region, anchor_mode=Anchor(mode="ANY"))]
    # execution
    query = Query(query_steps=query_plan)
    result = query.build(contacts)
    # test
    assert result.compute().shape[0] == 4
    assert sorted(result.compute().read_name.tolist()) == sorted(
        ["read1", "read2", "read3", "read4"]
    )


@pytest.mark.parametrize(
    "contact_fixture",
    [
        "example_2d_contacts_pandas",
        "example_2d_contacts_dask",
        "example_2d_contacts_duckdb",
    ],
)
def test_all_anchor_regions_returns_correct_contacts_multi_region(
    contact_fixture, multi_region, request
):
    """Test that all anchor regions returns correct contacts"""
    # setup
    contacts = request.getfixturevalue(contact_fixture)
    query_plan = [Overlap(regions=multi_region, anchor_mode=Anchor(mode="ALL"))]
    # execution
    query = Query(query_steps=query_plan)
    result = query.build(contacts)
    # test
    assert result.compute().shape[0] == 1
    assert sorted(result.compute().read_name.tolist()) == sorted(["read2"])


@pytest.mark.parametrize(
    "contact_fixture",
    [
        "example_2d_contacts_pandas",
        "example_2d_contacts_dask",
        "example_2d_contacts_duckdb",
    ],
)
def test_contacts_duplicated_for_multiple_overlapping_regions(
    contact_fixture, multi_region_2, request
):
    """
    This test verifies that when multiple overlapping regions are specified as anchor regions,
    the query returns duplicated contacts for each overlapping region.
    """
    # setup
    contacts = request.getfixturevalue(contact_fixture)
    query_plan = [Overlap(regions=multi_region_2, anchor_mode=Anchor(mode="ALL"))]
    # execution
    query = Query(query_steps=query_plan)
    result = query.build(contacts)
    # test
    assert result.compute().shape[0] == 2
    assert sorted(result.compute().read_name.tolist()) == sorted(["read2", "read2"])


@pytest.mark.parametrize(
    "contact_fixture,anchors,expected_reads",
    [
        (source_data, anchors, expected_reads)
        for source_data, anchors, expected_reads in zip(
            [
                "example_2d_contacts_pandas",
                "example_2d_contacts_dask",
                "example_2d_contacts_duckdb",
            ]
            * 3,
            [[1]] * 3 + [[2]] * 3 + [[1, 2]] * 3,
            [["read1", "read2", "read3"]] * 3
            + [["read2", "read4"]] * 3
            + [["read2"]] * 3,
        )
    ],
)
def test_specific_anchor_regions_returns_correct_contacts_multi_region(
    contact_fixture, anchors, expected_reads, multi_region, request
):
    """Test that specific anchor regions returns correct contacts"""
    # setup
    contacts = request.getfixturevalue(contact_fixture)
    query_plan = [
        Overlap(regions=multi_region, anchor_mode=Anchor(mode="ALL", anchors=anchors))
    ]
    # execution
    query = Query(query_steps=query_plan)
    result = query.build(contacts)
    # test
    assert result.compute().shape[0] == len(expected_reads)
    assert sorted(result.compute().read_name.tolist()) == sorted(expected_reads)


# validation problems


@pytest.mark.parametrize(
    "contact_fixture",
    [
        "example_2d_contacts_pandas",
        "example_2d_contacts_dask",
        "example_2d_contacts_duckdb",
    ],
)
def test_specific_anchor_region_not_in_contacts_raises_error(
    contact_fixture, single_region, request
):
    """Test that specific anchor region not in contacts raises error"""
    # setup
    contacts = request.getfixturevalue(contact_fixture)
    query_plan = [
        Overlap(regions=single_region, anchor_mode=Anchor(mode="ALL", anchors=[3]))
    ]
    with pytest.raises(ValueError):
        query = Query(query_steps=query_plan)
        query.build(contacts)
