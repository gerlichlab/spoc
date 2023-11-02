"""Tests for the contacts module."""

# pylint: disable=redefined-outer-name
import pytest
import pandas as pd
import pandera as pa
import numpy as np
import dask.dataframe as dd

from spoc import contacts, fragments

# pytlint: disable=unused-import
from ..fixtures.symmetry import (
    unlabelled_contacts_2d,
    labelled_binary_contacts_2d_sorted,
)


@pytest.fixture
def triplet_expander():
    """expander for triplets"""
    return fragments.FragmentExpander(number_fragments=3, contains_metadata=False)


@pytest.fixture
def triplet_expander_labelled():
    """expander for triplets"""
    return fragments.FragmentExpander(number_fragments=3, contains_metadata=True)


@pytest.fixture
def contact_manipulator():
    """manipulator for triplest"""
    return contacts.ContactManipulator()


@pytest.fixture
def bad_df():
    """bad df for testing"""
    return pd.DataFrame({"be": ["bop"]})


@pytest.fixture
def labelled_df():
    """Dataframe representing a labelled fragment file"""
    return pd.DataFrame(
        {
            "chrom": ["chr1"] * 6,
            "start": [1, 2, 3, 4, 5, 6],
            "end": [4, 5, 6, 7, 8, 9],
            "strand": [True] * 6,
            "read_name": ["dummy"] * 4 + ["dummy2"] * 2,
            "read_start": [1, 2, 3, 4, 5, 6],
            "read_end": [4, 5, 6, 7, 8, 9],
            "read_length": [1] * 6,
            "mapping_quality": [1, 2, 3, 4, 5, 6],
            "align_score": [1, 2, 3, 4, 5, 6],
            "align_base_qscore": [1, 2, 3, 4, 5, 6],
            "pass_filter": [True] * 6,
            "metadata": [
                "SisterA",
                "SisterB",
                "SisterA",
                "SisterB",
                "SisterA",
                "SisterB",
            ],
        }
    )


@pytest.fixture
def unlabelled_df():
    """Dataframe representing an unlabelled fragment file"""
    return pd.DataFrame(
        {
            "chrom": ["chr1"] * 6,
            "start": [1, 2, 3, 4, 5, 6],
            "end": [4, 5, 6, 7, 8, 9],
            "strand": [True] * 6,
            "read_name": ["dummy"] * 4 + ["dummy2"] * 2,
            "read_start": [1, 2, 3, 4, 5, 6],
            "read_end": [4, 5, 6, 7, 8, 9],
            "read_length": [1] * 6,
            "mapping_quality": [1, 2, 3, 4, 5, 6],
            "align_score": [1, 2, 3, 4, 5, 6],
            "align_base_qscore": [1, 2, 3, 4, 5, 6],
            "pass_filter": [True] * 6,
        }
    )


@pytest.fixture
def labelled_fragments(labelled_df):
    """labelled fragments"""
    return fragments.Fragments(labelled_df)


@pytest.fixture
def labelled_fragments_dask(labelled_df):
    """labelled fragments from a dask dataframe"""
    return fragments.Fragments(dd.from_pandas(labelled_df, npartitions=1))


@pytest.fixture
def unlabelled_fragments(unlabelled_df):
    """unlabelled fragments"""
    return fragments.Fragments(unlabelled_df)


@pytest.fixture
def unlabelled_fragments_dask(unlabelled_df):
    """unlabelled fragments from a dask dataframe"""
    return fragments.Fragments(dd.from_pandas(unlabelled_df, npartitions=1))


@pytest.mark.parametrize(
    "fragments, expander",
    [
        ("labelled_fragments", "triplet_expander_labelled"),
        ("labelled_fragments_dask", "triplet_expander_labelled"),
        ("unlabelled_fragments", "triplet_expander"),
        ("unlabelled_fragments_dask", "triplet_expander"),
    ],
)
def test_expander_drops_reads_w_too_little_fragments(expander, fragments, request):
    """Tests whether expander drops reads with too little fragments"""
    triplet_expander = request.getfixturevalue(expander)
    result = triplet_expander.expand(request.getfixturevalue(fragments)).data
    if isinstance(result, dd.DataFrame):
        result = result.compute()
    assert len(set(result.read_name)) == 1
    assert result.read_name[0] == "dummy"


@pytest.mark.parametrize(
    "fragments, expander",
    [
        ("labelled_fragments", "triplet_expander_labelled"),
        ("labelled_fragments_dask", "triplet_expander_labelled"),
        ("unlabelled_fragments", "triplet_expander"),
        ("unlabelled_fragments_dask", "triplet_expander"),
    ],
)
def test_expander_returns_correct_number_of_contacts(expander, fragments, request):
    """Tests whether expander returns correct number of contacts"""
    triplet_expander = request.getfixturevalue(expander)
    result = triplet_expander.expand(request.getfixturevalue(fragments)).data
    assert len(result) == 4


@pytest.mark.parametrize("fragments", ["labelled_fragments", "labelled_fragments_dask"])
def test_expander_returns_correct_contacts_labelled(
    triplet_expander_labelled, fragments, request
):
    """Tests whether expander returns correct contacts for labelled fragments"""
    df = request.getfixturevalue(fragments)
    result = triplet_expander_labelled.expand(df).data
    if isinstance(result, dd.DataFrame):
        result = result.compute()
    assert np.array_equal(result["start_1"].values, np.array([1, 1, 1, 2]))
    assert np.array_equal(result["end_1"].values, np.array([4, 4, 4, 5]))
    assert np.array_equal(result["start_2"].values, np.array([2, 2, 3, 3]))
    assert np.array_equal(result["end_2"].values, np.array([5, 5, 6, 6]))
    assert np.array_equal(result["start_3"].values, np.array([3, 4, 4, 4]))
    assert np.array_equal(result["end_3"].values, np.array([6, 7, 7, 7]))
    assert np.array_equal(
        result["metadata_1"].values,
        np.array(["SisterA", "SisterA", "SisterA", "SisterB"]),
    )
    assert np.array_equal(
        result["metadata_2"].values,
        np.array(["SisterB", "SisterB", "SisterA", "SisterA"]),
    )
    assert np.array_equal(
        result["metadata_3"].values,
        np.array(["SisterA", "SisterB", "SisterB", "SisterB"]),
    )


@pytest.mark.parametrize(
    "fragments", ["unlabelled_fragments", "unlabelled_fragments_dask"]
)
def test_expander_returns_correct_contacts_unlabelled(
    triplet_expander, fragments, request
):
    """Tests whether expander returns correct contacts for unlabelled fragments"""
    df = request.getfixturevalue(fragments)
    result = triplet_expander.expand(df).data
    if isinstance(result, dd.DataFrame):
        result = result.compute()
    assert np.array_equal(result["start_1"].values, np.array([1, 1, 1, 2]))
    assert np.array_equal(result["end_1"].values, np.array([4, 4, 4, 5]))
    assert np.array_equal(result["start_2"].values, np.array([2, 2, 3, 3]))
    assert np.array_equal(result["end_2"].values, np.array([5, 5, 6, 6]))
    assert np.array_equal(result["start_3"].values, np.array([3, 4, 4, 4]))
    assert np.array_equal(result["end_3"].values, np.array([6, 7, 7, 7]))
    assert "metadata_1" not in result.columns


def test_contacts_constructor_rejects_wrong_df(bad_df):
    """Tests whether contacts constructor rejects wrong df"""
    with pytest.raises(pa.errors.SchemaError):
        contacts.Contacts(bad_df, number_fragments=3)


def test_merge_works_for_good_pandas_df(
    triplet_expander, contact_manipulator, labelled_fragments
):
    """Tests whether merge works for good pandas df"""
    contacts = triplet_expander.expand(labelled_fragments)
    result = contact_manipulator.merge_contacts([contacts, contacts]).data
    assert result.shape[0] == 8
    assert result.shape[1] == contacts.data.shape[1]


def test_merge_works_for_good_dask_df(
    triplet_expander, contact_manipulator, labelled_fragments
):
    """Tests whether merge works for good dask df"""
    cont = triplet_expander.expand(labelled_fragments)
    contacts_dask = contacts.Contacts(
        dd.from_pandas(cont.data, npartitions=1), number_fragments=3
    )
    result = contact_manipulator.merge_contacts(
        [contacts_dask, contacts_dask]
    ).data.compute()
    assert result.shape[0] == 8
    assert result.shape[1] == cont.data.shape[1]


def test_merge_fails_for_pandas_dask_mixed(
    triplet_expander, contact_manipulator, labelled_fragments
):
    """Tests whether merge fails for pandas dask mixed"""
    with pytest.raises(ValueError):
        contacts_pandas = triplet_expander.expand(labelled_fragments)
        contacts_dask = contacts.Contacts(
            dd.from_pandas(contacts_pandas.data, npartitions=1), number_fragments=3
        )
        contact_manipulator.merge_contacts([contacts_pandas, contacts_dask])


def test_subset_metadata_fails_if_not_labelled(unlabelled_contacts_2d):
    """Tests whether subset fails if the datafrane is not laelled"""
    contact_manipulator = contacts.ContactManipulator()
    unlab_contacts = contacts.Contacts(unlabelled_contacts_2d)
    with pytest.raises(AssertionError):
        contact_manipulator.subset_on_metadata(unlab_contacts, ["A", "B"])


def test_subset_metadata_fails_if_pattern_longer_than_number_fragments(
    labelled_binary_contacts_2d_sorted,
):
    """Tests whether subset fails if the pattern is longer than number of fragments"""
    contact_manipulator = contacts.ContactManipulator()
    lab_contacts = contacts.Contacts(labelled_binary_contacts_2d_sorted)
    with pytest.raises(AssertionError):
        contact_manipulator.subset_on_metadata(lab_contacts, ["A", "B", "A"])


def test_subset_metadata_fails_if_pattern_contains_strings_not_in_metadata(
    labelled_binary_contacts_2d_sorted,
):
    """Tests whether subset fails if the pattern contains strings not in metadata"""
    contact_manipulator = contacts.ContactManipulator()
    lab_contacts = contacts.Contacts(labelled_binary_contacts_2d_sorted)
    with pytest.raises(AssertionError):
        contact_manipulator.subset_on_metadata(lab_contacts, ["A", "C"])


def test_subset_metadata_creates_correct_subset(labelled_binary_contacts_2d_sorted):
    """Tests whether subset creates the correct subset"""
    contact_manipulator = contacts.ContactManipulator()
    lab_contacts = contacts.Contacts(labelled_binary_contacts_2d_sorted)
    result = contact_manipulator.subset_on_metadata(lab_contacts, ["A", "B"])
    assert len(result.data) == 2
    assert result.data["metadata_1"].unique() == ["A"]
    assert result.data["metadata_2"].unique() == ["B"]
    assert result.metadata_combi == ["A", "B"]


# TODO: merge rejects labelled and unlabelled contacts

# TODO: merge fails for contacts of different oder
