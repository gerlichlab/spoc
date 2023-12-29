"""Tests for the contacts module."""

# pylint: disable=redefined-outer-name
import pytest
import pandas as pd

from spoc import contacts


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


def test_subset_metadata_fails_if_not_labelled(unlabelled_contacts_2d):
    """Tests whether subset fails if the datafrane is not laelled"""
    contact_manipulator = contacts.ContactManipulator()
    unlab_contacts = contacts.Contacts(unlabelled_contacts_2d)
    with pytest.raises(ValueError):
        contact_manipulator.subset_on_metadata(unlab_contacts, ["A", "B"])


def test_subset_metadata_fails_if_pattern_longer_than_number_fragments(
    labelled_binary_contacts_2d_sorted,
):
    """Tests whether subset fails if the pattern is longer than number of fragments"""
    contact_manipulator = contacts.ContactManipulator()
    lab_contacts = contacts.Contacts(labelled_binary_contacts_2d_sorted)
    with pytest.raises(ValueError):
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
