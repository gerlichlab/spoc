"""Tests for dealing with symmetry flipping for labelled and unlabelled contacts."""
# pylint: disable=redefined-outer-name
# pylint: disable=unused-import
import dask.dataframe as dd
import numpy as np
import pandas as pd
import pandera as pa
import pytest

from spoc.contacts import ContactManipulator
from spoc.contacts import Contacts


@pytest.mark.parametrize(
    "unflipped, flipped",
    [
        ("unlabelled_contacts_2d", "unlabelled_contacts_2d_flipped"),
        ("unlabelled_contacts_3d", "unlabelled_contacts_3d_flipped"),
    ],
)
def test_unlabelled_contacts_flipped_correctly(unflipped, flipped, request):
    """Test that unlabelled contacts are flipped correctly."""
    unflipped, flipped = request.getfixturevalue(unflipped), request.getfixturevalue(
        flipped
    )
    contacts = Contacts(unflipped)
    flipped_contacts = ContactManipulator().flip_symmetric_contacts(contacts)
    pd.testing.assert_frame_equal(flipped_contacts.data, flipped)


@pytest.mark.parametrize(
    "unflipped, flipped",
    [
        ("unlabelled_contacts_2d", "unlabelled_contacts_2d_flipped"),
        ("unlabelled_contacts_3d", "unlabelled_contacts_3d_flipped"),
    ],
)
def test_unlabelled_contacts_flipped_correctly_dask(unflipped, flipped, request):
    """Test that unlabelled contacts are flipped correctly."""
    unflipped, flipped = dd.from_pandas(
        request.getfixturevalue(unflipped), npartitions=1
    ), request.getfixturevalue(flipped)
    contacts = Contacts(unflipped)
    flipped_contacts = ContactManipulator().flip_symmetric_contacts(contacts)
    pd.testing.assert_frame_equal(
        flipped_contacts.data.compute().reset_index(drop=True),
        flipped.reset_index(drop=True),
    )


@pytest.mark.parametrize(
    "unsorted, sorted_contacts",
    [
        ("labelled_binary_contacts_2d", "labelled_binary_contacts_2d_sorted"),
        ("labelled_binary_contacts_3d", "labelled_binary_contacts_3d_sorted"),
    ],
)
def test_labelled_contacts_are_sorted_correctly(unsorted, sorted_contacts, request):
    """Test that labelled contacts are sorted correctly."""
    unsorted, sorted_contacts = request.getfixturevalue(
        unsorted
    ), request.getfixturevalue(sorted_contacts)
    contacts = Contacts(unsorted)
    result = ContactManipulator().sort_labels(contacts)
    pd.testing.assert_frame_equal(result.data, sorted_contacts)
    assert result.label_sorted


@pytest.mark.parametrize(
    "unsorted, sorted_contacts",
    [
        ("labelled_binary_contacts_2d", "labelled_binary_contacts_2d_sorted"),
        ("labelled_binary_contacts_3d", "labelled_binary_contacts_3d_sorted"),
    ],
)
def test_labelled_contacts_are_sorted_correctly_dask(
    unsorted, sorted_contacts, request
):
    """Test that labelled contacts are sorted correctly with
    underlying dask dataframe."""
    unsorted, sorted_contacts = dd.from_pandas(
        request.getfixturevalue(unsorted), npartitions=1
    ), request.getfixturevalue(sorted_contacts)
    contacts = Contacts(unsorted)
    result = ContactManipulator().sort_labels(contacts)
    pd.testing.assert_frame_equal(
        result.data.compute().reset_index(drop=True),
        sorted_contacts.reset_index(drop=True),
    )
    assert result.label_sorted


@pytest.mark.parametrize(
    "unequated, equated",
    [
        ("binary_contacts_not_equated_2d", "binary_contacts_equated_2d"),
        ("binary_contacts_not_equated_3d", "binary_contacts_equated_3d"),
        ("binary_contacts_not_equated_4d", "binary_contacts_equated_4d"),
    ],
)
def test_equate_binary_labels(unequated, equated, request):
    """Test that binary labels are equated correctly."""
    unequated, equated = request.getfixturevalue(unequated), request.getfixturevalue(
        equated
    )
    contacts = Contacts(unequated, label_sorted=True)
    result = ContactManipulator().equate_binary_labels(contacts)
    pd.testing.assert_frame_equal(result.data, equated)


@pytest.mark.parametrize(
    "unequated, equated",
    [
        ("binary_contacts_not_equated_2d", "binary_contacts_equated_2d"),
        ("binary_contacts_not_equated_3d", "binary_contacts_equated_3d"),
        ("binary_contacts_not_equated_4d", "binary_contacts_equated_4d"),
    ],
)
def test_equate_binary_labels_dask(unequated, equated, request):
    """Test that binary labels are equated correctly with
    underlying dask dataframe."""
    unequated, equated = dd.from_pandas(
        request.getfixturevalue(unequated), npartitions=1
    ), request.getfixturevalue(equated)
    contacts = Contacts(unequated, label_sorted=True)
    result = ContactManipulator().equate_binary_labels(contacts)
    pd.testing.assert_frame_equal(
        result.data.compute().reset_index(drop=True), equated.reset_index(drop=True)
    )


@pytest.mark.parametrize(
    "unflipped, flipped",
    [
        (
            "labelled_binary_contacts_2d_unflipped",
            "labelled_binary_contacts_2d_flipped",
        ),
        (
            "labelled_binary_contacts_3d_unflipped_example2",
            "labelled_binary_contacts_3d_flipped_example2",
        ),
        (
            "labelled_binary_contacts_3d_unflipped",
            "labelled_binary_contacts_3d_flipped",
        ),
    ],
)
def test_flip_labelled_contacts(unflipped, flipped, request):
    """Test that labelled contacts are flipped correctly."""
    unflipped, flipped = request.getfixturevalue(unflipped), request.getfixturevalue(
        flipped
    )
    contacts = Contacts(unflipped, label_sorted=True)
    result = ContactManipulator().flip_symmetric_contacts(contacts)
    pd.testing.assert_frame_equal(
        result.data.reset_index(drop=True), flipped.reset_index(drop=True)
    )


@pytest.mark.parametrize(
    "unflipped, flipped",
    [
        (
            "labelled_binary_contacts_2d_unflipped",
            "labelled_binary_contacts_2d_flipped",
        ),
        (
            "labelled_binary_contacts_3d_unflipped_example2",
            "labelled_binary_contacts_3d_flipped_example2",
        ),
        (
            "labelled_binary_contacts_3d_unflipped",
            "labelled_binary_contacts_3d_flipped",
        ),
    ],
)
def test_flip_labelled_contacts_dask(unflipped, flipped, request):
    """Test that labelled contacts are flipped correctly with
    underlying dask dataframe."""
    unflipped, flipped = dd.from_pandas(
        request.getfixturevalue(unflipped), npartitions=1
    ), request.getfixturevalue(flipped)
    contacts = Contacts(unflipped, label_sorted=True)
    result = ContactManipulator().flip_symmetric_contacts(contacts)
    pd.testing.assert_frame_equal(
        result.data.compute().reset_index(drop=True), flipped.reset_index(drop=True)
    )


@pytest.mark.parametrize(
    "unflipped, flipped",
    [
        (
            "unlabelled_contacts_diff_chrom_3d",
            "unlabelled_contacts_diff_chrom_3d_flipped",
        ),
        (
            "unlabelled_contacts_diff_chrom_2d",
            "unlabelled_contacts_diff_chrom_2d_flipped",
        ),
        (
            "unlabelled_contacts_diff_chrom_4d",
            "unlabelled_contacts_diff_chrom_4d_flipped",
        ),
        (
            "labelled_binary_contacts_diff_chrom_2d",
            "labelled_binary_contacts_diff_chrom_2d_flipped",
        ),
        (
            "labelled_binary_contacts_diff_chrom_3d",
            "labelled_binary_contacts_diff_chrom_3d_flipped",
        ),
    ],
)
def test_flip_unlabelled_contacts_different_chromosomes(unflipped, flipped, request):
    """Test that unlabelled contacts are flipped correctly with different chromosomes."""
    unflipped, flipped = request.getfixturevalue(unflipped), request.getfixturevalue(
        flipped
    )
    contacts = Contacts(unflipped)
    result = ContactManipulator().flip_symmetric_contacts(
        contacts, sort_chromosomes=True
    )
    pd.testing.assert_frame_equal(
        result.data.reset_index(drop=True).sort_index(axis=1),
        flipped.reset_index(drop=True).sort_index(axis=1),
    )
