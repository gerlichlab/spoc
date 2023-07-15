"""Tests for dealing with symmetry flipping for labelled and unlabelled contacts."""

import pytest
import pandas as pd
import pandera as pa
import numpy as np
import dask.dataframe as dd
from spoc.contacts import Contacts, ContactManipulator
from .fixtures.symmetry import (
            unlabelled_contacts_2d, 
            unlabelled_contacts_2d_flipped,
            unlabelled_contacts_3d,
            unlabelled_contacts_3d_flipped,
            labelled_binary_contacts_2d,
            labelled_binary_contacts_2d_sorted,
            labelled_binary_contacts_3d,
            labelled_binary_contacts_3d_sorted,
            binary_contacts_not_equated_2d,
            binary_contacts_not_equated_3d,
            binary_contacts_not_equated_4d,
            binary_contacts_equated_2d,
            binary_contacts_equated_3d,
            binary_contacts_equated_4d,
            labelled_binary_contacts_2d_unflipped,
            labelled_binary_contacts_2d_flipped,
            labelled_binary_contacts_3d_unflipped,
            labelled_binary_contacts_3d_unflipped_example2,
            labelled_binary_contacts_3d_flipped,
            labelled_binary_contacts_3d_flipped_example2
)



@pytest.mark.parametrize("unflipped, flipped",
                            [('unlabelled_contacts_2d', 'unlabelled_contacts_2d_flipped'),
                             ('unlabelled_contacts_3d', 'unlabelled_contacts_3d_flipped')])
def test_unlabelled_contacts_flipped_correctly(unflipped, flipped, request):
    unflipped, flipped = request.getfixturevalue(unflipped), request.getfixturevalue(flipped)
    contacts = Contacts(unflipped)
    flipped_contacts = ContactManipulator().flip_symmetric_contacts(contacts)
    pd.testing.assert_frame_equal(flipped_contacts.data, flipped)

@pytest.mark.parametrize("unflipped, flipped",
                            [('unlabelled_contacts_2d', 'unlabelled_contacts_2d_flipped'),
                             ('unlabelled_contacts_3d', 'unlabelled_contacts_3d_flipped')])
def test_unlabelled_contacts_flipped_correctly_dask(unflipped, flipped, request):
    unflipped, flipped = dd.from_pandas(request.getfixturevalue(unflipped), npartitions=1), request.getfixturevalue(flipped)
    contacts = Contacts(unflipped)
    flipped_contacts = ContactManipulator().flip_symmetric_contacts(contacts)
    pd.testing.assert_frame_equal(flipped_contacts.data.compute().reset_index(drop=True), flipped.reset_index(drop=True))


@pytest.mark.parametrize("unsorted, sorted_contacts",
                            [('labelled_binary_contacts_2d', 'labelled_binary_contacts_2d_sorted'),
                             ('labelled_binary_contacts_3d', 'labelled_binary_contacts_3d_sorted')])
def test_labelled_contacts_are_sorted_correctly(unsorted, sorted_contacts, request):
    unsorted, sorted_contacts = request.getfixturevalue(unsorted), request.getfixturevalue(sorted_contacts)
    contacts = Contacts(unsorted)
    result = ContactManipulator().sort_labels(contacts)
    pd.testing.assert_frame_equal(result.data, sorted_contacts)
    assert result.label_sorted

@pytest.mark.parametrize("unsorted, sorted_contacts",
                            [('labelled_binary_contacts_2d', 'labelled_binary_contacts_2d_sorted'),
                             ('labelled_binary_contacts_3d', 'labelled_binary_contacts_3d_sorted')])
def test_labelled_contacts_are_sorted_correctly_dask(unsorted, sorted_contacts, request):
    unsorted, sorted_contacts = dd.from_pandas(request.getfixturevalue(unsorted), npartitions=1), request.getfixturevalue(sorted_contacts)
    contacts = Contacts(unsorted)
    result = ContactManipulator().sort_labels(contacts)
    pd.testing.assert_frame_equal(result.data.compute().reset_index(drop=True), sorted_contacts.reset_index(drop=True))
    assert result.label_sorted


@pytest.mark.parametrize("unequated, equated",
                            [('binary_contacts_not_equated_2d', 'binary_contacts_equated_2d'),
                             ('binary_contacts_not_equated_3d', 'binary_contacts_equated_3d'),
                             ('binary_contacts_not_equated_4d', 'binary_contacts_equated_4d')])
def test_equate_binary_labels(unequated, equated, request):
    unequated, equated = request.getfixturevalue(unequated), request.getfixturevalue(equated)
    contacts = Contacts(unequated, label_sorted=True)
    result = ContactManipulator().equate_binary_labels(contacts)
    pd.testing.assert_frame_equal(result.data, equated)

@pytest.mark.parametrize("unequated, equated",
                            [('binary_contacts_not_equated_2d', 'binary_contacts_equated_2d'),
                             ('binary_contacts_not_equated_3d', 'binary_contacts_equated_3d'),
                             ('binary_contacts_not_equated_4d', 'binary_contacts_equated_4d')])
def test_equate_binary_labels_dask(unequated, equated, request):
    unequated, equated = dd.from_pandas(request.getfixturevalue(unequated), npartitions=1), request.getfixturevalue(equated)
    contacts = Contacts(unequated, label_sorted=True)
    result = ContactManipulator().equate_binary_labels(contacts)
    pd.testing.assert_frame_equal(result.data.compute().reset_index(drop=True), equated.reset_index(drop=True))


@pytest.mark.parametrize("unflipped, flipped",
                            [('labelled_binary_contacts_2d_unflipped', 'labelled_binary_contacts_2d_flipped'),
                             ('labelled_binary_contacts_3d_unflipped_example2', 'labelled_binary_contacts_3d_flipped_example2'),
                             ('labelled_binary_contacts_3d_unflipped', 'labelled_binary_contacts_3d_flipped')])
def test_flip_labelled_contacts(unflipped, flipped, request):
    unflipped, flipped = request.getfixturevalue(unflipped), request.getfixturevalue(flipped)
    contacts = Contacts(unflipped, label_sorted=True)
    result = ContactManipulator().flip_symmetric_contacts(contacts)
    pd.testing.assert_frame_equal(result.data.reset_index(drop=True), flipped.reset_index(drop=True))


@pytest.mark.parametrize("unflipped, flipped",
                            [('labelled_binary_contacts_2d_unflipped', 'labelled_binary_contacts_2d_flipped'),
                             ('labelled_binary_contacts_3d_unflipped_example2', 'labelled_binary_contacts_3d_flipped_example2'),
                             ('labelled_binary_contacts_3d_unflipped', 'labelled_binary_contacts_3d_flipped')])
def test_flip_labelled_contacts_dask(unflipped, flipped, request):
    unflipped, flipped = dd.from_pandas(request.getfixturevalue(unflipped), npartitions=1), request.getfixturevalue(flipped)
    contacts = Contacts(unflipped, label_sorted=True)
    result = ContactManipulator().flip_symmetric_contacts(contacts)
    pd.testing.assert_frame_equal(result.data.compute().reset_index(drop=True), flipped.reset_index(drop=True))
