"""Tests for dealing with symmetry flipping for labelled and unlabelled contacts."""

import pytest
import pandas as pd
import pandera as pa
import numpy as np
import dask.dataframe as dd
from spoc.contacts import Contacts, ContactManipulator

@pytest.fixture
def unlabelled_contacts_2d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 3000],
            "end_1": [200, 3000, 4000],
            "mapping_quality_1": [10, 10, 15],
            "align_score_1": [10, 10, 15],
            "align_base_qscore_1": [10, 10, 15],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 200, 300],
            "end_2": [2000, 300, 400],
            "mapping_quality_2": [10, 10, 10],
            "align_score_2": [10, 10, 10],
            "align_base_qscore_2": [10, 10, 10]
        }
    )

@pytest.fixture
def unlabelled_contacts_2d_dask(unlabelled_contacts_2d):
    return dd.from_pandas(unlabelled_contacts_2d, npartitions=2)

@pytest.fixture
def unlabelled_contacts_3d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 3000],
            "end_1": [200, 3000, 4000],
            "mapping_quality_1": [10, 10, 15],
            "align_score_1": [10, 10, 15],
            "align_base_qscore_1": [10, 10, 15],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 200, 300],
            "end_2": [2000, 300, 400],
            "mapping_quality_2": [10, 10, 10],
            "align_score_2": [10, 10, 10],
            "align_base_qscore_2": [10, 10, 10],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [250, 400, 100],
            "end_3": [300, 500, 200],
            "mapping_quality_3": [10, 10, 5],
            "align_score_3": [10, 10, 5],
            "align_base_qscore_3": [10, 10, 5]
        }
    )

@pytest.fixture
def unlabelled_contacts_3d_dask(unlabelled_contacts_3d):
    return dd.from_pandas(unlabelled_contacts_3d, npartitions=2)

@pytest.fixture
def unlabelled_contacts_2d_flipped():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 200, 300],
            "end_1": [200, 300, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 2000, 3000],
            "end_2": [2000, 3000, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15]
        }
    )


@pytest.fixture
def unlabelled_contacts_3d_flipped():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 200, 100],
            "end_1": [200, 300, 200],
            "mapping_quality_1": [10, 10, 5],
            "align_score_1": [10, 10, 5],
            "align_base_qscore_1": [10, 10, 5],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [250, 400, 300],
            "end_2": [300, 500, 400],
            "mapping_quality_2": [10, 10, 10],
            "align_score_2": [10, 10, 10],
            "align_base_qscore_2": [10, 10, 10],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [1000, 2000, 3000],
            "end_3": [2000, 3000, 4000],
            "mapping_quality_3": [10, 10, 15],
            "align_score_3": [10, 10, 15],
            "align_base_qscore_3": [10, 10, 15]
        }
    )

@pytest.fixture
def labelled_binary_contacts_2d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 200, 300],
            "end_1": [200, 300, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "meta_data_1": ["A", "B", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 2000, 3000],
            "end_2": [2000, 3000, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "meta_data_2": ["B", "A", "A"]
        }
    )

@pytest.fixture
def labelled_binary_contacts_2d_sorted():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 300],
            "end_1": [200, 3000, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "meta_data_1": ["A", "A", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 200, 3000],
            "end_2": [2000, 300, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "meta_data_2": ["B", "B", "A"]
        }
    )


@pytest.fixture
def labelled_binary_contacts_3d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 200, 300],
            "end_1": [200, 300, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "meta_data_1": ["A", "B", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 2000, 3000],
            "end_2": [2000, 3000, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "meta_data_2": ["B", "A", "A"],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [1000, 2000, 3000],
            "end_3": [2000, 3000, 4000],
            "mapping_quality_3": [10, 10, 15],
            "align_score_3": [10, 10, 15],
            "align_base_qscore_3": [10, 10, 15],
            "meta_data_3": ["B", "A", "A"]
        }
    )

@pytest.fixture
def labelled_binary_contacts_3d_sorted():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 300],
            "end_1": [200, 3000, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "meta_data_1": ["A", "A", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 2000, 3000],
            "end_2": [2000, 3000, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "meta_data_2": ["B", "A", "A"],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [1000, 200, 3000],
            "end_3": [2000, 300, 4000],
            "mapping_quality_3": [10, 10, 15],
            "align_score_3": [10, 10, 15],
            "align_base_qscore_3": [10, 10, 15],
            "meta_data_3": ["B", "B", "A"]
        }
    )

@pytest.fixture
def binary_contacts_not_equated_2d():
    return pd.DataFrame(
            {
                "read_name": ["read1", "read2", "read3"],
                "read_length": [100, 100, 100],
                "chrom_1": ["chr1", "chr1", "chr1"],
                "start_1": [100, 2000, 300],
                "end_1": [200, 3000, 400],
                "mapping_quality_1": [10, 10, 10],
                "align_score_1": [10, 10, 10],
                "align_base_qscore_1": [10, 10, 10],
                "meta_data_1": ["B", "A", "A"],
                "chrom_2": ["chr1", "chr1", "chr1"],
                "start_2": [1000, 200, 3000],
                "end_2": [2000, 300, 4000],
                "mapping_quality_2": [10, 10, 15],
                "align_score_2": [10, 10, 15],
                "align_base_qscore_2": [10, 10, 15],
                "meta_data_2": ["B", "B", "A"]
            }
        )

@pytest.fixture
def binary_contacts_not_equated_3d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 300],
            "end_1": [200, 3000, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "meta_data_1": ["A", "A", "B"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 2000, 3000],
            "end_2": [2000, 3000, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "meta_data_2": ["B", "A", "B"],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [1000, 200, 3000],
            "end_3": [2000, 300, 4000],
            "mapping_quality_3": [10, 10, 15],
            "align_score_3": [10, 10, 15],
            "align_base_qscore_3": [10, 10, 15],
            "meta_data_3": ["B", "B", "B"]
        }
    )

@pytest.fixture
def binary_contacts_not_equated_4d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 300],
            "end_1": [200, 3000, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "meta_data_1": ["A", "A", "B"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 2000, 3000],
            "end_2": [2000, 3000, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "meta_data_2": ["B", "A", "B"],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [1000, 200, 3000],
            "end_3": [2000, 300, 4000],
            "mapping_quality_3": [10, 10, 15],
            "align_score_3": [10, 10, 15],
            "align_base_qscore_3": [10, 10, 15],
            "meta_data_3": ["B", "B", "B"],
            "chrom_4": ["chr1", "chr1", "chr1"],
            "start_4": [1000, 200, 3000],
            "end_4": [2000, 300, 4000],
            "mapping_quality_4": [10, 10, 15],
            "align_score_4": [10, 10, 15],
            "align_base_qscore_4": [10, 10, 15],
            "meta_data_4": ["B", "B", "B"]
        }
    )

@pytest.fixture
def binary_contacts_equated_2d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 300],
            "end_1": [200, 3000, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "meta_data_1": ["A", "A", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 200, 3000],
            "end_2": [2000, 300, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "meta_data_2": ["A", "B", "A"]
        }
    )

@pytest.fixture
def binary_contacts_equated_3d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 300],
            "end_1": [200, 3000, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "meta_data_1": ["A", "A", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 2000, 3000],
            "end_2": [2000, 3000, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "meta_data_2": ["A", "A", "A"],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [1000, 200, 3000],
            "end_3": [2000, 300, 4000],
            "mapping_quality_3": [10, 10, 15],
            "align_score_3": [10, 10, 15],
            "align_base_qscore_3": [10, 10, 15],
            "meta_data_3": ["B", "B", "A"]
        }
    )

@pytest.fixture
def binary_contacts_equated_4d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 300],
            "end_1": [200, 3000, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "meta_data_1": ["A", "A", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 2000, 3000],
            "end_2": [2000, 3000, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "meta_data_2": ["A", "A", "A"],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [1000, 200, 3000],
            "end_3": [2000, 300, 4000],
            "mapping_quality_3": [10, 10, 15],
            "align_score_3": [10, 10, 15],
            "align_base_qscore_3": [10, 10, 15],
            "meta_data_3": ["A", "B", "A"],
            "chrom_4": ["chr1", "chr1", "chr1"],
            "start_4": [1000, 200, 3000],
            "end_4": [2000, 300, 4000],
            "mapping_quality_4": [10, 10, 15],
            "align_score_4": [10, 10, 15],
            "align_base_qscore_4": [10, 10, 15],
            "meta_data_4": ["B", "B", "A"]
        }
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


def test_equate_binary_lables_raises_error_for_unsorted_contacts(labelled_binary_contacts_2d):
    contacts = Contacts(labelled_binary_contacts_2d)
    with pytest.raises(AssertionError):
        ContactManipulator().equate_binary_labels(contacts)

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
