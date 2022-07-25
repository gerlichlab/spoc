from array import array
import pytest
import pandas as pd
import pandera as pa
import numpy as np
import dask.dataframe as dd

from spoc import contacts, models


@pytest.fixture
def triplet_expander():
    """expander for triplets"""
    return contacts.ContactExpander(number_fragments=3)


@pytest.fixture
def triplet_manipulator():
    """manipulator for triplest"""
    return contacts.ContactManipulator(number_fragments=3)


@pytest.fixture
def triplet_manipulator_dask():
    """manipulator for triplest"""
    return contacts.ContactManipulator(number_fragments=3, use_dask=True)


@pytest.fixture
def bad_df():
    return pd.DataFrame({"be": ["bop"]})


@pytest.fixture
def good_df():
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
            "is_labelled": [False, True, False, True, False, True],
            "sister_identity": [
                "SisterA",
                "SisterB",
                "SisterA",
                "SisterB",
                "SisterA",
                "SisterB",
            ],
        }
    )


def test_expander_rejects_df_w_wrong_structure(triplet_expander, bad_df):
    with pytest.raises(pa.errors.SchemaError):
        triplet_expander.expand(bad_df)


def test_expander_drops_reads_w_too_little_fragments(triplet_expander, good_df):
    result = triplet_expander.expand(good_df)
    assert len(set(result.read_name)) == 1
    assert result.read_name[0] == "dummy"


def test_expander_returns_correct_schema(triplet_expander, good_df):
    result = triplet_expander.expand(good_df)
    models.HigherOrderContactSchema(number_fragments=3).validate(result)


def test_expander_returns_correct_number_of_contacts(triplet_expander, good_df):
    result = triplet_expander.expand(good_df)
    assert len(result) == 4


def test_expander_returns_correct_contacts(triplet_expander, good_df):
    result = triplet_expander.expand(good_df)
    assert np.array_equal(result["start_1"].values, np.array([1, 1, 1, 2]))
    assert np.array_equal(result["end_1"].values, np.array([4, 4, 4, 5]))
    assert np.array_equal(result["start_2"].values, np.array([2, 2, 3, 3]))
    assert np.array_equal(result["end_2"].values, np.array([5, 5, 6, 6]))
    assert np.array_equal(result["start_3"].values, np.array([3, 4, 4, 4]))
    assert np.array_equal(result["end_3"].values, np.array([6, 7, 7, 7]))
    assert np.array_equal(
        result["sister_identity_1"].values,
        np.array(["SisterA", "SisterA", "SisterA", "SisterB"]),
    )
    assert np.array_equal(
        result["sister_identity_2"].values,
        np.array(["SisterB", "SisterB", "SisterA", "SisterA"]),
    )
    assert np.array_equal(
        result["sister_identity_3"].values,
        np.array(["SisterA", "SisterB", "SisterB", "SisterB"]),
    )


def test_merge_rejects_wrong_pandas_df(triplet_manipulator, bad_df):
    with pytest.raises(pa.errors.SchemaError):
        triplet_manipulator.merge_contacts([bad_df])


def test_merge_rejects_wrong_dask_df(triplet_manipulator_dask, bad_df):
    with pytest.raises(pa.errors.SchemaError):
        result = triplet_manipulator_dask.merge_contacts([bad_df])
        result.compute()


def test_merge_works_for_good_pandas_df(triplet_expander, triplet_manipulator, good_df):
    contacts = triplet_expander.expand(good_df)
    result = triplet_manipulator.merge_contacts([contacts, contacts])
    assert result.shape[0] == 8
    assert result.shape[1] == contacts.shape[1]


def test_merge_works_for_good_dask_df(
    triplet_expander, triplet_manipulator_dask, good_df
):
    contacts = triplet_expander.expand(good_df)
    contacts_dask = dd.from_pandas(contacts, npartitions=1)
    result = triplet_manipulator_dask.merge_contacts(
        [contacts_dask, contacts_dask]
    ).compute()
    assert result.shape[0] == 8
    assert result.shape[1] == contacts.shape[1]
