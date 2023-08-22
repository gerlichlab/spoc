"""Tests for label functionality"""

import pytest
import pandas as pd
import pandera as pa
import numpy as np
import dask.dataframe as dd

from spoc import fragments


@pytest.fixture
def empty_annotator():
    """sets up an annotator with empty
    label library"""
    return fragments.FragmentAnnotator(dict())


@pytest.fixture
def annotator_with_entries():
    """Annotator with entries"""
    return fragments.FragmentAnnotator(
        {"dummy_chr1_1_4": True, "dummy_chr1_2_5": False}
    )


@pytest.fixture
def bad_df():
    return pd.DataFrame({"be": ["bop"]})

@pytest.fixture
def unlabelled_df():
    return pd.DataFrame(
        {
            "chrom": pd.Series(["chr1"] * 3, dtype="category"),
            "start": [1, 2, 3],
            "end": [4, 5, 6],
            "strand": [True] * 3,
            "read_name": ["dummy"] * 3,
            "read_start": [1, 2, 3],
            "read_end": [4, 5, 6],
            "read_length": [1] * 3,
            "mapping_quality": [1, 2, 3],
            "align_score": [1, 2, 3],
            "align_base_qscore": [1, 2, 3],
            "pass_filter": [True] * 3,
        }
    )

@pytest.fixture
def unlabelled_fragments(unlabelled_df):
    return fragments.Fragments(unlabelled_df)


@pytest.fixture
def unlabelled_fragments_dask(unlabelled_df):
    return fragments.Fragments(dd.from_pandas(unlabelled_df, npartitions=1))

@pytest.fixture
def labelled_df():
    return pd.DataFrame(
        {
            "chrom": pd.Series(["chr1"] * 3, dtype="category"),
            "start": [1, 2, 3],
            "end": [4, 5, 6],
            "strand": [True] * 3,
            "read_name": ["dummy"] * 3,
            "read_start": [1, 2, 3],
            "read_end": [4, 5, 6],
            "read_length": [1] * 3,
            "mapping_quality": [1, 2, 3],
            "align_score": [1, 2, 3],
            "align_base_qscore": [1, 2, 3],
            "pass_filter": [True] * 3,
            "metadata": ["SisterA"] * 3,
        }
    )


def test_fragment_constructor_rejects_df_w_wrong_structure(bad_df):
    with pytest.raises(pa.errors.SchemaError):
        fragments.Fragments(bad_df)


def test_fragments_constructor_accepts_unlabelled_fragments(unlabelled_df):
    frag = fragments.Fragments(unlabelled_df)
    assert not frag.contains_metadata


def test_fragments_constructor_accepts_labelled_fragments(labelled_df):
    frag = fragments.Fragments(labelled_df)
    assert frag.contains_metadata

@pytest.mark.parametrize("fragments", ["unlabelled_fragments", "unlabelled_fragments_dask"])
def test_annotator_drops_unknown_fragments(
    annotator_with_entries, fragments, request
):
    labelled_fragments = annotator_with_entries.annotate_fragments(
        request.getfixturevalue(fragments)
    )
    # check length
    assert len(labelled_fragments.data) == 2

@pytest.mark.parametrize("fragments", ["unlabelled_fragments", "unlabelled_fragments_dask"])
def test_annotator_produces_correct_schema(
    annotator_with_entries, fragments, request
):
    labelled_fragments = annotator_with_entries.annotate_fragments(
        request.getfixturevalue(fragments)
    )
    # check schema (fragment constructor checks it)
    assert labelled_fragments.contains_metadata

@pytest.mark.parametrize("fragments", ["unlabelled_fragments", "unlabelled_fragments_dask"])
def test_annotator_calls_sisters_correctly(
    annotator_with_entries, fragments, request
):
    labelled_fragments = annotator_with_entries.annotate_fragments(
        request.getfixturevalue(fragments)
    )
    # check values
    expected = pd.Series(["SisterB", "SisterA"])
    if isinstance(request.getfixturevalue(fragments), dd.DataFrame):
        assert np.array_equal(labelled_fragments.data.metadata.values.compute(), expected.values)
    else:
        assert np.array_equal(labelled_fragments.data.metadata.values, expected)

@pytest.mark.parametrize("fragments", ["unlabelled_fragments", "unlabelled_fragments_dask"])
def test_annotator_maintains_dataframe_type(
    annotator_with_entries, fragments, request
):
    labelled_fragments = annotator_with_entries.annotate_fragments(
        request.getfixturevalue(fragments)
    )
    # check values
    assert isinstance(labelled_fragments.data, type(request.getfixturevalue(fragments).data))