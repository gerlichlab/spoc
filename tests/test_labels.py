"""Tests for label functionality"""

import pytest
import pandas as pd
import pandera as pa
import numpy as np

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
def unlabelled_fragments():
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
def labelled_fragments():
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
            "meta_data": ["SisterA"] * 3,
        }
    )


def test_fragment_constructor_rejects_df_w_wrong_structure(bad_df):
    with pytest.raises(pa.errors.SchemaError):
        fragments.Fragments(bad_df)


def test_fragments_constructor_accepts_unlabelled_fragments(unlabelled_fragments):
    frag = fragments.Fragments(unlabelled_fragments)
    assert not frag.contains_meta_data


def test_fragments_constructor_accepts_labelled_fragments(labelled_fragments):
    frag = fragments.Fragments(labelled_fragments)
    assert frag.contains_meta_data


def test_annotator_drops_unknown_fragments(
    annotator_with_entries, unlabelled_fragments
):
    labelled_fragments = annotator_with_entries.annotate_fragments(
        fragments.Fragments(unlabelled_fragments)
    )
    # check length
    assert len(labelled_fragments.data) == 2


def test_annotator_produces_correct_schema(
    annotator_with_entries, unlabelled_fragments
):
    labelled_fragments = annotator_with_entries.annotate_fragments(
        fragments.Fragments(unlabelled_fragments)
    )
    # check schema (fragment constructor checks it)
    assert labelled_fragments.contains_meta_data


def test_annotator_calls_sisters_correctly(
    annotator_with_entries, unlabelled_fragments
):
    labelled_fragments = annotator_with_entries.annotate_fragments(
        fragments.Fragments(unlabelled_fragments)
    )
    # check values
    expected = pd.Series(["SisterB", "SisterA"])
    np.array_equal(labelled_fragments.data.meta_data.values, expected.values)
