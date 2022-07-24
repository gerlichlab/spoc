"""Tests for label functionality"""

import pytest
import pandas as pd
import pandera as pa
import numpy as np

from spoc import labels, models


@pytest.fixture
def empty_annotator():
    """sets up an annotator with empty
    label library"""
    return labels.FragmentAnnotator(dict())


@pytest.fixture
def annotator_with_entries():
    """Annotator with entries"""
    return labels.FragmentAnnotator({"dummy_chr1_1_4": True, "dummy_chr1_2_5": False})


@pytest.fixture
def bad_df():
    return pd.DataFrame({"be": ["bop"]})


@pytest.fixture
def good_df():
    return pd.DataFrame(
        {
            "chrom": ["chr1"] * 3,
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


def test_annotator_rejects_df_w_wrong_structure(empty_annotator, bad_df):
    with pytest.raises(pa.errors.SchemaError):
        empty_annotator.annotate_fragments(bad_df)


def test_annotator_drops_unknown_fragments(annotator_with_entries, good_df):
    labelled_fragments = annotator_with_entries.annotate_fragments(good_df)
    # check length
    assert len(labelled_fragments) == 2


def test_annotator_produces_correct_schema(annotator_with_entries, good_df):
    labelled_fragments = annotator_with_entries.annotate_fragments(good_df)
    # check schema
    models.annotated_fragment_schema.validate(labelled_fragments)


def test_annotator_calls_labels_correctly(annotator_with_entries, good_df):
    labelled_fragments = annotator_with_entries.annotate_fragments(good_df)
    # check values
    expected = pd.Series([True, False])
    np.array_equal(labelled_fragments.is_labelled.values, expected.values)


def test_annotator_calls_sisters_correctly(annotator_with_entries, good_df):
    labelled_fragments = annotator_with_entries.annotate_fragments(good_df)
    # check values
    expected = pd.Series(["SisterB", "SisterA"])
    np.array_equal(labelled_fragments.sister_identity.values, expected.values)
