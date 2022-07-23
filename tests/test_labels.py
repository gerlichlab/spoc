"""Tests for label functionality"""

import pytest
import pandas as pd
from pandas.testing import assert_series_equal
import pandera as pa
import numpy as np

from spoc import labels, models


@pytest.fixture
def empty_annotator():
    """sets up an annotator with empty
    label library"""
    return labels.LabelAnnotator(dict())

@pytest.fixture
def annotator_with_entries():
    """Annotator with entries"""
    return labels.LabelAnnotator(
        {
            "dummy_chr1_1_4": True,
            "dummy_chr1_2_5": False
        }
    )

@pytest.fixture
def bad_df():
    return pd.DataFrame({
        "be": ["bop"]
    })

@pytest.fixture
def good_df():
    return pd.DataFrame(
        {
            "chrom": ["chr1"]*3,
            "start": [1, 2, 3],
            "end": [4, 5, 6],
            "strand": [True]*3,
            "read_name": ["dummy"]*3,
            "mapping_quality": [1,2,3],
            "align_score": [1,2,3],
            "align_base_qscore": [1,2,3]
        }
    )


def test_annotator_rejects_df_w_wrong_structure(empty_annotator, bad_df):
    with pytest.raises(pa.errors.SchemaError):
        empty_annotator.annotate_fragments(bad_df)


def test_annotator_gives_none_for_unknown_fragments(empty_annotator, good_df):
    labelled_fragments = empty_annotator.annotate_fragments(good_df, drop_uninformative=False)
    assert len(labelled_fragments) == 3
    assert np.all(pd.isna(labelled_fragments.is_labelled))

def test_annotator_filters_unknown_fragments(empty_annotator, good_df):
    labelled_fragments = empty_annotator.annotate_fragments(good_df, drop_uninformative=True)
    assert len(labelled_fragments) == 0

def test_annotator_annotates_fragments_correctly(annotator_with_entries, good_df):
    labelled_fragments = annotator_with_entries.annotate_fragments(good_df)
    # check length
    assert len(labelled_fragments) == 2
    # check schema
    models.labelled_fragment_schema.validate(labelled_fragments)
    # check values
    expected = pd.Series([True, False])
    np.array_equal(labelled_fragments.is_labelled.values, expected.values)