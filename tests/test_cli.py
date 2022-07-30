"""Tests for CLI of spoc"""

import pytest
import pandas as pd
from pandas.testing import assert_frame_equal
import numpy as np
import shutil
import os
from click.testing import CliRunner

from spoc import cli, models

@pytest.fixture
def good_annotated_porec_file():
    # setup
    os.mkdir("tmp")
    yield "tests/test_files/good_porec.lab.parquet"
    # teardown
    shutil.rmtree("tmp")

@pytest.fixture
def label_library_path():
    return "tests/test_files/ll1.pickle"

@pytest.fixture
def good_triplet_files():
    # setup
    os.mkdir("tmp")
    yield ["tests/test_files/good_contacts.triplets.parquet", "tests/test_files/good_contacts2.triplets.parquet"]
    # teardown
    shutil.rmtree("tmp")

@pytest.fixture
def good_porec_file():
    # setup
    os.mkdir("tmp")
    yield "tests/test_files/good_porec.parquet"
    # teardown
    shutil.rmtree("tmp")


def test_expand_triplets_works(good_annotated_porec_file):
    """Happy path of cli expand"""
    runner = CliRunner()
    output_path = "tmp/test_output1.parquet"
    runner.invoke(cli.expand, [good_annotated_porec_file, output_path])
    # check that file was created
    assert os.path.isfile(output_path)
    # check content of file
    result = pd.read_parquet(output_path)
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

def test_annotate_fragments_works(good_porec_file, label_library_path):
    """tests happy path for annotated fragments"""
    runner = CliRunner()
    output_path = "tmp/test_output2.parquet"
    runner.invoke(cli.annotate, [good_porec_file, label_library_path , output_path])
    # check that file was created
    assert os.path.isfile(output_path)
    # check content of file
    labelled_fragments = pd.read_parquet(output_path)
    assert len(labelled_fragments) == 2
    models.annotated_fragment_schema.validate(labelled_fragments)
    expected = pd.Series([True, False])
    np.array_equal(labelled_fragments.is_labelled.values, expected.values)
    expected = pd.Series(["SisterB", "SisterA"])
    np.array_equal(labelled_fragments.sister_identity.values, expected.values)


def test_merge_contacts_works(good_triplet_files):
    """happy path for merging contacts"""
    runner = CliRunner()
    output_path = "tmp/test_output3.parquet"
    result = runner.invoke(cli.contacts, [good_triplet_files[0], good_triplet_files[1], f"-o{output_path}"])
    # check content of file
    labelled_fragments = pd.read_parquet(output_path)
    assert len(labelled_fragments) == 8
    first_half = labelled_fragments.iloc[:4, :].reset_index(drop=True)
    second_half = labelled_fragments.iloc[4:, :].reset_index(drop=True)
    assert_frame_equal(first_half, second_half)


