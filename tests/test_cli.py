"""Tests for CLI of spoc"""

from cmath import exp
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
    yield [
        "tests/test_files/good_contacts.triplets.parquet",
        "tests/test_files/good_contacts2.triplets.parquet",
    ]
    # teardown
    shutil.rmtree("tmp")


@pytest.fixture
def good_triplet_file_for_pixels():
    # setup
    os.mkdir("tmp")
    yield "tests/test_files/good_contacts3.triplets.parquet"
    # teardown
    shutil.rmtree("tmp")


@pytest.fixture
def chromosome_sizes():
    return "tests/test_files/hg19.chrom.sizes"


@pytest.fixture
def good_porec_file():
    # setup
    os.mkdir("tmp")
    yield "tests/test_files/good_porec.parquet"
    # teardown
    shutil.rmtree("tmp")


@pytest.fixture
def expected_pixels_w_sister_sorting():
    return pd.DataFrame(
        {
            "chrom": ["chr1"] * 3,
            "start_1": [500_000, 10_000_000, 10_000_000],
            "start_2": [600_000, 6_000_000, 25_000_000],
            "start_3": [100_000, 25_000_000, 6_000_000],
            "contact_count": [1, 1, 1],
        }
    )


@pytest.fixture
def expected_pixels_wo_sister_sorting():
    return pd.DataFrame(
        {
            "chrom_1": ["chr1"] * 3,
            "start_1": [100_000, 5_000_000, 10_000_000],
            "chrom_2": ["chr1"] * 3,
            "start_2": [500_000, 7_000_000, 25_000_000],
            "chrom_3": ["chr1", "chr4", "chr1"],
            "start_3": [600_000, 2_000_000, 6_000_000],
            "contact_count": [1, 1, 2],
        }
    )


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
    runner.invoke(cli.annotate, [good_porec_file, label_library_path, output_path])
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
    result = runner.invoke(
        cli.contacts, [good_triplet_files[0], good_triplet_files[1], f"-o{output_path}"]
    )
    # check content of file
    labelled_fragments = pd.read_parquet(output_path)
    assert len(labelled_fragments) == 8
    first_half = labelled_fragments.iloc[:4, :].reset_index(drop=True)
    second_half = labelled_fragments.iloc[4:, :].reset_index(drop=True)
    assert_frame_equal(first_half, second_half)


def test_bin_contacts_w_sister_sorting(
    good_triplet_file_for_pixels, chromosome_sizes, expected_pixels_w_sister_sorting
):
    """happy path for binning contacts with sister sorting"""
    runner = CliRunner()
    output_path = "tmp/test_output4.parquet"
    result = runner.invoke(
        cli.bin_contacts,
        [good_triplet_file_for_pixels, chromosome_sizes, output_path, "-s", "-c"],
    )
    # check content of file
    pixels = pd.read_parquet(output_path)
    np.array_equal(pixels.values, expected_pixels_w_sister_sorting.values)


def test_bin_contacts_wo_sister_sorting(
    good_triplet_file_for_pixels, chromosome_sizes, expected_pixels_wo_sister_sorting
):
    """happy path for binning contacts without sister sorting"""
    runner = CliRunner()
    output_path = "tmp/test_output5.parquet"
    result = runner.invoke(
        cli.bin_contacts, [good_triplet_file_for_pixels, chromosome_sizes, output_path]
    )
    # check content of file
    pixels = pd.read_parquet(output_path)
    np.array_equal(pixels.values, expected_pixels_wo_sister_sorting.values)
