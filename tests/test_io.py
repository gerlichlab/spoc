"""This file tests the io module"""
import tempfile
import pytest
from itertools import product
import os
import json
import shutil
from pathlib import Path
import pandas as pd
from spoc.contacts import Contacts
from spoc.io import FileManager
from spoc.file_parameter_models import PixelParameters
import dask.dataframe as dd
from .fixtures.symmetry import (
        unlabelled_contacts_2d
)

CONTACT_PARAMETERS = (
    [2],
    [['A', 'B'], ['B', 'C']],
    [True, False],
    [True, False],
    [True, False],
)

def _create_tmp_dir():
    # check if tmp dir exists
    if not os.path.exists("tmp"):
        os.mkdir("tmp")
    else:
        # if it does, clear it
        shutil.rmtree("tmp")
        os.mkdir("tmp")

@pytest.fixture
def pixels_order_2():
    return pd.DataFrame(
        {
            "chrom": ["chr1"] * 6,
            "start_1": [1, 2, 3, 4, 5, 6],
            "start_2": [1, 2, 3, 4, 5, 6],
            "count": [1, 2, 3, 4, 5, 6],
        }
    )

@pytest.fixture
def pixels_order_3():
    return pd.DataFrame(
        {
            "chrom": ["chr1"] * 6,
            "start_1": [1, 2, 3, 4, 5, 6],
            "start_2": [1, 2, 3, 4, 5, 6],
            "start_3": [1, 2, 3, 4, 5, 6],
            "count": [1, 2, 3, 4, 5, 6],
        }
    )

@pytest.fixture
def example_pixels_w_metadata(pixels_order_2, pixels_order_3):
    # setup
    _create_tmp_dir()
    # create pixels directory
    pixels_dir = "tmp/pixels_test.parquet"
    os.mkdir(pixels_dir)
    expected_parameters = [
        PixelParameters(number_fragments=2, binsize=1000),
        PixelParameters(number_fragments=3, binsize=10_000, metadata_combi=['A', 'B', 'B'],
                        label_sorted=True, binary_labels_equal=True, symmetry_flipped=True),
        PixelParameters(number_fragments=2, binsize=100_000)
    ]
    paths = [
        Path("tmp/pixels_test.parquet/test1.parquet"),
        Path("tmp/pixels_test.parquet/test2.parquet"),
        Path("tmp/pixels_test.parquet/test3.parquet"),
    ]
    dataframes = [
        pixels_order_2,
        pixels_order_3,
        pixels_order_2
    ]
    # create pixels files
    for path, df in zip(paths, dataframes):
        df.to_parquet(path)
    # create metadata json file
    metadata = {
        "test1.parquet": expected_parameters[0].dict(),
        "test2.parquet": expected_parameters[1].dict(),
        "test3.parquet": expected_parameters[2].dict(),
    }
    with open(pixels_dir + '/metadata.json', 'w') as f:
        json.dump(metadata, f)
    yield pixels_dir, expected_parameters, paths, dataframes
    # teardown
    shutil.rmtree("tmp")





@pytest.mark.parametrize('params',
                         product(*CONTACT_PARAMETERS)
                        )
def test_write_read_contacts_global_parameters_w_metadata_pandas(unlabelled_contacts_2d, params):
    """Test writing and reading contacts metadata with pandas"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + '/' + '_'.join([str(x) for x in params]) + '.parquet'
        contacts = Contacts(unlabelled_contacts_2d, *params)
        FileManager().write_multiway_contacts(file_name, contacts)
        # read contacts
        contacts_read = FileManager().load_contacts(file_name)
        # check whether parameters are equal
        assert contacts.get_global_parameters() == contacts_read.get_global_parameters()

@pytest.mark.parametrize('params',
                         product(*CONTACT_PARAMETERS)
                        )
def test_write_read_contacts_global_parameters_w_metadata_dask(unlabelled_contacts_2d, params):
    """Test writing and reading contacts metadata """
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + '/' + '_'.join([str(x) for x in params]) + '.parquet'
        contacts = Contacts(dd.from_pandas(unlabelled_contacts_2d, npartitions=2), *params)
        FileManager().write_multiway_contacts(file_name, contacts)
        # read contacts
        contacts_read = FileManager(use_dask=True).load_contacts(file_name)
        # check whether parameters are equal
        assert contacts.get_global_parameters() == contacts_read.get_global_parameters()

@pytest.mark.parametrize('params',
                         product(*CONTACT_PARAMETERS)
                        )
def test_write_read_contacts_global_parameters_w_metadata_pandas_to_dask(unlabelled_contacts_2d, params):
    """Test writing with pandas and reading with dask"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + '/' + '_'.join([str(x) for x in params]) + '.parquet'
        contacts = Contacts(unlabelled_contacts_2d, *params)
        FileManager().write_multiway_contacts(file_name, contacts)
        # read contacts
        contacts_read = FileManager(use_dask=True).load_contacts(file_name)
        # check whether parameters are equal
        assert contacts.get_global_parameters() == contacts_read.get_global_parameters()

@pytest.mark.parametrize('params',
                         product(*CONTACT_PARAMETERS)
                        )
def test_write_read_contacts_global_parameters_w_metadata_dask_to_pandas(unlabelled_contacts_2d, params):
    """Test writing with pandas and reading with dask"""
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + '/' + '_'.join([str(x) for x in params]) + '.parquet'
        contacts = Contacts(dd.from_pandas(unlabelled_contacts_2d, npartitions=2), *params)
        FileManager().write_multiway_contacts(file_name, contacts)
        # read contacts
        contacts_read = FileManager(use_dask=False).load_contacts(file_name)
        # check whether parameters are equal
        assert contacts.get_global_parameters() == contacts_read.get_global_parameters()


def test_read_pixels_metadata_json(example_pixels_w_metadata):
    """Test reading pixels metadata json file"""
    pixels_dir, expected_parameters, _, _ = example_pixels_w_metadata
    # read metadata
    available_pixels = FileManager().list_pixels(pixels_dir)
    # check whether parameters are equal
    assert len(available_pixels) == len(expected_parameters)
    assert all(actual == expected for actual, expected in zip(available_pixels, expected_parameters))

def test_read_pixels_metadata_json_fails_gracefully():
    """Test reading pixels metadata json file"""
    # read metadata
    with pytest.raises(ValueError) as e:
        FileManager().list_pixels("bad_path")
        assert e.value == "Metadata file not found at bad_path/metadata.json"

def test_read_pixels_as_path(example_pixels_w_metadata):
    """Test reading pixels metadata json file"""
    pixels_dir, expected_parameters, paths, _ = example_pixels_w_metadata
    # read metadata
    for path, expected in zip(paths, expected_parameters):
        pixels = FileManager().load_pixels(pixels_dir, expected, load_dataframe=False)
        assert pixels.path == path
        assert pixels.get_global_parameters() == expected

def test_read_pixels_as_pandas_df(example_pixels_w_metadata):
    """Test reading pixels metadata json file"""
    pixels_dir, expected_parameters, paths, dataframes = example_pixels_w_metadata
    # read metadata
    for path, expected, df in zip(paths, expected_parameters, dataframes):
        pixels = FileManager(use_dask=False).load_pixels(pixels_dir, expected, load_dataframe=True)
        assert pixels.get_global_parameters() == expected
        assert pixels.data.equals(df)

def test_read_pixels_as_dask_df(example_pixels_w_metadata):
    """Test reading pixels metadata json file"""
    pixels_dir, expected_parameters, paths, dataframes = example_pixels_w_metadata
    # read metadata
    for path, expected, df in zip(paths, expected_parameters, dataframes):
        pixels = FileManager(use_dask=True).load_pixels(pixels_dir, expected, load_dataframe=True)
        assert pixels.get_global_parameters() == expected
        assert pixels.data.compute().equals(df)
