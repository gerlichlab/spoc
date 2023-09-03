"""This file tests the io module"""
import tempfile
import pytest
from itertools import product
from spoc.contacts import Contacts
from spoc.io import FileManager
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