"""This file tests the io module for contacts"""
# pylint: disable=redefined-outer-name
import tempfile
import os
import json
import shutil
from pathlib import Path
import pytest
import dask.dataframe as dd
from spoc.contacts import Contacts
from spoc.io import FileManager
from spoc.models.dataframe_models import DataMode
from spoc.models.file_parameter_models import ContactsParameters
from ..fixtures.symmetry import unlabelled_contacts_2d, labelled_binary_contacts_2d


def _create_tmp_dir():
    # check if tmp dir exists
    if not os.path.exists("tmp"):
        os.mkdir("tmp")
    else:
        # if it does, clear it
        shutil.rmtree("tmp")
        os.mkdir("tmp")


@pytest.fixture
def example_contacts_w_metadata(unlabelled_contacts_2d, labelled_binary_contacts_2d):
    # setup
    _create_tmp_dir()
    # create contacts directory
    contacts_dir = "tmp/contacts_test.parquet"
    os.mkdir(contacts_dir)
    expected_parameters = [
        ContactsParameters(number_fragments=2),
        ContactsParameters(number_fragments=2, metadata_combi=["A", "B"]),
        ContactsParameters(
            number_fragments=2, metadata_combi=["A", "B"], label_sorted=True
        ),
        ContactsParameters(
            number_fragments=2,
            metadata_combi=["A", "B"],
            label_sorted=True,
            symmetry_flipped=True,
        ),
    ]
    paths = [
        Path("tmp/contacts_test.parquet/test1.parquet"),
        Path("tmp/contacts_test.parquet/test2.parquet"),
        Path("tmp/contacts_test.parquet/test3.parquet"),
        Path("tmp/contacts_test.parquet/test4.parquet"),
    ]
    dataframes = [
        unlabelled_contacts_2d,
        labelled_binary_contacts_2d,
        labelled_binary_contacts_2d,
        labelled_binary_contacts_2d,
    ]
    # create pixels files
    for path, df in zip(paths, dataframes):
        df.to_parquet(path)
    # create metadata json file
    metadata = {
        "test1.parquet": expected_parameters[0].dict(),
        "test2.parquet": expected_parameters[1].dict(),
        "test3.parquet": expected_parameters[2].dict(),
        "test4.parquet": expected_parameters[3].dict(),
    }
    with open(contacts_dir + "/metadata.json", "w") as f:
        json.dump(metadata, f)
    yield contacts_dir, expected_parameters, paths, dataframes
    # teardown
    shutil.rmtree("tmp")


def test_read_contacts_metadata_json(example_contacts_w_metadata):
    """Test reading pixels metadata json file"""
    contacts_dir, expected_parameters, _, _ = example_contacts_w_metadata
    # read metadata
    available_contacts = FileManager().list_contacts(contacts_dir)
    # check whether parameters are equal
    assert len(available_contacts) == len(expected_parameters)
    assert all(
        actual == expected
        for actual, expected in zip(available_contacts, expected_parameters)
    )


def test_read_contacts_metadata_json_fails_gracefully():
    """Test reading contacts metadata json file"""
    # read metadata
    with pytest.raises(ValueError) as e:
        FileManager().list_contacts("bad_path")
        assert e.value == "Metadata file not found at bad_path/metadata.json"


def test_read_contacts_as_pandas_df(example_contacts_w_metadata):
    """Test reading contacts as pandas dataframe"""
    contacts_dir, expected_parameters, paths, dataframes = example_contacts_w_metadata
    # read metadata
    for path, expected, df in zip(paths, expected_parameters, dataframes):
        contacts = FileManager().load_contacts(contacts_dir, expected)
        assert contacts.get_global_parameters() == expected
        assert contacts.data.equals(df)


def test_read_contacts_as_dask_df(example_contacts_w_metadata):
    """Test reading contacts as pandas dataframe"""
    contacts_dir, expected_parameters, paths, dataframes = example_contacts_w_metadata
    # read metadata
    for path, expected, df in zip(paths, expected_parameters, dataframes):
        contacts = FileManager(DataMode.DASK).load_contacts(contacts_dir, expected)
        assert contacts.get_global_parameters() == expected
        assert contacts.data.compute().equals(df)


@pytest.mark.parametrize(
    "df, params",
    [
        ("unlabelled_contacts_2d", ContactsParameters(number_fragments=2)),
        (
            "labelled_binary_contacts_2d",
            ContactsParameters(
                number_fragments=2, metadata_combi=["A", "B"], symmetry_flipped=True
            ),
        ),
        (
            "labelled_binary_contacts_2d",
            ContactsParameters(
                number_fragments=2, metadata_combi=["A", "B"], label_sorted=True
            ),
        ),
    ],
)
def test_write_pandas_contacts_to_new_file(df, params, request):
    df = request.getfixturevalue(df)
    contacts = Contacts(df, **params.dict())
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + "/" + "test.parquet"
        FileManager().write_contacts(file_name, contacts)
        # check metadata
        metadata = FileManager().list_contacts(file_name)
        assert len(metadata) == 1
        assert metadata[0] == contacts.get_global_parameters()
        # read contacts
        contacts_read = FileManager().load_contacts(file_name, metadata[0])
        # check whether parameters are equal
        assert contacts.get_global_parameters() == contacts_read.get_global_parameters()
        assert contacts.data.equals(contacts_read.data)


@pytest.mark.parametrize(
    "df, params",
    [
        ("unlabelled_contacts_2d", ContactsParameters(number_fragments=2)),
        (
            "labelled_binary_contacts_2d",
            ContactsParameters(
                number_fragments=2, metadata_combi=["A", "B"], symmetry_flipped=True
            ),
        ),
        (
            "labelled_binary_contacts_2d",
            ContactsParameters(
                number_fragments=2, metadata_combi=["A", "B"], label_sorted=True
            ),
        ),
    ],
)
def test_write_dask_contacts_to_new_file(df, params, request):
    df = request.getfixturevalue(df)
    dask_df = dd.from_pandas(df, npartitions=2)
    contacts = Contacts(dask_df, **params.dict())
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + "/" + "test.parquet"
        FileManager().write_contacts(file_name, contacts)
        # check metadata
        metadata = FileManager().list_contacts(file_name)
        assert len(metadata) == 1
        assert metadata[0] == contacts.get_global_parameters()
        # read contacts
        contacts_read = FileManager().load_contacts(file_name, metadata[0])
        # check whether parameters are equal
        assert contacts.get_global_parameters() == contacts_read.get_global_parameters()
        assert contacts.data.compute().equals(contacts_read.data)


@pytest.mark.parametrize(
    "df1,df2,params",
    [
        (
            "unlabelled_contacts_2d",
            "labelled_binary_contacts_2d",
            [
                ContactsParameters(number_fragments=2),
                ContactsParameters(number_fragments=2, metadata_combi=["A", "B"]),
            ],
        ),
        (
            "labelled_binary_contacts_2d",
            "labelled_binary_contacts_2d",
            [
                ContactsParameters(number_fragments=2, metadata_combi=["A", "B"]),
                ContactsParameters(
                    number_fragments=2, metadata_combi=["A", "B"], label_sorted=True
                ),
            ],
        ),
        (
            "labelled_binary_contacts_2d",
            "labelled_binary_contacts_2d",
            [
                ContactsParameters(number_fragments=2, metadata_combi=["B", "A"]),
                ContactsParameters(
                    number_fragments=2, metadata_combi=["A", "B"], symmetry_flipped=True
                ),
            ],
        ),
    ],
)
def test_add_pandas_contacts_to_existing_file(df1, df2, params, request):
    df1, df2 = request.getfixturevalue(df1), request.getfixturevalue(df2)
    params_1, params_2 = params
    contacts1 = Contacts(df1, **params_1.dict())
    contacts2 = Contacts(df2, **params_2.dict())
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + "/" + "test.parquet"
        FileManager().write_contacts(file_name, contacts1)
        FileManager().write_contacts(file_name, contacts2)
        # check metadata
        metadata = FileManager().list_contacts(file_name)
        assert len(metadata) == 2
        # read contacts
        for contacts in [contacts1, contacts2]:
            contacts_read = FileManager().load_contacts(
                file_name, contacts.get_global_parameters()
            )
            # check whether parameters are equal
            assert (
                contacts.get_global_parameters()
                == contacts_read.get_global_parameters()
            )
            assert contacts.data.equals(contacts_read.data)


@pytest.mark.parametrize(
    "df1,df2,params",
    [
        (
            "labelled_binary_contacts_2d",
            "labelled_binary_contacts_2d",
            [
                ContactsParameters(number_fragments=2, metadata_combi=["A", "B"]),
                ContactsParameters(number_fragments=2, metadata_combi=["A", "B"]),
            ],
        ),
        (
            "labelled_binary_contacts_2d",
            "labelled_binary_contacts_2d",
            [
                ContactsParameters(
                    number_fragments=2, metadata_combi=["B", "A"], symmetry_flipped=True
                ),
                ContactsParameters(
                    number_fragments=2, metadata_combi=["B", "A"], symmetry_flipped=True
                ),
            ],
        ),
    ],
)
def test_adding_contacts_to_existing_file_twice_overwrites_contacts(
    df1, df2, params, request
):
    df1, df2 = request.getfixturevalue(df1), request.getfixturevalue(df2)
    params_1, params_2 = params
    contacts1 = Contacts(df1, **params_1.dict())
    contacts2 = Contacts(df2, **params_2.dict())
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + "/" + "test.parquet"
        FileManager().write_contacts(file_name, contacts1)
        FileManager().write_contacts(file_name, contacts2)
        # check metadata
        metadata = FileManager().list_contacts(file_name)
        assert len(metadata) == 1
        # read contacts
        contacts_read = FileManager().load_contacts(
            file_name, contacts1.get_global_parameters()
        )
        # check whether parameters are equal
        assert (
            contacts1.get_global_parameters() == contacts_read.get_global_parameters()
        )
        assert contacts1.data.equals(contacts_read.data)


@pytest.mark.parametrize(
    "df, params",
    [
        ("unlabelled_contacts_2d", ContactsParameters(number_fragments=2)),
    ],
)
def test_load_contacts_from_uri_fails_without_required_parameters(df, params, request):
    """Test loading contacts from uri fails without required parameters"""
    df = request.getfixturevalue(df)
    contacts = Contacts(df, **params.dict())
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + "/" + "test.parquet"
        FileManager().write_contacts(file_name, contacts)
        # try loading without required parameters
        with pytest.raises(ValueError) as e:
            Contacts.from_uri(file_name)


@pytest.mark.parametrize(
    "df, params",
    [
        ("unlabelled_contacts_2d", ContactsParameters(number_fragments=2)),
        (
            "labelled_binary_contacts_2d",
            ContactsParameters(number_fragments=2, metadata_combi=["A", "B"]),
        ),
        (
            "labelled_binary_contacts_2d",
            ContactsParameters(
                number_fragments=2, metadata_combi=["A", "B"], label_sorted=True
            ),
        ),
    ],
)
def test_load_contacts_from_uri_succeeds_exact_match(df, params, request):
    """Test loading contacts from uri succeeds with all required parameters"""
    df = request.getfixturevalue(df)
    contacts = Contacts(df, **params.dict())
    # get meata data parameter
    if params.dict()["metadata_combi"] is None:
        params.metadata_combi = "None"
    else:
        params.metadata_combi = str("".join(params.dict()["metadata_combi"]))
    uri = (
        str(params.dict()["number_fragments"])
        + "::"
        + str(params.dict()["metadata_combi"])
        + "::"
        + str(params.dict()["binary_labels_equal"])
        + "::"
        + str(params.dict()["symmetry_flipped"])
        + "::"
        + str(params.dict()["label_sorted"])
    )
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + "/" + "test.parquet"
        FileManager().write_contacts(file_name, contacts)
        # load contacts
        contacts_read = Contacts.from_uri(file_name + "::" + uri)
        assert contacts.get_global_parameters() == contacts_read.get_global_parameters()


@pytest.mark.parametrize(
    "df, params",
    [
        ("unlabelled_contacts_2d", ContactsParameters(number_fragments=2)),
        (
            "labelled_binary_contacts_2d",
            ContactsParameters(
                number_fragments=2, metadata_combi=["A", "B"], symmetry_flipped=True
            ),
        ),
        (
            "labelled_binary_contacts_2d",
            ContactsParameters(
                number_fragments=2, metadata_combi=["B", "A"], label_sorted=True
            ),
        ),
    ],
)
def test_load_contacts_from_uri_succeeds_partial_match(df, params, request):
    """Test loading contacts from uri succeeds with sufficient required parameters"""
    df = request.getfixturevalue(df)
    contacts = Contacts(df, **params.dict())
    # get meata data parameter
    if params.dict()["metadata_combi"] is None:
        params.metadata_combi = "None"
    else:
        params.metadata_combi = str("".join(params.dict()["metadata_combi"]))
    uri = (
        str(params.dict()["number_fragments"])
        + "::"
        + str(params.dict()["metadata_combi"])
    )
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + "/" + "test.parquet"
        FileManager().write_contacts(file_name, contacts)
        # load contacts
        contacts_read = Contacts.from_uri(file_name + "::" + uri)
        assert contacts.get_global_parameters() == contacts_read.get_global_parameters()


@pytest.mark.parametrize(
    "df, params",
    [
        (
            "labelled_binary_contacts_2d",
            [
                ContactsParameters(
                    number_fragments=2, metadata_combi=["A", "B"], symmetry_flipped=True
                ),
                ContactsParameters(
                    number_fragments=2,
                    metadata_combi=["A", "B"],
                    symmetry_flipped=False,
                ),
            ],
        ),
    ],
)
def test_load_contacts_from_uri_fails_with_ambiguous_specification(df, params, request):
    """Test loading contacts from uri fails with uri is ambiguous"""
    df = request.getfixturevalue(df)
    contacts = Contacts(df, **params[0].dict())
    contacts2 = Contacts(df, **params[1].dict())
    uri = (
        str(params[0].dict()["number_fragments"])
        + "::"
        + str(params[0].dict()["metadata_combi"])
    )
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + "/" + "test.parquet"
        FileManager().write_contacts(file_name, contacts)
        FileManager().write_contacts(file_name, contacts2)
        # load contacts
        with pytest.raises(ValueError) as e:
            Contacts.from_uri(file_name + "::" + uri)
