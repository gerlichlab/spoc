"""This file tests the io module for pixels"""
# pylint: disable=redefined-outer-name
import tempfile
import os
import json
import shutil
from pathlib import Path
import pytest
import pandas as pd
import dask.dataframe as dd
from spoc.io import FileManager
from spoc.models.file_parameter_models import PixelParameters
from spoc.pixels import Pixels


CONTACT_PARAMETERS = (
    [2],
    [("A", "B"), ("B", "C")],
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
def df_order_2():
    return pd.DataFrame(
        {
            "chrom": ["chr1"] * 6,
            "start_1": [1, 2, 3, 4, 5, 6],
            "start_2": [1, 2, 3, 4, 5, 6],
            "count": [1, 2, 3, 4, 5, 6],
        }
    )


@pytest.fixture
def df_order_3():
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
def example_pixels_w_metadata(df_order_2, df_order_3):
    # setup
    _create_tmp_dir()
    # create pixels directory
    pixels_dir = "tmp/pixels_test.parquet"
    os.mkdir(pixels_dir)
    expected_parameters = [
        PixelParameters(number_fragments=2, binsize=1000),
        PixelParameters(
            number_fragments=3,
            binsize=10_000,
            metadata_combi=("A", "B", "B"),
            label_sorted=True,
            binary_labels_equal=True,
            symmetry_flipped=True,
        ),
        PixelParameters(number_fragments=2, binsize=100_000),
    ]
    paths = [
        Path("tmp/pixels_test.parquet/test1.parquet"),
        Path("tmp/pixels_test.parquet/test2.parquet"),
        Path("tmp/pixels_test.parquet/test3.parquet"),
    ]
    dataframes = [df_order_2, df_order_3, df_order_2]
    # create pixels files
    for path, df in zip(paths, dataframes):
        df.to_parquet(path)
    # create metadata json file
    metadata = {
        "test1.parquet": expected_parameters[0].dict(),
        "test2.parquet": expected_parameters[1].dict(),
        "test3.parquet": expected_parameters[2].dict(),
    }
    with open(pixels_dir + "/metadata.json", "w") as f:
        json.dump(metadata, f)
    yield pixels_dir, expected_parameters, paths, dataframes
    # teardown
    shutil.rmtree("tmp")


def test_read_pixels_metadata_json(example_pixels_w_metadata):
    """Test reading pixels metadata json file"""
    pixels_dir, expected_parameters, _, _ = example_pixels_w_metadata
    # read metadata
    available_pixels = FileManager().list_pixels(pixels_dir)
    # check whether parameters are equal
    assert len(available_pixels) == len(expected_parameters)
    assert all(
        actual == expected
        for actual, expected in zip(available_pixels, expected_parameters)
    )


def test_read_pixels_metadata_json_fails_gracefully():
    """Test reading pixels metadata json file"""
    # read metadata
    with pytest.raises(ValueError) as e:
        FileManager().list_pixels("bad_path")
        assert e.value == "Metadata file not found at bad_path/metadata.json"


def test_read_pixels_as_pandas_df(example_pixels_w_metadata):
    """Test reading pixels metadata json file"""
    pixels_dir, expected_parameters, paths, dataframes = example_pixels_w_metadata
    # read metadata
    for path, expected, df in zip(paths, expected_parameters, dataframes):
        pixels = FileManager(use_dask=False).load_pixels(
            pixels_dir, expected
        )
        assert pixels.get_global_parameters() == expected
        assert pixels.data.equals(df)


def test_read_pixels_as_dask_df(example_pixels_w_metadata):
    """Test reading pixels metadata json file"""
    pixels_dir, expected_parameters, paths, dataframes = example_pixels_w_metadata
    # read metadata
    for path, expected, df in zip(paths, expected_parameters, dataframes):
        pixels = FileManager(use_dask=True).load_pixels(
            pixels_dir, expected
        )
        assert pixels.get_global_parameters() == expected
        assert pixels.data.compute().equals(df)


@pytest.mark.parametrize(
    "df, params",
    [
        ("df_order_2", PixelParameters(number_fragments=2, binsize=1000)),
        (
            "df_order_3",
            PixelParameters(
                number_fragments=3, binsize=10_000, metadata_combi=("A", "B", "B")
            ),
        ),
        ("df_order_2", PixelParameters(number_fragments=2, binsize=100_000)),
    ],
)
def test_write_pandas_pixels_to_new_file(df, params, request):
    df = request.getfixturevalue(df)
    pixels = Pixels(df, **params.dict())
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + "/" + "test.parquet"
        FileManager().write_pixels(file_name, pixels)
        # check metadata
        metadata = FileManager().list_pixels(file_name)
        assert len(metadata) == 1
        assert metadata[0] == pixels.get_global_parameters()
        # read pixels
        pixels_read = FileManager().load_pixels(file_name, metadata[0])
        # check whether parameters are equal
        assert pixels.get_global_parameters() == pixels_read.get_global_parameters()
        assert pixels.data.equals(pixels_read.data)


@pytest.mark.parametrize(
    "df, params",
    [
        ("df_order_2", PixelParameters(number_fragments=2, binsize=1000)),
        (
            "df_order_3",
            PixelParameters(
                number_fragments=3, binsize=10_000, metadata_combi=("A", "B", "B")
            ),
        ),
        ("df_order_2", PixelParameters(number_fragments=2, binsize=100_000)),
    ],
)
def test_write_dask_pixels_to_new_file(df, params, request):
    df = request.getfixturevalue(df)
    dask_df = dd.from_pandas(df, npartitions=2)
    pixels = Pixels(dask_df, **params.dict())
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + "/" + "test.parquet"
        FileManager().write_pixels(file_name, pixels)
        # check metadata
        metadata = FileManager().list_pixels(file_name)
        assert len(metadata) == 1
        # read pixels
        pixels_read = FileManager().load_pixels(file_name, metadata[0])
        # check whether parameters are equal
        assert pixels.get_global_parameters() == pixels_read.get_global_parameters()
        assert pixels.data.compute().equals(pixels_read.data)


@pytest.mark.parametrize(
    "df1,df2,params",
    [
        (
            "df_order_2",
            "df_order_3",
            [
                PixelParameters(number_fragments=2, binsize=1000),
                PixelParameters(number_fragments=3, binsize=10_000),
            ],
        ),
        (
            "df_order_3",
            "df_order_2",
            [
                PixelParameters(
                    number_fragments=3, binsize=10_000, metadata_combi=("A", "B", "B")
                ),
                PixelParameters(number_fragments=2, binsize=100),
            ],
        ),
        (
            "df_order_2",
            "df_order_3",
            [
                PixelParameters(number_fragments=2, binsize=100_000),
                PixelParameters(
                    number_fragments=3, binsize=10_000, metadata_combi=("A", "B", "B")
                ),
            ],
        ),
    ],
)
def test_add_pandas_pixels_to_existing_file(df1, df2, params, request):
    df1, df2 = request.getfixturevalue(df1), request.getfixturevalue(df2)
    params_1, params_2 = params
    pixels1 = Pixels(df1, **params_1.dict())
    pixels2 = Pixels(df2, **params_2.dict())
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + "/" + "test.parquet"
        FileManager().write_pixels(file_name, pixels1)
        FileManager().write_pixels(file_name, pixels2)
        # check metadata
        metadata = FileManager().list_pixels(file_name)
        assert len(metadata) == 2
        # read pixels
        for pixels in [pixels1, pixels2]:
            pixels_read = FileManager().load_pixels(
                file_name, pixels.get_global_parameters()
            )
            # check whether parameters are equal
            assert pixels.get_global_parameters() == pixels_read.get_global_parameters()
            assert pixels.data.equals(pixels_read.data)


@pytest.mark.parametrize(
    "df, params",
    [
        ("df_order_2", PixelParameters(number_fragments=2, binsize=1000)),
    ],
)
def test_load_pixels_from_uri_fails_without_required_parameters(df, params, request):
    """Test loading pixels from uri fails without required parameters"""
    df = request.getfixturevalue(df)
    pixels = Pixels(df, **params.dict())
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + "/" + "test.parquet"
        FileManager().write_pixels(file_name, pixels)
        # try loading without required parameters
        with pytest.raises(ValueError) as e:
            Pixels.from_uri(file_name)


@pytest.mark.parametrize(
    "df, params",
    [
        ("df_order_2", PixelParameters(number_fragments=2, binsize=1000)),
        (
            "df_order_2",
            PixelParameters(
                number_fragments=2, binsize=1000, metadata_combi=("A", "B")
            ),
        ),
        (
            "df_order_2",
            PixelParameters(
                number_fragments=2,
                binsize=1000,
                metadata_combi=("A", "B"),
                label_sorted=True,
            ),
        ),
    ],
)
def test_load_pixels_from_uri_succeeds_exact_match(df, params, request):
    """Test loading pixels from uri succeeds with all required parameters"""
    df = request.getfixturevalue(df)
    pixels = Pixels(df, **params.dict())
    # get meata data parameter
    if params.dict()["metadata_combi"] is None:
        params.metadata_combi = "None"
    else:
        params.metadata_combi = str("".join(params.dict()["metadata_combi"]))
    uri = (
        str(params.dict()["number_fragments"])
        + "::"
        + str(params.dict()["binsize"])
        + "::"
        + str(params.dict()["metadata_combi"])
        + "::"
        + str(params.dict()["binary_labels_equal"])
        + "::"
        + str(params.dict()["symmetry_flipped"])
        + "::"
        + str(params.dict()["label_sorted"])
        + "::"
        + str(params.dict()["same_chromosome"])
    )
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + "/" + "test.parquet"
        FileManager().write_pixels(file_name, pixels)
        # load pixels
        pixels_read = Pixels.from_uri(file_name + "::" + uri)
        assert pixels.get_global_parameters() == pixels_read.get_global_parameters()


@pytest.mark.parametrize(
    "df, params",
    [
        ("df_order_2", PixelParameters(number_fragments=2, binsize=1000)),
        (
            "df_order_2",
            PixelParameters(
                number_fragments=2, binsize=1000, metadata_combi=("A", "B")
            ),
        ),
        (
            "df_order_2",
            PixelParameters(
                number_fragments=2,
                binsize=1000,
                metadata_combi=["A", "B"],
                label_sorted=True,
            ),
        ),
    ],
)
def test_load_pixels_from_uri_succeeds_partial_match(df, params, request):
    """Test loading pixels from uri succeeds with sufficient required parameters"""
    df = request.getfixturevalue(df)
    pixels = Pixels(df, **params.dict())
    # get meata data parameter
    if params.dict()["metadata_combi"] is None:
        params.metadata_combi = "None"
    else:
        params.metadata_combi = str("".join(params.dict()["metadata_combi"]))
    uri = (
        str(params.dict()["number_fragments"])
        + "::"
        + str(params.dict()["binsize"])
        + "::"
        + str(params.dict()["metadata_combi"])
    )
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + "/" + "test.parquet"
        FileManager().write_pixels(file_name, pixels)
        # load pixels
        pixels_read = Pixels.from_uri(file_name + "::" + uri)
        assert pixels.get_global_parameters() == pixels_read.get_global_parameters()


@pytest.mark.parametrize(
    "df, params",
    [
        (
            "df_order_2",
            [
                PixelParameters(
                    number_fragments=2, binsize=1000, metadata_combi=("A", "B")
                ),
                PixelParameters(
                    number_fragments=2,
                    binsize=1000,
                    metadata_combi=["A", "B"],
                    label_sorted=True,
                ),
            ],
        )
    ],
)
def test_load_pixels_from_uri_fails_with_ambiguous_specification(df, params, request):
    """Test loading pixels from uri fails with uri is ambiguous"""
    df = request.getfixturevalue(df)
    pixels = Pixels(df, **params[0].dict())
    pixels2 = Pixels(df, **params[1].dict())
    uri = (
        str(params[0].dict()["number_fragments"])
        + "::"
        + str(params[0].dict()["binsize"])
    )
    with tempfile.TemporaryDirectory() as tmpdirname:
        file_name = tmpdirname + "/" + "test.parquet"
        FileManager().write_pixels(file_name, pixels)
        FileManager().write_pixels(file_name, pixels2)
        # load pixels
        with pytest.raises(ValueError) as e:
            Pixels.from_uri(file_name + "::" + uri)
