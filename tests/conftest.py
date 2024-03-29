"""Fixtures for testing symmetry.py"""
# pylint: disable=redefined-outer-name
from __future__ import annotations

import dask.dataframe as dd
import pandas as pd
import pytest


@pytest.fixture
def unlabelled_contacts_2d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 3000],
            "end_1": [200, 3000, 4000],
            "mapping_quality_1": [10, 10, 15],
            "align_score_1": [10, 10, 15],
            "align_base_qscore_1": [10, 10, 15],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 200, 300],
            "end_2": [2000, 300, 400],
            "mapping_quality_2": [10, 10, 10],
            "align_score_2": [10, 10, 10],
            "align_base_qscore_2": [10, 10, 10],
        },
    )


@pytest.fixture
def unlabelled_contacts_2d_dask(unlabelled_contacts_2d):
    return dd.from_pandas(unlabelled_contacts_2d, npartitions=2)


@pytest.fixture
def unlabelled_contacts_3d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 3000],
            "end_1": [200, 3000, 4000],
            "mapping_quality_1": [10, 10, 15],
            "align_score_1": [10, 10, 15],
            "align_base_qscore_1": [10, 10, 15],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 200, 300],
            "end_2": [2000, 300, 400],
            "mapping_quality_2": [10, 10, 10],
            "align_score_2": [10, 10, 10],
            "align_base_qscore_2": [10, 10, 10],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [250, 400, 100],
            "end_3": [300, 500, 200],
            "mapping_quality_3": [10, 10, 5],
            "align_score_3": [10, 10, 5],
            "align_base_qscore_3": [10, 10, 5],
        },
    )


@pytest.fixture
def unlabelled_contacts_3d_dask(unlabelled_contacts_3d):
    return dd.from_pandas(unlabelled_contacts_3d, npartitions=2)


@pytest.fixture
def unlabelled_contacts_2d_flipped():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 200, 300],
            "end_1": [200, 300, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 2000, 3000],
            "end_2": [2000, 3000, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
        },
    )


@pytest.fixture
def unlabelled_contacts_3d_flipped():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 200, 100],
            "end_1": [200, 300, 200],
            "mapping_quality_1": [10, 10, 5],
            "align_score_1": [10, 10, 5],
            "align_base_qscore_1": [10, 10, 5],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [250, 400, 300],
            "end_2": [300, 500, 400],
            "mapping_quality_2": [10, 10, 10],
            "align_score_2": [10, 10, 10],
            "align_base_qscore_2": [10, 10, 10],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [1000, 2000, 3000],
            "end_3": [2000, 3000, 4000],
            "mapping_quality_3": [10, 10, 15],
            "align_score_3": [10, 10, 15],
            "align_base_qscore_3": [10, 10, 15],
        },
    )


@pytest.fixture
def labelled_binary_contacts_2d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 200, 300],
            "end_1": [200, 300, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "metadata_1": ["A", "B", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 2000, 3000],
            "end_2": [2000, 3000, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "metadata_2": ["B", "A", "A"],
        },
    )


@pytest.fixture
def labelled_binary_contacts_2d_sorted():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 300],
            "end_1": [200, 3000, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "metadata_1": ["A", "A", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 200, 3000],
            "end_2": [2000, 300, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "metadata_2": ["B", "B", "A"],
        },
    )


@pytest.fixture
def labelled_binary_contacts_3d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 200, 300],
            "end_1": [200, 300, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "metadata_1": ["A", "B", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 2000, 3000],
            "end_2": [2000, 3000, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "metadata_2": ["B", "A", "A"],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [1000, 2000, 3000],
            "end_3": [2000, 3000, 4000],
            "mapping_quality_3": [10, 10, 15],
            "align_score_3": [10, 10, 15],
            "align_base_qscore_3": [10, 10, 15],
            "metadata_3": ["B", "A", "A"],
        },
    )


@pytest.fixture
def labelled_binary_contacts_3d_sorted():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 300],
            "end_1": [200, 3000, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "metadata_1": ["A", "A", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 2000, 3000],
            "end_2": [2000, 3000, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "metadata_2": ["B", "A", "A"],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [1000, 200, 3000],
            "end_3": [2000, 300, 4000],
            "mapping_quality_3": [10, 10, 15],
            "align_score_3": [10, 10, 15],
            "align_base_qscore_3": [10, 10, 15],
            "metadata_3": ["B", "B", "A"],
        },
    )


@pytest.fixture
def binary_contacts_not_equated_2d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 300],
            "end_1": [200, 3000, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "metadata_1": ["B", "A", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 200, 3000],
            "end_2": [2000, 300, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "metadata_2": ["B", "B", "A"],
        },
    )


@pytest.fixture
def binary_contacts_not_equated_3d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 300],
            "end_1": [200, 3000, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "metadata_1": ["A", "A", "B"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 2000, 3000],
            "end_2": [2000, 3000, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "metadata_2": ["B", "A", "B"],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [1000, 200, 3000],
            "end_3": [2000, 300, 4000],
            "mapping_quality_3": [10, 10, 15],
            "align_score_3": [10, 10, 15],
            "align_base_qscore_3": [10, 10, 15],
            "metadata_3": ["B", "B", "B"],
        },
    )


@pytest.fixture
def binary_contacts_not_equated_4d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 300],
            "end_1": [200, 3000, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "metadata_1": ["A", "A", "B"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 2000, 3000],
            "end_2": [2000, 3000, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "metadata_2": ["B", "A", "B"],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [1000, 200, 3000],
            "end_3": [2000, 300, 4000],
            "mapping_quality_3": [10, 10, 15],
            "align_score_3": [10, 10, 15],
            "align_base_qscore_3": [10, 10, 15],
            "metadata_3": ["B", "B", "B"],
            "chrom_4": ["chr1", "chr1", "chr1"],
            "start_4": [1000, 200, 3000],
            "end_4": [2000, 300, 4000],
            "mapping_quality_4": [10, 10, 15],
            "align_score_4": [10, 10, 15],
            "align_base_qscore_4": [10, 10, 15],
            "metadata_4": ["B", "B", "B"],
        },
    )


@pytest.fixture
def binary_contacts_equated_2d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 300],
            "end_1": [200, 3000, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "metadata_1": ["A", "A", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 200, 3000],
            "end_2": [2000, 300, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "metadata_2": ["A", "B", "A"],
        },
    )


@pytest.fixture
def binary_contacts_equated_3d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 300],
            "end_1": [200, 3000, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "metadata_1": ["A", "A", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 2000, 3000],
            "end_2": [2000, 3000, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "metadata_2": ["A", "A", "A"],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [1000, 200, 3000],
            "end_3": [2000, 300, 4000],
            "mapping_quality_3": [10, 10, 15],
            "align_score_3": [10, 10, 15],
            "align_base_qscore_3": [10, 10, 15],
            "metadata_3": ["B", "B", "A"],
        },
    )


@pytest.fixture
def binary_contacts_equated_4d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 300],
            "end_1": [200, 3000, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "metadata_1": ["A", "A", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 2000, 3000],
            "end_2": [2000, 3000, 4000],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "metadata_2": ["A", "A", "A"],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [1000, 200, 3000],
            "end_3": [2000, 300, 4000],
            "mapping_quality_3": [10, 10, 15],
            "align_score_3": [10, 10, 15],
            "align_base_qscore_3": [10, 10, 15],
            "metadata_3": ["A", "B", "A"],
            "chrom_4": ["chr1", "chr1", "chr1"],
            "start_4": [1000, 200, 3000],
            "end_4": [2000, 300, 4000],
            "mapping_quality_4": [10, 10, 15],
            "align_score_4": [10, 10, 15],
            "align_base_qscore_4": [10, 10, 15],
            "metadata_4": ["B", "B", "A"],
        },
    )


@pytest.fixture
def labelled_binary_contacts_2d_unflipped():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 3000],
            "end_1": [200, 3000, 4000],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "metadata_1": ["A", "A", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 200, 300],
            "end_2": [2000, 300, 400],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "metadata_2": ["B", "B", "A"],
        },
    )


@pytest.fixture
def labelled_binary_contacts_3d_unflipped():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3", "read4"],
            "read_length": [100, 100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1", "chr1"],
            "start_1": [100, 200, 300, 410],
            "end_1": [200, 300, 400, 510],
            "mapping_quality_1": [10, 10, 10, 10],
            "align_score_1": [10, 10, 10, 10],
            "align_base_qscore_1": [10, 10, 10, 10],
            "metadata_1": ["A", "A", "B", "A"],
            "chrom_2": ["chr1", "chr1", "chr1", "chr1"],
            "start_2": [1000, 2000, 3000, 4000],
            "end_2": [2000, 3000, 4000, 5000],
            "mapping_quality_2": [10, 10, 15, 20],
            "align_score_2": [10, 10, 15, 20],
            "align_base_qscore_2": [10, 10, 15, 20],
            "metadata_2": ["B", "A", "B", "A"],
            "chrom_3": ["chr1", "chr1", "chr1", "chr1"],
            "start_3": [100, 200, 310, 400],
            "end_3": [200, 300, 410, 500],
            "mapping_quality_3": [10, 10, 15, 14],
            "align_score_3": [10, 10, 15, 14],
            "align_base_qscore_3": [10, 10, 15, 14],
            "metadata_3": ["B", "B", "B", "A"],
        },
    )


@pytest.fixture
def labelled_binary_contacts_3d_unflipped_example2():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 5000, 800],
            "end_1": [200, 5500, 900],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "metadata_1": ["A", "A", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [200, 2000, 3000],
            "end_2": [300, 3000, 3200],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "metadata_2": ["A", "A", "A"],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [1000, 200, 3000],
            "end_3": [2000, 300, 4000],
            "mapping_quality_3": [10, 10, 14],
            "align_score_3": [10, 10, 14],
            "align_base_qscore_3": [10, 10, 14],
            "metadata_3": ["B", "B", "A"],
        },
    )


@pytest.fixture
def labelled_binary_contacts_2d_flipped():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 300],
            "end_1": [200, 3000, 400],
            "mapping_quality_1": [10, 10, 15],
            "align_score_1": [10, 10, 15],
            "align_base_qscore_1": [10, 10, 15],
            "metadata_1": ["A", "A", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 200, 3000],
            "end_2": [2000, 300, 4000],
            "mapping_quality_2": [10, 10, 10],
            "align_score_2": [10, 10, 10],
            "align_base_qscore_2": [10, 10, 10],
            "metadata_2": ["B", "B", "A"],
        },
    )


@pytest.fixture
def labelled_binary_contacts_3d_flipped():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3", "read4"],
            "read_length": [100, 100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1", "chr1"],
            "start_1": [100, 200, 300, 400],
            "end_1": [200, 300, 400, 500],
            "mapping_quality_1": [10, 10, 10, 14],
            "align_score_1": [10, 10, 10, 14],
            "align_base_qscore_1": [10, 10, 10, 14],
            "metadata_1": ["A", "A", "B", "A"],
            "chrom_2": ["chr1", "chr1", "chr1", "chr1"],
            "start_2": [100, 2000, 310, 410],
            "end_2": [200, 3000, 410, 510],
            "mapping_quality_2": [10, 10, 15, 10],
            "align_score_2": [10, 10, 15, 10],
            "align_base_qscore_2": [10, 10, 15, 10],
            "metadata_2": ["B", "A", "B", "A"],
            "chrom_3": ["chr1", "chr1", "chr1", "chr1"],
            "start_3": [1000, 200, 3000, 4000],
            "end_3": [2000, 300, 4000, 5000],
            "mapping_quality_3": [10, 10, 15, 20],
            "align_score_3": [10, 10, 15, 20],
            "align_base_qscore_3": [10, 10, 15, 20],
            "metadata_3": ["B", "B", "B", "A"],
        },
    )


@pytest.fixture
def labelled_binary_contacts_3d_flipped_example2():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 800],
            "end_1": [200, 3000, 900],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "metadata_1": ["A", "A", "A"],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [200, 5000, 3000],
            "end_2": [300, 5500, 3200],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "metadata_2": ["A", "A", "A"],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [1000, 200, 3000],
            "end_3": [2000, 300, 4000],
            "mapping_quality_3": [10, 10, 14],
            "align_score_3": [10, 10, 14],
            "align_base_qscore_3": [10, 10, 14],
            "metadata_3": ["B", "B", "A"],
        },
    )


@pytest.fixture
def unlabelled_contacts_diff_chrom_2d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr2", "chr3"],
            "start_1": [100, 200, 100],
            "end_1": [200, 300, 200],
            "mapping_quality_1": [10, 10, 15],
            "align_score_1": [10, 10, 15],
            "align_base_qscore_1": [10, 10, 15],
            "chrom_2": ["chr5", "chr1", "chr2"],
            "start_2": [1000, 2000, 300],
            "end_2": [2000, 3000, 400],
            "mapping_quality_2": [10, 10, 10],
            "align_score_2": [10, 10, 10],
            "align_base_qscore_2": [10, 10, 10],
        },
    )


@pytest.fixture
def unlabelled_contacts_diff_chrom_3d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr3"],
            "start_1": [100, 2000, 100],
            "end_1": [200, 3000, 200],
            "mapping_quality_1": [10, 10, 15],
            "align_score_1": [10, 10, 15],
            "align_base_qscore_1": [10, 10, 15],
            "chrom_2": ["chr1", "chr1", "chr2"],
            "start_2": [1000, 200, 300],
            "end_2": [2000, 300, 400],
            "mapping_quality_2": [10, 10, 10],
            "align_score_2": [10, 10, 10],
            "align_base_qscore_2": [10, 10, 10],
            "chrom_3": ["chr1", "chr1", "chr1"],
            "start_3": [250, 400, 500],
            "end_3": [300, 500, 600],
            "mapping_quality_3": [10, 10, 5],
            "align_score_3": [10, 10, 5],
            "align_base_qscore_3": [10, 10, 5],
        },
    )


@pytest.fixture
def unlabelled_contacts_diff_chrom_4d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 100],
            "end_1": [200, 3000, 200],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 200, 300],
            "end_2": [2000, 300, 400],
            "mapping_quality_2": [10, 10, 10],
            "align_score_2": [10, 10, 10],
            "align_base_qscore_2": [10, 10, 10],
            "chrom_3": ["chr4", "chr4", "chr4"],
            "start_3": [3000, 400, 300],
            "end_3": [4000, 500, 400],
            "mapping_quality_3": [10, 10, 10],
            "align_score_3": [10, 10, 10],
            "align_base_qscore_3": [10, 10, 10],
            "chrom_4": ["chr4", "chr4", "chr4"],
            "start_4": [250, 200, 500],
            "end_4": [350, 300, 600],
            "mapping_quality_4": [10, 10, 10],
            "align_score_4": [10, 10, 10],
            "align_base_qscore_4": [10, 10, 10],
        },
    )


@pytest.fixture
def unlabelled_contacts_diff_chrom_4d_flipped():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 200, 100],
            "end_1": [200, 300, 200],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "chrom_2": ["chr1", "chr1", "chr1"],
            "start_2": [1000, 2000, 300],
            "end_2": [2000, 3000, 400],
            "mapping_quality_2": [10, 10, 10],
            "align_score_2": [10, 10, 10],
            "align_base_qscore_2": [10, 10, 10],
            "chrom_3": ["chr4", "chr4", "chr4"],
            "start_3": [250, 200, 300],
            "end_3": [350, 300, 400],
            "mapping_quality_3": [10, 10, 10],
            "align_score_3": [10, 10, 10],
            "align_base_qscore_3": [10, 10, 10],
            "chrom_4": ["chr4", "chr4", "chr4"],
            "start_4": [3000, 400, 500],
            "end_4": [4000, 500, 600],
            "mapping_quality_4": [10, 10, 10],
            "align_score_4": [10, 10, 10],
            "align_base_qscore_4": [10, 10, 10],
        },
    )


@pytest.fixture
def unlabelled_contacts_diff_chrom_3d_flipped():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 200, 500],
            "end_1": [200, 300, 600],
            "mapping_quality_1": [10, 10, 5],
            "align_score_1": [10, 10, 5],
            "align_base_qscore_1": [10, 10, 5],
            "chrom_2": ["chr1", "chr1", "chr2"],
            "start_2": [250, 400, 300],
            "end_2": [300, 500, 400],
            "mapping_quality_2": [10, 10, 10],
            "align_score_2": [10, 10, 10],
            "align_base_qscore_2": [10, 10, 10],
            "chrom_3": ["chr1", "chr1", "chr3"],
            "start_3": [1000, 2000, 100],
            "end_3": [2000, 3000, 200],
            "mapping_quality_3": [10, 10, 15],
            "align_score_3": [10, 10, 15],
            "align_base_qscore_3": [10, 10, 15],
        },
    )


@pytest.fixture
def unlabelled_contacts_diff_chrom_2d_flipped():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr2"],
            "start_1": [100, 2000, 300],
            "end_1": [200, 3000, 400],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "chrom_2": ["chr5", "chr2", "chr3"],
            "start_2": [1000, 200, 100],
            "end_2": [2000, 300, 200],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
        },
    )


@pytest.fixture
def labelled_binary_contacts_diff_chrom_2d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 3000],
            "end_1": [200, 3000, 4000],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "metadata_1": ["A", "A", "A"],
            "chrom_2": ["chr1", "chr1", "chr2"],
            "start_2": [1000, 200, 300],
            "end_2": [2000, 300, 400],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "metadata_2": ["B", "B", "A"],
        },
    )


@pytest.fixture
def labelled_binary_contacts_diff_chrom_2d_flipped():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3"],
            "read_length": [100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1"],
            "start_1": [100, 2000, 3000],
            "end_1": [200, 3000, 4000],
            "mapping_quality_1": [10, 10, 10],
            "align_score_1": [10, 10, 10],
            "align_base_qscore_1": [10, 10, 10],
            "metadata_1": ["A", "A", "A"],
            "chrom_2": ["chr1", "chr1", "chr2"],
            "start_2": [1000, 200, 300],
            "end_2": [2000, 300, 400],
            "mapping_quality_2": [10, 10, 15],
            "align_score_2": [10, 10, 15],
            "align_base_qscore_2": [10, 10, 15],
            "metadata_2": ["B", "B", "A"],
        },
    )


@pytest.fixture
def labelled_binary_contacts_diff_chrom_3d():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3", "read4"],
            "read_length": [100, 100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1", "chr1"],
            "start_1": [100, 200, 300, 410],
            "end_1": [200, 300, 400, 510],
            "mapping_quality_1": [10, 10, 10, 10],
            "align_score_1": [10, 10, 10, 10],
            "align_base_qscore_1": [10, 10, 10, 10],
            "metadata_1": ["A", "A", "B", "A"],
            "chrom_2": ["chr2", "chr1", "chr5", "chr1"],
            "start_2": [1000, 2000, 3000, 4000],
            "end_2": [2000, 3000, 4000, 5000],
            "mapping_quality_2": [10, 10, 15, 20],
            "align_score_2": [10, 10, 15, 20],
            "align_base_qscore_2": [10, 10, 15, 20],
            "metadata_2": ["B", "A", "B", "A"],
            "chrom_3": ["chr3", "chr1", "chr6", "chr1"],
            "start_3": [100, 200, 310, 400],
            "end_3": [200, 300, 410, 500],
            "mapping_quality_3": [10, 10, 15, 14],
            "align_score_3": [10, 10, 15, 14],
            "align_base_qscore_3": [10, 10, 15, 14],
            "metadata_3": ["B", "B", "B", "A"],
        },
    )


@pytest.fixture
def labelled_binary_contacts_diff_chrom_3d_flipped():
    return pd.DataFrame(
        {
            "read_name": ["read1", "read2", "read3", "read4"],
            "read_length": [100, 100, 100, 100],
            "chrom_1": ["chr1", "chr1", "chr1", "chr1"],
            "start_1": [100, 200, 300, 400],
            "end_1": [200, 300, 400, 500],
            "mapping_quality_1": [10, 10, 10, 14],
            "align_score_1": [10, 10, 10, 14],
            "align_base_qscore_1": [10, 10, 10, 14],
            "metadata_1": ["A", "A", "B", "A"],
            "chrom_2": ["chr2", "chr1", "chr5", "chr1"],
            "start_2": [1000, 2000, 3000, 410],
            "end_2": [2000, 3000, 4000, 510],
            "mapping_quality_2": [10, 10, 15, 10],
            "align_score_2": [10, 10, 15, 10],
            "align_base_qscore_2": [10, 10, 15, 10],
            "metadata_2": ["B", "A", "B", "A"],
            "chrom_3": ["chr3", "chr1", "chr6", "chr1"],
            "start_3": [100, 200, 310, 4000],
            "end_3": [200, 300, 410, 5000],
            "mapping_quality_3": [10, 10, 15, 20],
            "align_score_3": [10, 10, 15, 20],
            "align_base_qscore_3": [10, 10, 15, 20],
            "metadata_3": ["B", "B", "B", "A"],
        },
    )
