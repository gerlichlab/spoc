"""Tests for the pixels module"""
# pylint: disable=redefined-outer-name
import dask.dataframe as dd
import numpy as np
import pandas as pd
import pytest

from spoc import pixels
from spoc.contacts import Contacts


@pytest.fixture
def genomic_binner():
    """genomic binner for pixels"""
    return pixels.GenomicBinner(bin_size=100_000)


@pytest.fixture
def contacts_df():
    """A dataframe containing contacts of order 3"""
    return pd.DataFrame(
        {
            "read_name": ["a", "b", "c", "d"],
            "read_length": [1, 2, 3, 4],
            "chrom_1": ["chr1", "chr1", "chr1", "chr1"],
            "start_1": [100_010, 5_000_010, 10_000_010, 10_000_010],
            "end_1": [100_015, 5_000_050, 10_000_050, 10_000_050],
            "mapping_quality_1": [1] * 4,
            "align_score_1": [1] * 4,
            "align_base_qscore_1": [1] * 4,
            "metadata_1": ["SisterA", "SisterB", "SisterA", "SisterA"],
            "chrom_2": ["chr1", "chr1", "chr1", "chr1"],
            "start_2": [500_010, 7_000_050, 25_000_800, 25_001_000],
            "end_2": [500_050, 7_000_070, 25_000_900, 25_002_000],
            "mapping_quality_2": [1] * 4,
            "align_score_2": [1] * 4,
            "align_base_qscore_2": [1] * 4,
            "metadata_2": ["SisterB", "SisterA", "SisterB", "SisterA"],
            "chrom_3": ["chr1", "chr4", "chr1", "chr1"],
            "start_3": [600_100, 2_000_300, 6_000_050, 6_000_010],
            "end_3": [600_200, 2_000_400, 6_000_600, 6_000_700],
            "mapping_quality_3": [1] * 4,
            "align_score_3": [1] * 4,
            "align_base_qscore_3": [1] * 4,
            "metadata_3": ["SisterB", "SisterA", "SisterA", "SisterB"],
        }
    )


@pytest.fixture
def expected_pixels():
    """A dataframe containing the expected pixels after binning"""
    return pd.DataFrame(
        {
            "chrom": ["chr1"] * 2,
            "start_1": [100_000, 10_000_000],
            "start_2": [500_000, 25_000_000],
            "start_3": [600_000, 6_000_000],
            "contact_count": [1, 2],
        }
    )


@pytest.fixture
def expected_pixels_different_chromosomes():
    """A dataframe containing the expected pixels after binning with different chromosomes
    enabled"""
    return pd.DataFrame(
        {
            "chrom_1": ["chr1"] * 3,
            "start_1": [100_000, 10_000_000, 5_000_000],
            "chrom_2": ["chr1"] * 3,
            "start_2": [500_000, 25_000_000, 7_000_000],
            "chrom_3": ["chr1", "chr1", "chr4"],
            "start_3": [600_000, 6_000_000, 2_000_000],
            "contact_count": [1, 2, 1],
        }
    )


def test_genomic_binner_bins_correctly_same_chromosome_pandas(
    genomic_binner, contacts_df, expected_pixels
):
    """Test if genomic binner bins correctly with same chromosome enabled"""
    contacts = Contacts(contacts_df)
    result = genomic_binner.bin_contacts(contacts)
    assert np.array_equal(result.data.values, expected_pixels.values)
    assert result.number_fragments == 3
    assert result.binsize == 100_000
    assert not result.binary_labels_equal
    assert not result.symmetry_flipped
    assert result.metadata_combi is None


def test_genomic_binner_bins_correctly_w_different_chromosomes_pandas(
    genomic_binner, contacts_df, expected_pixels_different_chromosomes
):
    """Test if genomic binner bins correctly with same chromosome disabled"""
    contacts = Contacts(contacts_df)
    result = genomic_binner.bin_contacts(contacts, same_chromosome=False)
    assert np.array_equal(
        result.data.values, expected_pixels_different_chromosomes.values
    )
    assert result.number_fragments == 3
    assert result.binsize == 100_000
    assert not result.binary_labels_equal
    assert not result.symmetry_flipped
    assert result.metadata_combi is None


def test_genomic_binner_bins_correctly_same_chromosome_dask(
    genomic_binner, contacts_df, expected_pixels
):
    """Test if genomic binner bins correctly with same chromosome enabled
    using a dask dataframe"""
    contacts = Contacts(dd.from_pandas(contacts_df, chunksize=1000))
    result = genomic_binner.bin_contacts(contacts)
    assert np.array_equal(result.data.compute().values, expected_pixels.values)
    assert result.number_fragments == 3
    assert result.binsize == 100_000
    assert not result.binary_labels_equal
    assert not result.symmetry_flipped
    assert result.metadata_combi is None


def test_genomic_binner_bins_correctly_w_different_chromosome_dask(
    genomic_binner, contacts_df, expected_pixels_different_chromosomes
):
    """Test if genomic binner bins correctly with same chromosome disabled
    using a dask dataframe"""
    contacts = Contacts(dd.from_pandas(contacts_df, chunksize=1000))
    result = genomic_binner.bin_contacts(contacts, same_chromosome=False)
    assert np.array_equal(
        result.data.compute().values, expected_pixels_different_chromosomes.values
    )
    assert result.number_fragments == 3
    assert result.binsize == 100_000
    assert not result.binary_labels_equal
    assert not result.symmetry_flipped
    assert result.metadata_combi is None
