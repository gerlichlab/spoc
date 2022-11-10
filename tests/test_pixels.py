import pytest
import dask.dataframe as dd
import pandas as pd
import numpy as np
import pandera as pa
from spoc import pixels


@pytest.fixture
def chromosome_sizes():
    return pd.read_csv(
        "./tests/test_files/hg19.chrom.sizes",
        sep="\t",
        header=None,
        names=["chrom", "size"],
        index_col=["chrom"],
        squeeze=True,
    )


@pytest.fixture
def genomic_binner(chromosome_sizes):
    """genomic binner for pixels"""
    return pixels.GenomicBinner(
        bin_size=100_000, chrom_sizes=chromosome_sizes, sort_sisters=False
    )


@pytest.fixture
def genomic_binner_sister_sensitive_w_flipping(chromosome_sizes):
    """genomic binner for pixels"""
    return pixels.GenomicBinner(
        bin_size=100_000,
        chrom_sizes=chromosome_sizes,
        sort_sisters=True,
        same_chromosome=True,
        flip_contacts=True
    )

@pytest.fixture
def genomic_binner_sister_sensitive(chromosome_sizes):
    """genomic binner for pixels"""
    return pixels.GenomicBinner(
        bin_size=100_000,
        chrom_sizes=chromosome_sizes,
        sort_sisters=True,
        same_chromosome=True
    )


@pytest.fixture
def bad_contacts():
    return pd.DataFrame({"be": ["bop"]})


@pytest.fixture
def contacts():
    df = pd.DataFrame(
        {
            "read_name": ["a", "b", "c", "d"],
            "read_length": [1, 2, 3, 4],
            "chrom_1": ["chr1", "chr1", "chr1", "chr1"],
            "start_1": [100_010, 5_000_010, 10_000_010, 10_000_010],
            "end_1": [100_015, 5_000_050, 10_000_050, 10_000_050],
            "mapping_quality_1": [1] * 4,
            "align_score_1": [1] * 4,
            "align_base_qscore_1": [1] * 4,
            "is_labelled_1": [True] * 4,
            "sister_identity_1": ["SisterA", "SisterB", "SisterA", "SisterA"],
            "chrom_2": ["chr1", "chr1", "chr1", "chr1"],
            "start_2": [500_010, 7_000_050, 25_000_800, 25_001_000],
            "end_2": [500_050, 7_000_070, 25_000_900, 25_002_000],
            "mapping_quality_2": [1] * 4,
            "align_score_2": [1] * 4,
            "align_base_qscore_2": [1] * 4,
            "is_labelled_2": [True] * 4,
            "sister_identity_2": ["SisterB", "SisterA", "SisterB", "SisterA"],
            "chrom_3": ["chr1", "chr4", "chr1", "chr1"],
            "start_3": [600_100, 2_000_300, 6_000_050, 6_000_010],
            "end_3": [600_200, 2_000_400, 6_000_600, 6_000_700],
            "mapping_quality_3": [1] * 4,
            "align_score_3": [1] * 4,
            "align_base_qscore_3": [1] * 4,
            "is_labelled_3": [True] * 4,
            "sister_identity_3": ["SisterB", "SisterA", "SisterA", "SisterB"],
        }
    )
    return dd.from_pandas(df, chunksize=1000)

@pytest.fixture
def contacts_w_lower_triangular():
    df = pd.DataFrame(
        {
            "read_name": ["a", "b", "c", "d"]*2,
            "read_length": [1, 2, 3, 4]*2,
            "chrom_1": ["chr1", "chr1", "chr1", "chr1"]*2,
            "start_1": [100_010, 5_000_010, 10_000_010, 10_000_010, 500_010, 7_000_050, 25_000_800, 25_001_000],
            "end_1": [100_015, 5_000_050, 10_000_050, 10_000_050, 500_050, 7_000_070, 25_000_900, 25_002_000],
            "mapping_quality_1": [1] * 8,
            "align_score_1": [1] * 8,
            "align_base_qscore_1": [1] * 8,
            "is_labelled_1": [True] * 8,
            "sister_identity_1": ["SisterA", "SisterB", "SisterA", "SisterA", "SisterB", "SisterA", "SisterB", "SisterA"],
            "chrom_2": ["chr1", "chr1", "chr1", "chr1"]*2,
            "start_2": [500_010, 7_000_050, 25_000_800, 25_001_000, 100_010, 5_000_010, 10_000_010, 10_000_010],
            "end_2": [500_050, 7_000_070, 25_000_900, 25_002_000, 100_015, 5_000_050, 10_000_050, 10_000_050],
            "mapping_quality_2": [1] * 8,
            "align_score_2": [1] * 8,
            "align_base_qscore_2": [1] * 8,
            "is_labelled_2": [True] * 8,
            "sister_identity_2": ["SisterB", "SisterA", "SisterB", "SisterA", "SisterA", "SisterB", "SisterA", "SisterA"],
            "chrom_3": ["chr1", "chr4", "chr1", "chr1"]*2,
            "start_3": [600_100, 2_000_300, 6_000_050, 6_000_010]*2,
            "end_3": [600_200, 2_000_400, 6_000_600, 6_000_700]*2,
            "mapping_quality_3": [1] * 8,
            "align_score_3": [1] * 8,
            "align_base_qscore_3": [1] * 8,
            "is_labelled_3": [True] * 8,
            "sister_identity_3": ["SisterB", "SisterA", "SisterA", "SisterB"]*2,
        }
    )
    return dd.from_pandas(df, chunksize=1000)


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
def expected_pixels_w_sister_sorting_and_flipping():
    return pd.DataFrame(
        {
            "chrom": ["chr1"] * 3,
            "start_1": [500_000, 10_000_000, 10_000_000],
            "start_2": [600_000, 6_000_000, 25_000_000],
            "start_3": [100_000, 25_000_000, 6_000_000],
            "contact_count": [2, 2, 2],
        }
    )


def test_genomic_binner_rejects_bad(genomic_binner, bad_contacts):
    with pytest.raises(pa.errors.SchemaError):
        genomic_binner.bin_contacts(bad_contacts)


def test_genomic_binner_bins_correctly_wo_sister_sorting(
    genomic_binner, contacts, expected_pixels_wo_sister_sorting
):
    result = genomic_binner.bin_contacts(contacts)
    np.array_equal(result.values, expected_pixels_wo_sister_sorting.values)


def test_genomic_binner_sorts_sisters_correctly(
    genomic_binner_sister_sensitive, contacts, expected_pixels_w_sister_sorting
):
    result = genomic_binner_sister_sensitive.bin_contacts(contacts)
    np.array_equal(result.values, expected_pixels_w_sister_sorting.values)

def test_genomic_binner_flips_triplets_correctly(
    genomic_binner_sister_sensitive_w_flipping,  contacts_w_lower_triangular, expected_pixels_w_sister_sorting_and_flipping
):
    result = genomic_binner_sister_sensitive_w_flipping.bin_contacts(contacts_w_lower_triangular)
    np.array_equal(result.values, expected_pixels_w_sister_sorting_and_flipping.values)
