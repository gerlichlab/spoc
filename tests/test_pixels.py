import pytest
import dask.dataframe as dd
import pandas as pd
import pandera as pa
from spoc import contacts, models, pixels


@pytest.fixture
def chromosome_sizes():
    return pd.read_csv(
            "./tests/test_files/hg19.chrom.sizes",
            sep="\t", header=None, names=["chrom", "size"], index_col=["chrom"], squeeze=True
        )

@pytest.fixture
def genomic_binner(chromosome_sizes):
    """genomic binner for pixels"""
    return pixels.GenomicBinner(bin_size=100_000, chrom_sizes=chromosome_sizes, sort_sisters=False)

@pytest.fixture
def good_contacts():
    df = pd.DataFrame(
        {
            "read_name": ["a", "b", "c", "d"],
            "read_length": [1,2,3,4],
            "chrom_1": ["chr1", "chr1", "chr1", "chr1"],
            "start_1": [100_000, 5_000_000, 10_000_000, 10_000_000],
            "end_1": [100_000, 5_000_000, 10_000_000, 10_000_000],
            "mapping_quality_1": [1]*4,
            "align_score_1": [1]*4,
            "align_base_qscore_1": [1]*4,
            "is_labelled_1": [True]*4,
            "sister_identity_1": ["SisterA", "SisterB", "SisterA", "SisterA"],
            "chrom_2":["chr1", "chr1", "chr1", "chr1"],
            "start_2": [500_000, 7_000_000, 25_000_000, 25_000_000],
            "end_2": [500_000, 7_000_000, 25_000_000, 25_000_000],
            "mapping_quality_2": [1]*4,
            "align_score_2": [1]*4,
            "align_base_qscore_2": [1]*4,
            "is_labelled_2": [True]*4,
            "sister_identity_2": ["SisterB", "SisterA", "SisterA", "SisterA"],
            "chrom_3": ["chr2", "chr3", "chr4", "chr4"],
            "start_3": [600_000, 2_000_000, 6_000_000, 6_000_000],
            "end_3": [600_000, 2_000_000, 6_000_000, 6_000_000],
            "mapping_quality_3": [1]*4,
            "align_score_3": [1]*4,
            "align_base_qscore_3": [1]*4,
            "is_labelled_3": [True]*4,
            "sister_identity_3": ["SisterB", "SisterA", "SisterB", "SisterB"],
        }
    )
    return dd.from_pandas(df, chunksize=1000)


def test_genomic_binner_bins_correctly_wo_sister_sorting(genomic_binner, good_contacts):
    result = genomic_binner.bin_contacts(good_contacts)
    import pdb; pdb.set_trace()