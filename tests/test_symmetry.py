import pytest
import pandas as pd
import pandera as pa
import numpy as np
import dask.dataframe as dd
from spoc.contacts import Contacts

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
            "align_base_qscore_2": [10, 10, 10]
        }
    )

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
            "align_base_qscore_3": [10, 10, 5]
        }
    )

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
            "align_base_qscore_2": [10, 10, 15]
        }
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
            "align_base_qscore_3": [10, 10, 15]
        }
    )

@pytest.mark.parametrize("unflipped, flipped",
                            [('unlabelled_contacts_2d', 'unlabelled_contacts_2d_flipped'),
                             ('unlabelled_contacts_3d', 'unlabelled_contacts_3d_flipped')])
def test_unlabelled_contacts_flipped_correctly(unflipped, flipped, request):
    unflipped, flipped = request.getfixturevalue(unflipped), request.getfixturevalue(flipped)
    contacts = Contacts(unflipped)
    flipped_contacts = contacts.flip_symmetric_contacts()
    pd.testing.assert_frame_equal(flipped_contacts.data, flipped)

@pytest.mark.parametrize("unflipped, flipped",
                            [('unlabelled_contacts_2d', 'unlabelled_contacts_2d_flipped'),
                             ('unlabelled_contacts_3d', 'unlabelled_contacts_3d_flipped')])
def test_unlabelled_contacts_flipped_correctly_dask(unflipped, flipped, request):
    unflipped, flipped = dd.from_pandas(request.getfixturevalue(unflipped), npartitions=1), request.getfixturevalue(flipped)
    contacts = Contacts(unflipped)
    flipped_contacts = contacts.flip_symmetric_contacts()
    pd.testing.assert_frame_equal(flipped_contacts.data.compute(), flipped)