"""Dataframe models"""

import pandera as pa

fragment_schema = pa.DataFrameSchema({
    "chrom": pa.Column(str, checks=[
        pa.Check.str_startswith("chr")
    ]),
    "start": pa.Column(int, checks=pa.Check.gt(0)),
    "end": pa.Column(int, checks=pa.Check.gt(0)),
    "strand": pa.Column(bool),
    "read_name": pa.Column(str),
    "mapping_quality": pa.Column(int),
    "align_score": pa.Column(int),
    "align_base_qscore": pa.Column(int)
})

labelled_fragment_schema = pa.DataFrameSchema({
    "chrom": pa.Column(str, checks=[
        pa.Check.str_startswith("chr")
    ]),
    "start": pa.Column(int, checks=pa.Check.gt(0)),
    "end": pa.Column(int, checks=pa.Check.gt(0)),
    "strand": pa.Column(bool),
    "read_name": pa.Column(str),
    "mapping_quality": pa.Column(int),
    "align_score": pa.Column(int),
    "align_base_qscore": pa.Column(int),
    "is_labelled": pa.Column(bool),
    "sister_identity": pa.Column(str, checks=[pa.Check(lambda x: x.isin(["SisterA", "SisterB"]))])
})