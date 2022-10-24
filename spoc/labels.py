"""This part of spoc is responsible for dealing with label information.
This incudes annotation of labeled fragments, contact-type and sister identity."""

from typing import Dict, Union
import pandas as pd
import numpy as np
from .models import FragmentSchema


class FragmentAnnotator:
    """Responsible for annotating labels and sister identity of mapped read fragments"""

    def __init__(self, label_library: Dict[str, bool]) -> None:
        self._label_library = label_library

    def _is_read_labelled(self, read_code: str) -> Union[bool, None]:
        """If a read is in the label dict, return its labeling state.
        If it is not in there, return None."""
        return self._label_library.get(read_code, None)

    @staticmethod
    def _assign_sister(data_frame) -> pd.Series:
        """assigns sister identity for a given row"""
        return pd.Series(
            np.select(
                [data_frame.strand.astype(bool) != data_frame.is_labelled.astype(bool)],
                ["SisterA"],
                default="SisterB",
            )
        )

    def _assign_label_state(self, data_frame: pd.DataFrame) -> pd.Series:
        """helper method that annotates a fragment data frame"""
        read_codes = data_frame.read_name.str.cat(
            [data_frame.chrom, data_frame.start.astype("str"),
                data_frame.end.astype("str")], sep="_"
        )
        return pd.Series(read_codes.apply(self._is_read_labelled))

    def annotate_fragments(self, fragments: pd.DataFrame) -> pd.DataFrame:
        """Takes fragment dataframe and returns a copy of it with its labelling state in a separate
        column with name `is_labelled`. If drop_uninformative is true, drops fragments that
        are not in label library."""
        # check schema
        FragmentSchema.validate(fragments)
        # assign fragments
        return (
            fragments.assign(is_labelled=self._assign_label_state)
            .dropna(subset=["is_labelled"])
            .assign(sister_identity=self._assign_sister)
        )
