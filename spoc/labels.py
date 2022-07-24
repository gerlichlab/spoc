"""This part of spoc is responsible for dealing with label information. This incudes annotation of labeled fragments, contact-type and sister identity."""

from typing import Dict, Union
import pandas as pd
from .models import fragment_schema


class FragmentAnnotator:
    """Responsible for annotating labels and sister identity of mapped read fragments"""

    def __init__(self, label_library: Dict[str, bool]) -> None:
        self._label_library = label_library

    def _is_read_labelled(self, read_code: str) -> Union[bool, None]:
        """If a read is in the label dict, return its labeling state.
        If it is not in there, return None."""
        return self._label_library.get(read_code, None)

    def _assign_sister(self, df_row) -> str:
        """assigns sister identity for a given row"""
        if df_row.strand != df_row.is_labelled:
            return "SisterA"
        return "SisterB"

    def _assign_label_state(self, df: pd.DataFrame) -> pd.Series:
        """helper method that annotates a fragment data frame"""
        read_codes = df.read_name.str.cat(
            [df.chrom, df.start.astype("str"), df.end.astype("str")], sep="_"
        )
        return read_codes.apply(self._is_read_labelled)

    def annotate_fragments(self, fragments: pd.DataFrame) -> pd.DataFrame:
        """Takes fragment dataframe and returns a copy of it with its labelling state in a separate
        column with name `is_labelled`. If drop_uninformative is true, drops fragments that
        are not in label library."""
        # check schema
        fragment_schema.validate(fragments)
        # assign fragments
        fragments_w_label = fragments.assign(
            is_labelled=self._assign_label_state
        ).dropna(subset=["is_labelled"])
        fragments_w_label.loc[:, "is_labelled"] = fragments_w_label.astype(bool)
        fragments_w_label.loc[:, "sister_identity"] = fragments_w_label.apply(
            self._assign_sister, axis=1
        )
        return fragments_w_label


class ContactAnnotator:
    """Responsible for annotating contacts TODO: implement"""
