"""This part of spoc is responsible for dealing with label information. This incudes annotation of labeled fragments, contact-type and sister identity."""

from typing import Dict, Union
import pandas as pd
from .models import fragment_schema


class LabelAnnotator():
    """Responsible for annotating labels of read/contact data structures"""

    def __init__(self, label_library: Dict[str, bool]) -> None:
        self._label_library = label_library
    
    def _is_read_labelled(self, read_code: str) -> Union[bool, None]:
        """If a read is in the label dict, return its labeling state.
        If it is not in there, return None."""
        return self._label_library.get(read_code, pd.NA)

    def _annotate_fragment_frame(self, df: pd.DataFrame) -> pd.Series:
        """helper method that annotates a fragment data frame"""
        read_codes = df.read_name.str.cat([df.chrom, df.start.astype("str"), df.end.astype("str")], sep="_")
        return read_codes.apply(self._is_read_labelled)

    def annotate_fragments(self, fragments: pd.DataFrame, drop_uninformative=True) -> pd.DataFrame:
        """Takes fragment dataframe and returns a copy of it with its labelling state in a separate
        column with name `is_labelled`. If drop_uninformative is true, drops fragments that
        are not in label library."""
        # check schema
        fragment_schema.validate(fragments)
        # assign fragments
        if drop_uninformative:
            return fragments.assign(is_labelled=self._annotate_fragment_frame).dropna(subset=["is_labelled"])
        return fragments.assign(is_labelled=self._annotate_fragment_frame)