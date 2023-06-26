"""This part of spoc is responsible for dealing aligned fragments that
have not yet been converted to contacts. It deals with label information
as well as expanding fragments to contacts."""

from typing import Dict, Union
import pandas as pd
import dask.dataframe as dd
import numpy as np
from itertools import combinations
from .dataframe_models import FragmentSchema, ContactSchema
from .contacts import Contacts


class Fragments:
    """Genomic fragments that can be labelled or not"""

    def __init__(self, fragment_frame: Union[pd.DataFrame, dd.DataFrame]) -> None:
        self._data = FragmentSchema.validate(fragment_frame)
        self._contains_meta_data = True if "meta_data" in fragment_frame.columns else False

    @property
    def data(self):
        return self._data

    @property
    def contains_meta_data(self):
        return self._contains_meta_data


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
            [
                data_frame.chrom,
                data_frame.start.astype("str"),
                data_frame.end.astype("str"),
            ],
            sep="_",
        )
        return pd.Series(read_codes.apply(self._is_read_labelled))

    def annotate_fragments(self, fragments: Fragments) -> Fragments:
        """Takes fragment dataframe and returns a copy of it with its labelling state in a separate
        column with name `is_labelled`. If drop_uninformative is true, drops fragments that
        are not in label library."""
        return Fragments(
            fragments.data.assign(is_labelled=self._assign_label_state)
            .dropna(subset=["is_labelled"])
            .assign(meta_data=self._assign_sister)
            .drop("is_labelled", axis=1)
        )


# TODO: make generic such that label library can hold arbitrary information
class FragmentExpander:
    """Expands n-way fragments over sequencing reads
    to yield contacts."""

    def __init__(self, number_fragments: int) -> None:
        self._number_fragments = number_fragments

    @staticmethod
    def _add_suffix(row, suffix, is_labelled):
        """expands contact fields"""
        output = {}
        for key in ContactSchema.get_contact_fields(is_labelled):
            output[key + f"_{suffix}"] = getattr(row, key)
        return output

    def expand(self, fragments: Fragments) -> Contacts:
        """expand contacts n-ways"""
        # expand fragments
        result = []
        keep_segments = fragments.data.query("pass_filter == True")
        for (read_name, read_df) in keep_segments.groupby("read_name", as_index=False):
            if len(read_df) < self._number_fragments:
                continue

            rows = list(
                read_df.sort_values(["read_start"], ascending=True)
                .assign(pos_on_read=lambda x: np.arange(len(x)))
                .itertuples()
            )

            read_length = rows[0].read_length
            for alignments in combinations(rows, self._number_fragments):
                contact = {"read_name": read_name, "read_length": read_length}
                # add reads
                for index, align in enumerate(alignments, start=1):
                    contact.update(
                        self._add_suffix(align, index, fragments.contains_meta_data)
                    )
                result.append(contact)
        return Contacts(
            pd.DataFrame(result),
            number_fragments=self._number_fragments
        )
