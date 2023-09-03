"""This part of spoc is responsible for dealing aligned fragments that
have not yet been converted to contacts. It deals with label information
as well as expanding fragments to contacts."""

from typing import Dict, Union
import pandas as pd
import dask.dataframe as dd
import numpy as np
from itertools import combinations
from .dataframe_models import FragmentSchema, ContactSchema, DataFrame
from .contacts import Contacts


class Fragments:
    """Genomic fragments that can be labelled or not"""

    def __init__(self, fragment_frame: DataFrame) -> None:
        self._data = FragmentSchema.validate(fragment_frame)
        self._contains_metadata = True if "metadata" in fragment_frame.columns else False

    @property
    def data(self):
        return self._data

    @property
    def contains_metadata(self):
        return self._contains_metadata
    
    @property
    def is_dask(self):
        return isinstance(self._data, dd.DataFrame)


# TODO: make generic such that label library can hold arbitrary information
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
            .assign(metadata=self._assign_sister)
            .drop("is_labelled", axis=1)
        )



class FragmentExpander:
    """Expands n-way fragments over sequencing reads
    to yield contacts."""

    def __init__(self, number_fragments: int, contains_metadata: bool = True) -> None:
        self._number_fragments = number_fragments
        self._contains_metadata = contains_metadata
        self._schema = ContactSchema(number_fragments, contains_metadata)

    @staticmethod
    def _add_suffix(row, suffix:int, contains_metadata:bool) -> Dict:
        """expands contact fields"""
        output = {}
        for key in ContactSchema.get_contact_fields(contains_metadata):
            output[key + f"_{suffix}"] = getattr(row, key)
        return output

    def _get_expansion_output_structure(self) -> pd.DataFrame:
        """returns expansion output dataframe structure for dask"""
        return pd.DataFrame(columns=list(self._schema._schema.columns.keys()) + ['level_2']).set_index(["read_name", 'read_length', 'level_2'])

    def _expand_single_read(self, read_df: pd.DataFrame, contains_metadata:bool) -> pd.DataFrame:
        """Expands a single read"""
        if len(read_df) < self._number_fragments:
            return pd.DataFrame()

        rows = list(
            read_df.sort_values(["read_start"], ascending=True)
            .assign(pos_on_read=lambda x: np.arange(len(x)))
            .itertuples()
        )

        result = []
        for alignments in combinations(rows, self._number_fragments):
            contact = {}
            # add reads
            for index, align in enumerate(alignments, start=1):
                contact.update(
                    self._add_suffix(align, index, contains_metadata)
                )
            result.append(contact)
        return pd.DataFrame(result)

    def expand(self, fragments: Fragments) -> Contacts:
        """expand contacts n-ways"""
        # construct dataframe type specific kwargs
        if fragments.is_dask:
            kwargs = dict(meta=self._get_expansion_output_structure())
        else:
            kwargs = dict()
        # expand
        contact_df = fragments.data\
                            .groupby(["read_name", "read_length"])\
                            .apply(self._expand_single_read, 
                                    contains_metadata=fragments.contains_metadata,
                                    **kwargs)\
                            .reset_index()\
                            .drop("level_2", axis=1)
        #return contact_df
        return Contacts(
            contact_df,
            number_fragments=self._number_fragments
        )
