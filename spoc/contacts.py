"""Managing multi-way contacts."""

from itertools import combinations
from typing import List, Union
import pandas as pd
import dask.dataframe as dd
import numpy as np
from typing import Union, Optional
from .dataframe_models import AnnotatedFragmentSchema, HigherOrderContactSchema


class Contacts:
    """N-way genomic contacts"""

    def __init__(self, contact_frame: Optional[Union[pd.DataFrame, dd.DataFrame]] = None, number_fragments: int = 3) -> None:
        self._schema = HigherOrderContactSchema(number_fragments=3)
        if isinstance(contact_frame, pd.DataFrame):
            self.is_dask = False
        else:
            self.is_dask = True
        self.number_fragments = number_fragments
        self._data = self._schema.validate(contact_frame)
    
    @property
    def data(self):
        return self._data
    
    @data.setter
    def data(self, contact_frame):
        self._data = self._schema.validate(contact_frame)


class ContactExpander:
    """Expands n-way contacts over sequencing reads"""

    # TODO: Implement upper triangular flipping of Triplets

    def __init__(self, number_fragments: int) -> None:
        self._number_fragments = number_fragments

    @staticmethod
    def _add_suffix(row, suffix):
        """expands contact fields"""
        output = {}
        for key in HigherOrderContactSchema.contact_fields:
            output[key + f"_{suffix}"] = getattr(row, key)
        return output

    def expand(self, fragments: pd.DataFrame) -> Contacts:
        """expand contacts n-ways"""
        # check input
        AnnotatedFragmentSchema.validate(fragments)
        # expand fragments
        result = []
        keep_segments = fragments.query("pass_filter == True")
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
                    contact.update(self._add_suffix(align, index))
                result.append(contact)
        return Contacts(pd.DataFrame(result), number_fragments=self._number_fragments)


class ContactManipulator:
    """Responsible for performing operations on
    contact data such as merging, splitting and subsetting."""

    def merge_contacts(
        self, merge_list: List[Contacts]
    ) -> Contacts:
        """Merge contacts"""
        # validate that merge is possible
        assert len(set([i.number_fragments for i in merge_list])) == 1, "All contacts need to have the same order!"
        assert len(set([i.is_dask for i in merge_list])) == 1, "Mixture of dask and pandas dataframes is not supported!"
        
        number_fragments = merge_list[0].number_fragments
        if merge_list[0].is_dask:
            return Contacts(dd.concat([i.data for i in merge_list]), number_fragments=number_fragments)
        return  Contacts(pd.concat([i.data for i in merge_list]), number_fragments=number_fragments)