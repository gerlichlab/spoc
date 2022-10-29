"""Managing multi-way contacts."""

from itertools import combinations
from typing import List, Union
import pandas as pd
import dask.dataframe as dd
import numpy as np
from .models import AnnotatedFragmentSchema, HigherOrderContactSchema


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

    def expand(self, fragments: pd.DataFrame) -> pd.DataFrame:
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
        return pd.DataFrame(result)


class ContactManipulator:
    """Responsible for performing operations on
    contact data such as merging, splitting and subsetting."""

    def __init__(self, number_fragments: int, use_dask: bool = False) -> None:
        self._number_fragments = number_fragments
        self._schema = HigherOrderContactSchema(number_fragments=number_fragments)
        self._use_dask = use_dask

    def _merge_contacts_pandas(self, merge_list: List[pd.DataFrame]) -> pd.DataFrame:
        # validate schema
        for data_frame in merge_list:
            self._schema.validate(data_frame)
        return pd.concat(merge_list)

    def _merge_contacts_dask(self, merge_list: List[dd.DataFrame]) -> dd.DataFrame:
        # validate schema
        merge_with_check = []
        for data_frame in merge_list:
            merge_with_check.append(
                self._schema.validate(data_frame)
            )  # needed to add check to task graph
        return dd.concat(merge_with_check)

    def merge_contacts(
        self, merge_list: List[Union[pd.DataFrame, dd.DataFrame]]
    ) -> Union[pd.DataFrame, dd.DataFrame]:
        """Merge contacts"""
        if not self._use_dask:
            return self._merge_contacts_pandas(merge_list)
        return self._merge_contacts_dask(merge_list)
