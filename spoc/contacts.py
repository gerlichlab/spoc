"""Managing multi-way contacts."""

from __future__ import annotations # needed for self reference in type hints
from typing import List, Union
import pandas as pd
import dask.dataframe as dd
from typing import Union, Optional
from itertools import permutations, product
from spoc.dataframe_models import ContactSchema
import numpy as np

DataFrame = Union[pd.DataFrame, dd.DataFrame]


class Contacts:
    """N-way genomic contacts"""

    def __init__(
        self,
        contact_frame: DataFrame,
        number_fragments: Optional[int] = None,
        metadata_combi: Optional[List[str]] = None,
        label_sorted: bool = False,
    ) -> None:
        self.contains_meta_data = "meta_data_1" in contact_frame.columns # All contacts contain at least one fragment
        if number_fragments is None:
            self.number_fragments = self._guess_number_fragments(contact_frame)
        else:
            self.number_fragments = number_fragments
        self._schema = ContactSchema(
            number_fragments=self.number_fragments, contains_meta_data=self.contains_meta_data
        )
        if isinstance(contact_frame, pd.DataFrame):
            self.is_dask = False
        else:
            self.is_dask = True
        self._data = self._schema.validate(contact_frame)
        self.metadata_combi = metadata_combi
        self.label_sorted = label_sorted

    def _guess_number_fragments(self, contact_frame: DataFrame) -> int:
        """Guesses the number of fragments from the contact frame"""
        return max(int(i.split("_")[1]) for i in contact_frame.columns if "start" in i)


    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, contact_frame):
        self._data = self._schema.validate(contact_frame)


    def subset_on_metadata(self, metadata_combi:List[str]):
        """Filters on metadata and sorts label states"""
        # TODO

    def __repr__(self) -> str:
        return f"<Contacts | order: {self.number_fragments} | contains metadata: {self.contains_meta_data}>"


class ContactManipulator:
    """Responsible for performing operations on
    contact data such as merging, splitting and subsetting."""

    def merge_contacts(self, merge_list: List[Contacts]) -> Contacts:
        """Merge contacts"""
        # validate that merge is possible
        assert (
            len(set([i.number_fragments for i in merge_list])) == 1
        ), "All contacts need to have the same order!"
        assert (
            len(set([i.is_dask for i in merge_list])) == 1
        ), "Mixture of dask and pandas dataframes is not supported!"
        # TODO: assert all have same labelling state

        number_fragments = merge_list[0].number_fragments
        if merge_list[0].is_dask:
            return Contacts(
                dd.concat([i.data for i in merge_list]),
                number_fragments=number_fragments,
            )
        return Contacts(
            pd.concat([i.data for i in merge_list]), number_fragments=number_fragments
        )

    @staticmethod
    def _generate_rename_columns(order):
        columns = ["chrom", "start", "end", "mapping_quality", "align_score", "align_base_qscore", "meta_data"]
        rename_columns = {}
        for i in range(len(order)):
            for column in columns:
                current_name = f"{column}_{i+1}"
                new_name = f"{column}_{order.index(i+1) + 1}"
                rename_columns[current_name] = new_name
        return rename_columns

    def _flip_unlabelled_contacts(self, df: DataFrame) -> DataFrame:
        """Flips contacts"""
        fragment_order = max(int(i.split("_")[1]) for i in df.columns if "start" in i)
        subsets = []
        for perm in permutations(range(1, fragment_order+1)):
            query = "<=".join([f"start_{i}" for i in perm])
            subsets.append(df.query(query).rename(columns=self._generate_rename_columns(perm)))
        # determine which method to use for concatenation
        if isinstance(df, pd.DataFrame):
            result = pd.concat(subsets)
        else:
            result = dd.concat(subsets)
        return result

    def sort_labels(self, contacts:Contacts) -> Contacts:
        """Sorts labels in ascending, alphabetical order"""
        if not contacts.contains_meta_data:
            raise ValueError("Sorting labels for unlabelled contacts is not implemented.")
        # get label values. TODO: this should be a method of contacts
        label_values = set()
        for i in range(contacts.number_fragments):
            label_values.update(contacts.data[f"meta_data_{i+1}"].unique())
        # iterate over all permutations of label values
        subsets = []
        for perm in product(label_values, repeat=contacts.number_fragments):
            query = " and ".join([f"meta_data_{i+1} == '{j}'" for i, j in enumerate(perm)])
            desired_order = [i + 1 for i in np.argsort(perm)]
            subsets.append(contacts.data.query(query).rename(columns=self._generate_rename_columns(desired_order)))
        # determine which method to use for concatenation
        if contacts.is_dask:
            result = dd.concat(subsets).sort_index()
        else:
            result = pd.concat(subsets).sort_index()
        return Contacts(result, number_fragments=contacts.number_fragments, label_sorted=True)


    def subset_on_metadata(self, contacts:Contacts, metadata_combi: List[List[str]]) -> Contacts:
        """Subset contacts based on metadata and sort label states.
        If the metadata combination list contains more than one element, all
        combinations will be assumed to be equivalent and renamed to the first
        element in the list."""
        # TODO


    def flip_symmetric_contacts(self, contacts: Contacts) -> Contacts:
        """Flips contacts based on inherent symmetry"""
        if contacts.contains_meta_data and contacts.metadata_combi is None:
            raise ValueError("""Flipping symmetry is only supported for pure metadata combinations.
                             Either subset or pass to constructor.""")
        if contacts.contains_meta_data:
            raise NotImplementedError("Flipping symmetry for labelled contacts is not yet implemented.")
        else:
            result = self._flip_unlabelled_contacts(contacts.data)
        return Contacts(result, number_fragments=contacts.number_fragments)