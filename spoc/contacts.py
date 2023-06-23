"""Managing multi-way contacts."""

from __future__ import annotations # needed for self reference in type hints
from typing import List, Union
import pandas as pd
import dask.dataframe as dd
from typing import Union, Optional, Dict
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
        binary_labels_equal: bool = False,
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
        self.binary_labels_equal = binary_labels_equal

    def _guess_number_fragments(self, contact_frame: DataFrame) -> int:
        """Guesses the number of fragments from the contact frame"""
        return max(int(i.split("_")[1]) for i in contact_frame.columns if "start" in i)

    def get_label_values(self) -> List[str]:
        """Returns all label values"""
        if not self.contains_meta_data:
            raise ValueError("Contacts do not contain metadata!")
        output = set()
        for i in range(self.number_fragments):
            if self.is_dask:
                output.update(self.data[f"meta_data_{i+1}"].unique().compute())
            else:
                output.update(self.data[f"meta_data_{i+1}"].unique())
        return output

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
            result = pd.concat(subsets).sort_index()
        else:
            result = dd.concat(subsets).reset_index()\
                                        .sort_values("index")\
                                        .set_index("index")
        return result

    def sort_labels(self, contacts:Contacts) -> Contacts:
        """Sorts labels in ascending, alphabetical order"""
        if not contacts.contains_meta_data:
            raise ValueError("Sorting labels for unlabelled contacts is not implemented.")
        # get label values. TODO: this should be a method of contacts
        label_values = contacts.get_label_values()
        # iterate over all permutations of label values
        subsets = []
        for perm in product(label_values, repeat=contacts.number_fragments):
            query = " and ".join([f"meta_data_{i+1} == '{j}'" for i, j in enumerate(perm)])
            desired_order = [i + 1 for i in np.argsort(perm)]
            subsets.append(contacts.data.query(query).rename(columns=self._generate_rename_columns(desired_order)))
        # determine which method to use for concatenation
        if contacts.is_dask:
            # this is a bit of a hack to get the index sorted. Dask does not support index sorting
            result = dd.concat(subsets).reset_index()\
                                        .sort_values("index")\
                                        .set_index("index")
        else:
            result = pd.concat(subsets).sort_index()
        return Contacts(result, number_fragments=contacts.number_fragments, label_sorted=True)


    def _generate_binary_label_mapping(self, label_values:List[str], number_fragments: int) -> Dict[str, str]:
        sorted_labels = sorted(label_values)
        mapping = {}
        for i in range(number_fragments + 1):
            target = [sorted_labels[0]]*(number_fragments - i) + [sorted_labels[-1]]*(i)
            source = [sorted_labels[0]]*(i) + [sorted_labels[-1]]*(number_fragments - i)
            if i <= (number_fragments // 2):
                mapping[tuple(source)] = tuple(target)
            else:
                mapping[tuple(source)] = ()
        return mapping

    def equate_binary_labels(self, contacts:Contacts) -> Contacts:
        """Binary labels often only carry information about whether
        they happen between the same or different fragments. This
        method equates these labels be replacing all equivalent binary labels with
        the alphabetically first label.
        For example, if we have a contact between two fragments
        that are labelled A and B, the label is either AB or BA. For most
        applications, there is no difference between these two contacts and this
        method would replace both labels with AB.
        """
        assert contacts.contains_meta_data, "Contacts do not contain metadata!"
        assert contacts.label_sorted, "Contacts are not label sorted!"
        # get label values
        label_values = contacts.get_label_values()
        assert len(label_values) == 2, "Equate binary labels only works for binary labels!"
        # generate mapping diectionary
        mapping = self._generate_binary_label_mapping(label_values, contacts.number_fragments)
        subsets = []
        for source, target in mapping.items():
            query = " and ".join([f"meta_data_{i+1} == '{j}'" for i, j in enumerate(source)])
            subset = contacts.data.query(query)
            # assign target labels to dataframe
            for i, j in enumerate(target):
                subset[f"meta_data_{i+1}"] = j
            subsets.append(subset)
        # determine which method to use for concatenation
        if contacts.is_dask:
            # this is a bit of a hack to get the index sorted. Dask does not support index sorting
            result = dd.concat(subsets).reset_index()\
                                        .sort_values("index")\
                                        .set_index("index")
        else:
            result = pd.concat(subsets).sort_index()
        return Contacts(result, number_fragments=contacts.number_fragments, label_sorted=True,
                        binary_labels_equal=True)


    def subset_on_metadata(self, contacts:Contacts, metadata_combi: List[str]) -> Contacts:
        """Subset contacts based on metadata and sort label states."""
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