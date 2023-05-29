"""Managing multi-way contacts."""

from __future__ import annotations # needed for self reference in type hints
from typing import List, Union
import pandas as pd
import dask.dataframe as dd
from typing import Union, Optional
from spoc.dataframe_models import ContactSchema
from spoc.symmetry import LabelledSymmetryFlipper, UnlabelledSymmetryFlipper


class Contacts:
    """N-way genomic contacts"""

    def __init__(
        self,
        contact_frame: Union[pd.DataFrame, dd.DataFrame],
        number_fragments: Optional[int] = None,
        metadata_combi: Optional[List[str]] = None
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
        self._metadata_combi = metadata_combi

    def _guess_number_fragments(self, contact_frame: Union[pd.DataFrame, dd.DataFrame]) -> int:
        """Guesses the number of fragments from the contact frame"""
        return max(int(i.split("_")[1]) for i in contact_frame.columns if "start" in i)


    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, contact_frame):
        self._data = self._schema.validate(contact_frame)


    def subset_on_metadata(self, metadata_combi:List[str]):
        """Filters on metadata"""
        # TODO

    def flip_symmetric_contacts(self) -> Contacts:
        """Flips contacts based on inherent symmetry"""
        if self.contains_meta_data and self._metadata_combi is None:
            raise ValueError("""Flipping symmetry is only supported for pure metadata combinations.
                             Either subset or pass to constructor.""")
        if self.contains_meta_data:
            flipper = LabelledSymmetryFlipper(self._metadata_combi)
        else:
            flipper = UnlabelledSymmetryFlipper()
        result = flipper.flip_contacts(self.data)
        return Contacts(result, number_fragments=self.number_fragments)

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