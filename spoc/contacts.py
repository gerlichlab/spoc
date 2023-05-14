"""Managing multi-way contacts."""

from typing import List, Union
import pandas as pd
import dask.dataframe as dd
from typing import Union, Optional
from .dataframe_models import ContactSchema


class Contacts:
    """N-way genomic contacts"""

    def __init__(
        self,
        contact_frame: Optional[Union[pd.DataFrame, dd.DataFrame]] = None,
        number_fragments: int = 3,
        metadata_combi: Optional[List[str]] = None
    ) -> None:
        self.contains_meta_data = "meta_data_1" in contact_frame.columns # All contacts contain at least one fragment
        self._schema = ContactSchema(
            number_fragments=number_fragments, contains_meta_data=self.contains_meta_data
        )
        if isinstance(contact_frame, pd.DataFrame):
            self.is_dask = False
        else:
            self.is_dask = True
        self.number_fragments = number_fragments
        self._data = self._schema.validate(contact_frame)
        self._metadata_combi = metadata_combi

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, contact_frame):
        self._data = self._schema.validate(contact_frame)


    def subset_on_metadata(self, metadata_combi:List[str]):
        """Filters on metadata"""
        # TODO

    def flip_symmetric_contacts(self) -> None:
        """Flips contacts based on inherent symmetry. TODO"""
        if self._metadata_combi is None:
            raise ValueError("""Flipping symmetry is only supported for pure metadata combinations.
                             Either subset or pass to constructor.""")
        raise NotImplementedError

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