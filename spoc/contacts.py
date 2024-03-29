"""Managing multi-way contacts."""
from __future__ import annotations  # needed for self reference in type hints

from itertools import permutations
from itertools import product
from typing import Dict
from typing import List
from typing import Optional
from typing import Tuple

import dask.dataframe as dd
import duckdb
import numpy as np
import pandas as pd

from spoc.models.dataframe_models import ContactSchema
from spoc.models.dataframe_models import DataFrame
from spoc.models.dataframe_models import DataMode
from spoc.models.dataframe_models import GenomicDataSchema
from spoc.models.file_parameter_models import ContactsParameters


class Contacts:
    """N-way genomic contacts

    Args:
        contact_frame (DataFrame): DataFrame containing the contact data.
        number_fragments (int, optional): Number of fragments. Defaults to None.
        metadata_combi (List[str], optional): List of metadata combinations. Defaults to None.
        label_sorted (bool, optional): Whether the labels are sorted. Defaults to False.
        binary_labels_equal (bool, optional): Whether the binary labels are equal. Defaults to False.
        symmetry_flipped (bool, optional): Whether the symmetry is flipped. Defaults to False.

    Attributes:
        contains_metadata (bool): Whether the contact data contains metadata.
        number_fragments (int): Number of fragments.
        is_dask (bool): Whether the contact data is a Dask DataFrame.
        metadata_combi (List[str]): List of metadata combinations.
        label_sorted (bool): Whether the labels are sorted.
        binary_labels_equal (bool): Whether the binary labels are equal.
        symmetry_flipped (bool): Whether the symmetry is flipped.
    """

    def __init__(
        self,
        contact_frame: DataFrame,
        number_fragments: Optional[int] = None,
        metadata_combi: Optional[List[str]] = None,
        label_sorted: bool = False,
        binary_labels_equal: bool = False,
        symmetry_flipped: bool = False,
    ) -> None:
        self.contains_metadata = (
            "metadata_1" in contact_frame.columns
        )  # All contacts contain at least one fragment
        if number_fragments is None:
            self.number_fragments = self._guess_number_fragments(contact_frame)
        else:
            self.number_fragments = number_fragments
        self._schema = ContactSchema(
            number_fragments=self.number_fragments,
            contains_metadata=self.contains_metadata,
        )
        # TODO: make this work for duckdb pyrelation -> switch to mode
        if isinstance(contact_frame, pd.DataFrame):
            self.data_mode = DataMode.PANDAS
        elif isinstance(contact_frame, dd.DataFrame):
            self.data_mode = DataMode.DASK
        elif isinstance(contact_frame, duckdb.DuckDBPyRelation):
            self.data_mode = DataMode.DUCKDB
        else:
            raise ValueError("Unknown data mode!")
        self._data = self._schema.validate(contact_frame)
        self.metadata_combi = metadata_combi
        self.label_sorted = label_sorted
        self.binary_labels_equal = binary_labels_equal
        self.symmetry_flipped = symmetry_flipped

    @staticmethod
    def from_uri(uri, mode=DataMode.PANDAS):
        """Construct contacts from uri.
        Will match parameters based on the following order:

        PATH::number_fragments::metadata_combi::binary_labels_equal::symmetry_flipped::label_sorted

        Path, number_fragments are required. The rest is optional
        and will be tried to match to the available contacts. If no match is found, or there is no
         unique match, an error is raised.
        Mode can be one of pandas|dask, which corresponds to the type of the pixel source.
        """
        # import here to avoid circular imports
        from spoc.io import FileManager

        return FileManager(mode).load_contacts(uri)

    def get_global_parameters(self) -> ContactsParameters:
        """Returns global parameters"""
        return ContactsParameters(
            number_fragments=self.number_fragments,
            metadata_combi=self.metadata_combi,
            label_sorted=self.label_sorted,
            binary_labels_equal=self.binary_labels_equal,
            symmetry_flipped=self.symmetry_flipped,
        )

    def get_schema(self) -> GenomicDataSchema:
        """Returns the schema of the underlying data"""
        return self._schema

    def _guess_number_fragments(self, contact_frame: DataFrame) -> int:
        """Guesses the number of fragments from the contact frame"""
        return max(int(i.split("_")[1]) for i in contact_frame.columns if "start" in i)

    def get_label_values(self) -> List[str]:
        """Returns all label values"""
        # TODO: This could be put in global metadata of parquet file
        if not self.contains_metadata:
            raise ValueError("Contacts do not contain metadata!")
        output = set()
        for i in range(self.number_fragments):
            if self.data_mode == DataMode.DASK:
                output.update(self.data[f"metadata_{i+1}"].unique().compute())
            elif self.data_mode == DataMode.PANDAS:
                output.update(self.data[f"metadata_{i+1}"].unique())
            else:
                raise ValueError("Label values not supported for duckdb!")
        return list(output)

    def get_chromosome_values(self) -> List[str]:
        """Returns all chromosome values"""
        # TODO: This could be put in global metadata of parquet file
        output = set()
        for i in range(self.number_fragments):
            if self.data_mode == DataMode.DASK:
                output.update(self.data[f"chrom_{i+1}"].unique().compute())
            elif self.data_mode == DataMode.PANDAS:
                output.update(self.data[f"chrom_{i+1}"].unique())
            else:
                raise ValueError("Chromosome values not supported for duckdb!")
        return list(output)

    @property
    def data(self):
        """Returns the contact data"""
        return self._data

    @data.setter
    def data(self, contact_frame):
        """Sets the contact data

        Args:
            contact_frame (DataFrame): DataFrame containing the contact data.
        """
        self._data = self._schema.validate(contact_frame)

    def __repr__(self) -> str:
        return f"<Contacts | order: {self.number_fragments} | contains metadata: {self.contains_metadata}>"


class ContactManipulator:
    """Responsible for performing operations on
    contact data such as merging, splitting and subsetting."""

    def merge_contacts(self, merge_list: List[Contacts]) -> Contacts:
        """Merge contacts

        Args:
            merge_list (List[Contacts]): List of Contacts objects to merge.

        Returns:
            Contacts: Merged Contacts object.
        """
        # validate that merge is possible
        if len({i.number_fragments for i in merge_list}) != 1:
            raise ValueError("All contacts need to have the same order!")
        if len({i.data_mode for i in merge_list}) != 1:
            raise ValueError("Mixture of dataframes is not supported!")
        # TODO: assert all have same labelling state
        number_fragments = merge_list[0].number_fragments
        if merge_list[0].data_mode == DataMode.DASK:
            return Contacts(
                dd.concat([i.data for i in merge_list]),
                number_fragments=number_fragments,
            )
        elif merge_list[0].data_mode == DataMode.PANDAS:
            return Contacts(
                pd.concat([i.data for i in merge_list]),
                number_fragments=number_fragments,
            )
        else:
            raise ValueError("Merging duckdb relations is not supported!")

    @staticmethod
    def _generate_rename_columns(order, start_index=1):
        columns = [
            "chrom",
            "start",
            "end",
            "mapping_quality",
            "align_score",
            "align_base_qscore",
            "metadata",
        ]
        rename_columns = {}
        for i in range(len(order)):
            for column in columns:
                current_name = f"{column}_{i+start_index}"
                new_name = f"{column}_{order.index(i+start_index) + start_index}"
                rename_columns[current_name] = new_name
        return rename_columns

    @staticmethod
    def _get_label_combinations(labels, order):
        sorted_labels = sorted(labels)
        combinations = set(
            tuple(sorted(i)) for i in product(sorted_labels, repeat=order)
        )
        return combinations

    @staticmethod
    def _get_combination_splits(combination):
        splits = []
        for index, (i, j) in enumerate(zip(combination[:-1], combination[1:])):
            if i != j:
                splits.append(index + 2)
        return [1] + splits + [len(combination) + 1]

    def _flip_unlabelled_contacts(
        self,
        df: DataFrame,
        start_index: Optional[int] = None,
        end_index: Optional[int] = None,
    ) -> DataFrame:
        """Flips contacts"""
        fragment_order = max(int(i.split("_")[1]) for i in df.columns if "start" in i)
        if start_index is None:
            start_index = 1
        if end_index is None:
            end_index = fragment_order + 1
        subsets = []
        for perm in permutations(range(start_index, end_index)):
            query = "<=".join([f"start_{i}" for i in perm])
            subsets.append(
                df.query(query).rename(
                    columns=self._generate_rename_columns(perm, start_index)
                )
            )
        # determine which method to use for concatenation
        if isinstance(df, pd.DataFrame):
            result = pd.concat(subsets).sort_index()
            # this is needed if there are reads with equal start positions
            result = result.loc[~result.index.duplicated(keep="first")]
        else:
            result = (
                dd.concat(subsets)
                .reset_index()
                .sort_values("index")
                .drop_duplicates(subset=["index"])
                .set_index("index")
            )
        return result

    def _flip_labelled_contacts(
        self, df: DataFrame, label_values: List[str]
    ) -> DataFrame:
        """Flips labelled contacts"""
        fragment_order = max(int(i.split("_")[1]) for i in df.columns if "start" in i)
        label_combinations = self._get_label_combinations(label_values, fragment_order)
        subsets = []
        for combination in label_combinations:
            splits = self._get_combination_splits(combination)
            # separate out name constanc_columns
            query = " and ".join(
                [f"metadata_{i} == '{j}'" for i, j in enumerate(combination, 1)]
            )
            candidate_frame = df.query(query)
            if len(candidate_frame) == 0:
                continue
            constant_df, variable_df = candidate_frame[
                ["read_name", "read_length"]
            ], candidate_frame.drop(["read_name", "read_length"], axis=1)
            split_frames = [constant_df]
            for start, end in zip(splits, splits[1:]):
                # get all columns wiht nubmer between start and
                subset_columns = [
                    i
                    for i in variable_df.columns
                    if start <= int(i.split("_")[-1]) < end
                ]
                # if only columns is present, no need for flipping
                if start + 1 == end:
                    split_frame = variable_df[subset_columns]
                else:
                    split_frame = self._flip_unlabelled_contacts(
                        variable_df[subset_columns], start, end
                    )
                split_frames.append(split_frame)
            # concatenate split frames
            if isinstance(df, pd.DataFrame):
                subset = pd.concat(split_frames, axis=1)
            else:
                subset = dd.concat(split_frames, axis=1)
            subsets.append(subset)
        # determine which method to use for concatenation
        if isinstance(df, pd.DataFrame):
            result = pd.concat(subsets).sort_index()
            # this is needed if there are reads with equal start positions
            result = result.loc[~result.index.duplicated(keep="first")]
        else:
            result = (
                dd.concat(subsets)
                .reset_index()
                .sort_values("index")
                .drop_duplicates(subset=["index"])
                .set_index("index")
            )
        return result

    def sort_labels(self, contacts: Contacts) -> Contacts:
        """Sorts labels in ascending, alphabetical order

        Args:
            contacts (Contacts): Contacts object to sort.

        Returns:
            Contacts: Sorted Contacts object.
        """
        if not contacts.contains_metadata:
            raise ValueError(
                "Sorting labels for unlabelled contacts is not implemented."
            )
        # get label values.
        label_values = contacts.get_label_values()
        # iterate over all permutations of label values
        subsets = []
        for perm in product(label_values, repeat=contacts.number_fragments):
            query = " and ".join(
                [f"metadata_{i+1} == '{j}'" for i, j in enumerate(perm)]
            )
            desired_order = [i + 1 for i in np.argsort(perm)]
            subsets.append(
                contacts.data.query(query).rename(
                    columns=self._generate_rename_columns(desired_order)
                )
            )
        # determine which method to use for concatenation
        if contacts.data_mode == DataMode.DASK:
            # this is a bit of a hack to get the index sorted. Dask does not support index sorting
            result = (
                dd.concat(subsets).reset_index().sort_values("index").set_index("index")
            )
        elif contacts.data_mode == DataMode.PANDAS:
            result = pd.concat(subsets).sort_index()
        else:
            raise ValueError("Sorting labels for duckdb relations is not implemented.")
        return Contacts(
            result, number_fragments=contacts.number_fragments, label_sorted=True
        )

    def _sort_chromosomes(self, df: DataFrame, number_fragments: int) -> DataFrame:
        """Sorts chromosomes in ascending, alphabetical order"""
        # iterate over all permutations of chromosomes that exist
        subsets = []
        if isinstance(df, dd.DataFrame):
            chromosome_conbinations = (
                df[[f"chrom_{i}" for i in range(1, number_fragments + 1)]]
                .drop_duplicates()
                .compute()
                .values.tolist()
            )
        else:
            chromosome_conbinations = (
                df[[f"chrom_{i}" for i in range(1, number_fragments + 1)]]
                .drop_duplicates()
                .values.tolist()
            )
        for perm in chromosome_conbinations:
            query = " and ".join([f"chrom_{i+1} == '{j}'" for i, j in enumerate(perm)])
            desired_order = [i + 1 for i in np.argsort(perm, kind="stable")]
            sorted_frame = df.query(query).rename(
                columns=self._generate_rename_columns(desired_order)
            )
            # ensure correct column order
            subsets.append(sorted_frame)
        # determine which method to use for concatenation
        if isinstance(df, dd.DataFrame):
            # this is a bit of a hack to get the index sorted. Dask does not support index sorting
            result = (
                dd.concat(subsets).reset_index().sort_values("index").set_index("index")
            )
        else:
            result = pd.concat(subsets).sort_index()
        return result

    def _generate_binary_label_mapping(
        self, label_values: List[str], number_fragments: int
    ) -> Dict[Tuple[str, ...], Tuple[str, ...]]:
        sorted_labels = sorted(label_values)
        mapping = {}
        for i in range(number_fragments + 1):
            target = [sorted_labels[0]] * (number_fragments - i) + [
                sorted_labels[-1]
            ] * (i)
            source = [sorted_labels[0]] * (i) + [sorted_labels[-1]] * (
                number_fragments - i
            )
            if i <= (number_fragments // 2):
                mapping[tuple(source)] = tuple(target)
            else:
                mapping[tuple(source)] = ()
        return mapping

    def equate_binary_labels(self, contacts: Contacts) -> Contacts:
        """
        Equate binary labels.

        Binary labels often only carry information about whether
        they happen between the same or different fragments. This
        method equates these labels be replacing all equivalent binary labels with
        the alphabetically first label.
        For example, if we have a contact between two fragments
        that are labelled B and B, the label will be replaced with AA.

        Args:
            contacts (Contacts): Contacts object to equate binary labels.

        Returns:
            Contacts: Contacts object with equated binary labels.

        """
        assert contacts.contains_metadata, "Contacts do not contain metadata!"
        if not contacts.label_sorted:
            contacts = self.sort_labels(contacts)
        # get label values
        label_values = contacts.get_label_values()
        assert (
            len(label_values) == 2
        ), "Equate binary labels only works for binary labels!"
        # generate mapping diectionary
        mapping = self._generate_binary_label_mapping(
            label_values, contacts.number_fragments
        )
        subsets = []
        for source, target in mapping.items():
            query = " and ".join(
                [f"metadata_{i+1} == '{j}'" for i, j in enumerate(source)]
            )
            subset = contacts.data.query(query)
            # assign target labels to dataframe
            for i, j in enumerate(target):
                subset[f"metadata_{i+1}"] = j
            subsets.append(subset)
        # determine which method to use for concatenation
        if contacts.data_mode == DataMode.DASK:
            # this is a bit of a hack to get the index sorted. Dask does not support index sorting
            result = (
                dd.concat(subsets).reset_index().sort_values("index").set_index("index")
            )
        elif contacts.data_mode == DataMode.PANDAS:
            result = pd.concat(subsets).sort_index()
        else:
            raise ValueError(
                "Equate binary labels for duckdb relations is not implemented."
            )
        return Contacts(
            result,
            number_fragments=contacts.number_fragments,
            label_sorted=True,
            binary_labels_equal=True,
        )

    def subset_on_metadata(
        self, contacts: Contacts, metadata_combi: List[str]
    ) -> Contacts:
        """Subset contacts based on metadata

        Args:
            contacts (Contacts): Contacts object to subset.
            metadata_combi (List[str]): List of metadata combinations to subset on.

        Returns:
            Contacts: Subsetted Contacts object.

        """
        # check if metadata is present
        if not contacts.contains_metadata:
            raise ValueError("Contacts do not contain metadata!")
        # check if metadata_combi has the correct length
        if not len(metadata_combi) == contacts.number_fragments:
            raise ValueError("Metadata combination does not match number of fragments!")
        # get label values
        label_values = contacts.get_label_values()
        # check if metadata_combi is compatible with label values
        assert all(
            i in label_values for i in metadata_combi
        ), "Metadata combination is not compatible with label values!"
        # subset contacts
        query = " and ".join(
            [f"metadata_{i+1} == '{j}'" for i, j in enumerate(metadata_combi)]
        )
        result = contacts.data.query(query)
        return Contacts(
            result,
            number_fragments=contacts.number_fragments,
            metadata_combi=metadata_combi,
            label_sorted=contacts.label_sorted,
            binary_labels_equal=contacts.binary_labels_equal,
            symmetry_flipped=contacts.symmetry_flipped,
        )

    def flip_symmetric_contacts(
        self, contacts: Contacts, sort_chromosomes: bool = False
    ) -> Contacts:
        """Flips contacts based on inherent symmetry

        Args:
            contacts (Contacts): Contacts object to flip symmetric contacts.
            sort_chromosomes (bool, optional): Whether to sort chromosomes. Defaults to False.

        Returns:
            Contacts: Contacts object with flipped symmetric contacts.

        """
        if contacts.contains_metadata:
            if not contacts.label_sorted:
                contacts = self.sort_labels(contacts)
            label_values = contacts.get_label_values()
            result = self._flip_labelled_contacts(contacts.data, label_values)
            if sort_chromosomes:
                result = self._sort_chromosomes(result, contacts.number_fragments)
            return Contacts(
                result,
                number_fragments=contacts.number_fragments,
                label_sorted=True,
                binary_labels_equal=contacts.binary_labels_equal,
                symmetry_flipped=True,
            )
        result = self._flip_unlabelled_contacts(contacts.data)
        if sort_chromosomes:
            result = self._sort_chromosomes(result, contacts.number_fragments)
        return Contacts(
            result,
            number_fragments=contacts.number_fragments,
            symmetry_flipped=True,
        )
