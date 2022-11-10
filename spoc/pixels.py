"""This part of spoc is responsible for binned,
higher order contacts in the form of 'genomic pixels'"""
from pathlib import Path
import pandas as pd
import dask.dataframe as dd
import bioframe as bf
import pyranges as pr
from spoc.models import HigherOrderContactSchema


class PersistedPixels:
    """Wrapper for path to pixel parquet file.
    Responsible for validation. This breaks the IO vs domain logic layer,
    but improves performance of pileup."""

    def __init__(self, path: str):
        # check whether path exists
        if not Path(path).exists():
            raise ValueError(f"Path: {path} does not exist!")
        self._path = path

    @property
    def path(self):
        return self._path


class GenomicBinner:
    """Bins higher order contacts into genomic bins of fixed size.
    Is capable of sorting genomic bins along columns based on sister chromatid
    identity"""

    def __init__(
        self,
        bin_size: int,
        chrom_sizes: pd.Series,
        same_chromosome: bool = False,
        contact_order: int = 3,
        sort_sisters: bool = True,
        flip_contacts: bool = False
    ) -> None:
        self._bins = self._create_bins(chrom_sizes, bin_size)
        self._same_chromosome = same_chromosome
        if contact_order != 3:
            raise NotImplementedError(
                "Only contacts of order 3 are currently supported!"
            )
        self._contact_order = contact_order
        self._input_schema = HigherOrderContactSchema(contact_order)
        self._sort_sisters = sort_sisters
        self._flip_contacts = flip_contacts

    @staticmethod
    def _create_bins(chrom_sizes: pd.Series, bin_size: int) -> pr.PyRanges:
        """Creates genomic bins of size bin_size"""
        bins_df = bf.binnify(chrom_sizes, bin_size)
        bins_df.index.name = "bin_id"
        return pr.PyRanges(
            bins_df.reset_index().rename(
                columns={"start": "Start", "end": "End", "chrom": "Chromosome"}
            )
        )

    @staticmethod
    def _create_pyranges_bin_query(query_df: pd.DataFrame) -> pr.PyRanges:
        query_df.columns = ["Chromosome", "Start", "contact_index"]
        return pr.PyRanges(query_df.assign(End=lambda df: df.Start + 1))

    def _intersect_with_genomic_bins(
        self, suffix_index: int, query_position: pr.PyRanges
    ) -> pd.DataFrame:
        suffix = f"_{suffix_index}"
        return (
            query_position.join(self._bins, how="right", suffix=suffix)
            .df[["contact_index", "Chromosome", f"Start{suffix}"]]
            .rename(
                columns={
                    "Chromosome": f"chrom{suffix}",
                    f"Start{suffix}": f"start{suffix}",
                }
            )
            .query("contact_index != -1")
        )

    def _assign_bins(self, data_frame: pd.DataFrame) -> pd.DataFrame:
        # capture empty dataframe
        if data_frame.empty:
            return pd.DataFrame(
                columns=[
                    "contact_index",
                    "chrom_1",
                    "start_1",
                    "chrom_2",
                    "start_2",
                    "chrom_3",
                    "start_3",
                ],
                dtype=int,
            )
        # create contact_index
        data_frame.loc[:, "contact_index"] = range(len(data_frame))
        # create pyranges
        query_positions = map(
            self._create_pyranges_bin_query,
            [
                data_frame[["chrom_1", "pos_1", "contact_index"]],
                data_frame[["chrom_2", "pos_2", "contact_index"]],
                data_frame[["chrom_3", "pos_3", "contact_index"]],
            ],
        )
        # assign bins
        assigned_bins = list(
            map(self._intersect_with_genomic_bins, range(1, 4), query_positions)
        )
        # merge output
        merged_output = assigned_bins[0]
        for assigned_bin in assigned_bins[1:]:
            merged_output = merged_output.merge(
                assigned_bin, on="contact_index", how="inner"
            )
        return merged_output

    @staticmethod
    def _get_sister_order_selectors(contacts: dd.DataFrame):
        """Returns boolean arrays that indicate the order of sister contacts."""
        is_cis_cis_trans = (
            contacts.sister_identity_1 == contacts.sister_identity_2
        ) & (contacts.sister_identity_1 != contacts.sister_identity_3)
        is_cis_trans_cis = (
            contacts.sister_identity_1 == contacts.sister_identity_3
        ) & (contacts.sister_identity_1 != contacts.sister_identity_2)
        is_trans_cis_cis = (
            contacts.sister_identity_2 == contacts.sister_identity_3
        ) & (contacts.sister_identity_1 != contacts.sister_identity_2)
        return is_cis_cis_trans, is_cis_trans_cis, is_trans_cis_cis

    @staticmethod
    def _reset_contact_columns(contacts: dd.DataFrame) -> dd.DataFrame:
        contacts.columns = [
            "read_name",
            "chrom_1",
            "start_1",
            "end_1",
            "chrom_2",
            "start_2",
            "end_2",
            "chrom_3",
            "start_3",
            "end_3",
        ]
        return contacts

    def _sort_sister_contacts(self, contacts: dd.DataFrame) -> dd.DataFrame:
        """Sorts sister contacts such that they are SisterA_SisterA_SisterB"""
        (
            is_cis_cis_trans,
            is_cis_trans_cis,
            is_trans_cis_cis,
        ) = self._get_sister_order_selectors(contacts)
        # split dataframes into different sorting subsets
        sister_triplets_cis_cis_trans = contacts.loc[
            is_cis_cis_trans,
            [
                "read_name",
                "chrom_1",
                "start_1",
                "end_1",
                "chrom_2",
                "start_2",
                "end_2",
                "chrom_3",
                "start_3",
                "end_3",
            ],
        ]
        sister_triplets_cis_trans_cis = contacts.loc[
            is_cis_trans_cis,
            [
                "read_name",
                "chrom_1",
                "start_1",
                "end_1",
                "chrom_3",
                "start_3",
                "end_3",
                "chrom_2",
                "start_2",
                "end_2",
            ],
        ]
        sister_triplets_trans_cis_cis = contacts.loc[
            is_trans_cis_cis,
            [
                "read_name",
                "chrom_2",
                "start_2",
                "end_2",
                "chrom_3",
                "start_3",
                "end_3",
                "chrom_1",
                "start_1",
                "end_1",
            ],
        ]
        # concatenate and return
        return dd.concat(
            list(
                map(
                    self._reset_contact_columns,
                    [
                        sister_triplets_cis_cis_trans,
                        sister_triplets_cis_trans_cis,
                        sister_triplets_trans_cis_cis,
                    ],
                )
            )
        )

    def _flip_symmetric_contacts(self, contacts: dd.DataFrame) -> dd.DataFrame:
        """Flips contacts to upper triangular such that
        start_1 <= start_2. This is a Triplet specific implementation.
        TODO: Think about what this means for other contact order.
        Symmetry is not trivial for higher dimensional matrices (see 
        https://math.stackexchange.com/questions/615119/how-to-generalize-symmetry-for-higher-dimensional-arrays).
        Ideally, we would have a contact class that has a notion of symmetry, and flipping would be offloaded
        to that contact class. For example, non-sister sensitive triplets would have full symmetry, i.e. permutations
        of indices do not matter. However, sister-sensitive triplets have a more restricted notion of symmetry since
        only permutation of the first two indices leaves information unchanged. This method implements contact flipping
        based on the latter notion of symmetry.
        """
        # create boolean indexer
        is_lower_triangular = (contacts.chrom_1 == contacts.chrom_2) & (contacts.start_1 > contacts.start_2)
        # select and flip
        lower_flipped = contacts.loc[is_lower_triangular, :]\
                        .rename(
                            columns={"start_1": "start_2", "start_2": "start_1", "end_1": "end_2", "end_2": "end_1"}
                            )
        return dd.concat(
            [
                contacts.loc[~is_lower_triangular, :], # upper triangular
                lower_flipped
            ]
        )


    @staticmethod
    def _assign_midpoints(contacts: dd.DataFrame) -> dd.DataFrame:
        """Collapses start-end to a middle position"""
        return (
            contacts.assign(pos_1=lambda df: (df.start_1 + df.end_1) // 2)
            .assign(pos_2=lambda df: (df.start_2 + df.end_2) // 2)
            .assign(pos_3=lambda df: (df.start_3 + df.end_3) // 2)
            .drop(["start_1", "end_1", "start_2", "end_2", "start_3", "end_3"], axis=1)
        )

    def bin_contacts(self, contacts: dd.DataFrame) -> pd.DataFrame:
        """Bins genomic contacts"""
        # validate input header
        self._input_schema.validate_header(
            contacts
        )  # TODO: validate this in the contacts constructor
        # check if sister sorting should be performed
        if self._sort_sisters:
            contacts = self._sort_sister_contacts(
                contacts
            )  # TODO: move this into contacts class?
        if self._flip_contacts:
                contacts = self._flip_symmetric_contacts(contacts) # TODO: move this into contact class
        contacts_w_midpoints = self._assign_midpoints(
            contacts
        )  # TODO: move this into contacts class?
        # assign bins
        triplet_bins = contacts_w_midpoints.map_partitions(
            self._assign_bins,
            meta=pd.DataFrame(
                columns=[
                    "contact_index",
                    "chrom_1",
                    "start_1",
                    "chrom_2",
                    "start_2",
                    "chrom_3",
                    "start_3",
                ],
                dtype=int,
            ),
        )
        # compute pixels -> assumption is that pixels fit into memory after computing
        # TODO: think about how to achieve sorting without computing.
        # Maybe wait for https://github.com/dask/dask/issues/958
        pixels = (
            triplet_bins.groupby(
                ["chrom_1", "start_1", "chrom_2", "start_2", "chrom_3", "start_3"],
                observed=True,
            )
            .size()
            .reset_index()
            .rename(columns={0: "contact_count"})
        ).compute()
        if self._same_chromosome:
            pixels = (
                pixels.loc[
                    (pixels.chrom_1.astype(str) == pixels.chrom_2.astype(str))
                    & (pixels.chrom_2.astype(str) == pixels.chrom_3.astype(str))
                ]
                .drop(["chrom_2", "chrom_3"], axis=1)
                .rename(columns={"chrom_1": "chrom"})
            )
            return pixels.sort_values(
                ["chrom", "start_1", "start_2", "start_3"]
            ).reset_index(drop=True)
        return pixels.sort_values(
            ["chrom_1", "start_1", "chrom_2", "start_2", "chrom_3", "start_3"]
        ).reset_index(drop=True)


class PixelManipulator:
    """Has methods to manipulate pixels such as:
    - Coarsening
    - Balancing
    - Transferring weights
    """
