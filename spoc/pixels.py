"""This part of spoc is responsible for binned higher order contacts in the form of 'genomic pixels'"""
from typing import Union
import pandas as pd
import dask.dataframe as dd
from cooler.util import binnify
import pyranges as pr


class GenomicBinner:
    """Bins higher order contacts into genomic bins of fixed size.
    Is capable of sorting genomic bins along columns based on sister chromatid
    identity"""

    def __init__(
        self, bin_size: int, chrom_sizes: pd.DataFrame, use_dask: bool = True
    ) -> None:
        self._use_dask = use_dask
        self._bins = self._create_bins(chrom_sizes, bin_size)

    def _create_bins(self, chrom_sizes: pd.DataFrame, bin_size: int) -> pr.PyRanges:
        bins_df = binnify(chrom_sizes, bin_size)
        bins_df.index.name = "bin_id"
        return pr.PyRanges(
            bins_df.reset_index().rename(
                columns={"start": "Start", "end": "End", "chrom": "Chromosome"}
            )
        )

    def _assign_bins(self, df:pd.DataFrame, bins:pr.PyRanges) -> pd.DataFrame:
        # create contact_index
        df.loc[:, "contact_index"] = range(len(df))
        # create pyranges
        bin_1 = pr.PyRanges(
            df[["chrom_1", "pos_1", "contact_index"]].rename(columns={"chrom_1": "Chromosome", "pos_1": "Start"})\
                                    .assign(End=lambda df: df.Start + 1)
        )
        bin_2 = pr.PyRanges(
            df[["chrom_2", "pos_2", "contact_index"]].rename(columns={"chrom_2": "Chromosome", "pos_2": "Start"})\
                                    .assign(End=lambda df: df.Start + 1)
        )
        bin_3 = pr.PyRanges(
            df[["chrom_3", "pos_3", "contact_index"]].rename(columns={"chrom_3": "Chromosome", "pos_3": "Start"})\
                                    .assign(End=lambda df: df.Start + 1)
        )
        # assign bin id
        bin_1_id = bin_1.join(bins, how="right").df[["contact_index", "bin_id"]].rename(columns={"bin_id": "bin_1_id"}).query("contact_index != -1")
        bin_2_id = bin_2.join(bins, how="right").df[["contact_index", "bin_id"]].rename(columns={"bin_id": "bin_2_id"}).query("contact_index != -1")
        bin_3_id = bin_3.join(bins, how="right").df[["contact_index", "bin_id"]].rename(columns={"bin_id": "bin_3_id"}).query("contact_index != -1")
        return bin_1_id.merge(bin_2_id, on="contact_index", how="inner").merge(bin_3_id, on="contact_index", how="inner")

    def _get_sister_order_selectors(self, contacts:dd.DataFrame):
        """Returns boolean arrays that indicate the order of sister contacts."""
        is_cis_cis_trans = (
            (contacts.sister_identity_1 == contacts.sister_identity_2) & (contacts.sister_identity_1 != contacts.sister_identity_3)
        )
        is_cis_trans_cis = (
            (contacts.sister_identity_1 == contacts.sister_identity_3) & (contacts.sister_identity_1 != contacts.sister_identity_2)
        )
        is_trans_cis_cis = (
            (contacts.sister_identity_2 == contacts.sister_identity_3) & (contacts.sister_identity_1 != contacts.sister_identity_2)
        )
        return is_cis_cis_trans, is_cis_trans_cis, is_trans_cis_cis

    def _reset_contact_columns(self, contacts:dd.DataFrame):
        output = contacts.copy()
        output.columns =  ["read_name","chrom_1", "start_1", "end_1", "chrom_2", "start_2", "end_2", "chrom_3", "start_3", "end_3"]
        return output

    def _sort_sister_contacts(self, contacts):
        """Sorts sister contacts such that they are SisterA_SisterA_SisterB"""
        is_cis_cis_trans, is_cis_trans_cis, is_trans_cis_cis = self._get_sister_order_selectors(contacts)
        # split dataframes into different sorting subsets
        sister_triplets_cis_cis_trans = contacts.loc[is_cis_cis_trans, ["read_name","chrom_1", "start_1", "end_1", "chrom_2", "start_2", "end_2", "chrom_3", "start_3", "end_3"]]
        sister_triplets_cis_trans_cis = contacts.loc[is_cis_trans_cis, ["read_name","chrom_1", "start_1", "end_1", "chrom_3", "start_3", "end_3", "chrom_2", "start_2", "end_2"]]
        sister_triplets_trans_cis_cis = contacts.loc[is_trans_cis_cis, ["read_name", "chrom_2", "start_2", "end_2", "chrom_3", "start_3", "end_3","chrom_1", "start_1", "end_1"]]
        # rename dataframes to same columns
        return dd.concat(list(map(self._reset_contact_columns, [sister_triplets_cis_cis_trans, sister_triplets_cis_trans_cis, sister_triplets_trans_cis_cis])))

    def _assign_midpoints(self, contacts:dd.DataFrame):
        """Collapses start end to a middle position"""
        return contacts.assign(pos_1=lambda df: (df.start_1 + df.end_1)//2)\
                                                    .assign(pos_2=lambda df: (df.start_2 + df.end_2)//2)\
                                                    .assign(pos_3=lambda df: (df.start_3 + df.end_3)//2)\
                                                    .drop(["start_1", "end_1", "start_2", "end_2", "start_3", "end_3"], axis=1)


    def bin_contacts(self, contacts: dd.DataFrame, sort_sisters: bool = True) -> pd.DataFrame:
        """Bins genomic contacts"""
        if sort_sisters:
            contacts = self._sort_sister_contacts(contacts)
        contacts_w_midpoints = self._assign_midpoints(contacts)
        # assign bins
        triplet_bins = contacts_w_midpoints.map_partitions(self._assign_bins, bins=self._bins, meta=pd.DataFrame(columns=["contact_index", "bin_1_id", "bin_2_id", "bin_3_id"], dtype=int))
        return triplet_bins.groupby(["bin_1_id", "bin_2_id", "bin_3_id"]).size().reset_index().rename(columns={0: "contact_count"})



class PixelManipulator:
    """Has methods to manipulate pixels such as:
        - Coarsening
        - Balancing
        - Transferring weights
    """