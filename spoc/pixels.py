"""This part of spoc is responsible for binned,
higher order contacts in the form of 'genomic pixels'"""
from functools import partial
from typing import List
from typing import Optional

import dask.dataframe as dd
import duckdb
import pandas as pd

from spoc.contacts import Contacts
from spoc.models.dataframe_models import DataFrame
from spoc.models.dataframe_models import DataMode
from spoc.models.dataframe_models import GenomicDataSchema
from spoc.models.dataframe_models import PixelSchema
from spoc.models.file_parameter_models import PixelParameters


class Pixels:
    """Genomic pixels of arbitrary order.
    Contain information about:
        - Bin size
        - Symmetry (whether contacts that where used to construct them where symmetric)
        - Order
        - Metadata combination (Whether the pixels represent a certain combination of metadata)
        - Whether binary labels are equal (e.g. whether AB pixels also represent BA pixels)


    Pixels can contain different data sources such as:
        - pandas dataframe
        - dask dataframe
        - path to a parquet file


    Args:
        pixel_source (Union[pd.DataFrame, dd.DataFrame, str]): The source of the pixel data.
        number_fragments (Optional[int], optional): The number of fragments. Defaults to None.
        binsize (Optional[int], optional): The bin size. Defaults to None.
        metadata_combi (Optional[List[str]], optional): The metadata combination. Defaults to None.
        label_sorted (bool, optional): Whether the labels are sorted. Defaults to False.
        binary_labels_equal (bool, optional): Whether binary labels are equal. Defaults to False.
        symmetry_flipped (bool, optional): Whether the pixels are symmetry flipped. Defaults to False.
        same_chromosome (bool, optional): Whether the pixels are on the same chromosome. Defaults to True.
    """

    def __init__(
        self,
        pixel_source: DataFrame,
        number_fragments: int,
        binsize: int,
        metadata_combi: Optional[List[str]] = None,
        label_sorted: bool = False,
        binary_labels_equal: bool = False,
        symmetry_flipped: bool = False,
        same_chromosome: bool = True,
    ):
        """Constructor for genomic pixels. pixel_source
        can be a pandas or dask dataframe or a path. Caveat is that
        if pixels are a path, source data is not validated."""
        self._schema = PixelSchema(
            number_fragments=number_fragments,
            same_chromosome=same_chromosome,
            binsize=binsize,
        )
        self._same_chromosome = same_chromosome
        self._number_fragments = number_fragments
        self._binsize = binsize
        self._binary_labels_equal = binary_labels_equal
        self._symmetry_flipped = symmetry_flipped
        self._metadata_combi = metadata_combi
        self._label_sorted = label_sorted
        # get data mode
        if isinstance(pixel_source, pd.DataFrame):
            self.data_mode = DataMode.PANDAS
        elif isinstance(pixel_source, dd.DataFrame):
            self.data_mode = DataMode.DASK
        elif isinstance(pixel_source, duckdb.DuckDBPyRelation):
            self.data_mode = DataMode.DUCKDB
        else:
            raise ValueError("Unknown data mode!")
        self._data = self._schema.validate(pixel_source)

    @staticmethod
    def from_uri(uri, mode=DataMode.PANDAS) -> "Pixels":
        """Construct pixels from uri.
        Will match parameters based on the following order:

        PATH::number_fragments::binsize::metadata_combi::binary_labels_equal::symmetry_flipped::label_sorted::same_chromosome

        Path, number_fragments and binsize are required. The rest is optional
        and will be tried to match to the available pixels. If no match is found, or there is no
         unique match, an error is raised.
        Mode can be one of pandas|dask|path, which corresponds to the type of the pixel source.


        Args:
            uri (str): The URI to construct the pixels from.
            mode (str, optional): The mode to use. Defaults to "path".

        Returns:
            Pixels: The constructed pixels.

        """
        # import here to avoid circular imports
        # pylint: disable=import-outside-toplevel
        from spoc.io import FileManager

        return FileManager(mode).load_pixels(uri)

    def get_global_parameters(self) -> PixelParameters:
        """Returns global parameters of pixels

        Returns:
            PixelParameters: The global parameters of the pixels.
        """
        return PixelParameters(
            number_fragments=self._number_fragments,
            binsize=self._binsize,
            metadata_combi=self._metadata_combi,
            label_sorted=self._label_sorted,
            binary_labels_equal=self._binary_labels_equal,
            symmetry_flipped=self._symmetry_flipped,
            same_chromosome=self._same_chromosome,
        )

    @property
    def data(self) -> DataFrame:
        """Returns pixels as dataframe

        Returns:
            DataFrame: The pixels as a dataframe.

        """
        return self._data

    @property
    def number_fragments(self) -> int:
        """Returns number of fragments in pixels

        Returns:
            int: The number of fragments in the pixels.
        """
        return self._number_fragments

    @property
    def binsize(self) -> int:
        """Returns binsize of pixels

        Returns:
            int: The binsize of the pixels.
        """
        return self._binsize

    @property
    def binary_labels_equal(self) -> bool:
        """Returns whether binary labels are equal

        Returns:
            bool: Whether binary labels are equal.
        """
        return self._binary_labels_equal

    @property
    def symmetry_flipped(self) -> bool:
        """Returns whether pixels are symmetry flipped

        Returns:
            bool: Whether pixels are symmetry flipped.
        """
        return self._symmetry_flipped

    @property
    def metadata_combi(self) -> Optional[List[str]]:
        """Returns metadata combination of pixels

        Returns:
            Optional[List[str]]: The metadata combination of the pixels.
        """
        return self._metadata_combi

    @property
    def same_chromosome(self) -> bool:
        """Returns whether pixels are on same chromosome


        Returns:
            bool: Whether pixels are on same chromosome.

        """
        return self._same_chromosome

    def get_schema(self) -> GenomicDataSchema:
        """Returns the schema of the underlying data"""
        return self._schema


class GenomicBinner:
    """Bins higher order contacts into genomic bins of fixed size.
    Is capable of sorting genomic bins along columns based on sister chromatid
    identity

    Args:
        bin_size (int): The size of the genomic bins.
    """

    def __init__(self, bin_size: int) -> None:
        self._bin_size = bin_size

    def _get_assigned_bin_output_structure(self, contact_order: int):
        columns = [f"chrom_{index}" for index in range(1, contact_order + 1)] + [
            f"start_{index}" for index in range(1, contact_order + 1)
        ]
        return pd.DataFrame(columns=columns).astype(int)

    def _assign_bins(
        self, data_frame: pd.DataFrame, contact_order: int
    ) -> pd.DataFrame:
        # capture empty dataframe
        if data_frame.empty:
            return self._get_assigned_bin_output_structure(contact_order)
        return data_frame.assign(
            **{
                f"start_{index}": (data_frame[f"pos_{index}"] // self._bin_size)
                * self._bin_size
                for index in range(1, contact_order + 1)
            }
        ).filter(regex="(chrom|start)")

    def _assign_midpoints(
        self, contacts: dd.DataFrame, contact_order: int
    ) -> dd.DataFrame:
        """Collapses start-end to a middle position"""
        return contacts.assign(
            **{
                f"pos_{index}": (contacts[f"start_{index}"] + contacts[f"end_{index}"])
                // 2
                for index in range(1, contact_order + 1)
            }
        ).drop(
            [
                c
                for index in range(1, contact_order + 1)
                for c in [f"start_{index}", f"end_{index}"]
            ],
            axis=1,
        )

    def bin_contacts(self, contacts: Contacts, same_chromosome: bool = True) -> Pixels:
        """Bins genomic contacts

        Args:
            contacts (Contacts): The genomic contacts to bin.
            same_chromosome (bool, optional): Whether to only retain pixels on the same chromosome. Defaults to True.

        Returns:
            Pixels: The binned genomic pixels.

        """
        contact_order = contacts.number_fragments
        contacts_w_midpoints = self._assign_midpoints(contacts.data, contact_order)
        if contacts.data_mode == DataMode.DASK:
            contact_bins = contacts_w_midpoints.map_partitions(
                partial(self._assign_bins, contact_order=contact_order),
                meta=self._get_assigned_bin_output_structure(contact_order),
            )
        elif contacts.data_mode == DataMode.PANDAS:
            contact_bins = self._assign_bins(contacts_w_midpoints, contact_order)
        else:
            raise ValueError(f"Data mode: {contacts.data_mode} not supported!")
        pixels = (
            contact_bins.groupby(
                [
                    c
                    for index in range(1, contact_order + 1)
                    for c in [f"chrom_{index}", f"start_{index}"]
                ],
                observed=True,
            )
            .size()
            .reset_index()
            .rename(columns={0: "count"})
        )
        # only retain pixels on same chromosome
        if same_chromosome:
            pixels = (
                pixels.loc[
                    (pixels.chrom_1.astype(str) == pixels.chrom_2.astype(str))
                    & (pixels.chrom_2.astype(str) == pixels.chrom_3.astype(str))
                ]
                .drop(
                    [f"chrom_{index}" for index in range(2, contact_order + 1)],
                    axis=1,
                )
                .rename(columns={"chrom_1": "chrom"})
            )
            # sort pixels
            pixels_sorted = pixels.sort_values(
                ["chrom"] + [f"start_{index}" for index in range(1, contact_order + 1)]
            ).reset_index(drop=True)
        else:
            pixels_sorted = pixels.sort_values(
                [f"chrom_{index}" for index in range(1, contact_order + 1)]
                + [f"start_{index}" for index in range(1, contact_order + 1)]
            ).reset_index(drop=True)
        # construct pixels and return
        return Pixels(
            pixels_sorted,
            same_chromosome=same_chromosome,
            number_fragments=contact_order,
            binsize=self._bin_size,
            binary_labels_equal=contacts.binary_labels_equal,
            symmetry_flipped=contacts.symmetry_flipped,
            metadata_combi=contacts.metadata_combi,
        )


class PixelManipulator:
    """Has methods to manipulate pixels such as:
    - Coarsening
    - Balancing
    - Transferring weights
    """
