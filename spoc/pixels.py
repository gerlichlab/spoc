"""This part of spoc is responsible for binned,
higher order contacts in the form of 'genomic pixels'"""
from pathlib import Path
from typing import Union, Optional, List
import pandas as pd
import dask.dataframe as dd
from spoc.dataframe_models import PixelSchema
from spoc.file_parameter_models import PixelParameters
from spoc.contacts import Contacts



class Pixels:
    """Genomic pixels of arbitrary order.
    Contain information about:
        - Bin size
        - Symmetry (whehter contacts that where used to construct them where symmetric)
        - Order
        - Metadata combination (Whether the pixels represent a certain combination of metadata)
        - Whether binary labels are equal (e.g. whether AB pixles also represent BA pixels)

    Pixels can contain different data sources such as:
        - pandas dataframe
        - dask dataframe
        - path to a parquet file
    """

    def __init__(
        self,
        pixel_source: Union[pd.DataFrame, dd.DataFrame, str],
        number_fragments: Optional[int] = None,
        binsize: Optional[int] = None,
        metadata_combi: Optional[List[str]] = None,
        label_sorted: bool = False,
        binary_labels_equal: bool = False,
        symmetry_flipped:bool = False,
        same_chromosome: bool = True
    ):
        """Constructor for genomic pixels. pixel_source
        can be a pandas or dask dataframe or a path. Caveate is that
        if pixels are a path, source data is not validated."""
        self._schema = PixelSchema(number_fragments=number_fragments, same_chromosome=same_chromosome)
        self._same_chromosome = same_chromosome
        self._number_fragments = number_fragments
        self._binsize = binsize
        self._binary_labels_equal = binary_labels_equal
        self._symmetry_flipped = symmetry_flipped
        self._metadata_combi = metadata_combi
        self._label_sorted = label_sorted
        if isinstance(pixel_source, pd.DataFrame) or isinstance(
            pixel_source, dd.DataFrame
        ):
            self._data = self._schema.validate(pixel_source)
            self._path = None
        else:
            # check whether path exists
            if not Path(pixel_source).exists():
                raise ValueError(f"Path: {pixel_source} does not exist!")
            self._path = Path(pixel_source)
            self._data = None

    def from_uri(uri, mode="path"):
        """Construct pixels from uri.
        Will match parameters based on the following order:

        PATH::number_fragments::binsize::metadata_combi::binary_labels_equal::symmetry_flipped::label_sorted::same_chromosome

        PAth, number_fragments and binsize are required. The rest is optional
        and will be tried to match to the available pixels. If no match is found, or there is no
         uniue match, an error is raised.
        Mode can be one of pandas|dask|path, which corresponds to the type of the pixel source.
        """
        # import here to avoid circular imports
        from spoc.io import FileManager
        # Define uir parameters
        PARAMETERS = ['number_fragments', 'binsize', 'metadata_combi', 'binary_labels_equal', 'symmetry_flipped', 'label_sorted', 'same_chromosome']
        # parse uri
        uri = uri.split("::")
        # validate uri
        if len(uri) < 3:
            raise ValueError(
                f"Uri: {uri} is not valid. Must contain at least Path, number_fragments and binsize"
            )
        params = {
            key:value for key, value in zip(PARAMETERS, uri[1:])
        }
        # rewrite metadata_combi parameter
        if 'metadata_combi' in params.keys() and params['metadata_combi'] != 'None':
            params['metadata_combi'] = str(list(params['metadata_combi']))
        # read mode
        if mode == "path":
            load_dataframe = False
            use_dask = False
        elif mode == "pandas":
            load_dataframe = True
            use_dask = False
        else:
            load_dataframe = True
            use_dask = True
        # get availabe pixels
        available_pixels = FileManager().list_pixels(uri[0])
        # filter pixels
        matched_pixels = [
            pixel for pixel in available_pixels 
                    if all( params[key] == str(pixel.dict()[key]) for key in params.keys())
        ]
        # check whether there is a unique match
        if len(matched_pixels) == 0:
            raise ValueError(
                f"No pixels found for uri: {uri}"
            )
        elif len(matched_pixels) > 1:
            raise ValueError(
                f"Multiple pixels found for uri: {uri}"
            )
        return FileManager(use_dask=use_dask).load_pixels(uri[0], matched_pixels[0], load_dataframe=load_dataframe)



    def get_global_parameters(self):
        """Returns global parameters of pixels"""
        return PixelParameters(
            number_fragments=self._number_fragments,
            binsize=self._binsize,
            metadata_combi=self._metadata_combi,
            label_sorted=self._label_sorted,
            binary_labels_equal=self._binary_labels_equal,
            symmetry_flipped=self._symmetry_flipped,
            same_chromosome=self._same_chromosome
        )

    @property
    def path(self):
        return self._path

    @property
    def data(self):
        return self._data

    @property
    def number_fragments(self):
        return self._number_fragments

    @property
    def binsize(self):
        return self._binsize

    @property
    def binary_labels_equal(self):
        return self._binary_labels_equal

    @property
    def symmetry_flipped(self):
        return self._symmetry_flipped

    @property
    def metadata_combi(self):
        return self._metadata_combi

    @property
    def same_chromosome(self):
        return self._same_chromosome



class GenomicBinner:
    """Bins higher order contacts into genomic bins of fixed size.
    Is capable of sorting genomic bins along columns based on sister chromatid
    identity"""

    def __init__(
        self,
        bin_size: int
    ) -> None:
        self._bin_size = bin_size
        self._contact_order = None

    def _get_assigned_bin_output_structure(self):
        columns = (
            [
                    f'chrom_{index}' for index in range(1 , self._contact_order + 1)
            ] +
            [
                    f'start_{index}' for index in range(1 , self._contact_order + 1)
            ]
        )
        return pd.DataFrame(columns=columns).astype(int)

    def _assign_bins(self, data_frame: pd.DataFrame) -> pd.DataFrame:
        # capture empty dataframe
        if data_frame.empty:
            return self._get_assigned_bin_output_structure()
        return data_frame.assign(
            **{
                f"start_{index}": (data_frame[f"pos_{index}"] // self._bin_size) * self._bin_size
                    for index in range(1, self._contact_order + 1)
            }
        ).filter(regex='(chrom|start)')

    def _assign_midpoints(self, contacts: dd.DataFrame) -> dd.DataFrame:
        """Collapses start-end to a middle position"""
        return (
            contacts
                .assign(
                    **{
                        f"pos_{index}": (contacts[f'start_{index}'] + contacts[f'end_{index}'])//2
                            for index in range(1, self._contact_order + 1)
                    }
                )
                .drop(
                    [
                        c for index in range(1, self._contact_order + 1) for c in [f'start_{index}', f'end_{index}']
                    ], axis=1
                )
        )

    def bin_contacts(self, contacts: Contacts, same_chromosome: bool =True) -> dd.DataFrame:
        """Bins genomic contacts"""
        self._contact_order = contacts.number_fragments
        contacts_w_midpoints = self._assign_midpoints(
            contacts.data
        )
        if contacts.is_dask:
            contact_bins = contacts_w_midpoints.map_partitions(
                self._assign_bins,
                meta=self._get_assigned_bin_output_structure()
            )
        else:
            contact_bins = self._assign_bins(contacts_w_midpoints)
        pixels = (
            contact_bins.groupby(
                [
                    c for index in range(1, self._contact_order + 1)
                        for c in [f'chrom_{index}', f'start_{index}']
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
                .drop([f'chrom_{index}' for index in range(2, self._contact_order + 1)], axis=1)
                .rename(columns={"chrom_1": "chrom"})
            )
            # sort pixels
            pixels_sorted = pixels.sort_values(
                ['chrom'] + [f'start_{index}' for index in range(1, self._contact_order + 1)]
            ).reset_index(drop=True)
        else:
            pixels_sorted = pixels.sort_values(
                [f'chrom_{index}' for index in range(1, self._contact_order + 1)] +
                [f'start_{index}' for index in range(1, self._contact_order + 1)]
            ).reset_index(drop=True)
        # construct pixels and return
        return Pixels(
            pixels_sorted,
            same_chromosome=same_chromosome,
            number_fragments=self._contact_order,
            binsize=self._bin_size,
            binary_labels_equal=contacts.binary_labels_equal,
            symmetry_flipped=contacts.symmetry_flipped,
            metadata_combi=contacts.metadata_combi
        )


class PixelManipulator:
    """Has methods to manipulate pixels such as:
    - Coarsening
    - Balancing
    - Transferring weights
    """
