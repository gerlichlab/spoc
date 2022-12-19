"""Snipping strategies for Snipper that implement specific snipping functionality"""
from abc import ABC, abstractmethod
from functools import partial
from enum import Enum
from typing import Union
from multiprocess.pool import ThreadPool
import pandas as pd
import numpy as np
from sparse import COO
import duckdb
import bioframe
from spoc.pixels import PersistedPixels


class SnippingValues(Enum):
    """Which values the snippet should consist of."""

    ICCF = 0
    OBSEXP = 1


class SnippingStrategy(ABC):
    """Defines interface for snipping strategy"""

    def __init__(
        self,
        bin_size: int,
        half_window_size: int,
        snipping_value: Union[str, SnippingValues],
        **kwargs,
    ):
        """Defines the values that should be snipped"""
        if isinstance(snipping_value, str):
            # check whether string refers to a snipping strategy
            snipping_value = SnippingValues[snipping_value.upper()]
        self._snipping_value = snipping_value
        self._half_window_size = half_window_size
        if snipping_value == SnippingValues.OBSEXP:
            self._n_random_regions: int = kwargs.get("n_random_regions", 5000)
            self._genome: str = kwargs.get("genome", "hg19")
        self._bin_size = bin_size

    @staticmethod
    def _get_random_coordinates(n_coordinates: int, length: int, genome: str):
        """Number of coordinates will not be exactly returned, due to rounding
        when distributing to chromosomes."""
        chrom_sizes = bioframe.fetch_chromsizes(genome)
        chrom_fractions = chrom_sizes / chrom_sizes.sum()
        # accumulate output
        chrom_frames = []
        for chrom in chrom_sizes.index:
            max_coord = chrom_sizes[chrom]
            number_to_choose = int(chrom_fractions[chrom] * n_coordinates)
            starts = np.random.randint(0, max_coord, size=number_to_choose)
            ends = starts + length
            chrom_frames.append(
                pd.DataFrame({"chrom": chrom, "start": starts, "end": ends})
            )
        return pd.concat(chrom_frames)

    def __repr__(self) -> str:
        return f"<Triplet1DSnippingStrategy>"

    @abstractmethod
    def get_params():
        raise NotImplementedError

    @abstractmethod
    def snip(
        self,
        pixels: Union[pd.DataFrame, PersistedPixels],
        snip_positions: pd.DataFrame,
        threads: int = 2,
        override_snipping_value: bool = False,
    ) -> pd.DataFrame:
        """Do snipping"""
        raise NotImplementedError


class TripletCCT1DSnippingStrategy(SnippingStrategy):
    """Implements snipping of 2D-regions based on a set of 1D-regions.
    The 1D-regions are taken from the last column of pixels, which is assumed
    to contain contacts on the second sister chromatid. The first two pixels
    are assumed to contain contacts form the first sister chromatid."""

    CIS_SNIPPING_QUERY = """
        SELECT
            t.position_id,
            FLOOR((p.start_1 - t.pos)/{bin_size}::float) as offset_1,
            FLOOR((p.start_2 - t.pos)/{bin_size}::float) as offset_2,
            SUM(p.contact_count) as contacts
        FROM {source_table} as p
        INNER JOIN chunk as t ON t.chrom = p.chrom
            and
                abs(FLOOR((p.start_1 - t.pos)/{bin_size}::float))::int <= {pixel_offset}
            and
                abs(FLOOR((p.start_2 - t.pos)/{bin_size}::float))::int <= {pixel_offset}
            and
                abs(FLOOR((p.start_3 - (t.pos + {relative_offset}))/{bin_size}::float))::int <= {cis_selector_offset}
        GROUP BY 1,2,3
    """

    def __init__(
        self,
        bin_size: int,
        half_window_size: int,
        snipping_value: SnippingValues,
        **kwargs,
    ):
        """Override constructor to add additional parameters."""
        super().__init__(bin_size, half_window_size, snipping_value, **kwargs)
        self._position_slack: int = kwargs.get(
            "position_slack", self._bin_size
        )  # default is binsize
        self._relative_offset: int = kwargs.get("relative_offset", 0)
        self._expected = None

    def get_params(self) -> dict:
        return {
            "bin_size": self._bin_size,
            "half_window_size": self._half_window_size,
            "snipping_value": self._snipping_value,
            "position_slack": self._position_slack,
            "relative_offset": self._relative_offset,
        }

    def _align_positions_to_bins(self, trans_positions: pd.DataFrame):
        """Adds index and round positions to bins"""
        if "pos" not in trans_positions.columns:
            trans_positions = trans_positions.assign(
                pos=lambda df_: (df_.start + df_.end) // 2
            )
        return trans_positions.assign(
            position_id=lambda df_: range(len(df_)),
            pos=lambda df_: (df_.pos // self._bin_size) * self._bin_size,
        )

    def _get_array_coordinates_from_offset(self, offset: pd.Series):
        """Transform offsets to start from 0 to be used as array index"""
        return (offset + (self._half_window_size // self._bin_size)).astype(int).values

    def _reduce_snipping_frame(self, snips: pd.DataFrame):
        """Takes concat result of snipping and reduces
        it along the region dimension"""
        output_size = 2 * (self._half_window_size // self._bin_size) + 1
        return (
            COO(
                (
                    snips.position_id.values,
                    self._get_array_coordinates_from_offset(snips.offset_1),
                    self._get_array_coordinates_from_offset(snips.offset_2),
                ),
                snips.contacts.values,
                # TODO: fix shape such that if ids are missing from the end of the input, input shape will not be affected
                shape=(np.max(snips.position_id) + 1, output_size, output_size),
            )
            .mean(axis=0) # this reduces the result along the region id dimension
            .todense()
        )

    def _get_genomic_extent(self):
        output_size = 2 * (self._half_window_size // self._bin_size) + 1
        return [
            f"{i * self._bin_size//1000 - self._half_window_size//1000} kb"
            for i in range(output_size)
        ]

    def _create_expected_for_cis_snipping(
        self,
        pixels: Union[pd.DataFrame, PersistedPixels],
        threads: int = 2,
    ):
        if self._expected is None:
            random_regions = self._get_random_coordinates(
                self._n_random_regions, length=100, genome=self._genome
            )
            self._expected = self.snip(
                pixels,
                random_regions,
                threads=threads,
                override_snipping_value=True,
            )
            return self._expected
        return self._expected

    def snip(
        self,
        pixels: Union[pd.DataFrame, PersistedPixels],
        snip_positions: pd.DataFrame,
        threads: int = 2,
        override_snipping_value: bool = False,
    ):
        """Snips cis sister windows based on supplied trans positions."""
        # dispatch
        with ThreadPool(processes=threads) as pool:
            result = pool.map(
                partial(
                    self._snip_cis_windows,
                    pixels=pixels,
                ),
                np.array_split(self._align_positions_to_bins(snip_positions), threads),
            )
        # reduce along positions
        dense_matrix = self._reduce_snipping_frame(pd.concat(result))
        # check whether expected is needed
        if (self._snipping_value == SnippingValues.OBSEXP) and (
            not override_snipping_value
        ):
            expected = self._create_expected_for_cis_snipping(pixels, threads)
            output = dense_matrix / expected
        else:
            output = dense_matrix
        genomic_extent = self._get_genomic_extent()
        return pd.DataFrame(output, columns=genomic_extent, index=genomic_extent)

    def _snip_cis_windows(
        self, chunk: pd.DataFrame, pixels: Union[pd.DataFrame, PersistedPixels]
    ):
        # convert parameters into pixel units
        pixel_offset = self._half_window_size // self._bin_size
        cis_selector_offset = self._position_slack // self._bin_size
        # create local connection. No need to close it, is closed when reference to it goes out of scope
        local_connection = duckdb.connect()
        local_connection.register("chunk", chunk)
        # register pixels if needed
        if isinstance(pixels, PersistedPixels):
            source_table = f"read_parquet('{pixels.path}')"
        else:
            source_table = "pixel_frame"
            local_connection.register(source_table, pixels)
        # do snipping
        return local_connection.execute(
            self.CIS_SNIPPING_QUERY.format(
                source_table=source_table,
                pixel_offset=pixel_offset,
                cis_selector_offset=cis_selector_offset,
                bin_size=self._bin_size,
                relative_offset=self._relative_offset,
            )
        ).df()
