"""Snipping strategies for Snipper that implement specific snipping functionality"""
from abc import ABC, abstractmethod
from functools import lru_cache, partial
from multiprocess.pool import ThreadPool
from uuid import uuid4
from enum import Enum
from typing import Union
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix
import duckdb
import bioframe
import dask.dataframe as dd

class SnippingValues(Enum):
    """Which values the snippet should consist of."""
    ICCF = 0
    OBSEXP = 1

class SnippingStrategy(ABC):
    """Defines interface for snipping strategy"""
    
    def __init__(self, bin_size:int, snipping_value: SnippingValues, **kwargs):
        """Defines the values that should be snipped"""
        self._snipping_value = snipping_value
        if snipping_value == SnippingValues.OBSEXP:
            self._n_random_regions = kwargs['n_random_regions']
            self._genome = kwargs['genome']
        self._bin_size = bin_size

    @staticmethod
    def _get_random_coordinates(n_coordinates: int, length: int, genome: str):
        """Number of coordinates will not be exactly returned, due to rounding
        when distributing to chromosomes."""
        chrom_sizes = bioframe.fetch_chromsizes(genome)
        chrom_fractions = chrom_sizes/chrom_sizes.sum()
        # accumulate output
        chrom_frames = []
        for chrom in chrom_sizes.index:
            max_coord = chrom_sizes[chrom]
            number_to_choose = int(chrom_fractions[chrom] * n_coordinates)
            starts = np.random.randint(0, max_coord, size=number_to_choose)
            ends = starts + length
            chrom_frames.append(pd.DataFrame({
                "chrom": chrom,
                "start": starts,
                "end":ends 
            }))
        return pd.concat(chrom_frames)


    @abstractmethod
    def snip(self, *args, **kwargs) -> pd.DataFrame:
        """Do snipping"""


class Triplet1DSnippingStrategy(SnippingStrategy):
    """Implements snipping of 2D-regions based on a set of 1D-regions.
    The 1D-regions are taken from the last column of pixels, which is assumed
    to contain contacts on the second sister chromatid. The first two pixels
    are assumed to contain contacts form the first sister chromatid."""

    CIS_SNIPPING_QUERY = """
        WITH individual_regions as (
            SELECT
                t.position_id,
                FLOOR((p.start_1 - t.pos)/{bin_size}::float) as offset_1,
                FLOOR((p.start_2 - t.pos)/{bin_size}::float) as offset_2,
                SUM(p.contact_count) as contacts
            FROM pixel_frame as p
            Inner JOIN {chunk_name} as t ON t.chrom = p.chrom
                and
                    abs(FLOOR((p.start_1 - t.pos)/{bin_size}::float))::int <= {pixel_offset}
                and
                    abs(FLOOR((p.start_2 - t.pos)/{bin_size}::float))::int <= {pixel_offset}
                and
                    abs(FLOOR((p.start_3 - (t.pos + {relative_offset}))/{bin_size}::float))::int <= {cis_selector_offset}
            GROUP BY 1,2,3
        )

        SELECT
            ind.offset_1 as offset_1,
            ind.offset_2 as offset_2,
            AVG(ind.contacts) as contacts
        FROM individual_regions as ind
        GROUP BY 1, 2
    """

    def _prepare_trans_positions(self, trans_positions: pd.DataFrame):
        """Adds index and round positions to bins"""
        if "pos" not in trans_positions.columns:
            trans_positions = trans_positions.assign(
                pos=lambda df_: (df_.start + df_.end)//2
            )
        return trans_positions.assign(
            position_id=lambda df_: range(len(df_)),
            pos=lambda df_: (df_.pos // self._bin_size) * self._bin_size,
        )

    @lru_cache(maxsize=100)
    def _create_expected_for_cis_snipping(
        self,
        pixels: pd.DataFrame,
        half_window_size: int,
        position_slack: int,
        relative_offset: int,
        threads: int = 2,
        
    ):
        random_regions = self._get_random_coordinates(self._n_random_regions, length=100, genome=self._genome)
        return self.snip(pixels, random_regions, half_window_size, position_slack, relative_offset, threads=threads, override_snipping_value=True)

    def snip(
        self,
        pixels: pd.DataFrame,
        trans_positions: pd.DataFrame,
        half_window_size: int,
        position_slack: Union[int, None] = None,
        relative_offset: int = 0,
        threads: int = 2,
        override_snipping_value: bool = False
    ):
        """Snips cis sister windows based on supplied trans positions."""
        # default is to take central pixel
        if position_slack is None:
            position_slack = self._bin_size
        # dispatch
        with ThreadPool(processes=threads) as pool:
            result = pool.map(
                partial(
                    self._snip_cis_window,
                    pixels=pixels,
                    half_window_size=half_window_size,
                    position_slack=position_slack,
                    relative_offset=relative_offset,
                ),
                np.array_split(self._prepare_trans_positions(trans_positions), threads),
            )
        # reduce result
        reduced = (
            pd.concat(result)
            .groupby(["offset_1", "offset_2"])
            .mean()
            .reset_index()
            .assign(
                array_index_1=lambda df_: df_.offset_1
                + (half_window_size // self._bin_size),
                array_index_2=lambda df_: df_.offset_2
                + (half_window_size // self._bin_size),
            )
        )
        # construct dense matrix
        output_size = 2 * (half_window_size // self._bin_size) + 1
        dense_matrix = coo_matrix(
            (
                reduced.contacts.values,
                (reduced.array_index_1.values, reduced.array_index_2.values),
            ),
            shape=(output_size, output_size),
        ).toarray()
        # check whether expected is needed
        if (self._snipping_value == SnippingValues.OBSEXP) and (not override_snipping_value):
            expected  = self._create_expected_for_cis_snipping(pixels, half_window_size, position_slack, relative_offset, threads)
            output = dense_matrix/expected
        else:
            output = dense_matrix
        genomic_extent = [f"{i * self._bin_size//1000 - half_window_size//1000} kb" for i in range(output_size)]
        return pd.DataFrame(output, columns=genomic_extent, index=genomic_extent)

    def _snip_cis_window(
        self,
        chunk: pd.DataFrame,
        pixels: pd.DataFrame,
        half_window_size: int,
        position_slack: int,
        relative_offset: int,
    ):
        # convert parameters into pixel units
        pixel_offset = half_window_size // self._bin_size
        cis_selector_offset = position_slack // self._bin_size
        # create local connection. No need to close it, is closed when reference to it goes out of scope
        local_connection = duckdb.connect()
        chunk_name = f"local_chunk_{uuid4().hex}"
        local_connection.register(chunk_name, chunk)
        # register pixels
        local_connection.register('pixel_frame', pixels)
        # do snipping
        return local_connection.execute(
            self.CIS_SNIPPING_QUERY.format(
                pixel_offset=pixel_offset,
                cis_selector_offset=cis_selector_offset,
                bin_size=self._bin_size,
                chunk_name=chunk_name,
                relative_offset=relative_offset,
            )
        ).df()