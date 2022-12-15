"""Snipping strategies for Snipper that implement specific snipping functionality"""
from abc import ABC, abstractmethod, abstractproperty
from enum import Enum
from functools import partial
from typing import Dict, List, Tuple, Union

import bioframe
import duckdb
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as ss
from multiprocess.pool import ThreadPool
from sparse import COO

from ..pixels import PersistedPixels


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


class TripletCCTSnippingStrategy(SnippingStrategy):
    SNIPPING_QUERY = """
        SELECT poi.position_id AS position_id,
               CEIL((triplets.start_1 - poi.pos) / {binsize}::float)::int AS bin_1,
               CEIL((triplets.start_2 - poi.pos) / {binsize}::float)::int AS bin_2,
               CEIL((triplets.start_3 - poi.pos) / {binsize}::float)::int AS bin_3,
               SUM(triplets.contact_count)::int AS contacts
        FROM {triplets} AS triplets
        INNER JOIN poi ON poi.chrom = triplets.chrom
        WHERE poi.snip_start <= triplets.start_1 AND triplets.start_1 < poi.snip_end AND
              poi.snip_start <= triplets.start_2 AND triplets.start_2 < poi.snip_end AND
              poi.snip_start <= triplets.start_3 AND triplets.start_3 < poi.snip_end
        GROUP BY position_id, bin_1, bin_2, bin_3
    """
    
    REDUCTION_QUERY = abstractproperty()
    
    @abstractmethod
    def _create_snip_regions(self, snip_positions: pd.DataFrame) -> pd.DataFrame:
        pass
    
    @abstractmethod
    def _get_snips(self, pixels, snip_regions, threads, **kwargs) -> pd.DataFrame:
        pass
    
    @abstractmethod
    def _flip_coords_by_strand(self, snip_regions: pd.DataFrame, snips: pd.DataFrame, strand: str):
        pass
    
    @staticmethod
    def _convert_coords(coords: pd.Series, strand: str) -> pd.Series:
        if strand == "+":
            return coords
        elif strand == "-":
            return -1 * coords
        else:
            raise ValueError

    @abstractmethod
    def _aggregate_snips(self, snips: pd.DataFrame):
        pass

    @abstractmethod
    def _get_genomic_coordinates(self) -> List[str]:
        pass


class TripletCCT2DSnippingStrategy(TripletCCTSnippingStrategy):
    REDUCTION_QUERY = """
        SELECT position_id, bin_3, SUM(contacts)::int AS contacts
        FROM full_snips
        WHERE {offset_11} <= bin_1 AND bin_1 <= {offset_12} AND
            {offset_21} <= bin_2 AND bin_2 <= {offset_22}
        GROUP BY position_id, bin_3
    """
    
    def __init__(self, bin_size, half_window_size, snipping_value, offsets, **kwargs):
        super().__init__(bin_size, half_window_size, snipping_value, **kwargs)
        self._offsets = self._parse_offsets(offsets)
        self._expected = None
    
    def _parse_offsets(self, supplied_offsets) -> Dict[str, int]:
        if isinstance(supplied_offsets, str):
            raise NotImplementedError
        elif isinstance(supplied_offsets, tuple):
            offsets = dict.fromkeys(('offset_11', 'offset_12', 'offset_21', 'offset_22'))
            if len(supplied_offsets) != 2:
                raise ValueError
            for i, offsets_set in enumerate(supplied_offsets):
                if len(offsets_set) != 2:
                    raise ValueError
                offset_set_1, offset_set_2 = offsets_set
                if not isinstance(offset_set_1, int) or not isinstance(offset_set_2, int):
                    raise ValueError
                offsets[f"offset_{i + 1}1"] = offset_set_1
                offsets[f"offset_{i + 1}2"] = offset_set_2
        else:
            raise ValueError
        return offsets
        
    def get_params(self):
        return None

    def snip(self, pixels, snip_positions, strand=None, threads=1, ci=.95) -> pd.DataFrame:
        snip_regions = self._create_snip_regions(snip_positions)
        snips = self._get_snips(pixels, snip_regions, threads)
        if strand is not None:
            self._flip_coords_by_strand(snip_regions, snips, strand)
        avg_profile, ci_margin = self._aggregate_snips(snips, ci=ci)
        return pd.DataFrame({'avg': avg_profile, 'margin': ci_margin}, index=self._get_genomic_coordinates())  # Maybe it would be better to put genomic coordinates as another column instead of index?
    
    def _create_snip_regions(self, snip_positions):
        if "pos" not in snip_positions.columns:
            snip_regions = snip_positions.assign(pos=lambda df_: (df_['start'] + df_['end']) // 2)
        else:
            snip_regions = snip_positions.copy()
        snip_regions['snip_start'] = snip_regions['pos'] - self._half_window_size - self._bin_size + 1
        snip_regions['snip_end'] = snip_regions['pos'] + self._half_window_size + 1
        snip_regions['position_id'] = range(len(snip_regions))
        return snip_regions
    
    def _get_snips(self, pixels, snip_regions, threads):
        con = duckdb.connect()
        con.execute(f"PRAGMA threads={threads};")
        con.register('poi', snip_regions)
        if isinstance(pixels, PersistedPixels):
            triplets = f"read_parquet('{pixels.path}')"
        else:
            triplets = "pixels"
            con.register(triplets, pixels)
        snips_rel = con.query(self.SNIPPING_QUERY.format(triplets=triplets, binsize=self._bin_size))
        reduced_rel = snips_rel.query('full_snips', self.REDUCTION_QUERY.format(offset_11=self._offsets['offset_11'] // self._bin_size,
                                                                                offset_12=self._offsets['offset_12'] // self._bin_size,
                                                                                offset_21=self._offsets['offset_21'] // self._bin_size,
                                                                                offset_22=self._offsets['offset_22'] // self._bin_size))
        return reduced_rel.df()
    
    def _flip_coords_by_strand(self, snip_regions: pd.DataFrame, snips: pd.DataFrame, strand: str):
        if not strand in snip_regions.columns:
            raise ValueError
        regions_strands = dict(snip_regions[['position_id', strand]].to_records(index=False))
        corrected_coords = snips.groupby('position_id', group_keys=False).apply(lambda df: self._convert_coords(df['bin_3'], regions_strands[df.name]))
        snips['bin_3'] = corrected_coords

    def _aggregate_snips(self, snips: pd.DataFrame, ci: float = .95) -> Tuple[pd.DataFrame, pd.DataFrame]:
        snip_size = 2 * (self._half_window_size // self._bin_size) + 1
        n_positions = snips['position_id'].nunique()
        bin_1_snip_size = abs(self._offsets['offset_12'] // self._bin_size - self._offsets['offset_11'] // self._bin_size) + 1
        bin_2_snip_size = abs(self._offsets['offset_22'] // self._bin_size - self._offsets['offset_21'] // self._bin_size) + 1
        snip_groupsize = bin_1_snip_size * bin_2_snip_size
        
        half_window_bins = self._half_window_size // self._bin_size
        
        sparse_matrix = COO((snips['position_id'].values, snips['bin_3'].values + half_window_bins),
                            data=snips['contacts'].values,
                            shape=(n_positions, snip_size))

        avg_profile = sparse_matrix.mean(axis=0).todense() / snip_groupsize
        var_profile = sparse_matrix.var(axis=0).todense() / snip_groupsize ** 2
        t_val = ss.t(n_positions - 1).isf((1 - ci) / 2)
        ci_margin = t_val * (var_profile / n_positions) ** 0.5
        return avg_profile, ci_margin
    
    def _get_genomic_coordinates(self):
        output_size = 2 * (self._half_window_size // self._bin_size) + 1
        return [f"{i * self._bin_size // 1000 - self._half_window_size // 1000} kb"
                for i in range(output_size)]


def plot_1d_profile(profile: pd.DataFrame, line_kwargs: Dict = {}, fill_kwargs: Dict = {}):
    plt.plot(profile.index, profile['avg'], color='black', **line_kwargs)
    plt.fill_between(profile.index, profile['avg'], profile['avg'] + profile['margin'], color='black', alpha=.3, **fill_kwargs)
    plt.fill_between(profile.index, profile['avg'], profile['avg'] - profile['margin'], color='black', alpha=.3, **fill_kwargs)
    plt.xticks(rotation=60)
