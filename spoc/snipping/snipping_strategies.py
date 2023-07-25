"""Snipping strategies for Snipper that implement specific snipping functionality"""
from abc import ABC, abstractmethod, abstractproperty
from enum import Enum
from functools import partial
from multiprocessing.pool import ThreadPool
from typing import Any, Dict, List, Optional, Tuple, Union

import bioframe as bf
import duckdb
import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy.stats as ss
from sparse import COO

from ..pixels import PersistedPixels

PixelsData = Union[pd.DataFrame, PersistedPixels]
IntArray = npt.NDArray[np.int_]
FloatArray = npt.NDArray[np.float_]


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
        chrom_sizes = bf.fetch_chromsizes(genome)
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
        cls_name = self.__class__.__name__
        return f"{cls_name}({', '.join((f'{key}={value}' for key, value in self.get_params().items()))})"

    @abstractmethod
    def get_params() -> Dict[str, Any]:
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


class TripletCCTSnippingStrategy(SnippingStrategy):
    SNIPPING_QUERY = """
        SELECT poi.position_id AS position_id,
               poi.strand_sign AS strand_sign,
               CEIL((triplets.start_1 - poi.pos) / {binsize}::float)::int AS bin_1,
               CEIL((triplets.start_2 - poi.pos) / {binsize}::float)::int AS bin_2,
               CEIL((triplets.start_3 - poi.pos) / {binsize}::float)::int AS bin_3,
               SUM(triplets.contact_count) AS contacts
        FROM {triplets} AS triplets
        INNER JOIN poi ON poi.chrom = triplets.chrom
        WHERE poi.snip_start <= triplets.start_1 AND triplets.start_1 < poi.snip_end AND
              poi.snip_start <= triplets.start_2 AND triplets.start_2 < poi.snip_end AND
              poi.snip_start <= triplets.start_3 AND triplets.start_3 < poi.snip_end
        GROUP BY position_id, strand_sign, bin_1, bin_2, bin_3
    """

    STRAND_FLIP_QUERY = """
        SELECT position_id,
               CASE strand_sign WHEN 1 THEN bin_1 WHEN -1 THEN -1 * bin_2 END AS bin_1,
               CASE strand_sign WHEN 1 THEN bin_2 WHEN -1 THEN -1 * bin_1 END AS bin_2,
               CASE strand_sign WHEN 1 then bin_3 WHEN -1 THEN -1 * bin_3 END AS bin_3,
               contacts
        FROM snips
    """
    
    REDUCTION_QUERY = abstractproperty()
    
    @abstractmethod
    def _get_snips(self, pixels: PixelsData, snip_regions: pd.DataFrame, threads: int, **kwargs) -> pd.DataFrame:
        pass

    @abstractmethod
    def _aggregate_snips(self, snips: pd.DataFrame):
        pass

    def _get_genomic_coordinates(self) -> List[str]:    
        output_size = 2 * (self._half_window_size // self._bin_size) + 1
        return [f"{i * self._bin_size // 1000 - self._half_window_size // 1000} kb"
                for i in range(output_size)]

    def _create_snip_regions(self,
                             snip_positions: pd.DataFrame,
                             strand: Optional[str] = None,
                             blacklist: Optional[pd.DataFrame] = None,
                             filterfield: str = 'pos') -> pd.DataFrame:
        if "pos" not in snip_positions.columns:
            snip_regions = snip_positions.assign(pos=lambda df_: (df_['start'] + df_['end']) // 2)
        else:
            snip_regions = snip_positions.copy()

        snip_regions['snip_start'] = snip_regions['pos'] - self._half_window_size - self._bin_size + 1
        snip_regions['snip_end'] = snip_regions['pos'] + self._half_window_size + 1

        if strand is None:
            snip_regions['strand_sign'] = 1
        elif strand not in snip_regions.columns:
            raise ValueError
        else:
            snip_regions['strand_sign'] = snip_regions[strand].map({"+": 1, "-": -1})
        

        if blacklist is not None:
            if filterfield == 'pos':
                snip_regions = snip_regions.assign(start=lambda df: df['pos'], end=lambda df: df['start'] + 1)
                snip_regions = bf.setdiff(snip_regions, blacklist).drop(['start', 'end'], axis=1)
            elif filterfield == 'snip':
                snip_regions = bf.setdiff(snip_regions, blacklist, cols1=('chrom', 'snip_start', 'snip_end'))
            else:
                raise ValueError('filterfield value must be one of "pos", "snip"')

        snip_regions['position_id'] = range(len(snip_regions))  # This is used later for averaging. After filtering by blacklist, there are less regions, which is accounted in averaging.
        return snip_regions


class TripletCCT1DSnippingStrategy(TripletCCTSnippingStrategy):
    REDUCTION_QUERY = """
        SELECT position_id, bin_1, bin_2, SUM(contacts) as contacts
        FROM full_snips
        WHERE {offset_1} <= bin_3 AND  bin_3 <= {offset_2}
        GROUP BY position_id, bin_1, bin_2
    """
    
    def __init__(self, bin_size: int, half_window_size: int, snipping_value: SnippingValues, offset: Union[int, Tuple[int, int]], **kwargs):
        super().__init__(bin_size, half_window_size, snipping_value, **kwargs)
        self._offset = self._parse_offset(offset)
    
    def _parse_offset(self, supplied_offset: Union[int, Tuple[int, int]]) -> Dict[str, int]:
        if isinstance(supplied_offset, int) or isinstance(supplied_offset, np.integer):
            supplied_offset = int(supplied_offset)
            offsets = {'offset_1': supplied_offset, 'offset_2': supplied_offset}
        elif isinstance(supplied_offset, tuple):
            if len(supplied_offset) != 2:
                raise ValueError(f"tuple offset must be of length 2.")
            first_coord, second_coord = supplied_offset
            if not isinstance(first_coord, int) or not isinstance(second_coord, int):
                raise ValueError(f"Both offsets in a tuple offset must be integers.")
            offset_1 = min(first_coord, second_coord)
            offset_2 = max(first_coord, second_coord)
            offsets = {'offset_1': offset_1, 'offset_2': offset_2}
        else:
            raise ValueError('offset must be either int or a tuple of 2 ints.')
        return offsets
    
    def get_params(self) -> Dict[str, Any]:
        offset = set(self._offset.values())
        if len(offset) == 1:
            (offset,) = offset
        else:
            offset = tuple(offset)
        return {"bin_size": self._bin_size,
                "half_window_size": self._half_window_size,
                "snipping_value": self._snipping_value,
                "offset": offset}
    
    def snip(self,
             pixels: PixelsData,
             snip_positions: pd.DataFrame,
             strand: Optional[str] = None,
             blacklist: Optional[pd.DataFrame] = None,
             filterfield: str = 'pos',
             threads: int = 1,
             symmetrize: bool = True,
             mem_limit: Optional[int] = None) -> pd.DataFrame:
        """Mem limit in GB"""
        snip_regions = self._create_snip_regions(snip_positions, strand, blacklist=blacklist, filterfield=filterfield)
        n_regions = len(snip_regions)  # This is used later for averaging. After filtering by blacklist, there are less regions, which is accounted in averaging.
        snips = self._get_snips(pixels, snip_regions, threads, mem_limit)
        avg_obs_matrix = self._aggregate_snips(snips, n_regions, symmetrize=symmetrize)
        if self._snipping_value == SnippingValues.ICCF:
            avg_matrix = avg_obs_matrix
        elif self._snipping_value == SnippingValues.OBSEXP:
            avg_exp_matrix = self._create_expected_snips(pixels, threads, blacklist, filterfield)
            avg_matrix = avg_obs_matrix / avg_exp_matrix
        else:
            raise NotImplementedError
        genomic_coordinates = self._get_genomic_coordinates()
        return pd.DataFrame(avg_matrix, columns=genomic_coordinates, index=genomic_coordinates)

    def _get_snips(self, pixels: PixelsData, snip_regions: pd.DataFrame, threads: int, mem_limit: Optional[int] = None) -> pd.DataFrame:
        if isinstance(pixels, PersistedPixels):
            if mem_limit is not None:
                mem_limit = mem_limit / threads
            with ThreadPool(threads) as p:
                snip_results = p.map(partial(self._execute_snipping, pixels=pixels, threads=1, mem_limit=mem_limit),
                                     np.array_split(snip_regions, threads))
                snips = pd.concat(snip_results, ignore_index=True)
        elif isinstance(pixels, pd.DataFrame):
            snips = self._execute_snipping(snip_regions, pixels, threads, mem_limit)
        else:
            raise NotImplementedError("pixels must be an instance of either PersistedPixels or pandas.DataFrame.")
        return snips

    def _execute_snipping(self, snip_regions: pd.DataFrame, pixels: PixelsData, threads: int = 1, mem_limit: Optional[int] = None) -> pd.DataFrame:
        con = duckdb.connect()
        con.execute(f"PRAGMA threads={threads};")
        if mem_limit is not None:
            con.execute(f"PRAGMA memory_limit='{mem_limit}GB';")
        con.register('poi', snip_regions)
        if isinstance(pixels, PersistedPixels):
            triplets = f"read_parquet('{pixels.path}')"
        else:
            triplets = "pixels"
            con.register(triplets, pixels)
        snips_rel = con.query(self.SNIPPING_QUERY.format(triplets=triplets, binsize=self._bin_size))
        strand_rel = snips_rel.query('snips', self.STRAND_FLIP_QUERY)
        reduced_rel = strand_rel.query('full_snips', self.REDUCTION_QUERY.format(offset_1=self._offset['offset_1'] // self._bin_size,
                                                                                 offset_2=self._offset['offset_2'] // self._bin_size))
        return reduced_rel.df()
    
    def _aggregate_snips(self, snips: pd.DataFrame, n_regions: int, symmetrize: bool = True) -> FloatArray:
        snip_size = 2 * (self._half_window_size // self._bin_size) + 1
        bin_3_snip_size = abs(self._offset['offset_2'] // self._bin_size - self._offset['offset_1'] // self._bin_size) + 1
        snip_groupsize = bin_3_snip_size

        half_window_bins = self._half_window_size // self._bin_size

        sparse_matrix = COO((snips['position_id'].values,
                             snips['bin_1'].values + half_window_bins,
                             snips['bin_2'].values + half_window_bins),
                            data=snips['contacts'].values,
                            shape=(n_regions, snip_size, snip_size))
        avg_matrix = sparse_matrix.mean(axis=0).todense() / snip_groupsize
        if symmetrize:
            avg_matrix = avg_matrix + avg_matrix.T - np.diag(np.diag(avg_matrix))
        return avg_matrix
    
    def _create_expected_snips(self, pixels: PixelsData, threads: int, blacklist: Optional[pd.DataFrame] = None, filterfield: str = 'pos') -> FloatArray:
        n_random = self._n_random_regions
        len_random = 100
        genome = self._genome
        random_intervals = self._get_random_coordinates(n_random, len_random, genome)
        random_snip_regions = self._create_snip_regions(random_intervals, blacklist=blacklist, filterfield=filterfield)
        n_random_snips = len(random_snip_regions)  # This is used later for averaging. After filtering by blacklist, there are less regions, which is accounted in averaging.
        random_snips = self._get_snips(pixels, random_snip_regions, threads)
        avg_random_matrix = self._aggregate_snips(random_snips, n_random_snips)
        return avg_random_matrix


class TripletCCT2DSnippingStrategy(TripletCCTSnippingStrategy):
    REDUCTION_QUERY = """
        SELECT position_id, bin_3, SUM(contacts) AS contacts
        FROM full_snips
        WHERE {offset_11} <= bin_1 AND bin_1 <= {offset_12} AND
            {offset_21} <= bin_2 AND bin_2 <= {offset_22}
        GROUP BY position_id, bin_3
    """
    _selection_order: Dict[str, int] = {'u': 1, 'c': 2, 'd': 3}
    
    def __init__(self, bin_size: int, half_window_size: int, snipping_value: SnippingValues, offsets: Union[str, Tuple[Tuple[int, int], Tuple[int, int]]], **kwargs):
        super().__init__(bin_size, half_window_size, snipping_value, **kwargs)
        self._offsets = self._parse_offsets(offsets)
    
    def _parse_offsets(self, supplied_offsets: Union[str, Tuple[Tuple[int, int], Tuple[int, int]]]) -> Dict[str, int]:
        if isinstance(supplied_offsets, str):
            if len(supplied_offsets) != 2:
                raise ValueError(f"The length of supplied offsets {supplied_offsets} is not equal to 2.")
            first_coord, second_coord = supplied_offsets
            if first_coord not in self._selection_order or second_coord not in self._selection_order:
                raise ValueError(f"Valid characters for str offsets are {' '.join(self._selection_order.keys())}")
            row_selection = max(first_coord, second_coord, key=lambda x: self._selection_order[x])  # careful here: it depends on the sorting order of the data
            col_selection = min(first_coord, second_coord, key=lambda x: self._selection_order[x])
            offset_11, offset_12 = self._get_offset_by_selection(row_selection)
            offset_21, offset_22 = self._get_offset_by_selection(col_selection)
            offsets = {"offset_11": offset_11, "offset_12": offset_12, "offset_21": offset_21, "offset_22": offset_22}
        elif isinstance(supplied_offsets, tuple):
            offsets = dict.fromkeys(('offset_11', 'offset_12', 'offset_21', 'offset_22'))
            if len(supplied_offsets) != 2:
                raise ValueError(f"The length of supplied offsets {supplied_offsets} is not equal to 2.")
            for i, offsets_set in enumerate(supplied_offsets):
                if len(offsets_set) != 2:
                    raise ValueError(f"The length of the pair {i} of offsets {offsets_set} is not equal to 2.")
                offset_set_1, offset_set_2 = offsets_set
                if not isinstance(offset_set_1, int) or not isinstance(offset_set_2, int):
                    raise ValueError(f"One of the offsets {offset_set_1}, {offset_set_2} are not instances of int.")
                offsets[f"offset_{i + 1}1"] = offset_set_1
                offsets[f"offset_{i + 1}2"] = offset_set_2
        else:
            raise ValueError(f"Supplied offsets {supplied_offsets} are not an instance of neither str nor tuple.")
        return offsets  # We'd better have it as a dataclass
    
    def _get_offset_by_selection(self, selection: str) -> Tuple[int, int]:
        pad_factor = self._half_window_size
        small_pad_from_zero = int(round(pad_factor / 3, -int(np.log10(self._bin_size))))
        if selection == 'u':
            return (-pad_factor, -small_pad_from_zero)
        elif selection == 'c':
            return (-small_pad_from_zero, small_pad_from_zero)
        elif selection == 'd':
            return (small_pad_from_zero, pad_factor)
        else:
            raise ValueError('selection is not one of u, c, d')
        
    def get_params(self):
        return {"bin_size": self._bin_size,
                "half_window_size": self._half_window_size,
                "snipping_value": self._snipping_value,
                "offsets": ((self._offsets['offset_11'], self._offsets['offset_12']),
                            (self._offsets['offset_21'], self._offsets['offset_22']))}

    def snip(self,
             pixels: PixelsData,
             snip_positions: pd.DataFrame,
             strand: Optional[str] = None,
             blacklist: Optional[pd.DataFrame] = None,
             filterfield: str = 'pos',
             threads: int = 1,
             ci: float = .95,
             mem_limit: Optional[int] = None) -> pd.DataFrame:
        snip_regions = self._create_snip_regions(snip_positions, strand, blacklist=blacklist, filterfield=filterfield)
        n_regions = len(snip_regions)  # This is used later for averaging. After filtering by blacklist, there are less regions, which is accounted in averaging.
        snips = self._get_snips(pixels, snip_regions, threads, mem_limit)
        avg_obs_profile, var_obs_profile = self._aggregate_snips(snips, n_regions)
        if self._snipping_value == SnippingValues.ICCF:
            avg_profile = avg_obs_profile
            ci_margin = self._ci_mean_1_sample(var_obs_profile, n_regions, ci=ci)
        elif self._snipping_value == SnippingValues.OBSEXP:
            avg_exp_profile, var_exp_profile = self._create_expected_snips(pixels, threads, blacklist, filterfield=filterfield)
            avg_profile = avg_obs_profile / avg_exp_profile
            n_exp = self._n_random_regions
            ci_margin = self._ci_ratio_mean_sample(avg_obs_profile, var_obs_profile, n_regions, avg_exp_profile, var_exp_profile, n_exp, ci)
        else:
            raise NotImplementedError
        return pd.DataFrame({'avg': avg_profile, 'margin': ci_margin}, index=self._get_genomic_coordinates())  # Maybe it would be better to put genomic coordinates as another column instead of index?
    
    def _get_snips(self, pixels: PixelsData, snip_regions: pd.DataFrame, threads: int, mem_limit: Optional[int] = None) -> pd.DataFrame:
        if isinstance(pixels, PersistedPixels):
            if mem_limit is not None:
                mem_limit = mem_limit / threads
            with ThreadPool(threads) as p:
                snip_results = p.map(partial(self._execute_snipping, pixels=pixels, threads=1, mem_limit=mem_limit),
                                     np.array_split(snip_regions, threads))
                snips = pd.concat(snip_results, ignore_index=True)
        elif isinstance(pixels, pd.DataFrame):
            snips = self._execute_snipping(snip_regions, pixels, threads, mem_limit)
        else:
            raise NotImplementedError("pixels must be an instance of either PersistedPixels or pandas.DataFrame.")
        return snips
    
    def _execute_snipping(self, snip_regions: pd.DataFrame, pixels: PixelsData, threads: int = 1, mem_limit: Optional[int] = None) -> pd.DataFrame:
        con = duckdb.connect()
        con.execute(f"PRAGMA threads={threads};")
        if mem_limit is not None:
            con.execute(f"PRAGMA memory_limit='{mem_limit}GB';")
        con.register('poi', snip_regions)
        if isinstance(pixels, PersistedPixels):
            triplets = f"read_parquet('{pixels.path}')"
        else:
            triplets = "pixels"
            con.register(triplets, pixels)
        snips_rel = con.query(self.SNIPPING_QUERY.format(triplets=triplets, binsize=self._bin_size))
        strand_rel = snips_rel.query('snips', self.STRAND_FLIP_QUERY)
        reduced_rel = strand_rel.query('full_snips', self.REDUCTION_QUERY.format(offset_11=self._offsets['offset_11'] // self._bin_size,
                                                                                 offset_12=self._offsets['offset_12'] // self._bin_size,
                                                                                 offset_21=self._offsets['offset_21'] // self._bin_size,
                                                                                 offset_22=self._offsets['offset_22'] // self._bin_size))
        return reduced_rel.df()

    def _aggregate_snips(self, snips: pd.DataFrame, n_regions: int) -> Tuple[FloatArray, FloatArray]:
        snip_size = 2 * (self._half_window_size // self._bin_size) + 1
        bin_1_snip_size = abs(self._offsets['offset_12'] // self._bin_size - self._offsets['offset_11'] // self._bin_size) + 1
        bin_2_snip_size = abs(self._offsets['offset_22'] // self._bin_size - self._offsets['offset_21'] // self._bin_size) + 1
        snip_groupsize = bin_1_snip_size * bin_2_snip_size
        
        half_window_bins = self._half_window_size // self._bin_size

        sparse_matrix = COO((snips['position_id'].values, snips['bin_3'].values + half_window_bins),
                            data=snips['contacts'].values,
                            shape=(n_regions, snip_size))

        avg_profile = sparse_matrix.mean(axis=0).todense() / snip_groupsize
        var_profile = sparse_matrix.var(axis=0).todense() / snip_groupsize ** 2
        return avg_profile, var_profile
    
    def _create_expected_snips(self, pixels: PixelsData, threads: int, blacklist: Optional[pd.DataFrame] = None, filterfield: str = 'pos') -> Tuple[FloatArray, FloatArray]:
        n_random = self._n_random_regions
        len_random = 100
        genome = self._genome
        random_intervals = self._get_random_coordinates(n_random, len_random, genome)
        random_snip_regions = self._create_snip_regions(random_intervals, blacklist=blacklist, filterfield=filterfield)
        n_random_snips = len(random_snip_regions)  # This is used later for averaging. After filtering by blacklist, there are less regions, which is accounted in averaging.
        random_snips = self._get_snips(pixels, random_snip_regions, threads)
        avg_random_profile, var_random_profile = self._aggregate_snips(random_snips, n_random_snips)
        return avg_random_profile, var_random_profile
    
    @staticmethod
    def _ci_mean_1_sample(var_profile: FloatArray, n_obs: int, ci: float) -> FloatArray:
        t_value = ss.t(n_obs - 1).isf((1 - ci) / 2)
        ci_margin = t_value * (var_profile / n_obs) ** 0.5
        return ci_margin
    
    @staticmethod
    def _ci_mean_2_sample(var_obs: FloatArray, n_obs: int, var_exp: FloatArray, n_exp: int, ci: float) -> FloatArray:
        t_value = ss.t(n_obs + n_exp - 2).isf((1 - ci) / 2)
        ci_margin = t_value * np.sqrt(var_obs / n_obs + var_exp / n_exp)
        return ci_margin
    
    @staticmethod
    def _ci_ratio_mean_sample(avg_obs: FloatArray, var_obs: FloatArray, n_obs: int, avg_exp: FloatArray, var_exp: FloatArray, n_exp: int, ci: float) -> FloatArray:
        t_value = ss.t(n_obs + n_exp - 2).isf((1 - ci) / 2)
        std_ratio = np.abs(avg_obs / avg_exp) * np.sqrt(var_obs / avg_obs ** 2 / n_obs + var_exp / avg_exp ** 2 / n_exp)
        ci_margin = t_value * std_ratio
        return ci_margin
