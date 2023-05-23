"""Specific implementations of symmetry flipping"""
from abc import ABC, abstractmethod
from itertools import permutations
from typing import Union
import pandas as pd
import dask.dataframe as dd

DataFrame = Union[pd.DataFrame, dd.DataFrame]


class SymmetryFlipper(ABC):

    @abstractmethod
    def flip_contacts(self, df:DataFrame) -> DataFrame:
        """Flips contacts based on symmetry"""


class UnlabelledSymmetryFlipper(SymmetryFlipper):
    """Flips contacts based on symmetry when
    all contacts are equal.
     Logic -> all orders of contacts are equal
    """

    @staticmethod
    def _generate_rename_columns(order):
        columns = ["chrom", "start", "end", "mapping_quality", "align_score", "align_base_qscore"]
        rename_columns = {}
        for i in range(len(order)):
            for column in columns:
                current_name = f"{column}_{i+1}"
                new_name = f"{column}_{order[i]}"
                rename_columns[current_name] = new_name
        return rename_columns

    def flip_contacts(self, df: DataFrame) -> DataFrame:
        """Flips contacts"""
        fragment_order = max(int(i.split("_")[1]) for i in df.columns if "start" in i)
        subsets = []
        for perm in permutations(range(1, fragment_order+1)):
            query = "<=".join([f"start_{i}" for i in perm])
            subsets.append(df.query(query).rename(columns=self._generate_rename_columns(perm)))
        # determine which method to use for concatenation
        if isinstance(df, pd.DataFrame):
            result = pd.concat(subsets)
        else:
            result = dd.concat(subsets)
        return result
        


class LabelledSymmetryFlipper(SymmetryFlipper):
    """Flips contacts based on symmetry when
    all contacts are equal.
    Logic -> all orders of a certain marking state are equal
    e.g. for sisters ABA is the same as AAB and BAA, but not the same as AAA or BBB
    There is a parameter that decides whether ABB is the same as BAA (need to think of naming)
    For metadata with three different states, the last parameter is not important.
    ACC is different from ABB always.
    Sorts labels based on state so ABA will be returned as AAB. Default is alphabetical, 
    option to pass custom order.
    """

    def flip_contacts(self, df: DataFrame) -> DataFrame:
        """Flips contacts"""
        # TODO