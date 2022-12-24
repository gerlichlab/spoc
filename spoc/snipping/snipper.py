"""Snipping class that orchestrates IO und snipping strategies to implement snipping behavior"""
from typing import List, Union, Dict, Optional
import pandas as pd
from spoc.snipping.snipping_strategies import SnippingStrategy
from spoc.pixels import PersistedPixels


class Snipper:
    def __init__(self, strategies: List[SnippingStrategy]) -> None:
        self._strategies = strategies

    def snip(
        self,
        pixels: Union[str, pd.DataFrame],
        snip_positions: pd.DataFrame,
        threads: int = 2,
        **kwargs: Dict,
    ) -> List[pd.DataFrame]:
        if isinstance(pixels, str):
            pixels = PersistedPixels(pixels)
        return [
            strategy.snip(pixels, snip_positions, threads=threads, **kwargs)
            for strategy in self._strategies
        ]
