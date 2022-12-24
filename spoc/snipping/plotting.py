from typing import Dict, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as ss
import seaborn as sns
from matplotlib.colors import LogNorm


def plot_1d_profile(profile: pd.DataFrame, color: Optional[str] = None, line_kwargs: Dict = {}, fill_kwargs: Dict = {}):
    plt.plot(profile.index, profile['avg'], color=color, **line_kwargs)
    plt.fill_between(profile.index, profile['avg'], profile['avg'] + profile['margin'], alpha=.3, color=color, **fill_kwargs)
    plt.fill_between(profile.index, profile['avg'], profile['avg'] - profile['margin'], alpha=.3, color=color, **fill_kwargs)
    plt.xticks(rotation=60)


def get_data_lims(mtx: pd.DataFrame) -> Tuple[float, float]:
    data_min = mtx.min().min()
    data_max = mtx.max().max()
    return data_min, data_max


def get_plot_lims(data_min: float, data_max: float) -> Tuple[float, float]:
    log2_highest_dev = max(abs(np.log2(data_max)), abs(np.log2(data_min)))
    vmax = 2 ** log2_highest_dev
    vmin = 2 ** (-log2_highest_dev)
    return vmin, vmax


class MidPointLogNorm(LogNorm):
    # https://stackoverflow.com/a/48632237 with modifications
    def __init__(self, vmin=None, vmax=None, center=None, base=2, clip=False):
        LogNorm.__init__(self, vmin=vmin, vmax=vmax, clip=clip)
        self.center = center
        self.base = base
    def __call__(self, value, clip=None):
        result, is_scalar = self.process_value(value)
        x = [np.log(self.vmin) / np.log(self.base), np.log(self.center) / np.log(self.base), np.log(self.vmax) / np.log(self.base)]
        y = [0, 0.5, 1]
        return np.ma.array(np.interp(np.log(value) / np.log(self.base), x, y), mask=result.mask, copy=False)

    
def plot_oe_mtx(mtx: pd.DataFrame, cmap: str = 'RdBu_r'):
    data_min, data_max = get_data_lims(mtx)
    vmin, vmax = get_plot_lims(data_min, data_max)
    g = sns.heatmap(mtx, cmap=cmap, norm=MidPointLogNorm(vmin, vmax, center=1), square=True, cbar_kws=dict(format='%.2e'))
    cbar = g.collections[0].colorbar
    cbar.ax.set_ylim(data_min, data_max)
    major_ticks = cbar.get_ticks()
    major_ticks = np.append(major_ticks, [data_min, data_max])
    major_ticks = [tick for tick in major_ticks if data_min <= tick <= data_max]
    cbar.set_ticks(major_ticks)
    minor_ticks = cbar.get_ticks(minor=True)
    bounded_minor_ticks = [tick for tick in minor_ticks if data_min <= tick <= data_max]
    minor_tick_intermode = ss.mode(np.diff(bounded_minor_ticks), keepdims=False).mode
    adjusted_minor_ticks = [tick for tick in bounded_minor_ticks if abs(tick - data_min) > minor_tick_intermode and abs(tick - data_max) > minor_tick_intermode]
    cbar.set_ticks(adjusted_minor_ticks, minor=True)
    return g
