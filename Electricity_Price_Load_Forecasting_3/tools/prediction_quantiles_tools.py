"""
Utility functions for quantile prediction
"""

# Author: Alessandro Brusaferri
# License: Apache-2.0 license


import numpy as np
from typing import List
import matplotlib.pyplot as plt
import pandas as pd

def build_alpha_quantiles_map(target_alpha: List, target_quantiles: List):
    """
    Build the map between PIs coverage levels and related quantiles
    """
    alpha_q = {'med': target_quantiles.index(0.5)}
    for alpha in target_alpha:
        alpha_q[alpha] = {
            'l': target_quantiles.index(alpha / 2),
            'u': target_quantiles.index(1 - alpha / 2),
        }
    return alpha_q


def fix_quantile_crossing(preds: np.array):
    """
    Fix crossing in the predicted quantiles by means of post-hoc sorting
    """
    return np.sort(preds, axis=-1)


def plot_quantiles(results: pd.DataFrame, target: str):
    """
    Plot predicted quantiles
    """
    title = target
    idx = results[target].index
    fig1, ax1 = plt.subplots()
    for i in results.columns.to_list():
        ax1.plot(idx, results[i], linestyle="-", color='steelblue', linewidth=0.9)

    ax1.plot(idx, results[target], '-', color='firebrick', label='$y_{true}$')
    ax1.grid()
    ax1.legend()
    ax1.set_ylabel("Predicted quantiles")
    ax1.set_title(title)
    fig1.show()
