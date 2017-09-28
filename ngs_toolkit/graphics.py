#!/usr/bin/env python


import matplotlib.pyplot as plt
import pandas as pd


# Set settings
pd.set_option("date_dayfirst", True)


def barmap(x, figsize=None, square=False, row_colors=None, z_score=None):
    """
    Plot a heatmap-style grid with barplots.

    :param pandas.DataFrame x: DataFrame with numerical values to plot. If DataFrame, indexes will be used as labels.
    :param tuple figsize: Size in inches (width, height) of figure to produce.
    :param bool square: Whether resulting figure should be square.
    :param list row_colors: Iterable of colors to use for each row.
    :param int z_score: Whether input matrix `x` should be Z-score transformed row-wise (0) or column-wise (1).
    :raises AssertionError: if length of `row_colors` does not match size of provided Y axis from matrix `x`.
    """
    y_size, x_size = x.shape

    # Check provided row_colors match provided matrix
    if row_colors is not None:
        assert len(row_colors) == y_size, "Length of row_colors does not match size of provided Y axis."

    # Z-score transform
    if z_score is not None:
        assert z_score in [1, 0], "z_score must be one of 0 (row-wise) or 1 (column-wise)."
        from scipy.stats import zscore
        rows, cols = x.index, x.columns
        x = zscore(x, axis=z_score)
        x = pd.DataFrame(x, index=rows, columns=cols)

    # Get figure dimentions if not given
    if figsize is None:
        square_size = 0.2  # in inches
        figsize = (x_size * square_size, y_size * square_size)
        if square:
            figsize = (figsize[0], figsize[0])
    fig, axis = plt.subplots(y_size, 1, figsize=figsize, gridspec_kw={"wspace": 0, "hspace": 0})

    # Plot row-by-row
    for i, row in enumerate(range(y_size)):
        color = row_colors[i] if row_colors is not None else None
        axis[i].bar(left=range(x_size), height=x.iloc[i, :], width=0.95, align="center", color=[color] * x_size if row_colors is not None else None)
        # axis[i].set_yticklabels(axis[i].get_yticklabels(), visible=False)
        if i != y_size - 1:
            axis[i].set_xticklabels(axis[i].get_xticklabels(), visible=False)
        axis[i].margins(0.0)

    axis[-1].set_xticks(range(x_size))
    axis[-1].set_xticklabels(x.columns, rotation='vertical', ha="center", va="center")

    for i, ax in enumerate(axis):
        m = max(ax.get_yticks())
        # ax.set_yticks([m])
        # ax.set_yticklabels([str(m)], rotation='horizontal', ha="center", va="center", visible=True)
        ax.set_ylabel(str(x.index[i]), rotation='horizontal', ha="center", va="center", visible=True)

    return fig
