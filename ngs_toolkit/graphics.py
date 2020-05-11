#!/usr/bin/env python


import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ngs_toolkit import _CONFIG, _LOGGER


def savefig(fig, file_name, **kwargs):
    from ngs_toolkit.utils import record_analysis_output
    from ngs_toolkit.utils import get_timestamp

    if isinstance(fig, sns.axisgrid.Grid):
        fig = fig.fig
    default_kwargs = _CONFIG["preferences"]["graphics"]["figure_saving"]
    default_kwargs.update(kwargs)

    if _CONFIG["preferences"]["report"]["timestamp_figures"]:
        s = file_name.split(".")
        end = s[-1]
        f = ".".join(s[:-1])
        file_name = ".".join([f, get_timestamp(), end])

    fig.savefig(file_name, **default_kwargs)
    if _CONFIG["preferences"]["graphics"]["close_saved_figures"]:
        plt.close(fig)

    if _CONFIG["preferences"]["report"]["record_figures"]:
        record_analysis_output(file_name)


def clustermap_varieties(
        df, output_dir, output_prefix,
        steps=['base', 'z_score', 'sorted'],
        rasterized=True, labels=False,
        quantity="Expression", **kwargs):
    """
    CLustered heatmaps for various data transformations/orders.
    Specifically it will plot heatmaps with raw values, zscores,
    rows and columns clustered or ordered with provided indexes.

    Parameters
    ----------
    df : :obj:`pandas.DataFrame`
        DataFrame to plot.
    output_dir : :obj:`str`
        Directory to save plots.
    output_prefix : :obj:`str`
        Plot prefix.
    steps : :obj:`list`, optional
        Types of plots to produce. Defaults to all possible kinds.
    rasterized : :obj:`bool`, optional
        Whether to reasterize heatmaps. The default is True.
    labels : {:obj:`None`, :obj:`str`}, optional
        Whether row or column labels should be plotted.
        Pass a boolean to control both, or "rows" or "columns" to control each.
    quantity : :obj:`str`, optional
        A label for the numerical quantity. The default is "Expression".
    **kwargs : :obj:`dict`, optional
        Additional keyword arguments are passed to :func:`~seaborn.clustermap`.
    """
    import scipy
    from ngs_toolkit.graphics import (
        clustermap_fix_label_orientation as fixlabels, savefig)

    v = np.absolute(scipy.stats.zscore(df, axis=1)).flatten().max()
    v += v / 10.0

    common_kwargs = dict(
        rasterized=rasterized,
        xticklabels=True if labels in [True, "cols"] else False,
        yticklabels=True if labels in [True, "rows"] else False)
    c_kwargs = dict(cbar_kws={"label": quantity}, robust=True)
    cz_kwargs = dict(
        cbar_kws={"label": quantity + " Z-score"},
        z_score=0, cmap="RdBu_r", vmin=-v, vmax=v, center=0)
    s_kwargs = dict(**c_kwargs, row_cluster=False, col_cluster=False)
    sz_kwargs = dict(**cz_kwargs, row_cluster=False, col_cluster=False)
    sz_kwargs = dict(**cz_kwargs, row_cluster=False, col_cluster=False)

    b = "base" in steps
    s = "sorted" in steps
    z = "zscore" in steps
    t = "threshold" in steps
    pars = list()
    if b:
        pars += [("", c_kwargs)]
    if s:
        pars += [(".sorted", s_kwargs)]
    if z:
        pars += [(".z_score", cz_kwargs)]
    if s and z:
        pars += [(".z_score.sorted", sz_kwargs)]
    if t:
        pars += [""]

    for label, params in pars:
        file = os.path.join(output_dir, output_prefix + label + ".svg")
        savefig(fixlabels(sns.clustermap(df, **common_kwargs, **params, **kwargs)), file)


def barmap(x, figsize=None, square=False, row_colors=None, z_score=None, ylims=None):
    """
    Plot a heatmap-style grid with barplots.

    Parameters
    ----------
    x : :obj:`pandas.DataFrame`
        DataFrame with numerical values to plot.
        If DataFrame, indexes will be used as labels.
    figsize : :obj:`tuple`
        Size in inches (width, height) of figure to produce.
    square: :obj:`bool`
        Whether resulting figure should be square.
    row_colors : :obj:`list`
        Iterable of colors to use for each row.
    z_score : :obj:`int`
        Whether input matrix `x` should be Z-score transformed row-wise (0) or column-wise (1).

    Returns
    ----------
    :class:`matplotlib.pyplot.Figure`
        Figure object

    Raises
    ----------
    :obj:`ValueError`
        If length of ``row_colors`` does not match size of provided Y axis from matrix ``x``.
    """
    from scipy.stats import zscore

    y_size, x_size = x.shape

    # Check provided row_colors match provided matrix
    if row_colors is not None:
        msg = "Length of row_colors does not match size of provided Y axis."
        if not len(row_colors) == y_size:
            _LOGGER.error(msg)
            raise ValueError(msg)

    # Z-score transform
    if z_score is not None:
        if z_score not in [1, 0]:
            _LOGGER.error("z_score must be one of 0 (row-wise) or 1 (column-wise).")
            raise ValueError(msg)
        rows, cols = x.index, x.columns
        x = zscore(x, axis=z_score)
        x = pd.DataFrame(x, index=rows, columns=cols)

    # Get figure dimensions if not given
    if figsize is None:
        square_size = 0.2  # in inches
        figsize = (x_size * square_size, y_size * square_size)
        if square:
            figsize = (figsize[0], figsize[0])
    fig, axis = plt.subplots(
        y_size, 1, figsize=figsize, gridspec_kw={"wspace": 0, "hspace": 0}
    )

    # Plot row-by-row
    for i, row in enumerate(range(y_size)):
        color = row_colors[i] if row_colors is not None else None
        axis[i].bar(
            left=range(x_size),
            height=x.iloc[i, :],
            width=0.95,
            align="center",
            color=[color] * x_size if row_colors is not None else None,
        )
        # axis[i].set_yticklabels(axis[i].get_yticklabels(), visible=False)
        axis[i].axhline(0, linestyle="--", color="black", linewidth=0.5)
        if i != y_size - 1:
            axis[i].set_xticklabels(axis[i].get_xticklabels(), visible=False)
        axis[i].margins(0.0)
        if ylims is not None:
            axis[i].set_ylim(ylims)

    axis[-1].set_xticks(range(x_size))
    axis[-1].set_xticklabels(
        x.columns, rotation="vertical", ha="center", va="top", fontsize=figsize[1]
    )

    for i, ax in enumerate(axis):
        # m = max(ax.get_yticks())
        # ax.set_yticks([m])
        # ax.set_yticklabels([str(m)], rotation='horizontal', ha="center",
        #                    va="center", visible=True, fontsize=figsize[1])
        ax.set_ylabel(
            str(x.index[i]),
            rotation="horizontal",
            ha="right",
            va="center",
            visible=True,
            fontsize=figsize[1],
        )

    return fig


# def bubbleplot(x, y, radius_var=None, data=None, hue=None, ax=None, orient="vertical"):
#     """
#     """

#     if data is not None:
#         x = data[x]
#         y = data[y]

#     if ax is None:
#         fig, ax = plt.subplots(1)

#     # For the case when x is categorical
#     cat_num = len(x.unique())


def radar_plot(
        data,
        subplot_var="patient_id",
        group_var="timepoint",
        radial_vars=["NaiveBcell", "SwitchedBcell", "UnswitchedBcell"],
        cmap="inferno",
        scale_to_max=True,
):
    """
    Heavy inspiration from here: https://matplotlib.org/examples/api/radar_chart.html

    Parameters
    ----------
    data : :obj:`pandas.DataFrame`
    subplot_var : :obj:`str`
    group_var : :obj:`str`
    radial_vars : :obj:`list`
    cmap : :obj:`str`
        Matplotlib colormap to use.
    scale_to_max: :obj:`bool`
        Whether values will be scaled
    """
    from matplotlib.projections.polar import PolarAxes
    from matplotlib.path import Path
    from matplotlib.projections import register_projection
    from matplotlib.spines import Spine

    def radar_factory(num_vars, frame="circle"):
        """Create a radar chart with `num_vars` axes.

        This function creates a RadarAxes projection and registers it.

        Parameters
        ----------
        num_vars : :obj:`int`
            Number of variables for radar chart.
        frame : {'circle' | 'polygon'}
            Shape of frame surrounding axes.

        """
        # calculate evenly-spaced axis angles
        theta = np.linspace(0, 2 * np.pi, num_vars, endpoint=False)
        # rotate theta such that the first axis is at the top
        theta += np.pi / 2

        def draw_poly_patch(self):
            verts = unit_poly_verts(theta)
            return plt.Polygon(verts, closed=True, edgecolor="k")

        def draw_circle_patch(self):
            # unit circle centered on (0.5, 0.5)
            return plt.Circle((0.5, 0.5), 0.5)

        patch_dict = {"polygon": draw_poly_patch, "circle": draw_circle_patch}
        if frame not in patch_dict:
            raise ValueError("unknown value for `frame`: %s" % frame)

        class RadarAxes(PolarAxes):
            name = "radar"
            # use 1 line segment to connect specified points
            RESOLUTION = 1
            # define draw_frame method
            draw_patch = patch_dict[frame]

            def fill(self, *args, **kwargs):
                """Override fill so that line is closed by default"""
                closed = kwargs.pop("closed", True)
                return super(RadarAxes, self).fill(closed=closed, *args, **kwargs)

            def plot(self, *args, **kwargs):
                """Override plot so that line is closed by default"""
                lines = super(RadarAxes, self).plot(*args, **kwargs)
                for line in lines:
                    self._close_line(line)

            @staticmethod
            def _close_line(line):
                x, y = line.get_data()
                # FIXME: markers at x[0], y[0] get doubled-up
                if x[0] != x[-1]:
                    x = np.concatenate((x, [x[0]]))
                    y = np.concatenate((y, [y[0]]))
                    line.set_data(x, y)

            def set_varlabels(self, labels):
                self.set_thetagrids(np.degrees(theta), labels)

            def _gen_axes_patch(self):
                return self.draw_patch()

            def _gen_axes_spines(self):
                if frame == "circle":
                    return PolarAxes._gen_axes_spines(self)
                # The following is a hack to get the spines (i.e. the axes frame)
                # to draw correctly for a polygon frame.

                # spine_type must be 'left', 'right', 'top', 'bottom', or `circle`.
                spine_type = "circle"
                verts = unit_poly_verts(theta)
                # close off polygon by repeating first vertex
                verts.append(verts[0])
                path = Path(verts)

                spine = Spine(self, spine_type, path)
                spine.set_transform(self.transAxes)
                return {"polar": spine}

        register_projection(RadarAxes)
        return theta

    def unit_poly_verts(theta):
        """Return vertices of polygon for subplot axes.

        This polygon is circumscribed by a unit circle centered at (0.5, 0.5)
        """
        x0, y0, r = [0.5] * 3
        verts = [(r * np.cos(t) + x0, r * np.sin(t) + y0) for t in theta]
        return verts

    theta = radar_factory(len(radial_vars), frame="polygon")

    n_row = n_col = int(np.ceil(np.sqrt(len(data[subplot_var].unique()))))

    fig, axes = plt.subplots(
        figsize=(9, 9), nrows=n_row, ncols=n_col, subplot_kw=dict(projection="radar")
    )
    fig.subplots_adjust(wspace=0.25, hspace=0.20, top=0.85, bottom=0.05)

    colors = plt.get_cmap(cmap, len(data[group_var].unique()))
    # Plot the four cases from the example data on separate axes
    for ax, subplot_title in zip(axes.flatten(), data[subplot_var].unique()):

        case_data = data[data[subplot_var] == subplot_title]

        if scale_to_max:
            case_data.loc[:, radial_vars] = (
                case_data.loc[:, radial_vars] / case_data.loc[:, radial_vars].max()
            )

        ax.set_rgrids([0.2, 0.4, 0.6, 0.8])
        ax.set_title(
            subplot_title,
            weight="bold",
            size="medium",
            position=(0.5, 1.1),
            horizontalalignment="center",
            verticalalignment="center",
        )
        for group in case_data[group_var].unique():
            d = case_data.loc[case_data[group_var] == group, radial_vars].squeeze()
            color = colors(np.log2(1 + group) / np.log2(1 + 256))
            ax.plot(theta, d, label=group, color=color)
            ax.fill(theta, d, label=group, alpha=0.25, facecolor=color)
            # ax.set_ylim(0.7, 0.85)
        ax.set_varlabels(radial_vars)

        # add legend relative to top-left plot
        ax.legend(loc=(0.9, 0.95), labelspacing=0.1, fontsize="small")

    return fig


def add_colorbar_to_axis(collection, label=None, position="right", size="5%", pad=0.05):
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(collection.axes)
    cax = divider.append_axes(position, size=size, pad=pad)
    plt.colorbar(mappable=collection, cax=cax, label=label, alpha=1)


def clustermap_fix_label_orientation(grid, fontsize="xx-small", **kwargs):
    grid.ax_heatmap.set_xticklabels(
        grid.ax_heatmap.get_xticklabels(), rotation=90, fontsize=fontsize, **kwargs
    )
    grid.ax_heatmap.set_yticklabels(
        grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize=fontsize, **kwargs
    )
    return grid


def clustermap_rasterize_heatmap(grid):
    q = [
        x
        for x in grid.ax_heatmap.get_children()
        if isinstance(x, matplotlib.collections.QuadMesh)
    ][0]
    q.set_rasterized(True)


def plot_projection(
        df,
        color_dataframe,
        dims,
        output_file,
        attributes_to_plot,
        plot_max_dims=8,
        rasterized=False,
        plot_group_centroids=True,
        axis_ticklabels=True,
        axis_ticklabels_name="PC",
        axis_lines=True,
        legends=False,
        always_legend=False,
):
    """
    Plot a low dimentionality projection of samples.

    Parameters
    ----------
    df : :obj:`pandas.DataFrame`
        Dataframe with sample projections.

    color_dataframe : :obj:`pandas.DataFrame`
        Dataframe of RGB tuples for sample i in attribute j.

    dims : :obj:`int`
        Number of dimensions to plot

    output_file : :obj:`str`
        Path to figure output file

    attributes_to_plot : :obj:`list`
        List of levels in df.index to plot

    plot_max_dims : number, optional
        Maximum number of dimensions to plot.
        Defaults to 8.

    plot_group_centroids: :obj:`bool`, optional
        Whether centroids of each sample group should be plotted alongside samples.
        Will be square shaped.
        Defaults to True.

    axis_ticklabels: :obj:`bool`, optional
        Whether axis ticks and tick labels should be plotted.
        Defaults to False.

    axis_lines: :obj:`bool`, optional
        Whether (0, 0) dashed lines should be plotted.
        Defaults to True.

    legends: :obj:`bool`, optional
        Whether legends for group colours should be plotted.
        Defaults to False.

    always_legend: :obj:`bool`, optional
        Whether legends for group colours should be plotted in every figure panel.
        If False, will plot just on first/last figure panel.
        Defaults to False.
    """
    from collections import OrderedDict

    # Make sure dataframes are index-sorted
    df = df.sort_index()
    color_dataframe = color_dataframe.sort_index()

    n_attr = len(attributes_to_plot)
    fig, axis = plt.subplots(dims, n_attr, figsize=(4 * n_attr, 4 * dims))
    if n_attr == dims == 1:
        axis = np.array([[axis]])
    elif (n_attr == 1) and (dims > 1):
        axis = axis.reshape((dims, 1))
    elif (n_attr > 1) and (dims == 1):
        axis = axis.reshape((1, n_attr))
    for pc in range(dims):
        for i, attr in enumerate(attributes_to_plot):

            # TODO: make plotting faster specially for many samples
            for j, sample in enumerate(df.index):
                sample = pd.Series(sample, index=df.index.names)
                label = df.index.get_level_values(attr)[j]
                axis[pc, i].scatter(
                    df.loc[sample["sample_name"], pc],
                    df.loc[sample["sample_name"], pc + 1],
                    s=30,
                    color=color_dataframe.loc[sample["sample_name"], attr],
                    alpha=0.75,
                    label=label,
                    rasterized=rasterized,
                )

            # Plot groups
            values = df.index.get_level_values(attr)
            numeric = (
                (not any([isinstance(x, str) for x in values])) and not
                all((values.dropna() == True) | (values.dropna() == False)))
            if plot_group_centroids and not numeric:
                try:
                    df2 = df.groupby(attr).mean()
                except IndexError:
                    _LOGGER.error("Could not plot centroids of projection!")
                    continue
                # get the color of each attribute group
                cd = color_dataframe.loc[:, attr]
                cd.name = None
                cd.index = df.index.get_level_values(attr)
                cd = cd.apply(
                    lambda x: tuple(x) if isinstance(x, list) else x
                )  # fix for deduplicating lists
                cd = cd.reset_index().drop_duplicates().set_index(attr).drop_duplicates()

                for j, group in enumerate(df2.index):
                    if group not in df2.index:
                        continue
                    axis[pc, i].scatter(
                        df2.loc[group, pc],
                        df2.loc[group, pc + 1],
                        marker="s",
                        s=50,
                        color=cd.loc[group].squeeze(),
                        alpha=0.95,
                        label=group,
                        rasterized=rasterized,
                    )
                    axis[pc, i].text(
                        df2.loc[group, pc],
                        df2.loc[group, pc + 1],
                        group,
                        color=cd.loc[group].squeeze(),
                        alpha=0.95,
                    )

            # Graphics
            if pc == 0:
                axis[pc, i].set_title(attr)
            axis[pc, i].set_xlabel("{} {}".format(axis_ticklabels_name, pc + 1))
            if i == 0:
                axis[pc, i].set_ylabel("{} {}".format(axis_ticklabels_name, pc + 2))
            if not axis_ticklabels:
                axis[pc, i].set_xticklabels(
                    axis[pc, i].get_xticklabels(), visible=False
                )
                axis[pc, i].set_yticklabels(
                    axis[pc, i].get_yticklabels(), visible=False
                )
            if axis_lines:
                axis[pc, i].axhline(0, linestyle="--", color="black", alpha=0.3)
                axis[pc, i].axvline(0, linestyle="--", color="black", alpha=0.3)

            if legends:
                # Unique legend labels
                handles, labels = axis[pc, i].get_legend_handles_labels()
                by_label = OrderedDict(zip(labels, handles))
                if (
                    any([isinstance(c, str) for c in by_label.keys()])
                    and len(by_label) <= plot_max_dims
                ):
                    # if not any([re.match("^\d", c) for c in by_label.keys()]):
                    if always_legend:
                        axis[pc, i].legend(by_label.values(), by_label.keys())
                    else:
                        if pc == (dims - 1):
                            axis[pc, i].legend(by_label.values(), by_label.keys())
    savefig(fig, output_file)


def plot_region_context_enrichment(
    enr,
    output_dir="results",
    output_prefix="region_type_enrichment",
    across_attribute=None,
    pvalue=0.05,
    top_n=5,
):
    """
    Plot results of ATACSeqAnalysis.region_context_enrichment.

    Parameters
    ----------
    enr : :obj:`pandas.DataFrame`
        Results of region_context_enrichment.

    output_dir : :obj:`str`, optional
        Directory to save plots to.
        Defaults to "results".

    optional output_prefix : :obj:`str`, optional
        Prefix to use when saveing plots.
        Defaults to "region_type_enrichment"

    across_attribute : :obj:`str`, optional
        Column in enrichment matrix to plot results across e.g. "comparison_name"
        when results matrix contains the result of various comparisons.
        Defaults to None (not used).

    pvalue : float, optional
        Value at which to plot a line marking the significant level.
        Defaults to 0.05.

    top_n : :obj:`int`, optional
        Number of features to label in volcano plot.
        Defaults to 5.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # replace inf p-values with
    x, y, z = "log2_odds_ratio", "-log10(p-value)", "region_type"
    enr = enr.sort_values(x)

    if across_attribute is not None:
        level = enr.loc[:, across_attribute].unique()
        n = len(level)
    else:
        level = [""]
        across_attribute = "intercept"
        enr.loc[:, across_attribute] = ""
        n = 1
    row = col = int(np.ceil(np.sqrt(n)))

    # barplots
    for var in [x, y]:
        fig, axis = plt.subplots(row, col, figsize=(col * 3, row * 3))
        if row == col == 1:
            axis = np.array([axis])
        axis = axis.flatten()
        for i, value in enumerate(level):
            axis[i].set_title(value)
            sns.barplot(
                data=enr.loc[enr[across_attribute] == value, :].reset_index(),
                x=var,
                y="region",
                hue=z,
                orient="horizontal",
                dodge=False,
                ax=axis[i],
            )
            if (enr[var] < 0).any():
                m = enr[var].abs().max()
                m += m * 0.1
                axis[i].set_xlim((-m, m))
                axis[i].axvline(0, linestyle="--", color="grey")
        savefig(
            fig, os.path.join(output_dir, output_prefix + ".{}.barplot.svg".format(var))
        )

    # volcano plot
    fig, axis = plt.subplots(row, col, figsize=(col * 3, row * 3))
    if row == col == 1:
        axis = np.array([axis])
    axis = axis.flatten()
    for i, value in enumerate(level):
        spec = enr.loc[enr[across_attribute] == value, :]
        axis[i].set_title(value)
        axis[i].axvline(0, linestyle="--", color="grey")
        axis[i].axhline(-np.log10(pvalue), linestyle="--", color="grey")
        sns.scatterplot(data=spec.reset_index(), x=x, y=y, hue=z, ax=axis[i])
        # label top points
        for n, s in spec.nlargest(top_n, y).iterrows():
            axis[i].text(s[x], s[y], s=n)
        # equal x axis
        ll = spec.loc[:, x].abs().max()
        ll += ll * 0.1
        axis[i].set_xlim((-ll, ll))
    savefig(
        fig,
        os.path.join(
            output_dir, output_prefix + ".volcano_plot.top_{}_labeled.svg".format(top_n)
        ),
    )


def plot_comparison_correlations(
    diff, output_dir, output_prefix="comparison_correlations"
):
    """
    Plot pairwise log fold changes for various comparisons.


    Parameters
    ----------
    diff : :obj:`pandas.DataFrame`
        Dataframe with differential results

    output_dir : :obj:`str`
        Output directory for plots.

    output_prefix : :obj:`str`, optional
        Prefix for plots.
        Defaults to "comparison_correlations"
    """
    import scipy

    comps = diff["comparison_name"].unique()
    n = len(comps)
    rows = cols = int(np.ceil(np.sqrt(n)))
    v = diff["log2FoldChange"].abs().max()
    # vmax = (-np.log10(diff['pvalue'])).max()

    fig, axis = plt.subplots(
        rows, cols, figsize=(cols * 3, rows * 3), tight_layout=True
    )
    for i, (base, comps) in enumerate(comps):
        a = diff.loc[diff["comparison_name"] == base, "log2FoldChange"]
        sign = (a > 0).astype(int).replace(0, -1)
        c = sign * (-np.log10(diff.loc[diff["comparison_name"] == base, "pvalue"]))
        vmax = np.sqrt(c.abs().max())
        for j, comp in enumerate(comps):
            b = diff.loc[diff["comparison_name"] == comp, "log2FoldChange"]
            axis[i, j].set_xlabel(base)
            axis[i, j].set_ylabel(comp)
            # equal limits and x = y line
            axis[i, j].set_xlim((-v, v))
            axis[i, j].set_ylim((-v, v))
            axis[i, j].plot((-v, -v), (v, v), linestyle="--", alpha=0.75, color="black")
            # quadrant lines
            axis[i, j].axhline(0, linestyle="--", alpha=0.5, color="black")
            axis[i, j].axvline(0, linestyle="--", alpha=0.5, color="black")

            # fit and plot regression line
            d = a.to_frame(name=base).join(b.to_frame(name=comp)).dropna()
            slope, intercept, r, p, error = scipy.stats.linregress(d[base], d[comp])
            x = np.linspace(-v, v, 1000)
            y = x * slope + intercept
            axis[i, j].plot(x, y, color="orange")
            axis[i, j].text(
                -7,
                7,
                s="y = {:.2f}x + {:.2f}\n".format(slope, intercept)
                + "r = {:.3f}\n".format(r)
                + "p = {:.3E}\n".format(p),
                fontsize=7,
                va="top",
            )

            # lastly the actual scatter
            axis[i, j].scatter(
                a,
                b,
                s=2,
                alpha=0.5,
                rasterized=True,
                c=c,
                cmap="RdBu_r",
                vmin=-vmax,
                vmax=vmax,
            )  # color="grey")

        for ax in axis[i, (j + 1):]:
            ax.set_visible(False)
    savefig(fig, os.path.join(output_dir, output_prefix + ".lmfit.svg"))

    fig, axis = plt.subplots(
        rows, cols, figsize=(cols * 3, rows * 3), tight_layout=True
    )
    for i, (base, comps) in enumerate(comps):
        a = diff.loc[diff["comparison_name"] == base, "log2FoldChange"]
        sign = (a > 0).astype(int).replace(0, -1)
        c = sign * (-np.log10(diff.loc[diff["comparison_name"] == base, "pvalue"]))
        vmax = np.sqrt(c.abs().max())
        for j, comp in enumerate(comps):
            b = diff.loc[diff["comparison_name"] == comp, "log2FoldChange"]
            axis[i, j].set_xlabel(base)
            axis[i, j].set_ylabel(comp)
            # equal limits and x = y line
            axis[i, j].set_xlim((-v, v))
            axis[i, j].set_ylim((-v, v))
            axis[i, j].plot((-v, -v), (v, v), linestyle="--", alpha=0.75, color="black")
            # quadrant lines
            axis[i, j].axhline(0, linestyle="--", alpha=0.5, color="black")
            axis[i, j].axvline(0, linestyle="--", alpha=0.5, color="black")

            # fit and plot regression line
            d = a.to_frame(name=base).join(b.to_frame(name=comp)).dropna()
            slope, intercept, r, p, error = scipy.stats.linregress(d[base], d[comp])
            x = np.linspace(-v, v, 1000)
            y = x * slope + intercept
            axis[i, j].plot(x, y, color="orange")
            axis[i, j].text(
                -7,
                7,
                s="y = {:.2f}x + {:.2f}\n".format(slope, intercept)
                + "r = {:.3f}\n".format(r)
                + "p = {:.3E}\n".format(p),
                fontsize=7,
                va="top",
            )

            # lastly the actual scatter
            axis[i, j].scatter(a, b, s=2, alpha=0.5, rasterized=True, color="grey")

        for ax in axis[i, (j + 1) :]:
            ax.set_visible(False)
    savefig(fig, os.path.join(output_dir, output_prefix + ".lmfit.no_color.svg"))
