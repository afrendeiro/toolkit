#!/usr/bin/env python


import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ngs_toolkit import _CONFIG
from ngs_toolkit import _LOGGER


def barmap(x, figsize=None, square=False, row_colors=None, z_score=None, ylims=None):
    """
    Plot a heatmap-style grid with barplots.

    :param pandas.DataFrame x: DataFrame with numerical values to plot. If DataFrame,
                               indexes will be used as labels.
    :param tuple figsize: Size in inches (width, height) of figure to produce.
    :param bool square: Whether resulting figure should be square.
    :param list row_colors: Iterable of colors to use for each row.
    :param int z_score: Whether input matrix `x` should be Z-score transformed row-wise (0)
                        or column-wise (1).
    :raises AssertionError: if length of `row_colors` does not match size of provided Y axis
                            from matrix `x`.
    """
    y_size, x_size = x.shape

    # Check provided row_colors match provided matrix
    if row_colors is not None:
        msg = "Length of row_colors does not match size of provided Y axis."
        if not len(row_colors) == y_size:
            _LOGGER.error(msg)
            raise AssertionError(msg)

    # Z-score transform
    if z_score is not None:
        if z_score not in [1, 0]:
            _LOGGER.error("z_score must be one of 0 (row-wise) or 1 (column-wise).")
            raise AssertionError(msg)
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
        axis[i].bar(
            left=range(x_size), height=x.iloc[i, :],
            width=0.95, align="center", color=[color] * x_size if row_colors is not None else None)
        # axis[i].set_yticklabels(axis[i].get_yticklabels(), visible=False)
        axis[i].axhline(0, linestyle="--", color="black", linewidth=0.5)
        if i != y_size - 1:
            axis[i].set_xticklabels(axis[i].get_xticklabels(), visible=False)
        axis[i].margins(0.0)
        if ylims is not None:
            axis[i].set_ylim(ylims)

    axis[-1].set_xticks(range(x_size))
    axis[-1].set_xticklabels(
        x.columns, rotation='vertical', ha="center", va="top", fontsize=figsize[1])

    for i, ax in enumerate(axis):
        # m = max(ax.get_yticks())
        # ax.set_yticks([m])
        # ax.set_yticklabels([str(m)], rotation='horizontal', ha="center",
        #                    va="center", visible=True, fontsize=figsize[1])
        ax.set_ylabel(
            str(x.index[i]), rotation='horizontal', ha="right", va="center",
            visible=True, fontsize=figsize[1])

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
        subplot_var="patient_id", group_var="timepoint",
        radial_vars=["NaiveBcell", "SwitchedBcell", "UnswitchedBcell"],
        cmap="inferno", scale_to_max=True):
    """
    
    :param pandas.DataFrame data:
    :param str subplot_var:
    :param str group_var:
    :param list radial_vars:
    :param cmap str: Matplotlib colormap to use.
    :param bool scale_to_max: Whether values will be scaled

    Heavy inspiration from here: https://matplotlib.org/examples/api/radar_chart.html
    """
    from matplotlib.path import Path
    from matplotlib.spines import Spine
    from matplotlib.projections.polar import PolarAxes
    from matplotlib.projections import register_projection

    def radar_factory(num_vars, frame='circle'):
        """Create a radar chart with `num_vars` axes.

        This function creates a RadarAxes projection and registers it.

        Parameters
        ----------
        num_vars : int
            Number of variables for radar chart.
        frame : {'circle' | 'polygon'}
            Shape of frame surrounding axes.

        """
        # calculate evenly-spaced axis angles
        theta = np.linspace(0, 2* np.pi, num_vars, endpoint=False)
        # rotate theta such that the first axis is at the top
        theta += np.pi / 2

        def draw_poly_patch(self):
            verts = unit_poly_verts(theta)
            return plt.Polygon(verts, closed=True, edgecolor='k')

        def draw_circle_patch(self):
            # unit circle centered on (0.5, 0.5)
            return plt.Circle((0.5, 0.5), 0.5)

        patch_dict = {'polygon': draw_poly_patch, 'circle': draw_circle_patch}
        if frame not in patch_dict:
            raise ValueError('unknown value for `frame`: %s' % frame)

        class RadarAxes(PolarAxes):
            name = 'radar'
            # use 1 line segment to connect specified points
            RESOLUTION = 1
            # define draw_frame method
            draw_patch = patch_dict[frame]

            def fill(self, *args, **kwargs):
                """Override fill so that line is closed by default"""
                closed = kwargs.pop('closed', True)
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
                if frame == 'circle':
                    return PolarAxes._gen_axes_spines(self)
                # The following is a hack to get the spines (i.e. the axes frame)
                # to draw correctly for a polygon frame.

                # spine_type must be 'left', 'right', 'top', 'bottom', or `circle`.
                spine_type = 'circle'
                verts = unit_poly_verts(theta)
                # close off polygon by repeating first vertex
                verts.append(verts[0])
                path = Path(verts)

                spine = Spine(self, spine_type, path)
                spine.set_transform(self.transAxes)
                return {'polar': spine}

        register_projection(RadarAxes)
        return theta

    def unit_poly_verts(theta):
        """Return vertices of polygon for subplot axes.

        This polygon is circumscribed by a unit circle centered at (0.5, 0.5)
        """
        x0, y0, r = [0.5] * 3
        verts = [(r*np.cos(t) + x0, r*np.sin(t) + y0) for t in theta]
        return verts

    theta = radar_factory(len(radial_vars), frame='polygon')

    n_row = n_col = int(np.ceil(np.sqrt(len(data[subplot_var].unique()))))

    fig, axes = plt.subplots(figsize=(9, 9), nrows=n_row, ncols=n_col,
                             subplot_kw=dict(projection='radar'))
    fig.subplots_adjust(wspace=0.25, hspace=0.20, top=0.85, bottom=0.05)

    colors = plt.get_cmap(cmap, len(data[group_var].unique()))
    # Plot the four cases from the example data on separate axes
    for ax, subplot_title in zip(axes.flatten(), data[subplot_var].unique()):

        case_data = data[data[subplot_var] == subplot_title]

        if scale_to_max:
            case_data.loc[:, radial_vars] = (
                case_data.loc[:, radial_vars] / case_data.loc[:, radial_vars].max())

        ax.set_rgrids([0.2, 0.4, 0.6, 0.8])
        ax.set_title(subplot_title, weight='bold', size='medium', position=(0.5, 1.1),
                     horizontalalignment='center', verticalalignment='center')
        for group in case_data[group_var].unique():
            d = case_data.loc[case_data[group_var] == group, radial_vars].squeeze()
            color = colors(np.log2(1 + group) / np.log2(1 + 256))
            ax.plot(theta, d, label=group, color=color)
            ax.fill(theta, d, label=group, alpha=0.25, facecolor=color)
            # ax.set_ylim(0.7, 0.85)
        ax.set_varlabels(radial_vars)

        # add legend relative to top-left plot
        ax.legend(loc=(0.9, .95),
                  labelspacing=0.1, fontsize='small')

    return fig


def add_colorbar_to_axis(collection, label=None, position="right", size="5%", pad=0.05):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(collection.axes)
    cax = divider.append_axes(position, size=size, pad=pad)
    plt.colorbar(mappable=collection, cax=cax, label=label, alpha=1)


def clustermap_fix_label_orientation(grid, fontsize="xx-small", **kwargs):
    grid.ax_heatmap.set_xticklabels(
        grid.ax_heatmap.get_xticklabels(), rotation=90, fontsize=fontsize, **kwargs)
    grid.ax_heatmap.set_yticklabels(
        grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize=fontsize, **kwargs)


def clustermap_rasterize_heatmap(grid):
    import matplotlib
    q = [x for x in grid.ax_heatmap.get_children()
         if isinstance(x, matplotlib.collections.QuadMesh)][0]
    q.set_rasterized(True)


def savefig(fig, file_name, **kwargs):
    if isinstance(fig, sns.axisgrid.Grid):
        fig = fig.fig
    default_kwargs = _CONFIG['graphics']['settings']['figure_saving']
    default_kwargs.update(kwargs)
    fig.savefig(file_name, **default_kwargs)
    if _CONFIG['preferences']['graphics']['close_saved_figures']:
        plt.close(fig)


def plot_projection(
        df,
        color_dataframe,
        dims,
        output_file,
        attributes_to_plot,
        plot_max_dims=8,
        rasterized=False, plot_group_centroids=True,
        axis_ticklabels=True, axis_lines=True,
        legends=False, always_legend=False):
    """
    Plot a low dimentionality projection of samples.

    :param df: Dataframe with sample projections
    :type df: pandas.DataFrame
    :param color_dataframe: Dataframe of RGB tuples for sample i in attribute j.
    :type color_dataframe: pandas.DataFrame
    :param dims: Number of dimentions to plot
    :type dims: int
    :param output_file: Path to figure output file
    :type output_file: str
    :param attributes_to_plot: List of levels in df.index to plot
    :type attributes_to_plot: list
    :param plot_max_dims: Maximum number of sample attributes to plot for each factor in plot legend.
                          Defaults to 20.
    :type plot_max_dims: int, optional
    :param plot_max_dims: Maximum number of principal components to plot. This only affects plotting.
                         All dims will be calculated.
                         Defaults to 8.
    :type plot_max_dims: number, optional
    :param plot_group_centroids: Whether centroids of each sample group should be plotted alongside
                                 samples. Will be square shaped.
                                 Defaults to True.
    :type plot_group_centroids: bool, optional
    :param axis_ticklabels: Whether MDS and PCA axis ticks and ticklabels should be plotted.
                            Defaults to False.
    :type axis_ticklabels: bool, optional
    :param axis_lines: Whether (0, 0) dashed lines should be plotted in MDS and PCA.
                       Defaults to True.
    :type axis_lines: bool, optional
    :param legends: Whether legends for group colours should be plotted in MDS and PCA.
                    Defaults to False.
    :type legends: bool, optional
    :param always_legend: Whether legends for group colours should be plotted in every figure
                          panel in MDS and PCA.
                          If False, will plot just on first/last figure panel.
                          Defaults to False.
    :type always_legend: bool, optional
    """
    from collections import OrderedDict

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
            for j, sample in enumerate(df.index):
                sample = pd.Series(sample, index=df.index.names)
                label = df.index.get_level_values(attr)[j]
                axis[pc, i].scatter(
                    df.loc[sample['sample_name'], pc],
                    df.loc[sample['sample_name'], pc + 1],
                    s=30, color=color_dataframe.loc[attr, sample['sample_name']],
                    alpha=0.75, label=label, rasterized=rasterized)

            # Plot groups
            if plot_group_centroids:
                df2 = df.groupby(attr).mean()
                # get the color of each attribute group
                cd = color_dataframe.loc[attr, :]
                cd.name = None
                cd.index = df.index.get_level_values(attr)
                cd = cd.apply(lambda x: tuple(x) if isinstance(x, list) else x)  # fix for deduplicating lists
                cd = cd.reset_index().drop_duplicates().set_index(attr)
                for j, group in enumerate(df2.index):
                    axis[pc, i].scatter(
                        df2.loc[group, pc],
                        df2.loc[group, pc + 1],
                        marker="s", s=50, color=cd.loc[group].squeeze(),
                        alpha=0.95, label=group, rasterized=rasterized)
                    axis[pc, i].text(
                        df2.loc[group, pc],
                        df2.loc[group, pc + 1], group,
                        color=cd.loc[group].squeeze(), alpha=0.95)

            # Graphics
            axis[pc, i].set_title(attr)
            axis[pc, i].set_xlabel("PC {}".format(pc + 1))
            axis[pc, i].set_ylabel("PC {}".format(pc + 2))
            if not axis_ticklabels:
                axis[pc, i].set_xticklabels(axis[pc, i].get_xticklabels(), visible=False)
                axis[pc, i].set_yticklabels(axis[pc, i].get_yticklabels(), visible=False)
            if axis_lines:
                axis[pc, i].axhline(0, linestyle="--", color="black", alpha=0.3)
                axis[pc, i].axvline(0, linestyle="--", color="black", alpha=0.3)

            if legends:
                # Unique legend labels
                handles, labels = axis[pc, i].get_legend_handles_labels()
                by_label = OrderedDict(zip(labels, handles))
                if any([isinstance(c, str) for c in by_label.keys()]) and len(by_label) <= plot_max_dims:
                    # if not any([re.match("^\d", c) for c in by_label.keys()]):
                    if always_legend:
                        axis[pc, i].legend(by_label.values(), by_label.keys())
                    else:
                        if pc == (dims - 1):
                            axis[pc, i].legend(
                                by_label.values(), by_label.keys())
    savefig(fig, output_file)


def plot_region_structure_results(
        enr,
        output_dir, output_prefix="region_type_enrichment",
        across_attribute=None,
        pvalue=0.05, top_n=5):
    """
    Plot results of ATACSeqAnalysis.region_context_enrichment.

    Attributes:

    :param enr:
        Results of region_context_enrichment.

    :param str output_dir:
        Directory to save plots to.

    :param str, optional output_prefix:
        Prefix to use when saveing plots.
        Defaults to "region_type_enrichment"

    :param str,optional across_attribute:
        Column in enrichment matrix to plot results across e.g. "comparison_name"
        when results matrix contains the result of various comparisons.
        Defaults to None (not used).

    :param float,optional pvalue:
        Value at which to plot a line marking the significant level.
        Defaults to 0.05.

    :param optional top_n:
        Number of features to label in volcano plot.
        Defaults to 5.

    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # replace inf p-values with
    x, y, z = 'log2_odds_ratio', '-log10(p-value)', "region_type"
    enr = enr.sort_values(x)

    if across_attribute is not None:
        level = enr.loc[:, across_attribute].unique()
        n = len(level)
    else:
        level = [""]
        across_attribute = 'intercept'
        enr.loc[:, across_attribute] = ""
        n = 1
    row = col = int(np.ceil(np.sqrt(n)))

    # barplots
    for i, var in enumerate([x, y]):
        fig, axis = plt.subplots(row, col, figsize=(col * 3, row * 3))
        if row == col == 1:
            axis = np.array([axis])
        axis = axis.flatten()
        for i, value in enumerate(level):
            axis[i].set_title(value)
            sns.barplot(
                data=enr.loc[enr[across_attribute] == value, :].reset_index(),
                x=var, y="region", hue=z,
                orient="horizontal", dodge=False, ax=axis[i])
            if (enr[var] < 0).any():
                m = enr[var].abs().max()
                m += m * 0.1
                axis[i].set_xlim((-m, m))
                axis[i].axvline(0, linestyle="--", color="grey")
        savefig(fig, os.path.join(output_dir, output_prefix + ".{}.barplot.svg".format(var)))

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
        sns.scatterplot(
            data=spec.reset_index(),
            x=x, y=y, hue=z, ax=axis[i])
        # label top points
        for s in spec.sort_values(y).tail(top_n).index:
            axis[i].text(spec.loc[s, x], spec.loc[s, y], s=s)
        # equal x axis
        ll = spec.loc[:, x].abs().max()
        ll += ll * 0.1
        axis[i].set_xlim((-ll, ll))
    savefig(fig, os.path.join(output_dir, output_prefix + ".volcano_plot.top_{}_labeled.svg".format(top_n)))
