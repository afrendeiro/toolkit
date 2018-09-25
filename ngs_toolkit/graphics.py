#!/usr/bin/env python


import matplotlib.pyplot as plt
import pandas as pd


def barmap(x, figsize=None, square=False, row_colors=None, z_score=None, ylims=None):
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
        axis[i].axhline(0, linestyle="--", color="black", linewidth=0.5)
        if i != y_size - 1:
            axis[i].set_xticklabels(axis[i].get_xticklabels(), visible=False)
        axis[i].margins(0.0)
        if ylims is not None:
            axis[i].set_ylim(ylims)

    axis[-1].set_xticks(range(x_size))
    axis[-1].set_xticklabels(x.columns, rotation='vertical', ha="center", va="top", fontsize=figsize[1])

    for i, ax in enumerate(axis):
        m = max(ax.get_yticks())
        # ax.set_yticks([m])
        # ax.set_yticklabels([str(m)], rotation='horizontal', ha="center", va="center", visible=True, fontsize=figsize[1])
        ax.set_ylabel(str(x.index[i]), rotation='horizontal', ha="right", va="center", visible=True, fontsize=figsize[1])

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

    import numpy as np
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

            def _close_line(self, line):
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
            case_data.loc[:, radial_vars] = case_data.loc[:, radial_vars] / case_data.loc[:, radial_vars].max()

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
        legend = ax.legend(loc=(0.9, .95),
                        labelspacing=0.1, fontsize='small')

    return fig


def add_colorbar_to_axis(points, axis):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(axis)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(mappable=cs, cax=cax)
