#!/usr/bin/env python

import os

import numpy as np
import pandas as pd
from tqdm import tqdm


def pickle_me(function):
    """
    Decorator for some methods of Analysis class.
    Important: Pickled function cannot have positional arguments!
    """
    def wrapper(obj, timestamp=False, *args, **kwargs):
        import pickle
        function(obj, *args, **kwargs)
        if timestamp:
            import time
            import datetime
            ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d-%H%M%S')
            p = obj.pickle_file.replace(".pickle", ".{}.pickle".format(ts))
        else:
            p = obj.pickle_file
        pickle.dump(obj, open(p, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    return wrapper


class Analysis(object):
    """
    Generic class holding functions and data from a typical NGS analysis.

    Other modules implement classes inheriting from this that in general contain data type-specific
    functions (e.g. `ngs_toolkit.atacseq.ATACSeqAnalysis` has a `get_consensus_sites` function to genrate a peak consensus map).

    Objects of this type can be used to store data (e.g. dataframes), variables (e.g. paths to files or configurations) and are
    easily serializable (saved to a file as an object) for rapid loading and cross-environment portability.
    See the `ngs_toolkit.general.Analysis.to_pickle`, `ngs_toolkit.general.Analysis.from_pickle` and `ngs_toolkit.general.Analysis.update` functions for this.

    :param name: Name of the analysis. Defaults to "analysis".
    :type name: str, optional
    :param samples: List of `looper.models.Sample` objects that this analysi is tied to. Defaults to None.
    :type samples: list, optional
    :param prj: A `looper.models.Project` object that this analysis is tied to. Defaults to None.
    :type prj: looper.models.Project, optional
    :param data_dir: Directory containing processed data that will be input to the analysis.
                     This is in principle not required. Defaults to "data".
    :type data_dir: str, optional
    :param results_dir: Directory to contain outputs produced by the analysis. Defaults to "results".
    :type results_dir: str, optional
    :param pickle_file: A pickle file to serialize the object. Defaults to "`name`.pickle".
    :type pickle_file: str, optional
    :param from_pickle: Whether the analysis should be loaded from an existing serialized analysis object in `pickle_file`. Defaults to False.
    :type from_pickle: bool, optional

    Additional keyword arguments will simply be stored as object attributes.
    """

    def __init__(
            self,
            name="analysis",
            samples=None,
            prj=None,
            data_dir="data",
            results_dir="results",
            pickle_file=None,
            from_pickle=False,
            **kwargs):
        # parse kwargs with default
        self.name = name
        self.data_dir = data_dir
        self.results_dir = results_dir
        self.samples = samples
        self.prj = prj
        if pickle_file is None:
            pickle_file = os.path.join(results_dir, "{}.pickle".format(self.name))
        self.pickle_file = pickle_file

        # If no samples are provided but the prj object has samples,
        # set them as well
        if self.samples is None:
            if self.prj is not None:
                if hasattr(self.prj, "samples"):
                    self.samples = prj.samples

        for directory in [self.data_dir, self.results_dir]:
            if not os.path.exists(directory):
                os.makedirs(directory)

        # parse remaining kwargs
        self.__dict__.update(kwargs)

        # reload itself if required
        if from_pickle:
            self.update()

    @pickle_me
    def to_pickle(self):
        """
        Serialize object (i.e. save to disk).
        """
        pass

    def from_pickle(self, pickle_file=None):
        """
        Load object from pickle file.

        :param pickle_file: Pickle file to load. By default this is the object's attribute `pickle_file`.
        :type pickle_file: str, optional
        """
        import cPickle as pickle
        if pickle_file is None:
            pickle_file = self.pickle_file
        return pickle.load(open(pickle_file, 'rb'))

    def update(self, pickle_file=None):
        """
        Update all of the object's attributes with the attributes from a serialized object (i.e. object stored in a file) object.

        :param pickle_file: Pickle file to load. By default this is the object's attribute `pickle_file`.
        :type pickle_file: str, optional
        """
        self.__dict__.update(self.from_pickle(pickle_file=pickle_file).__dict__)

    def get_level_colors(
            self, index=None, matrix="accessibility", levels=None,
            pallete="Paired", cmap="RdBu_r", nan_color=(0.662745, 0.662745, 0.662745, 0.5),
            # TODO: implement dataframe return
            as_dataframe=False):
        """
        Get tuples of floats representing a colour for a sample in a given variable in a dataframe's index
        (particularly useful with MultiIndex dataframes).

        By default, will act on the columns and its levels of an `matrix` dataframe of self. Other `index` and `levels` can
        be passed for customization.

        Will try to guess if each variable is categorical or numerical and return either colours from a colour `pallete`
        or a `cmap`, respectively with null values set to `nan_color` (a 4-value tuple with floats).

        :param index: Pandas Index to use. If not provided (default == None), this will be the column Index of the provided `matrix`.
        :type index: pandas.Index, optional
        :param matrix: Name of analysis attribute containing a dataframe with pandas.MultiIndex columns to use.
                       Defaults to "accessibility".
        :type matrix: str, optional
        :param levels: Levels of multiindex to restrict to. Defaults to all in index.
        :type levels: list, optional
        :param pallete: Name of matplotlib color palete to use with categorical levels.
                        See matplotlib.org/examples/color/colormaps_reference.html. Defaults to "Paired".
        :type pallete: str, optional
        :param cmap: Name of matplotlib color palete to use with numerical levels.
                     See matplotlib.org/examples/color/colormaps_reference.html. Defaults to "RdBu_r".
        :type cmap: str, optional
        :param nan_color: Color for missing (i.e. NA) values. Defaults to (0.662745, 0.662745, 0.662745, 0.5) == grey.
        :type nan_color: tuple, optional
        :param as_dataframe: Whether a dataframe should be return. Defaults to False. Not implemented yet.
        :type as_dataframe: bool, optional
        :returns: List of list tuples (matrix) of shape (level, sample) with rgb values of each of the variable.
        :rtype: {list}
        """
        import numpy as np
        import matplotlib
        import matplotlib.pyplot as plt
        from collections import Counter

        if index is None:
            index = getattr(self, matrix).columns

        if levels is not None:
            drop = [l.name for l in index.levels if l.name not in levels]
            index = index.droplevel(drop)

        # Handle special case of single level
        if type(index) is pd.core.indexes.base.Index:
            index = pd.MultiIndex.from_arrays([index.values], names=[index.name])

        _cmap = plt.get_cmap(cmap)
        _pallete = plt.get_cmap(pallete)

        colors = list()
        for level in index.levels:
            # For empty levels (all values nan), return nan colour
            if len(level) == 0:
                colors.append([nan_color] * len(index))
                print(level.name, "NaN")
                continue
            # determine the type of data in each level
            # TODO: check this works in all cases
            most_common = Counter([type(x) for x in level]).most_common()[0][0]
            print(level.name, most_common)

            # Add either colors based on categories or numerical scale
            if most_common in [int, float, np.float32, np.float64, np.int32, np.int64]:
                values = index.get_level_values(level.name)
                # Create a range of either 0-100 if only positive values are found
                # or symmetrically from the maximum absolute value found
                if not any(values.dropna() < 0):
                    norm = matplotlib.colors.Normalize(vmin=0, vmax=100)
                else:
                    r = max(abs(values.min()), abs(values.max()))
                    norm = matplotlib.colors.Normalize(vmin=-r, vmax=r)

                col = _cmap(norm(values))
                # replace color for nan cases
                col[np.where(index.get_level_values(
                    level.name).to_series().isnull().tolist())] = nan_color
                colors.append(col.tolist())
            else:
                n = len(set(index.get_level_values(level.name)))
                # get n equidistant colors
                p = [_pallete(1. * i / n) for i in range(n)]
                color_dict = dict(zip(list(set(index.get_level_values(level.name))), p))
                # color for nan cases
                color_dict[np.nan] = nan_color
                col = [color_dict[x] for x in index.get_level_values(level.name)]
                colors.append(col)

        return colors


def count_reads_in_intervals(bam, intervals):
    """
    Count total number of reads in a iterable holding strings
    representing genomic intervals of the type chrom:start-end.

    :param bam: BAM file.
    :type bam: str
    :param intervals: List of strings with genomic coordinates in format "chrom:start-end".
    :type intervals: list
    :returns: Dict of read counts for each interval.
    :rtype: dict
    """
    import pysam
    counts = dict()

    bam = pysam.AlignmentFile(bam, mode='rb')

    chroms = ["chr" + str(x) for x in range(1, 23)] + ["chrX", "chrY"]

    for interval in intervals:
        if interval.split(":")[0] not in chroms:
            continue
        counts[interval] = bam.count(region=interval.split("|")[0])
    bam.close()

    return counts


def normalize_quantiles_r(array):
    """
    Quantile normalization with a R implementation.
    Requires the "rpy2" library and the R library "preprocessCore".

    Install in following way:
    # install package
    # R
    # source('http://bioconductor.org/biocLite.R')
    # biocLite('preprocessCore')

    :param array: Numeric array to normalize.
    :type array: numpy.array
    :returns: Normalized numeric array.
    :rtype: numpy.array
    """
    import numpy as np
    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()

    robjects.r('require("preprocessCore")')
    normq = robjects.r('normalize.quantiles')
    return np.array(normq(array))


def normalize_quantiles_p(df_input):
    """
    Quantile normalization with a ure Python implementation.
    Code from https://github.com/ShawnLYU/Quantile_Normalize.
    Not fully reviewed, use at own risk!

    :param df_input: Dataframe to normalize.
    :type df_input: pandas.DataFrame
    :returns: Normalized numeric array.
    :rtype: numpy.array
    """
    df = df_input.copy()
    # compute rank
    dic = {}
    for col in df:
        dic.update({col: sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis=1).tolist()
    # sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df


def unsupervised_analysis(
        analysis,
        data_type="ATAC-seq",
        quant_matrix=None,
        samples=None,
        attributes_to_plot=["cell_type"],
        plot_prefix=None,
        test_pc_association=True,
        display_corr_values=False,
        plot_max_attr=20,
        plot_max_pcs=8,
        plot_group_centroids=True,
        axis_ticklabels=False,
        axis_lines=True,
        legends=False,
        always_legend=False,
        prettier_sample_names=True,
        rasterized=False,
        dpi=300,
        output_dir="{results_dir}/unsupervised"):
    """
    Apply unsupervised clustering and dimensionality reduction methods (MDS, PCA) on numeric matrix.
    Colours and labels samples by their attributes as given in `attributes_to_plot`.

    For PCA analysis, if `test_pc_association` is `True`, will compute association of PCs with sample
    attributes given in `attributes_to_plot`. For numeric attributes, the Pearson correlation will be computed
    and for categoriacal, a pairwise Kruskal-Wallis H-test (ANOVA).

    :param analysis: Analysis object to perform analysis for.
    :type analysis: ngs_toolkit.general.Analysis
    :param data_type: Data type. One of "ATAC-seq" or "RNA-seq". Defaults to "ATAC-seq".
    :type data_type: str, optional
    :param quant_matrix: Name of analysis attribute contatining the numeric dataframe to perform analysis on.
                         Defaults to "accessibility" if data_type is ATAC-seq and "expression_annotated" if data_type is RNA-seq.
                         This matrix should have a pandas.MultiIndex as column index.
    :type quant_matrix: str, optional
    :param samples: List of sample objects to restrict analysis to. Defaults to None.
    :type samples: list, optional
    :param attributes_to_plot: List of attributes shared between sample groups should be plotted. Defaults to ["cell_type"].
    :type attributes_to_plot: list, optional
    :param plot_prefix: Prefix for output files.
                        Defaults to "all_sites" if data_type is ATAC-seq and "all_genes" if data_type is RNA-seq.
    :type plot_prefix: str, optional
    :param test_pc_association: Whether a test of association of principal components and variables in `attributes_to_plot` should be conducted. Defaults to True.
    :type test_pc_association: bool, optional
    :param display_corr_values: Whether values in heatmap of sample correlations should be displayed overlaid on top of colours. Defaults to False.
    :type display_corr_values: bool, optional
    :param plot_max_attr: Maximum number of sample attributes to plot for each factor in plot legend. Defaults to 20.
    :type plot_max_attr: int, optional
    :param plot_max_pcs: Maximum number of principal components to plot. This only affects plotting. All PCs will be calculated. Defaults to 8.
    :type plot_max_pcs: number, optional
    :param plot_group_centroids: Whether centroids of each sample group should be plotted alongside samples. Will be square shaped. Defaults to True.
    :type plot_group_centroids: bool, optional
    :param axis_ticklabels: Whether MDS and PCA axis ticks and ticklabels should be plotted. Defaults to False.
    :type axis_ticklabels: bool, optional
    :param axis_lines: Whether (0, 0) dashed lines should be plotted in MDS and PCA. Defaults to True.
    :type axis_lines: bool, optional
    :param legends: Whether legends for group colours should be plotted in MDS and PCA. Defaults to False.
    :type legends: bool, optional
    :param always_legend: Whether legends for group colours should be plotted in every figure panel in MDS and PCA.
                          If False, will plot just on first/last figure panel. Defaults to False.
    :type always_legend: bool, optional
    :param prettier_sample_names: Whether it should attempt to prettify sample names by removing the data type from plots. Defaults to True.
    :type prettier_sample_names: bool, optional
    :param rasterized: Whether elements with many objects should be rasterized. Defaults to False.
    :type rasterized: bool, optional
    :param dpi: Definition of rasterized image in dots per inch. Defaults to 300dpi.
    :type dpi: number, optional
    :param output_dir: Directory for generated files and plots. Defaults to "{results_dir}/unsupervised".
    :type output_dir: str, optional

    :returns: None
    """
    from sklearn.decomposition import PCA
    from sklearn.manifold import MDS
    from collections import OrderedDict
    import re
    import itertools
    from scipy.stats import kruskal
    from scipy.stats import pearsonr
    import matplotlib.pyplot as plt
    import seaborn as sns

    if data_type == "ATAC-seq":
        if plot_prefix is None:
            plot_prefix = "all_sites"
        if quant_matrix is None:
            quant_matrix = "accessibility"
    elif data_type == "RNA-seq":
        if plot_prefix is None:
            plot_prefix = "all_genes"
        if quant_matrix is None:
            quant_matrix = "expression_annotated"
    else:
        raise ValueError("Data types can only be 'ATAC-seq' or 'RNA-seq'.")

    if "{results_dir}" in output_dir:
        output_dir = output_dir.format(results_dir=analysis.results_dir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    matrix = getattr(analysis, quant_matrix)

    if type(matrix.columns) is not pd.core.indexes.multi.MultiIndex:
        raise TypeError("Provided quantification matrix must have columns with MultiIndex.")

    if samples is None:
        samples = [
            s for s in analysis.samples if s.name in matrix.columns.get_level_values("sample_name")]

    # remove attributes with all NaNs
    attributes_to_plot = [attr for attr in attributes_to_plot if not pd.isnull(
        matrix.columns.get_level_values(attr)).all()]

    # This will always be a matrix for all samples
    color_dataframe = pd.DataFrame(
        analysis.get_level_colors(index=matrix.columns, levels=attributes_to_plot),
        index=attributes_to_plot, columns=matrix.columns.get_level_values("sample_name"))
    # will be filtered now by the requested samples if needed
    color_dataframe = color_dataframe[[s.name for s in samples]]

    # All regions, matching samples (provided samples in matrix)
    X = matrix.loc[:, matrix.columns.get_level_values(
        "sample_name").isin([s.name for s in samples])]

    # TODO: Re-implement to accomodate multiindex
    # if prettier_sample_names:
    #     X.columns = (
    #         color_dataframe.columns
    #         .str.replace("ATAC-seq_", "")
    #         .str.replace("RNA-seq_", "")
    #         .str.replace("ChIP-seq_", ""))

    # Pairwise correlations
    g = sns.clustermap(
        # yticklabels=sample_display_names,
        X.astype(float).corr(), xticklabels=False, yticklabels=True, annot=display_corr_values,
        cmap="Spectral_r", figsize=(0.2 * X.shape[1], 0.2 * X.shape[1]), cbar_kws={"label": "Pearson correlation"}, row_colors=color_dataframe.values.tolist())
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize='xx-small')
    g.ax_heatmap.set_xlabel(None, visible=False)
    g.ax_heatmap.set_ylabel(None, visible=False)
    g.fig.savefig(os.path.join(output_dir, "{}.{}.corr.clustermap.svg".format(
        analysis.name, plot_prefix)), bbox_inches='tight', dpi=dpi)

    # MDS
    mds = MDS(n_jobs=-1)
    x_new = mds.fit_transform(X.T)
    # transform again
    xx = pd.DataFrame(x_new, index=X.columns, columns=list(range(x_new.shape[1])))
    # xx = x.apply(lambda j: (j - j.mean()) / j.std(), axis=0)

    fig, axis = plt.subplots(1, len(attributes_to_plot),
                             figsize=(4 * len(attributes_to_plot), 4 * 1))
    axis = axis.flatten()
    for i, attr in enumerate(attributes_to_plot):
        for j, sample in enumerate(xx.index):
            sample = pd.Series(sample, index=X.columns.names)
            try:
                label = getattr(sample, attributes_to_plot[i])
            except AttributeError:
                label = np.nan
            axis[i].scatter(
                xx.loc[sample['sample_name'], 0],
                xx.loc[sample['sample_name'], 1],
                s=50, color=color_dataframe.loc[attr, sample['sample_name']], alpha=0.75, label=label, rasterized=rasterized)

        # Plot groups
        if plot_group_centroids:
            xx2 = xx.groupby(attr).mean()
            # get the color of each attribute group
            cd = color_dataframe.loc[attr]
            cd.name = None
            cd.index = X.columns.get_level_values(attr)
            cd = cd.reset_index().drop_duplicates().set_index(attr)
            for j, group in enumerate(xx2.index):
                axis[i].scatter(
                    xx2.loc[group, 0],
                    xx2.loc[group, 1],
                    marker="s", s=50, color=cd.loc[group].squeeze(), alpha=0.95, label=group, rasterized=rasterized)
                axis[i].text(
                    xx2.loc[group, 0],
                    xx2.loc[group, 1], group,
                    color=cd.loc[group].squeeze(), alpha=0.95)

        # Graphics
        axis[i].set_title(attributes_to_plot[i])
        axis[i].set_xlabel("MDS 1")
        axis[i].set_ylabel("MDS 2")
        if not axis_ticklabels:
            axis[i].set_xticklabels(axis[i].get_xticklabels(), visible=False)
            axis[i].set_yticklabels(axis[i].get_yticklabels(), visible=False)
        if axis_lines:
            axis[i].axhline(0, linestyle="--", color="black", alpha=0.3)
            axis[i].axvline(0, linestyle="--", color="black", alpha=0.3)

        if legends:
            # Unique legend labels
            handles, labels = axis[i].get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))
            if any([type(c) in [str, unicode] for c in by_label.keys()]) and len(by_label) <= 20:
                # if not any([re.match("^\d", c) for c in by_label.keys()]):
                axis[i].legend(by_label.values(), by_label.keys())
    fig.savefig(os.path.join(output_dir, "{}.{}.mds.svg".format(
        analysis.name, plot_prefix)), bbox_inches='tight', dpi=dpi)

    # PCA
    # TODO: Restrict PCA analysis to plot_max_pcs.
    pca = PCA()
    x_new = pca.fit_transform(X.T)
    # transform again
    xx = pd.DataFrame(x_new, index=X.columns, columns=list(range(x_new.shape[0])))

    # plot % explained variance per PC
    fig, axis = plt.subplots(1)
    axis.plot(
        range(1, len(pca.explained_variance_) + 1),  # all PCs
        (pca.explained_variance_ / pca.explained_variance_.sum()) * 100, 'o-')  # % of total variance
    axis.axvline(len(attributes_to_plot), linestyle='--')
    axis.set_xlabel("PC")
    axis.set_ylabel("% variance")
    sns.despine(fig)
    fig.savefig(os.path.join(output_dir, "{}.{}.pca.explained_variance.svg".format(
        analysis.name, plot_prefix)), bbox_inches='tight', dpi=dpi)

    # Write % variance expained to disk
    pd.Series((pca.explained_variance_ / pca.explained_variance_.sum()) * 100, name="PC").to_csv(
        os.path.join(output_dir, "{}.{}.pca.explained_variance.csv".format(
            analysis.name, plot_prefix)))

    # plot pca
    pcs = min(xx.shape[0] - 1, plot_max_pcs)
    fig, axis = plt.subplots(pcs, len(attributes_to_plot), figsize=(
        4 * len(attributes_to_plot), 4 * pcs))
    for pc in range(pcs):
        for i, attr in enumerate(attributes_to_plot):
            for j, sample in enumerate(xx.index):
                sample = pd.Series(sample, index=X.columns.names)
                try:
                    label = getattr(samples[j], attr)
                except AttributeError:
                    label = np.nan
                axis[pc, i].scatter(
                    xx.loc[sample['sample_name'], :].loc[:, pc],
                    xx.loc[sample['sample_name'], :].loc[:, pc + 1],
                    s=30, color=color_dataframe.loc[attr, sample['sample_name']], alpha=0.75, label=label, rasterized=rasterized)

            # Plot groups
            if plot_group_centroids:
                xx2 = xx.groupby(attr).mean()
                # get the color of each attribute group
                cd = color_dataframe.loc[attr]
                cd.name = None
                cd.index = X.columns.get_level_values(attr)
                cd = cd.reset_index().drop_duplicates().set_index(attr)
                for j, group in enumerate(xx2.index):
                    axis[pc, i].scatter(
                        xx2.loc[group, pc],
                        xx2.loc[group, pc + 1],
                        marker="s", s=50, color=cd.loc[group].squeeze(), alpha=0.95, label=group, rasterized=rasterized)
                    axis[pc, i].text(
                        xx2.loc[group, pc],
                        xx2.loc[group, pc + 1], group,
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
                if any([type(c) in [str, unicode] for c in by_label.keys()]) and len(by_label) <= plot_max_attr:
                    # if not any([re.match("^\d", c) for c in by_label.keys()]):
                    if always_legend:
                        axis[pc, i].legend(by_label.values(), by_label.keys())
                    else:
                        if pc == pcs - 1:
                            axis[pc, i].legend(
                                by_label.values(), by_label.keys())
    fig.savefig(os.path.join(output_dir, "{}.{}.pca.svg".format(
        analysis.name, plot_prefix)), bbox_inches="tight")

    # Test association of PCs with attributes
    if not test_pc_association:
        return

    associations = list()
    for pc in range(pcs):
        for attr in attributes_to_plot:
            print("PC {}; Attribute {}.".format(pc + 1, attr))

            # Get all values of samples for this attr
            groups = xx.index.get_level_values(attr)

            # Determine if attr is categorical or continuous
            if all([type(i) in [str, bool] for i in groups]) or len(groups) == 2:
                variable_type = "categorical"
            elif all([type(i) in [int, float, np.int64, np.float64] for i in groups]):
                variable_type = "numerical"
            else:
                print("attr %s cannot be tested." % attr)
                associations.append([pc + 1, attr, variable_type, np.nan, np.nan, np.nan])
                continue

            if variable_type == "categorical":
                # It categorical, test pairwise combinations of attributes
                for group1, group2 in itertools.combinations(groups, 2):
                    g1_mask = xx.index.get_level_values(attr) == group1
                    g2_mask = xx.index.get_level_values(attr) == group2

                    g1_values = xx.loc[g1_mask, pc]
                    g2_values = xx.loc[g2_mask, pc]

                    # Test ANOVA (or Kruskal-Wallis H-test)
                    p = kruskal(g1_values, g2_values)[1]

                    # Append
                    associations.append([pc + 1, attr, variable_type, group1, group2, p])

            elif variable_type == "numerical":
                # It numerical, calculate pearson correlation
                pc_values = xx.loc[:, pc]
                trait_values = xx.index.get_level_values(attr)
                p = pearsonr(pc_values, trait_values)[1]

                associations.append([pc + 1, attr, variable_type, np.nan, np.nan, p])

    associations = pd.DataFrame(associations, columns=[
                                "pc", "attribute", "variable_type", "group_1", "group_2", "p_value"])

    # write
    associations.to_csv(os.path.join(output_dir, "{}.{}.pca.variable_principle_components_association.csv".format(
        analysis.name, plot_prefix)), index=False)

    # Plot
    # associations[associations['p_value'] < 0.05].drop(['group_1', 'group_2'], axis=1).drop_duplicates()
    # associations.drop(['group_1', 'group_2'], axis=1).drop_duplicates().pivot(index="pc", columns="attribute", values="p_value")
    pivot = associations.groupby(["pc", "attribute"]).min()['p_value'].reset_index().pivot(
        index="pc", columns="attribute", values="p_value").dropna(axis=1)

    # heatmap of -log p-values
    g = sns.clustermap(-np.log10(pivot), row_cluster=False, annot=True,
                       cbar_kws={"label": "-log10(p_value) of association"}, square=True, rasterized=rasterized)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
    g.fig.savefig(os.path.join(output_dir, "{}.{}.pca.variable_principle_components_association.svg".format(
        analysis.name, plot_prefix)), bbox_inches="tight", dpi=dpi)

    # heatmap of masked significant
    g = sns.clustermap((pivot < 0.05).astype(int), row_cluster=False, cbar_kws={
                       "label": "significant association"}, square=True, rasterized=rasterized)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
    g.fig.savefig(os.path.join(output_dir, "{}.{}.pca.variable_principle_components_association.masked.svg".format(
        analysis.name, plot_prefix)), bbox_inches="tight", dpi=dpi)


def deseq_analysis(
        count_matrix, experiment_matrix, comparison_table, formula,
        output_dir, output_prefix,
        overwrite=True, alpha=0.05):
    """
    Perform differential comparison analysis with DESeq2.

    .. note::
        Do not include hyphens ("-") in any of the samples or groups names!
        R freaks out with this.

    # TODO: fix hyphens in names issue

    :param count_matrix: Data frame of shape (samples, variables) with raw read counts.
    :type count_matrix: pandas.DataFrame
    :param experiment_matrix: Data frame with columns "sample_name" and any other variables used in the `formula`.
    :type experiment_matrix: pandas.DataFrame
    :param comparison_table: Data frame with columns "comparison_name", "sample_group" and sample_name".
    :type comparison_table: pandas.DataFrame
    :param formula: Formula to test in R/patsy notation. Usually something like "~ batch + group".
    :type formula: str
    :param output_dir: Output directory for produced files.
    :type output_dir: str
    :param output_prefix: Prefix to add to produced files.
    :type output_prefix: str
    :param overwrite: Whether files existing should be overwritten. Defaults to True.
    :type overwrite: bool, optional
    :param alpha: Significance level to reject null hypothesis.
                  This in practice has no effect as results for all features will be returned.
                  Defaults to 0.05.
    :type alpha: number, optional
    :returns: Data frame with results, statistics for each feature.
    :rtype: pandas.DataFrame
    """
    import pandas as pd
    from tqdm import tqdm
    from rpy2.robjects import numpy2ri, pandas2ri
    import rpy2.robjects as robjects
    numpy2ri.activate()
    pandas2ri.activate()

    def r2pandas_df(r_df):
        import numpy as np
        df = pd.DataFrame(np.asarray(r_df)).T
        df.columns = [str(x) for x in r_df.colnames]
        df.index = [str(x) for x in r_df.rownames]
        return df

    robjects.r('require("DESeq2")')
    _as_formula = robjects.r('as.formula')
    _DESeqDataSetFromMatrix = robjects.r('DESeqDataSetFromMatrix')
    _DESeq = robjects.r('DESeq')
    _results = robjects.r('results')
    _as_data_frame = robjects.r('as.data.frame')

    # order experiment and count matrices in same way
    experiment_matrix = experiment_matrix.set_index("sample_name").loc[count_matrix.columns, :]

    # save the matrices just in case
    count_matrix.to_csv(os.path.join(output_dir, output_prefix + ".count_matrix.tsv"), sep="\t")
    experiment_matrix.to_csv(os.path.join(output_dir, output_prefix +
                                          ".experiment_matrix.tsv"), sep="\t")
    comparison_table.to_csv(os.path.join(output_dir, output_prefix +
                                         ".comparison_table.tsv"), sep="\t")

    # Rename samples to avoid R errors with sample names containing symbols
    count_matrix.columns = ["S{}".format(i) for i in range(len(count_matrix.columns))]
    experiment_matrix.index = ["S{}".format(i) for i in range(len(experiment_matrix.index))]

    # Run DESeq analysis
    dds = _DESeqDataSetFromMatrix(
        countData=count_matrix.astype(int),
        colData=experiment_matrix,
        design=_as_formula(formula))
    dds = _DESeq(dds, parallel=True)
    # _save(dds, file=os.path.join(output_dir, output_prefix + ".deseq_dds_object.Rdata"))

    results = pd.DataFrame()
    for comp in tqdm(comparison_table["comparison_name"].drop_duplicates().sort_values()):
        out_file = os.path.join(output_dir, output_prefix + ".deseq_result.{}.csv".format(comp))
        if not overwrite and os.path.exists(out_file):
            continue
        print("Doing comparison '{}'".format(comp))
        a = comparison_table.loc[
            (comparison_table["comparison_name"] == comp) &
            (comparison_table["comparison_side"] >= 1),
            "sample_group"].drop_duplicates().squeeze()
        b = comparison_table.loc[
            (comparison_table["comparison_name"] == comp) &
            (comparison_table["comparison_side"] <= 0),
            "sample_group"].drop_duplicates().squeeze()

        contrast = ["sample_group" + a, "sample_group" + b]
        res = _as_data_frame(_results(dds, contrast=contrast, alpha=alpha,
                                      independentFiltering=False, parallel=True))

        # convert to pandas dataframe
        res2 = r2pandas_df(res)
        res2.loc[:, "comparison_name"] = comp

        # save
        res2.to_csv(out_file)
        # append
        results = results.append(res2.reset_index(), ignore_index=True)

    # save all
    results.to_csv(os.path.join(output_dir, output_prefix +
                                ".deseq_result.all_comparisons.csv"), index=False)

    # return
    return results


def deseq_results_to_bed_file(
        deseq_result_file, bed_file, sort=True, ascending=False, normalize=False,
        significant_only=False, alpha=0.05, abs_fold_change=1.):
    """
    Write BED file with fold changes from DESeq2 as score value.
    """

    df = pd.read_csv(deseq_result_file, index_col=0)

    assert "log2FoldChange" in df.columns.tolist()

    if sort is True:
        df = df.sort_values("log2FoldChange", ascending=ascending)

    if significant_only is True:
        df = df.loc[(df['padj'] < alpha) & (df['log2FoldChange'].abs() > abs_fold_change), :]

    # decompose index string (chrom:start-end) into columns
    df['chrom'] = map(lambda x: x[0], df.index.str.split(":"))
    r = pd.Series(map(lambda x: x[1], df.index.str.split(":")))
    df['start'] = map(lambda x: x[0], r.str.split("-"))
    df['end'] = map(lambda x: x[1], r.str.split("-"))
    df['name'] = df.index
    if normalize:
        from sklearn.preprocessing import MinMaxScaler
        MinMaxScaler(feature_range=(0, 1000)).fit_transform(df["log2FoldChange"])
    df['score'] = df["log2FoldChange"]

    df[["chrom", "start", "end", "name", "score"]].to_csv(
        bed_file, sep="\t", header=False, index=False)


def least_squares_fit(
        quant_matrix, design_matrix, test_model,
        null_model="~ 1", standardize_data=True,
        multiple_correction_method="fdr_bh"):
    """
    Fit a least squares model with only categorical predictors.
    Computes p-values by comparing the log likelihood ratio of the chosen model to a `null_model`.

    :param quant_matrix: A Data frame of shape (samples, variables).
    :type quant_matrix: pandas.DataFrame
    :param design_matrix: A Data frame of shape (samples, variables) with all the variables in `test_model`.
    :type design_matrix: pandas.DataFrame
    :param test_model: Model design to test in R/patsy notation.
    :type test_model: str
    :param null_model: Null model design in R/patsy notation. Defaults to "~ 1".
    :type null_model: str, optional
    :param standardize_data: Whether data should be standardized prior to fitting. Defaults to True.
    :type standardize_data: bool, optional
    :param multiple_correction_method: Method to use for multiple test correction.
                                       See statsmodels.sandbox.stats.multicomp.multipletests. Defaults to "fdr_bh".
    :type multiple_correction_method: str, optional
    :returns: Statistics of model fitting and comparison between models for each feature.
    :rtype: pandas.DataFrame
    """
    from sklearn.preprocessing import StandardScaler
    import patsy
    from scipy.linalg import lstsq
    from scipy import stats
    from statsmodels.sandbox.stats.multicomp import multipletests

    # # to test
    # quant_matrix = np.random.random(10000000).reshape(100, 100000)
    # P = np.concatenate([[0] * 50, [1] * 50])
    # Q = np.concatenate([[0] * 25, [1] * 25] + [[0] * 25, [1] * 25])
    # design_matrix = pd.DataFrame([P, Q], index=["P", "Q"]).T
    # quant_matrix = quant_matrix.T * (1 + (design_matrix.sum(axis=1) * 4).values)
    # quant_matrix = pd.DataFrame(quant_matrix.T)
    # test_model = "~ Q + P"
    # null_model = "~ Q"

    if standardize_data:
        norm = StandardScaler()
        quant_matrix = pd.DataFrame(
            norm.fit_transform(quant_matrix),
            index=quant_matrix.index, columns=quant_matrix.columns)

    A1 = patsy.dmatrix(test_model, design_matrix)
    betas1, residuals1, _, _ = lstsq(A1, quant_matrix)

    A0 = patsy.dmatrix(null_model, design_matrix)
    betas0, residuals0, _, _ = lstsq(A0, quant_matrix)

    results = pd.DataFrame(betas1.T, columns=A1.design_info.column_names,
                           index=quant_matrix.columns)

    # Calculate the log-likelihood ratios
    n = float(quant_matrix.shape[0])
    results['model_residuals'] = residuals1
    results['null_residuals'] = residuals0
    results['model_log_likelihood'] = (-n / 2.) * np.log(2 * np.pi) - \
        n / 2. * np.log(results['model_residuals'] / n) - n / 2.
    results['null_log_likelihood'] = (-n / 2.) * np.log(2 * np.pi) - \
        n / 2. * np.log(results['null_residuals'] / n) - n / 2.

    results['log_likelihood_ratio'] = results['model_log_likelihood'] - \
        results['null_log_likelihood']
    results['D_statistic'] = 2 * results['log_likelihood_ratio']
    results['p_value'] = stats.chi2.sf(
        results['log_likelihood_ratio'], df=betas1.shape[0] - betas0.shape[0])
    results['q_value'] = multipletests(results['p_value'], method=multiple_correction_method)[1]

    if not standardize_data:
        results["mean"] = quant_matrix.mean(axis=0)

    return results


def differential_from_bivariate_fit(
        comparison_table, matrix,
        output_dir, output_prefix,
        n_bins=250, multiple_correction_method="fdr_bh",
        plot=True, palette="colorblind", make_values_positive=False):
    """
    Perform differential analysis using a bivariate gaussian fit
    on the relationship between mean and fold-change for each comparison.

    :param pandas.DataFrame comparison_table: Dataframe with 'comparison_name', 'comparison_side' and 'sample_name', 'sample_group' columns.
    :param pandas.DataFrame matrix: Matrix of `n_features, n_samples` with normalized, log-transformed values to perform analysis on.
    :param str output_dir: Output directory
    :param str output_prefix: Prefix for outputs.
    :param int n_bins: Number of bins of mean values along which to standardize fold-changes.
    :param str multiple_correction_method: Multiple correction method from `statsmodels.sandbox.stats.multicomp.multipletests`.
    :param bool plot: Whether to generate plots.
    :param str palette: Color palette to use. This can be any matplotlib palette and is passed to `sns.color_palette`.
    """
    from scipy.stats import gaussian_kde
    from statsmodels.sandbox.stats.multicomp import multipletests
    import matplotlib.pyplot as plt

    comparisons = comparison_table['comparison_name'].drop_duplicates().sort_values()
    if plot:
        fig, axis = plt.subplots(
            2, len(comparisons),
            figsize=(4 * len(comparisons), 4 * 2), sharex=True, sharey='row')

    if make_values_positive:
        matrix = matrix + abs(matrix.min().min())

    results = pd.DataFrame()
    for i, comparison in enumerate(comparisons):
        print("Doing comparison '{}'".format(comparison))
        out_file = os.path.join(output_dir, output_prefix + ".fit_result.{}.csv".format(comparison))

        sa = comparison_table.loc[
            (comparison_table['comparison_name'] == comparison) &
            (comparison_table['comparison_side'] == 1),
            "sample_name"]
        ga = comparison_table.loc[
            (comparison_table['comparison_name'] == comparison) &
            (comparison_table['comparison_side'] == 1),
            "sample_group"].squeeze()
        sb = comparison_table.loc[
            (comparison_table['comparison_name'] == comparison) &
            (comparison_table['comparison_side'] == 0),
            "sample_name"]
        gb = comparison_table.loc[
            (comparison_table['comparison_name'] == comparison) &
            (comparison_table['comparison_side'] == 0),
            "sample_group"].squeeze()
        a = matrix.loc[:, sa].mean(axis=1)
        a.name = ga
        b = matrix.loc[:, sb].mean(axis=1)
        b.name = gb

        # assemble stats
        res = a.to_frame()
        res = res.join(b)
        res['global_mean'] = matrix.mean(axis=1)
        res['global_std'] = matrix.std(axis=1)
        res['comparison_mean'] = res.mean(axis=1)
        res['comparison_std'] = res.mean(axis=1)
        res['log2FoldChange'] = np.log2(res[ga] / res[gb])
        res['comparison_name'] = comparison

        # standardize fold change
        bounds = np.linspace(0, res['comparison_mean'].max(), n_bins)
        for (start, end) in zip(bounds[:-2], bounds[1: -1]):
            r = res.loc[(res['comparison_mean'] > start) &
                        (res['comparison_mean'] < end)].index
            v = res.loc[r, 'log2FoldChange']
            res.loc[r, 'norm_log2FoldChange'] = (v - np.nanmean(v)) / np.nanstd(v)

        # let's try a bivariate gaussian kernel
        # separately for positive and negative to avoid biases in center of mass
        kernel = gaussian_kde(res.loc[res['norm_log2FoldChange'] > 0, [
            "comparison_mean", "norm_log2FoldChange"]].T.values)
        res.loc[res['norm_log2FoldChange'] > 0, "density"] = kernel(
            res.loc[res['norm_log2FoldChange'] > 0, ["comparison_mean", "norm_log2FoldChange"]].T.values)
        kernel = gaussian_kde(res.loc[res['norm_log2FoldChange'] <= 0, [
            "comparison_mean", "norm_log2FoldChange"]].T.values)
        res.loc[res['norm_log2FoldChange'] <= 0, "density"] = kernel(
            res.loc[res['norm_log2FoldChange'] <= 0, ["comparison_mean", "norm_log2FoldChange"]].T.values)

        # Let's calculate something like an empirical p-value on the density
        res['pvalue'] = (res['density'] - res['density'].min()) / \
            (res['density'].max() - res['density'].min())
        res['padj'] = multipletests(res['pvalue'].fillna(1), method=multiple_correction_method)[1]

        res['direction'] = (res['norm_log2FoldChange'] >= 0).astype(int).replace(0, -1)
        res.to_csv(out_file)

        if plot:
            axis[0, i].scatter(res["comparison_mean"],
                               res["log2FoldChange"], alpha=0.2, s=5, color=sns.color_palette(palette)[0], rasterized=True)
            axis[0, i].axhline(0, color="black", linestyle="--")
            axis[1, i].scatter(res["comparison_mean"],
                               res["norm_log2FoldChange"], alpha=0.2, s=5, color=sns.color_palette(palette)[0], rasterized=True)
            diff = res.loc[(res['pvalue'] < 0.05) & (
                res['norm_log2FoldChange'].abs() >= 2), :].index
            axis[0, i].scatter(res.loc[diff, "comparison_mean"],
                               res.loc[diff, "log2FoldChange"], alpha=0.2, s=5, color=sns.color_palette(palette)[1], rasterized=True)
            axis[1, i].scatter(res.loc[diff, "comparison_mean"],
                               res.loc[diff, "norm_log2FoldChange"], alpha=0.2, s=5, color=sns.color_palette(palette)[1], rasterized=True)
            axis[1, i].axhline(0, color="black", linestyle="--")
            axis[0, i].set_title(comparison + "\n" + ga + " vs " + gb)
            axis[1, i].set_xlabel("Comparison mean")
            if i == 0:
                axis[0, i].set_ylabel("log2 fold-change")
                axis[1, i].set_ylabel("Norm(log2 fold-change)")

        results = results.append(res.reset_index(), ignore_index=True)

    # save figure
    fig.savefig(os.path.join(output_dir, output_prefix +
                             ".deseq_result.all_comparisons.scatter.svg"), dpi=300, bbox_inches="tight")

    # save all
    results = results.set_index("index")
    results.to_csv(os.path.join(output_dir, output_prefix +
                                ".deseq_result.all_comparisons.csv"), index=True)

    return results


# def independent_filtering(df, alpha=0.05, n_quantiles=100):
#     """
#     """
#     raise NotImplementedError
#     import numpy as np

#     req_columns = ["pvalue", "baseMean"]
#     assert all([x in df.columns for x in req_columns])

#     # compute quantiles accross mean and pvalue distributions
#     stats = pd.DataFrame()
#     p = (np.arange(n_quantiles) / float(n_quantiles)) * 100
#     p = np.append(p, 100.)
#     for start, end in zip(p, p[1:]):
#         m = np.log2(1 + df['baseMean'])
#         i = df.loc[
#             (m >= np.percentile(m, start)) &
#             (m <= np.percentile(m, end)), :].index
#         stats.loc[start, "n"] = i.shape[0]
#         stats.loc[start, "mean"] = df.loc[i, "baseMean"].mean()
#         stats.loc[start, "mean_p"] = df.loc[i, "pvalue"].mean()
#         stats.loc[start, "n_sig_p"] = (df.loc[i, "pvalue"] < alpha).sum()

#     # plot
#     fig, axis = plt.subplots(1, 2, figsize=(4 * 2, 4 * 1))
#     axis[0].scatter(stats.index, stats.loc[:, "n_sig_p"])
#     axis[1].scatter(stats.index, -np.log10(stats.loc[:, "mean_p"]))

#     # choose inflection point

#     return


# def fit_curve():
#     import numpy as np
#     import matplotlib.pyplot as plt
#     from scipy.optimize import curve_fit

#     def func(x, a, b, c):
#         return a * np.exp(-b * x) + c

#     plt.plot(stats.index, -np.log10(stats['mean_p']), 'b-', label='data')

#     popt, pcov = curve_fit(func, stats.index, -np.log10(stats['mean_p']))
#     plt.plot(stats.index, func(stats.index, *popt), 'r-', label='fit')

#     popt, pcov = curve_fit(func, stats.index, stats['mean_p'], bounds=(0, [3., 2., 1.]))
#     plt.plot(stats.index, func(stats.index, *popt), 'g--', label='fit-with-bounds')
#     plt.xlabel('x')
#     plt.ylabel('y')
#     plt.legend()
#     plt.show()


def differential_analysis(
        analysis,
        comparison_table,
        data_type="ATAC-seq",
        samples=None,
        covariates=None,
        output_dir="results/differential_analysis_{data_type}",
        output_prefix="differential_analysis",
        alpha=0.05,
        overwrite=True,
        distributed=False):
    """
    Perform differential regions/genes across samples that are associated with a certain trait.
    Currently the only implementation is with DESeq2.
    This implies the rpy2 library and the respective R library are installed.

    For other implementations of differential analysis see `ngs_toolkit.general.least_squares_fit`
    and `ngs_toolkit.general.differential_from_bivariate_fit`.

    :param analysis: A ngs_toolkit Analysis object.
    :type analysis: ngs_toolkit.general.Analysis
    :param comparison_table: A dataframe with 'comparison_name', 'comparison_side' and 'sample_name', 'sample_group' columns.
    :type comparison_table: pandas.DataFrame
    :param data_type: Type of data under analysis. One of "ATAC-seq" or "RNA-seq". Defaults to "ATAC-seq".
    :type data_type: str, optional
    :param samples: Samples to limit analysis to. If None, defaults to all samples in analysis object.
    :type samples: list, optional
    :param covariates: Additional variables to take into account in the model fitting. Defaults to None.
    :type covariates: list, optional
    :param output_dir: Output directory for analysis. Defaults to "results/differential_analysis_{data_type}".
                       If containing "{data_type}", will format string with variable.
    :type output_dir: str, optional
    :param output_prefix: Prefix for output files. Defaults to "differential_analysis".
    :type output_prefix: str, optional
    :param alpha: Significance level to use in diferential analysis.
                  Results for all features will be returned nonetheless. Defaults to 0.05.
    :type alpha: float, optional
    :param overwrite: Whether results should be overwritten in case they already exist. Defaults to True.
    :type overwrite: bool, optional
    :param distributed: Whether analysis should be distributed in a computing cluster for each comparison.
                        Currently, only a SLURM implementation is available.
                        If `True`, will not return results. Defaults to False.
    :type distributed: bool, optional

    :returns pandas.DataFrame: DataFrame with analysis results for all comparisons.
                               Will be `None` if `distributed` is `True`.
    """
    # Check comparisons
    # check comparison table has required columns
    req_attrs = ['comparison_name', 'comparison_side', 'sample_name', 'sample_group']
    if not all([x in comparison_table.columns for x in req_attrs]):
        raise AssertionError(
            "Given comparison table does not have all of '{}' columns.".format("', '".join(req_attrs)))
    # check all comparisons have samples in two sides
    if not all(comparison_table.groupby("comparison_name")["comparison_side"].nunique() == 2):
        raise AssertionError("All comparisons must have samples in each side of the comparison.")
    # check if any comparison and sample group has samples disagreeing in their side
    if not all(comparison_table.groupby(['comparison_name', "sample_group"])['comparison_side'].nunique() == 1):
        raise AssertionError(
            "Samples in same comparison and group must agree on their side of the comparison.")

    # Handle samples under analysis
    if samples is None:
        samples = analysis.samples
    samples = [s for s in samples if s.name in comparison_table['sample_name'].tolist()]

    # Make output dir
    if "{data_type}" in output_dir:
        output_dir = output_dir.format(data_type=data_type)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get matrix and samples
    if data_type == "ATAC-seq":
        count_matrix = analysis.coverage
    elif data_type == "RNA-seq":
        count_matrix = analysis.expression_matrix_counts
    else:
        raise AssertionError("Given data type does not match 'ATAC-seq' or 'RNA-seq'.")

    if samples is None:
        samples = analysis.samples
    samples = [s for s in samples if
               (s.name in comparison_table['sample_name'].tolist()) &
               (s.name in count_matrix.columns)]
    count_matrix = count_matrix[[s.name for s in samples]]

    # Get experiment matrix
    # by getting any other relevant covariates as required
    if covariates is not None:
        sample_table = pd.DataFrame([s.as_series() for s in samples])
        # check all covariates are in the samples and none is null
        if not all([x in sample_table.columns for x in covariates]):
            raise AssertionError(
                "Not all of the specified covariate variables are in the selected samples.")
        if sample_table[covariates].isnull().any().any():
            raise AssertionError(
                "None of the selected samples can have a Null value in the specified covariate variables.")

        # add covariates to comparison table
        comparison_table = comparison_table.set_index("sample_name").join(
            sample_table.set_index("sample_name")[covariates]).reset_index()

    # Make table for DESeq2
    experiment_matrix = comparison_table[
        ["sample_name", "sample_group"] + (covariates if covariates is not None else [])
    ].drop_duplicates()

    # Make formula for DESeq2
    formula = "~ {}sample_group".format(
        " + ".join(covariates) + " + " if covariates is not None else "")

    # Run DESeq2 analysis
    if not distributed:
        results = deseq_analysis(
            count_matrix, experiment_matrix, comparison_table,
            formula, output_dir, output_prefix, alpha=alpha, overwrite=overwrite)
        return results

    else:
        from pypiper.ngstk import NGSTk
        import textwrap
        tk = NGSTk()
        for comparison_name in comparison_table['comparison_name'].drop_duplicates():
            # make directory for comparison input/output
            out = os.path.join(os.path.abspath(output_dir), comparison_name)
            if not os.path.exists(out):
                os.makedirs(out)

            comp = comparison_table[comparison_table['comparison_name'] == comparison_name]
            comp.to_csv(os.path.join(out, "comparison_table.csv"), index=False)

            exp = experiment_matrix[experiment_matrix['sample_name'].isin(
                comp['sample_name'].tolist())]
            exp.to_csv(os.path.join(out, "experiment_matrix.csv"), index=False)

            count = count_matrix[comp['sample_name'].drop_duplicates()]
            count.to_csv(os.path.join(out, "count_matrix.csv"), index=True)

            job_name = "deseq_job.{}".format(comparison_name)
            log_file = os.path.join(out, job_name + ".log")
            job_file = os.path.join(out, job_name + ".sh")

            cmd = tk.slurm_header(
                job_name=job_name, output=log_file,
                cpus_per_task=2, mem_per_cpu=16000)
            # TODO: add DESeq2 script to toolkit and make path configurable
            cmd += "python ~/deseq_parallel.py"
            cmd += " --output_prefix {}".format(output_prefix)
            cmd += " --formula '{}'".format(formula)
            cmd += " --alpha {}".format(alpha)
            if overwrite:
                cmd += " --overwrite"
            cmd += "  {}\n".format(out)
            cmd += tk.slurm_footer()

            with open(job_file, "w") as handle:
                handle.write(textwrap.dedent(cmd))

            tk.slurm_submit_job(job_file)


def collect_differential_analysis(
        comparison_table,
        data_type="ATAC-seq",
        output_dir="results/differential_analysis_{data_type}",
        output_prefix="differential_analysis",
        permissive=True,
        overwrite=False):
    """
    Collect results from DESeq2 differential analysis.
    Particularly useful when runing `ngs_toolkit.general.differential_analysis` with
    `distributed` as `True`.

    :param comparison_table: A dataframe with 'comparison_name', 'comparison_side' and 'sample_name', 'sample_group' columns.
    :type comparison_table: pandas.DataFrame
    :param data_type: Type of data under analysis. One of "ATAC-seq" or "RNA-seq". Defaults to "ATAC-seq".
    :type data_type: str, optional
    :param output_dir: Output directory for analysis. Defaults to "results/differential_analysis_{data_type}".
                       If containing "{data_type}", will format string with variable.
    :type output_dir: str, optional
    :param output_prefix: Prefix for output files. Defaults to "differential_analysis".
    :type output_prefix: str, optional
    :param permissive: Whether non-existing files should be skipped or an error be thrown. Defaults to True.
    :type permissive: bool, optional
    :param overwrite: Whether results should be overwritten in case they already exist. Defaults to True.
    :type overwrite: bool, optional

    :returns pandas.DataFrame: DataFrame with analysis results for all comparisons.
                               Will be `None` if `overwrite` is `False` and a results file already exists.
    """
    from tqdm import tqdm

    # Make output dir
    if "{data_type}" in output_dir:
        output_dir = output_dir.format(data_type=data_type)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    results_file = os.path.join(output_dir, output_prefix + ".deseq_result.all_comparisons.csv")
    if not overwrite and os.path.exists(results_file):
        print("Differential analysis results '{}' already exist and argument `overwrite` is True.".format(
            results_file))
        return None

    comps = comparison_table["comparison_name"].drop_duplicates().sort_values()
    results = pd.DataFrame()
    for comp in tqdm(comps, desc="comparisons", total=len(comps)):
        out_file = os.path.join(output_dir, comp, output_prefix +
                                ".deseq_result.{}.csv".format(comp))
        # print("Collecting comparison '{}'".format(comp))
        # read
        try:
            res2 = pd.read_csv(out_file, index_col=0)
        except IOError as e:
            if permissive:
                print("Results file for comparison '{}' do not exist. Skipping.".format(comp))
                continue
            else:
                raise e
        # append
        results = results.append(res2.reset_index(), ignore_index=True)

    # save all
    results.to_csv(results_file, index=False)
    return results


def differential_overlap(
        differential,
        total,
        data_type="ATAC-seq",
        output_dir="results/differential_analysis_{data_type}",
        output_prefix="differential_analysis"):
    """
    Visualize intersection of sets of differential regions/genes.

    :param pandas.DataFrame differential: DataFrame containing result of comparisons filtered for features considered as differential.
    """
    import numpy as np
    import itertools
    import matplotlib.pyplot as plt
    import matplotlib
    import seaborn as sns
    from tqdm import tqdm
    from scipy.stats import fisher_exact
    from statsmodels.sandbox.stats.multicomp import multipletests

    if "{data_type}" in output_dir:
        output_dir = output_dir.format(data_type=data_type)

    if data_type == "ATAC-seq":
        unit = "region"
    elif data_type == "RNA-seq":
        unit = "gene"

    if "direction" not in differential.columns:
        differential["direction"] = differential["log2FoldChange"].apply(
            lambda x: "up" if x > 0 else "down")

    differential.index.name = "index"
    differential["intersect"] = 1
    piv = pd.pivot_table(differential.reset_index(), index='index', columns=[
                         'comparison_name', 'direction'], values='intersect', fill_value=0)

    intersections = pd.DataFrame(
        columns=["group1", "group2", "dir1", "dir2", "size1", "size2", "intersection", "union"])
    perms = list(itertools.permutations(piv.T.groupby(
        level=['comparison_name', 'direction']).groups.items(), 2))
    for ((k1, dir1), i1), ((k2, dir2), i2) in tqdm(perms, total=len(perms)):
        i1 = set(piv[i1][piv[i1] == 1].dropna().index)
        i2 = set(piv[i2][piv[i2] == 1].dropna().index)
        intersections = intersections.append(
            pd.Series(
                [k1, k2, dir1, dir2, len(i1), len(i2), len(i1.intersection(i2)), len(i1.union(i2))],
                index=["group1", "group2", "dir1", "dir2",
                       "size1", "size2", "intersection", "union"]
            ),
            ignore_index=True
        )
    # convert to %
    intersections['intersection'] = intersections['intersection'].astype(float)
    intersections['perc_1'] = intersections['intersection'] / intersections['size1'] * 100.
    intersections['perc_2'] = intersections['intersection'] / intersections['size2'] * 100.
    intersections['intersection_max_perc'] = intersections[['perc_1', 'perc_2']].max(axis=1)

    # calculate p-value from Fisher's exact test
    intersections['a'] = total - intersections[['size1', 'size2', 'intersection']].sum(axis=1)
    intersections['b'] = intersections['size1'] - intersections['intersection']
    intersections['c'] = intersections['size2'] - intersections['intersection']
    intersections['d'] = intersections['intersection']

    for i, row in intersections[['d', 'b', 'c', 'a']].astype(int).iterrows():
        odds, p = fisher_exact(
            row
            .values
            .reshape((2, 2)),
            alternative="greater")
        intersections.loc[i, 'odds_ratio'] = odds
        intersections.loc[i, 'p_value'] = p
    # intersections['q_value'] = intersections['p_value'] * intersections.shape[0]
    intersections['q_value'] = multipletests(intersections['p_value'])[1]
    intersections['log_p_value'] = (-np.log10(intersections['p_value'])
                                    ).fillna(0).replace(np.inf, 300)
    intersections['log_q_value'] = (-np.log10(intersections['q_value'])
                                    ).fillna(0).replace(np.inf, 300)

    # save
    intersections.to_csv(os.path.join(output_dir, output_prefix +
                                      ".differential_overlap.csv"), index=False)
    intersections = pd.read_csv(os.path.join(
        output_dir, output_prefix + ".differential_overlap.csv"))

    for metric, label, description, fill_value in [
            ("intersection", "intersection", "total in intersection", 0),
            ("intersection_max_perc", "percentage_overlap", "max of intersection %", 0),
            ("log_p_value", "significance", "p-value", 0)]:
        print(metric)
        # make pivot tables
        piv_up = pd.pivot_table(
            intersections[(intersections['dir1'] == "up") & (intersections['dir2'] == "up")],
            index="group1", columns="group2", values=metric).fillna(fill_value)
        piv_down = pd.pivot_table(
            intersections[(intersections['dir1'] == "down") & (intersections['dir2'] == "down")],
            index="group1", columns="group2", values=metric).fillna(fill_value)
        if metric == "intersection":
            piv_up = np.log10(1 + piv_up)
            piv_down = np.log10(1 + piv_down)
        np.fill_diagonal(piv_up.values, np.nan)
        np.fill_diagonal(piv_down.values, np.nan)

        # heatmaps
        if metric == "intersection_max_perc":
            extra = {"vmin": 0, "vmax": 100}
        else:
            extra = {}
        fig, axis = plt.subplots(1, 2, figsize=(8 * 2, 8), subplot_kw={"aspect": 'equal'})
        sns.heatmap(piv_down, square=True, cmap="Blues", cbar_kws={
                    "label": "Concordant {}s ({})".format(unit, description)}, ax=axis[0], **extra)
        sns.heatmap(piv_up, square=True, cmap="Reds", cbar_kws={
                    "label": "Concordant {}s ({})".format(unit, description)}, ax=axis[1], **extra)
        axis[0].set_title("Downregulated {}s".format(unit))
        axis[1].set_title("Upregulated {}s".format(unit))
        axis[0].set_xticklabels(axis[0].get_xticklabels(), rotation=90, ha="center")
        axis[1].set_xticklabels(axis[1].get_xticklabels(), rotation=90, ha="center")
        axis[0].set_yticklabels(axis[0].get_yticklabels(), rotation=0, ha="right")
        axis[1].set_yticklabels(axis[1].get_yticklabels(), rotation=0, ha="right")
        fig.savefig(os.path.join(output_dir, output_prefix +
                                 ".differential_overlap.{}.up_down_split.svg".format(label)), bbox_inches="tight")

        # combined heatmap
        # with upregulated {}s in upper square matrix and downredulated in down square
        piv_combined = pd.DataFrame(np.triu(piv_up), index=piv_up.index,
                                    columns=piv_up.columns).replace(0, np.nan)
        piv_combined.update(pd.DataFrame(np.tril(-piv_down), index=piv_down.index,
                                         columns=piv_down.columns).replace(0, np.nan))
        piv_combined = piv_combined.fillna(fill_value)
        if metric == "intersection":
            piv_combined = np.log10(1 + piv_combined)
        np.fill_diagonal(piv_combined.values, np.nan)

        if metric == "intersection_max_perc":
            extra = {"vmin": -150, "vmax": 150}
        else:
            extra = {}
        fig, axis = plt.subplots(1, figsize=(8, 8))
        sns.heatmap(piv_combined, square=True, cmap="RdBu_r", center=0, cbar_kws={
                    "label": "Concordant {}s ({})".format(unit, description)}, ax=axis, **extra)
        axis.set_xticklabels(axis.get_xticklabels(), rotation=90, ha="center")
        axis.set_yticklabels(axis.get_yticklabels(), rotation=0, ha="right")
        fig.savefig(os.path.join(output_dir, output_prefix +
                                 ".differential_overlap.{}.up_down_together.svg".format(label)), bbox_inches="tight")

        # Rank plots
        if metric == "log_pvalue":
            r = pd.melt(piv_combined.reset_index(), id_vars=[
                        'group1'], var_name="group2", value_name="agreement")
            r = r.dropna().sort_values('agreement')
            r = r.iloc[range(0, r.shape[0], 2)]
            r['rank'] = r['agreement'].rank(ascending=False)

            fig, axis = plt.subplots(1, 3, figsize=(3 * 4, 4))
            axis[0].scatter(r['rank'], r['agreement'])
            axis[0].axhline(0, linestyle="--", color="black")
            axis[1].scatter(r['rank'].tail(10), r['agreement'].tail(10))
            axis[2].scatter(r['rank'].head(10), r['agreement'].head(10))
            for i, row in r.tail(10).iterrows():
                axis[1].text(
                    row['rank'], row['agreement'],
                    s=row['group1'] + "\n" + row['group2'],
                    fontsize=5)
            for i, row in r.head(10).iterrows():
                axis[2].text(
                    row['rank'], row['agreement'],
                    s=row['group1'] + "\n" + row['group2'],
                    fontsize=5)
            for ax in axis:
                ax.set_ylabel("Agreement (-log(p-value))")
                ax.set_xlabel("Rank")
            sns.despine(fig)
            fig.savefig(os.path.join(output_dir, output_prefix +
                                     ".differential_overlap.{}.agreement.rank.svg".format(label)), bbox_inches="tight")

        # Observe disagreement
        # (overlap of down-regulated with up-regulated and vice-versa)
        piv_up = pd.pivot_table(
            intersections[(intersections['dir1'] == "up") & (intersections['dir2'] == "down")],
            index="group1", columns="group2", values=metric)
        piv_down = pd.pivot_table(
            intersections[(intersections['dir1'] == "down") & (intersections['dir2'] == "up")],
            index="group1", columns="group2", values=metric)

        piv_disagree = pd.concat([piv_up, piv_down]).groupby(level=0).max()
        if metric == "intersection":
            piv_disagree = np.log10(1 + piv_disagree)
        np.fill_diagonal(piv_disagree.values, np.nan)

        fig, axis = plt.subplots(1, 2, figsize=(16, 8), subplot_kw={"aspect": 'equal'})
        sns.heatmap(piv_disagree, square=True, cmap="Greens", cbar_kws={
                    "label": "Discordant {}s ({})".format(unit, description)}, ax=axis[0])
        # sns.heatmap(np.log2(1 + piv_disagree), square=True, cmap="Greens", cbar_kws={"label": "Discordant {}s (log2)".format(unit)}, ax=axis[1][0])

        norm = matplotlib.colors.Normalize(vmin=0, vmax=piv_disagree.max().max())
        cmap = plt.get_cmap("Greens")
        log_norm = matplotlib.colors.Normalize(vmin=np.log2(
            1 + piv_disagree).min().min(), vmax=np.log2(1 + piv_disagree).max().max())
        for j, g2 in enumerate(piv_disagree.index):
            for i, g1 in enumerate(piv_disagree.columns):
                # print(len(piv_disagree.index) - (j + 0.5), len(piv_disagree.index) - (i + 0.5))
                axis[1].scatter(
                    len(piv_disagree.index) - (j + 0.5), len(piv_disagree.index) - (i + 0.5),
                    s=(100 ** (norm(piv_disagree.loc[g1, g2]))) - 1, color=cmap(norm(piv_disagree.loc[g1, g2])), marker="o")
                # axis[1][1].scatter(
                #     len(piv_disagree.index) - (j + 0.5), len(piv_disagree.index) - (i + 0.5),
                #     s=(100 ** (log_norm(np.log2(1 + piv_disagree).loc[g1, g2]))) - 1, color=cmap(log_norm(np.log2(1 + piv_disagree).loc[g1, g2])), marker="o")

        axis[0].set_xticklabels(axis[0].get_xticklabels(), rotation=90, ha="center")
        axis[0].set_yticklabels(axis[0].get_yticklabels(), rotation=0, ha="right")
        axis[1].set_xlim((0, len(piv_disagree.index)))
        axis[1].set_ylim((0, len(piv_disagree.columns)))
        fig.savefig(os.path.join(output_dir, output_prefix +
                                 ".differential_overlap.{}.disagreement.svg".format(label)), bbox_inches="tight")

        # Rank plots
        if metric == "log_pvalue":
            r = pd.melt(piv_disagree.reset_index(), id_vars=[
                        'group1'], var_name="group2", value_name="disagreement")
            r = r.dropna().sort_values('disagreement')
            r = r.iloc[range(0, r.shape[0], 2)]
            r['rank'] = r['disagreement'].rank(ascending=False)

            fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 4), subplot_kw={"aspect": 'equal'})
            axis[0].scatter(r['rank'], r['disagreement'])
            axis[1].scatter(r['rank'].tail(10), r['disagreement'].tail(10))
            for i, row in r.tail(10).iterrows():
                axis[1].text(
                    row['rank'], row['disagreement'],
                    s=row['group1'] + "\n" + row['group2'],
                    fontsize=5)
            for ax in axis:
                ax.set_ylabel("Disagreement (-log(p-value))")
                ax.set_xlabel("Rank")
            sns.despine(fig)
            fig.savefig(os.path.join(output_dir, output_prefix +
                                     ".differential_overlap.{}.disagreement.rank.svg".format(label)), bbox_inches="tight")


def plot_differential(
        analysis,
        results,
        comparison_table=None,
        samples=None,
        matrix=None,
        only_comparison_samples=False,
        data_type="ATAC-seq",
        alpha=0.05,
        corrected_p_value=True,
        fold_change=None,
        diff_based_on_rank=False,
        max_rank=1000,
        ranking_variable="padj",
        respect_stat_thresholds=True,
        output_dir="results/differential_analysis_{data_type}",
        output_prefix="differential_analysis",
        plot_each_comparison=True,
        mean_column="baseMean",
        log_fold_change_column="log2FoldChange",
        p_value_column="pvalue",
        adjusted_p_value_column="padj",
        comparison_column="comparison_name",
        rasterized=True,
        dpi=300,
        robust=False,
        feature_labels=False,
        group_wise_colours=False,
        group_variables=None):
    """
    Plot differential features (e.g. chromatin region, genes) discovered with supervised group comparisons
    by ``ngs_toolkit.general.differential_analysis``.
    This will plot number and direction of discovered features, scatter, MA and volcano plots for each comparison
    and joint heatmaps of log fold changes, normalized values or Z-scores of individual samples or groups in the differential features.


    :param analysis: Analysis object.
    :type analysis: ngs_toolkit.general.Analysis
    :param results: Data frame with differential analysis results. See ``ngs_toolkit.general.differential_analysis`` for more information.
    :type results: pandas.DataFrame
    :param comparison_table: Comparison table, defaults to None. If provided, group-wise plots will be produced.
    :type comparison_table: pandas.DataFrame, optional
    :param samples: List of sample objects to restrict analysis to. Defaults to all samples in analysis object.
    :type samples: list, optional
    :param matrix: Matrix of quantification to use for plotting feature values across samples/groups.
                   Defaults to either "accessibility" for ATAC-seq analysis or "expression" for RNA-seq.
    :type matrix: str, optional
    :param only_comparison_samples: Whether to use only samples present in the `comparison_table`. Defaults to False.
    :type only_comparison_samples: bool, optional
    :param data_type: The data type being analyzed. Currently supported is "ATAC-seq" or "RNA-seq". Defaults to "ATAC-seq".
    :type data_type: str, optional
    :param alpha: Significance level to consider a feature differential. Defaults to 0.05.
    :type alpha: float, optional
    :param corrected_p_value: Whether to use a corrected p-valueto consider a feature differential. Defaults to True.
    :type corrected_p_value: bool, optional
    :param fold_change: Effect size (fold change) to consider a feature differential. Considers absolute values. Defaults to None.
    :type fold_change: float, optional
    :param diff_based_on_rank: Whether a feature should be considered differential based on its rank. Defaults to False
    :type diff_based_on_rank: bool, optional
    :param max_rank: Rank to use when using `diff_based_on_rank`. Defaults to 1000.
    :type max_rank: int, optional
    :param ranking_variable: Which variable to use for ranking when using `diff_based_on_rank`. Defaults to "padj".
    :type ranking_variable: str, optional
    :param respect_stat_thresholds: Whether the statistical thresholds from `alpha` and `fold_change` should still be respected when using `diff_based_on_rank`. Defaults to True
    :type respect_stat_thresholds: bool, optional
    :param output_dir: Directory to create output files. Defaults to "results/differential_analysis_{data_type}"
    :type output_dir: str, optional
    :param output_prefix: Prefix to use when creating output files. Defaults to "differential_analysis".
    :type output_prefix: str, optional
    :param plot_each_comparison: Whether each comparison should be plotted in scatter, MA and volcano plots. Useful to turn off with many comparisons. Defaults to True.
    :type plot_each_comparison: bool, optional
    :param mean_column: Column  in `results` data frame containing values for mean values across samples. Defaults to "baseMean".
    :type mean_column: str, optional
    :param log_fold_change_column: Column  in `results` data frame containing values for log2FoldChange values across samples. Defaults to "log2FoldChange".
    :type log_fold_change_column: str, optional
    :param p_value_column: Column  in `results` data frame containing values for p-values across samples. Defaults to "pvalue".
    :type p_value_column: str, optional
    :param adjusted_p_value_column: Column  in `results` data frame containing values for adjusted p-values across samples. Defaults to "padj".
    :type adjusted_p_value_column: str, optional
    :param comparison_column: Column  in `results` data frame containing the name of the comparison. Defaults to "comparison_name".
    :type comparison_column: str, optional
    :param rasterized: Whether plots with many objects should be rasterized. Defaults to True.
    :type rasterized: bool, optional
    :param dpi: Rasterization resolution (dpi). Defaults to 300.
    :type dpi: number, optional
    :param robust: Whether heatmap color scale ranges should be robust (using quantiles) rather than extreme values. Useful for noisy/extreme data. Defaults to False.
    :type robust: bool, optional
    :param feature_labels: Whether features (regions/genes) should be labeled in heatmaps. Defaults to False.
    :type feature_labels: bool, optional
    :param group_wise_colours: Whether groups of samples should be coloured in heatmaps. Defaults to False.
    :type group_wise_colours: bool, optional
    :param group_variables: Which variables to colour if `group_wise_colours` if True. Defaults to None (must be given).
    :type group_variables: list, optional
    """
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns

    req_attrs = [mean_column, log_fold_change_column,
                 p_value_column, adjusted_p_value_column, comparison_column]
    if not all([x in results.columns for x in req_attrs]):
        raise AssertionError(
            "Results dataframe must have '{}' columns.".format(", ".join(req_attrs)))

    # Make output dir
    if "{data_type}" in output_dir:
        output_dir = output_dir.format(data_type=data_type)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get matrix and samples
    results = results.copy()
    if data_type == "ATAC-seq":
        if matrix is None:
            matrix = analysis.accessibility
        results.index.name = "region"
        var_name = "region"
        quantity = "Accessibility"
        unit = "RPM"
    elif data_type == "RNA-seq":
        if matrix is None:
            matrix = analysis.expression
        results.index.name = "gene_name"
        var_name = "gene"
        quantity = "Expression"
        unit = "TPM"
    else:
        raise AssertionError("Given data type does not match 'ATAC-seq' or 'RNA-seq'.")

    if samples is None:
        samples = analysis.samples
    samples = [s for s in samples if s.name in matrix.columns]
    if only_comparison_samples and comparison_table is not None:
        samples = [s for s in samples if s.name in comparison_table['sample_name'].tolist()]
    matrix = matrix[[s.name for s in samples]]

    # Handle group colouring
    if group_wise_colours:
        if type(group_variables) is None:
            raise AssertionError(
                "If `group_wise_colours` is True, a list of `group_variables` must be passed.")

        # This will always be a matrix for all samples
        color_dataframe = pd.DataFrame(
            analysis.get_level_colors(index=matrix.columns, levels=group_variables),
            index=group_variables, columns=matrix.columns.get_level_values("sample_name"))
        # will be filtered now by the requested samples if needed
        color_dataframe = color_dataframe[[s.name for s in samples]]
        color_dataframe = color_dataframe.loc[:, matrix.columns.get_level_values("sample_name")]

    # Extract significant based on p-value and fold-change
    if fold_change is not None:
        fc = (results[log_fold_change_column].abs() > fold_change)
    else:
        fc = [True] * results.shape[0]
    if corrected_p_value:
        p_var = adjusted_p_value_column
    else:
        p_var = p_value_column
    results.loc[(results[p_var] < alpha) & fc, "diff"] = True
    results.loc[:, "diff"] = results.loc[:, "diff"].fillna(False)
    # Declare significant based on top ranked features
    if diff_based_on_rank:
        for comparison in results[comparison_column].unique():
            if ranking_variable == log_fold_change_column:
                i = results.loc[results[comparison_column] == comparison,
                                ranking_variable].abs().sort_values().tail(max_rank).index
            else:
                i = results.loc[results[comparison_column] == comparison,
                                ranking_variable].sort_values().head(max_rank).index
            results.loc[(results[comparison_column] == comparison) &
                        results.index.isin(i), "diff_rank"] = True
        results.loc[:, "diff_rank"] = results.loc[:, "diff_rank"].fillna(False)
        if respect_stat_thresholds:
            results.loc[:, 'diff'] = (results.loc[:, 'diff'] == True) & (
                results.loc[:, 'diff_rank'] == True)
        else:
            results.loc[:, 'diff'] = results.loc[:, 'diff_rank']

    # Annotate direction of change
    results.loc[:, "direction"] = results.loc[:, log_fold_change_column].apply(
        lambda x: "up" if x >= 0 else "down")

    if results.loc[:, "diff"].sum() < 1:
        print("No significantly different regions found in any comparison.")
        return

    # PLOTS
    comparisons = sorted(results[comparison_column].drop_duplicates())
    n_side = int(np.ceil(np.sqrt(len(comparisons))))

    # P-value distributions
    fig, axis = plt.subplots(1, 1, figsize=(4, 4))
    sns.distplot(results[p_value_column].dropna(), kde=False, ax=axis)
    axis.set_xlabel("P-value")
    axis.set_ylabel("Frequency")
    sns.despine(fig)
    fig.savefig(os.path.join(output_dir, output_prefix +
                             ".pvalue_distribution.svg"), bbox_inches="tight")

    if plot_each_comparison:
        # per comparison
        g = sns.FacetGrid(data=results, col=comparison_column, col_wrap=n_side)
        g.map(sns.distplot, p_value_column, kde=False)
        for ax in g.axes:
            ax.set_xlabel("P-value")
            ax.set_ylabel("Frequency")
        sns.despine(g.fig)
        g.fig.savefig(os.path.join(output_dir, output_prefix +
                                   ".pvalue_distribution.per_comparison.svg"), bbox_inches="tight")

    # Number of differential vars
    n_vars = float(matrix.shape[0])
    total_diff = results.groupby([comparison_column])[
        'diff'].sum().sort_values(ascending=False).reset_index()
    split_diff = results.groupby([comparison_column, "direction"])[
        'diff'].sum().sort_values(ascending=False).reset_index()
    split_diff.loc[split_diff['direction'] == 'down', "diff"] *= -1
    split_diff['label'] = split_diff['comparison_name'] + ", " + split_diff['direction']
    total_diff['diff_perc'] = (total_diff['diff'] / n_vars) * 100
    split_diff['diff_perc'] = (split_diff['diff'] / n_vars) * 100

    fig, axis = plt.subplots(3, 2, figsize=(4 * 2, 4 * 3))
    sns.barplot(data=total_diff, x="diff", y="comparison_name", orient="h", ax=axis[0, 0])
    sns.barplot(data=total_diff, x="diff_perc", y="comparison_name", orient="h", ax=axis[0, 1])
    sns.barplot(data=split_diff, x="diff", y="label", orient="h", ax=axis[1, 0])
    sns.barplot(data=split_diff, x="diff_perc", y="label", orient="h", ax=axis[1, 1])
    sns.barplot(data=split_diff, x="diff", y="comparison_name",
                hue="direction", orient="h", ax=axis[2, 0])
    sns.barplot(data=split_diff, x="diff_perc", y="comparison_name",
                hue="direction", orient="h", ax=axis[2, 1])
    for ax in axis[:, 0]:
        ax.set_xlabel("N. diff")
    for ax in axis[:, 1]:
        ax.set_xlabel("N. diff (% of total)")
    for ax in axis[1:, :].flatten():
        ax.axvline(0, linestyle="--", color="black", alpha=0.6)
    m = split_diff['diff'].abs().max()
    axis[2, 0].set_xlim((-m, m))
    m = split_diff['diff_perc'].abs().max()
    axis[2, 1].set_xlim((-m, m))
    sns.despine(fig)
    fig.savefig(os.path.join(output_dir, output_prefix +
                             ".number_differential.directional.svg"), bbox_inches="tight")

    if plot_each_comparison:
        # Pairwise scatter plots
        # TODO: add different shadings of red for various levels of significance
        # TODO: do this for scatter, MA, volcano plots
        if comparison_table is not None:
            fig, axes = plt.subplots(n_side, n_side, figsize=(
                n_side * 4, n_side * 4), sharex=True, sharey=True)
            if n_side > 1 or n_side > 1:
                axes = iter(axes.flatten())
            else:
                axes = iter([axes])
            for comparison in comparisons:
                c = comparison_table[comparison_table[comparison_column] == comparison]
                a = c.loc[c['comparison_side'] >= 1, "sample_name"]
                b = c.loc[c['comparison_side'] <= 0, "sample_name"]

                a = matrix[[s.name for s in samples if s.name in a.tolist() and s.library ==
                            data_type]].mean(axis=1)
                b = matrix[[s.name for s in samples if s.name in b.tolist() and s.library ==
                            data_type]].mean(axis=1)

                # Hexbin plot
                ax = next(axes)
                ax.hexbin(b, a, alpha=0.85, cmap="Greys", color="black", edgecolors="white",
                          linewidths=0, bins='log', mincnt=1, rasterized=True)

                # Scatter for significant
                diff_vars = results[(results[comparison_column] == comparison)
                                    & (results["diff"] == True)].index
                if diff_vars.shape[0] > 0:
                    ax.scatter(b.loc[diff_vars], a.loc[diff_vars], alpha=0.1, color="red", s=2)
                ax.set_title(comparison)
                ax.set_xlabel("Down")
                ax.set_ylabel("Up")
                # x = y square
                lims = [np.min([ax.get_xlim(), ax.get_ylim()]),
                        np.max([ax.get_xlim(), ax.get_ylim()])]
                ax.plot(lims, lims, linestyle='--', alpha=0.5, zorder=0, color="black")
                ax.set_aspect('equal')
                ax.set_xlim(lims)
                ax.set_ylim(lims)
            for ax in axes:
                ax.set_visible(False)
            sns.despine(fig)
            fig.savefig(os.path.join(output_dir, output_prefix +
                                     ".scatter_plots.svg"), bbox_inches="tight", dpi=dpi)

        # Volcano plots
        fig, axes = plt.subplots(n_side, n_side, figsize=(
            n_side * 4, n_side * 4), sharex=False, sharey=False)
        if n_side > 1 or n_side > 1:
            axes = iter(axes.flatten())
        else:
            axes = iter([axes])
        for comparison in comparisons:
            t = results.loc[results[comparison_column] == comparison, :]

            # Hexbin plot
            ax = next(axes)
            ax.hexbin(
                t[log_fold_change_column], -np.log10(t[p_value_column]),
                alpha=0.85, cmap="Greys", color="black", edgecolors="white", linewidths=0, bins='log', mincnt=1, rasterized=True)

            # Scatter for significant
            diff_vars = t[t["diff"] == True].index
            if diff_vars.shape[0] > 0:
                ax.scatter(t.loc[diff_vars, log_fold_change_column], -
                           np.log10(t.loc[diff_vars, p_value_column]), alpha=0.1, color="red", s=2)
            ax.set_title(comparison)
            ax.set_xlabel("log2(Fold-change)")
            ax.set_ylabel("-log10(P-value)")
            ax.axvline(0, linestyle='--', alpha=0.5, zorder=0, color="black")
            l = np.max([abs(i) for i in ax.get_xlim()])
            ax.set_xlim(-l, l)

            # Add lines of significance
            ax.axhline(-np.log10(t.loc[t["diff"] == True, p_value_column].max()),
                       linestyle='--', alpha=0.5, zorder=0, color="black")
            if fold_change is not None:
                ax.axvline(-fold_change, linestyle='--', alpha=0.5, zorder=0, color="black")
                ax.axvline(fold_change, linestyle='--', alpha=0.5, zorder=0, color="black")
        for ax in axes:
            ax.set_visible(False)
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, output_prefix +
                                 ".volcano_plots.svg"), bbox_inches="tight", dpi=dpi)

        # MA plots
        fig, axes = plt.subplots(n_side, n_side, figsize=(
            n_side * 4, n_side * 4), sharex=False, sharey=False)
        if n_side > 1 or n_side > 1:
            axes = iter(axes.flatten())
        else:
            axes = iter([axes])
        for comparison in comparisons:
            t = results.loc[results[comparison_column] == comparison, :]

            # Hexbin plot
            ax = next(axes)
            ax.hexbin(
                np.log10(t[mean_column]), t[log_fold_change_column],
                alpha=0.85, cmap="Greys", color="black", edgecolors="white", linewidths=0, bins='log', mincnt=1, rasterized=True)

            # Scatter for significant
            diff_vars = t[t["diff"] == True].index
            if diff_vars.shape[0] > 0:
                ax.scatter(np.log10(t.loc[diff_vars, mean_column]), t.loc[diff_vars,
                                                                          log_fold_change_column], alpha=0.1, color="red", s=2)
            ax.set_title(comparison)
            ax.set_xlabel("Amplitude")
            ax.set_ylabel("log2(Fold-change)")
            ax.axhline(0, linestyle='--', alpha=0.5, zorder=0, color="black")
            l = np.max([abs(i) for i in ax.get_ylim()])
            ax.set_ylim(-l, l)

            # Add lines of significance
            if fold_change is not None:
                ax.axhline(-fold_change, linestyle='--', alpha=0.5, zorder=0, color="black")
                ax.axhline(fold_change, linestyle='--', alpha=0.5, zorder=0, color="black")
        for ax in axes:
            ax.set_visible(False)
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, output_prefix +
                                 ".ma_plots.svg"), bbox_inches="tight", dpi=dpi)

    # Observe values of variables across all comparisons
    all_diff = results[results["diff"] == True].index.drop_duplicates()
    if type(matrix.columns) is pd.MultiIndex:
        sample_cols = matrix.columns.get_level_values("sample_name").tolist()
    else:
        sample_cols = matrix.columns.tolist()

    if comparison_table is not None:
        if results[comparison_column].drop_duplicates().shape[0] > 1:
            groups = pd.DataFrame()
            for sample_group in comparison_table["sample_group"].drop_duplicates():
                c = comparison_table.loc[comparison_table["sample_group"]
                                         == sample_group, "sample_name"].drop_duplicates()
                if c.shape[0] > 0:
                    groups.loc[:, sample_group] = matrix[[
                        d for d in c if d in sample_cols]].mean(axis=1)

            if groups.empty:
                # It seems comparisons were not done in a all-versus-all fashion
                for group in comparison_table["sample_group"].drop_duplicates():
                    c = comparison_table.loc[comparison_table["sample_group"]
                                             == group, "sample_name"].drop_duplicates()
                    if c.shape[0] > 0:
                        groups.loc[:, group] = matrix[c].mean(axis=1)

            # Select only differential regions from groups
            groups = groups.loc[all_diff, :]

            figsize = (max(5, 0.12 * groups.shape[1]), 5)
            # Heatmaps
            # Comparison level
            g = sns.clustermap(
                groups.corr(),
                xticklabels=False, yticklabels=True, cbar_kws={
                    "label": "Pearson correlation\non differential {}s".format(var_name)},
                cmap="BuGn", metric="correlation", rasterized=True, figsize=(figsize[0], figsize[0]))
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(),
                                         rotation=0, fontsize="xx-small")
            g.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.groups.clustermap.corr.svg".format(
                var_name)), bbox_inches="tight", dpi=dpi, metric="correlation")

            g = sns.clustermap(
                groups,
                xticklabels=True, yticklabels=feature_labels, cbar_kws={
                    "label": "{} of\ndifferential {}s".format(quantity, var_name)},
                cmap="BuGn", robust=robust, metric="correlation", rasterized=True, figsize=figsize)
            g.ax_heatmap.set_ylabel("Differential {}s (n = {})".format(var_name, groups.shape[0]))
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(),
                                         rotation=90, fontsize="xx-small")
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(),
                                         rotation=0, fontsize="xx-small")
            g.fig.savefig(os.path.join(output_dir, output_prefix +
                                       ".diff_{}.groups.clustermap.svg".format(var_name)), bbox_inches="tight", dpi=dpi)

            g = sns.clustermap(
                groups,
                xticklabels=True, yticklabels=feature_labels, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}s".format(quantity, var_name)},
                cmap="RdBu_r", robust=robust, metric="correlation", rasterized=True, figsize=figsize)
            g.ax_heatmap.set_ylabel("Differential {}s (n = {})".format(var_name, groups.shape[0]))
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(),
                                         rotation=90, fontsize="xx-small")
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(),
                                         rotation=0, fontsize="xx-small")
            g.fig.savefig(os.path.join(output_dir, output_prefix +
                                       ".diff_{}.groups.clustermap.z0.svg".format(var_name)), bbox_inches="tight", dpi=dpi)

    # Fold-changes and P-values
    # pivot table of genes vs comparisons
    fold_changes = pd.pivot_table(
        results.loc[all_diff, :].reset_index(),
        index=results.index.name, columns=comparison_column,
        values=log_fold_change_column).fillna(0)
    p_values = -np.log10(pd.pivot_table(
        results.loc[all_diff, :].reset_index(),
        index=results.index.name, columns=comparison_column,
        values=adjusted_p_value_column))

    # fold
    if fold_changes.shape[1] > 1:
        figsize = (max(5, 0.12 * fold_changes.shape[1]), 5)

        g = sns.clustermap(
            fold_changes.corr(),
            xticklabels=False, yticklabels=True, cbar_kws={
                "label": "Pearson correlation\non fold-changes"},
            cmap="BuGn", vmin=0, vmax=1, metric="correlation", rasterized=True, figsize=(figsize[0], figsize[0]))
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(),
                                     rotation=0, fontsize="xx-small")
        g.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.groups.fold_changes.clustermap.corr.svg".format(
            var_name)), bbox_inches="tight", dpi=dpi, metric="correlation")

        try:
            g = sns.clustermap(
                fold_changes.loc[all_diff, :],
                xticklabels=True, yticklabels=feature_labels, cbar_kws={
                    "label": "Fold-change of\ndifferential {}s".format(var_name)},
                cmap="RdBu_r", robust=robust, metric="correlation", rasterized=True, figsize=figsize)
            g.ax_heatmap.set_ylabel("Differential {}s (n = {})".format(
                var_name, fold_changes.loc[all_diff, :].shape[0]))
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(),
                                         rotation=90, fontsize="xx-small")
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(),
                                         rotation=0, fontsize="xx-small")
            g.fig.savefig(os.path.join(output_dir, output_prefix +
                                       ".diff_{}.groups.fold_changes.clustermap.svg".format(var_name)), bbox_inches="tight", dpi=dpi)
        except FloatingPointError as e:
            pass

    # Sample level
    if type(matrix.columns) is pd.core.indexes.multi.MultiIndex:
        matrix.columns = matrix.columns.get_level_values("sample_name")

    matrix2 = matrix.loc[all_diff, :]
    figsize = (max(5, 0.12 * matrix2.shape[1]), 5)
    if group_wise_colours:
        extra = {"col_colors": color_dataframe.values}
    else:
        extra = {}

    g = sns.clustermap(
        matrix2.corr(),
        yticklabels=True, xticklabels=False,
        cbar_kws={"label": "Pearson correlation\non differential {}s".format(var_name)},
        cmap="BuGn", metric="correlation", figsize=(figsize[0], figsize[0]), rasterized=rasterized, robust=robust, **extra)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    g.fig.savefig(os.path.join(output_dir, output_prefix +
                               ".diff_{}.samples.clustermap.corr.svg".format(var_name)), bbox_inches="tight", dpi=dpi)

    g = sns.clustermap(
        matrix2,
        yticklabels=feature_labels, cbar_kws={
            "label": "{} of\ndifferential {}s".format(quantity, var_name)},
        xticklabels=True, vmin=0, cmap="RdBu_r", metric="correlation", figsize=figsize, rasterized=rasterized, robust=robust, **extra)
    g.ax_heatmap.set_ylabel("Differential {}s (n = {})".format(var_name, matrix2.shape[0]))
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    g.fig.savefig(os.path.join(output_dir, output_prefix +
                               ".diff_{}.samples.clustermap.svg".format(var_name)), bbox_inches="tight", dpi=dpi)

    g = sns.clustermap(
        matrix2,
        yticklabels=feature_labels, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}s".format(quantity, var_name)},
        xticklabels=True, cmap="RdBu_r", metric="correlation", figsize=figsize, rasterized=rasterized, robust=robust, **extra)
    g.ax_heatmap.set_ylabel("Differential {}s (n = {})".format(var_name, matrix2.shape[0]))
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
    g.fig.savefig(os.path.join(output_dir, output_prefix +
                               ".diff_{}.samples.clustermap.z0.svg".format(var_name)), bbox_inches="tight", dpi=dpi)


def lola(bed_files, universe_file, output_folder, genome="hg19", cpus=8):
    """
    Performs location overlap analysis (LOLA) on bedfiles with regions sets.
    """
    import rpy2.robjects as robj

    run = robj.r("""
        function(bedFiles, universeFile, outputFolder) {
            library("LOLA")

            userUniverse  <- LOLA::readBed(universeFile)

            dbPath1 = "/data/groups/lab_bock/shared/resources/regions/LOLACore/{genome}/"
            dbPath2 = "/data/groups/lab_bock/shared/resources/regions/customRegionDB/{genome}/"
            regionDB = loadRegionDB(c(dbPath1, dbPath2))

            if (typeof(bedFiles) == "character") {
                userSet <- LOLA::readBed(bedFiles)
                lolaResults = runLOLA(list(userSet), userUniverse, regionDB, cores={cpus})
                writeCombinedEnrichment(lolaResults, outFolder=outputFolder)
            } else if (typeof(bedFiles) == "double") {
                for (bedFile in bedFiles) {
                    userSet <- LOLA::readBed(bedFile)
                    lolaResults = runLOLA(list(userSet), userUniverse, regionDB, cores={cpus})
                    writeCombinedEnrichment(lolaResults, outFolder=outputFolder)
                }
            }
        }
    """.format(genome=genome, cpus=cpus))

    # convert the pandas dataframe to an R dataframe
    run(bed_files, universe_file, output_folder)


def bed_to_fasta(bed_file, fasta_file, genome="hg19", genome_2bit=None):
    import os
    import pandas as pd

    if genome_2bit is None:
        genome_2bit = "~/resources/genomes/{genome}/{genome}.2bit".format(genome=genome)

    # write name column
    bed = pd.read_csv(bed_file, sep='\t', header=None)
    bed['name'] = bed[0] + ":" + bed[1].astype(str) + "-" + bed[2].astype(str)
    bed[1] = bed[1].astype(int)
    bed[2] = bed[2].astype(int)
    bed.to_csv(bed_file + ".tmp.bed", sep='\t', header=None, index=False)

    # do enrichment
    cmd = "twoBitToFa {0} -bed={1} {2}".format(genome_2bit, bed_file + ".tmp.bed", fasta_file)

    os.system(cmd)
    # os.system("rm %s" % bed_file + ".tmp.bed")


def meme_ame(input_fasta, output_dir, background_fasta=None, organism="human"):
    import os

    dbs = {
        "human": "~/resources/motifs/motif_databases/HUMAN/HOCOMOCOv10.meme",
        "mouse": "~/resources/motifs/motif_databases/MOUSE/uniprobe_mouse.meme"
    }
    # shuffle input in no background is provided
    if background_fasta is None:
        shuffled = input_fasta + ".shuffled"
        cmd = """
        fasta-dinucleotide-shuffle -c 1 -f {0} > {1}
        """.format(input_fasta, shuffled)
        os.system(cmd)

    cmd = """
    ame --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 0.05 \\
    --control {0} -o {1} {2} {3}
    """.format(
        background_fasta if background_fasta is not None else shuffled,
        output_dir, input_fasta, dbs[organism])
    os.system(cmd)
    # os.system("rm %s" % shuffled)


def parse_ame(ame_dir):
    """
    Parse results of MEME-AME motif enrichment.

    :param ame_dir: Directory with MEME-AME results.
    :type ame_dir: str
    :returns: Data frame with enrichment statistics for each found TF motif.
    :rtype: pandas.DataFrame
    :raises: IOError
    """
    import os
    import pandas as pd

    with open(os.path.join(ame_dir, "ame.txt"), 'r') as handle:
        lines = handle.readlines()

    output = list()
    for line in lines:
        # skip header lines
        if line[0] not in [str(i) for i in range(10)]:
            continue

        # get motif string and the first half of it (simple name)
        motif = line.strip().split(" ")[5].split("_")[0]
        # get corrected p-value
        q_value = float(line.strip().split(" ")[-2])
        # append
        output.append((motif, q_value))

    return pd.Series(dict(output))


def homer_motifs(bed_file, output_dir, genome="hg19"):
    cmd = "findMotifsGenome.pl {bed} {genome}r {out_dir} \
    -size 1000 -h -p 2 -len 8,10,12,14 -noknown".format(
        bed=bed_file, genome=genome, out_dir=output_dir
    )
    os.system(cmd)


def parse_homer(homer_dir):
    """
    Parse results of HOMER findMotifs.pl de novo motif enrichment.

    :param homer_dir: Directory with HOMER results.
    :type homer_dir: str
    :returns: Data frame with enrichment statistics for each found TF motif.
    :rtype: pandas.DataFrame
    :raises: IOError
    """
    import glob
    import os
    import re
    import pandas as pd

    motif_htmls = sorted(glob.glob(os.path.join(homer_dir, "motif*.info.html")))

    if len(motif_htmls) < 1:
        raise IOError("Homer directory does not contain any discovered motifs.")

    output = pd.DataFrame()
    for motif_html in motif_htmls:

        motif = int(re.sub(".info.html", "", re.sub(
            os.path.join(homer_dir, "motif"), "", motif_html)))

        with open(motif_html, 'r') as handle:
            content = handle.read()

        # Parse table with motif info
        info_table = content[
            re.search("""<TABLE border="1" cellpading="0" cellspacing="0">""", content).end():
            re.search("</TABLE>", content).start()].strip()

        info_table = pd.DataFrame([x.split("</TD><TD>")
                                   for x in info_table.replace("<TR><TD>", "").split("</TD></TR>")])
        info_table.columns = ["description", "value"]
        info_table["description"] = info_table["description"].str.strip()
        info_table["motif"] = motif

        # Add most probable known motif name
        info_table["known_motif"] = content[
            re.search("<H4>", content).end():
            re.search("</H4>", content).start()]

        # append
        output = output.append(info_table, ignore_index=True)

    return output.sort_values("motif")


def enrichr(dataframe, gene_set_libraries=None, kind="genes"):
    """
    Use Enrichr on a list of genes (currently only genes supported through the API).
    """
    import json
    import requests
    import pandas as pd
    from tqdm import tqdm

    ENRICHR_ADD = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    ENRICHR_RETRIEVE = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'

    if gene_set_libraries is None:
        gene_set_libraries = [
            'GO_Biological_Process_2015',
            "ChEA_2015",
            "KEGG_2016",
            "ESCAPE",
            # "Epigenomics_Roadmap_HM_ChIP-seq",
            "ENCODE_TF_ChIP-seq_2015",
            "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",
            # "ENCODE_Histone_Modifications_2015",
            # "OMIM_Expanded",
            # "TF-LOF_Expression_from_GEO",
            "Single_Gene_Perturbations_from_GEO_down",
            "Single_Gene_Perturbations_from_GEO_up",
            # "Disease_Perturbations_from_GEO_down",
            # "Disease_Perturbations_from_GEO_up",
            "Drug_Perturbations_from_GEO_down",
            "Drug_Perturbations_from_GEO_up",
            "WikiPathways_2016",
            "Reactome_2016",
            "BioCarta_2016",
            "NCI-Nature_2016"
        ]

    results = pd.DataFrame()
    for gene_set_library in tqdm(gene_set_libraries):
        print("Using enricher on %s gene set library." % gene_set_library)

        if kind == "genes":
            # Build payload with bed file
            attr = "\n".join(dataframe["gene_name"].dropna().tolist())
        elif kind == "regions":
            # Build payload with bed file
            attr = "\n".join(dataframe[['chrom', 'start', 'end']].apply(
                lambda x: "\t".join([str(i) for i in x]), axis=1).tolist())

        payload = {
            'list': (None, attr),
            'description': (None, gene_set_library)
        }
        # Request adding gene set
        response = requests.post(ENRICHR_ADD, files=payload)
        if not response.ok:
            raise Exception('Error analyzing gene list')

        # Track gene set ID
        user_list_id = json.loads(response.text)['userListId']

        # Request enriched sets in gene set
        response = requests.get(
            ENRICHR_RETRIEVE + query_string % (user_list_id, gene_set_library)
        )
        if not response.ok:
            raise Exception('Error fetching enrichment results')

        # Get enriched sets in gene set
        res = json.loads(response.text)
        # If there's no enrichemnt, continue
        if len(res) < 0:
            continue

        # Put in dataframe
        res = pd.DataFrame([pd.Series(s) for s in res[gene_set_library]])
        if res.shape[0] == 0:
            continue
        if len(res.columns) == 7:
            res.columns = ["rank", "description", "p_value", "z_score",
                           "combined_score", "genes", "adjusted_p_value"]
        elif len(res.columns) == 9:
            res.columns = ["rank", "description", "p_value", "z_score", "combined_score",
                           "genes", "adjusted_p_value", "old_p_value", "old_adjusted_p_value"]

        # Remember gene set library used
        res["gene_set_library"] = gene_set_library

        # Append to master dataframe
        results = results.append(res, ignore_index=True)

    return results


def run_enrichment_jobs(
        analysis_name, results_dir, genome,
        background_bed="results/{PROJECT_NAME}_peak_set.bed"):
    """
    Submit enrichment jobs for a specifc analysis.
    """
    # LOLA
    cmds = ["""for F in `find {results_dir} -name "*_regions.bed"`; do
DIR=`dirname $F`
if [ ! -f ${{DIR}}/allEnrichments.txt ]; then
echo $DIR $F
sbatch -J lola.$F -o $F.lola.log -p shortq -c 8 --mem 24000 \
--wrap "Rscript ~/jobs/run_LOLA.R $F {background_bed} {GENOME}"
fi
done""".format(results_dir=results_dir, background_bed=background_bed, GENOME=genome)]

    # AME
    dbs = {
        "human": "~/resources/motifs/motif_databases/HUMAN/HOCOMOCOv10.meme",
        "mouse": "~/resources/motifs/motif_databases/MOUSE/uniprobe_mouse.meme"}
    omap = {"hg38": "human", "hg19": "human", "mm10": "mouse"}

    cmds += ["""for F in `find {results_dir} -name "*_regions.fa"`; do
DIR=`dirname $F`
if [ ! -f ${{DIR}}/ame.html ]; then
echo $DIR $F
sbatch -J "meme_ame.${{F}}" -o "${{F}}.meme_ame.log" -p shortq -c 1 --mem 4000 \
--wrap "fasta-dinucleotide-shuffle -c 1 -f "${{F}}" > "${{F}}".shuffled.fa; \
ame --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 0.05 \
--control "${{F}}".shuffled.fa -o "${{DIR}}" "${{F}}" {motifs}"
fi
done""".format(results_dir=results_dir, motifs=dbs[omap[genome]])]

    # HOMER
    cmds += ["""for F in `find {results_dir} -name "*_regions.bed"`; do
DIR=`dirname $F`
if [ ! -f ${{DIR}}/homerResults.html ]; then
echo $DIR $F
sbatch -J "homer.${{F}}" -o "${{F}}.homer.log" -p shortq -c 8 --mem 20000 \
--wrap "findMotifsGenome.pl ${{F}} {GENOME}r ${{DIR}} -size 1000 -h -p 2 -len 8,10,12,14 -noknown"
fi
done""".format(results_dir=results_dir, GENOME=genome)]

    # Enrichr
    cmds += ["""for F in `find {results_dir} -name "*.gene_symbols.txt"`; do
if [ ! -f ${{F/gene_symbols.txt/enrichr.csv}} ]; then
echo $F
sbatch -J enrichr.$F -o $F.enrichr.log -p shortq -c 1 --mem 4000 \
--wrap "python -u ~/jobs/run_Enrichr.py --input-file "$F" --output-file "${{F/gene_symbols.txt/enrichr.csv}}" "
fi
done""".format(results_dir=results_dir)]
    cmds += ["""for F in `find {results_dir} -name "*_genes.symbols.txt"`; do
if [ ! -f ${{F/symbols.txt/enrichr.csv}} ]; then
echo $F
sbatch -J enrichr.$F -o $F.enrichr.log -p shortq -c 1 --mem 4000 \
--wrap "python -u ~/jobs/run_Enrichr.py --input-file "$F" --output-file "${{F/symbols.txt/enrichr.csv}}" "
fi
done""".format(results_dir=results_dir)]

    for cmd in cmds:
        os.system(cmd)


def differential_enrichment(
    analysis,
    differential,
    data_type="ATAC-seq",
    output_dir="results/differential_analysis_{data_type}",
    output_prefix="differential_analysis",
    genome="hg19",
    directional=True,
    max_diff=1000,
    sort_var="pvalue",
    as_jobs=True
):
    """
    Perform various types of enrichment analysis given a dataframe
    of the results from differential analysis.
    Performs enrichment of gene sets (RNA-seq and ATAC-seq), genomic regions and TF motifs (ATAC-seq only).

    :param analysis: Analysis object.
    :type analysis: ngs_toolkit.general.Analysis
    :param differential: Data frame with differential results as produced by
                         ``ngs_toolkit.general.differential_analysis``. Must contain a "comparison_name" column.
    :type differential: pandas.DataFrame
    :param data_type: Data type. One of "ATAC-seq" and "RNA-seq". Defaults to "ATAC-seq".
    :type data_type: str, optional
    :param output_dir: Directory to create output files. Defaults to "results/differential_analysis_{data_type}".
    :type output_dir: str, optional
    :param output_prefix: Prefix to use when creating output files. Defaults to "differential_analysis".
    :type output_prefix: str, optional
    :param genome: Genome assembly of the analysis. Defaults to "hg19".
    :type genome: str, optional
    :param directional: Whether enrichments should be performed in a direction-dependent way
                        (up-regulated and down-regulated features separately).
                        This requires a column named "log2FoldChange" to exist. Defaults to True.
    :type directional: bool, optional
    :param max_diff: Number of maximum features to perform enrichment for ranked by variable in `max_diff`. Defaults to 1000.
    :type max_diff: number, optional
    :param sort_var: Variable to sort for when setting `max_diff`. Defaults to "pvalue".
    :type sort_var: str, optional
    :param as_jobs: One of "serial" or "job". Defaults to True.
    :type as_jobs: bool, optional
    """
    import pandas as pd
    from tqdm import tqdm

    serial = not as_jobs

    if data_type == "ATAC-seq":
        from ngs_toolkit.atacseq import characterize_regions_function
        matrix = analysis.coverage_annotated
        lola_enr = pd.DataFrame()
        meme_enr = pd.DataFrame()
        homer_enr = pd.DataFrame()
        pathway_enr = pd.DataFrame()
    elif data_type == "RNA-seq":
        from ngs_toolkit.general import enrichr
        matrix = analysis.expression
        pathway_enr = pd.DataFrame()
    else:
        raise AssertionError("`data_type` must match one of 'ATAC-seq' or 'RNA-seq'.")

    if "{data_type}" in output_dir:
        output_dir = output_dir.format(data_type=data_type)

    # Examine each region cluster
    max_diff = max_diff
    comps = differential['comparison_name'].drop_duplicates()
    for comp in tqdm(comps, total=len(comps)):
        if directional:
            # Separate in up/down-regulated genes
            params = [(np.less, 0, "down", "head"), (np.greater, 0, "up", "tail")]
        else:
            # get all genes
            params = [(np.less, np.inf, "all", "head")]

        for f, arg, direction, top in params:
            if directional:
                diff = differential.loc[
                    (differential["comparison_name"] == comp) &
                    (f(differential["log2FoldChange"], arg)), :].index
            else:
                diff = differential.loc[
                    (differential["comparison_name"] == comp), :].index

            # Handle extremes of regions
            if diff.shape[0] < 1:
                continue
            if diff.shape[0] > max_diff:
                if directional:
                    diff = (
                        getattr(
                            differential[
                                (differential["comparison_name"] == comp) &
                                (f(differential["log2FoldChange"], arg))]
                            [sort_var].sort_values(), top)
                        (max_diff).index)
                else:
                    diff = (
                        getattr(
                            differential[
                                (differential["comparison_name"] == comp)]
                            [sort_var].sort_values(), top)
                        (max_diff).index)

            # Add data_type specific info
            comparison_df = matrix.loc[diff, :]

            # Characterize
            comparison_dir = os.path.join(output_dir, "{}.{}".format(comp, direction))
            if not os.path.exists(comparison_dir):
                os.makedirs(comparison_dir)

            if data_type == "RNA-seq":
                print("Doing genes of comparison '{}', direction '{}'.".format(comp, direction))
                comparison_df.index.name = "gene_name"
                # write gene names to file
                (comparison_df
                 .reset_index()['gene_name']
                 .drop_duplicates()
                 .sort_values()
                 .to_csv(os.path.join(comparison_dir, output_prefix + ".gene_symbols.txt"), header=None, index=False))

                if serial:
                    if not os.path.exists(os.path.join(comparison_dir, output_prefix + ".enrichr.csv")):
                        enr = enrichr(comparison_df.reset_index())
                        enr.to_csv(os.path.join(comparison_dir,
                                                output_prefix + ".enrichr.csv"), index=False)
                    else:
                        enr = pd.read_csv(os.path.join(
                            comparison_dir, output_prefix + ".enrichr.csv"))
                        enr["comparison_name"] = comp
                        pathway_enr = pathway_enr.append(enr, ignore_index=True)
            else:
                print("Doing regions of comparison '{}', direction '{}'.".format(comp, direction))

                # do the suite of enrichment analysis
                characterize_regions_function(
                    analysis, comparison_df,
                    output_dir=comparison_dir, prefix=output_prefix, run=serial, genome=genome)

                # collect enrichments
                if serial:
                    # read/parse
                    meme_motifs = parse_ame(comparison_dir)
                    homer_motifs = parse_homer(os.path.join(comparison_dir, "homerResults"))
                    lola = pd.read_csv(os.path.join(comparison_dir, "allEnrichments.tsv"), sep="\t")
                    enr = pd.read_csv(os.path.join(comparison_dir, output_prefix +
                                                   "_regions.enrichr.csv"), index=False, encoding='utf-8')
                    # label
                    for d in [lola, enr, motifs]:
                        d["comparison_name"] = comp
                        d["direction"] = direction
                        d["label"] = "{}.{}".format(comp, direction)
                    # append
                    meme_enr = meme_enr.append(meme_motifs, ignore_index=True)
                    homer_enr = homer_enr.append(homer_motifs, ignore_index=True)
                    lola_enr = lola_enr.append(lola, ignore_index=True)
                    pathway_enr = pathway_enr.append(enr, ignore_index=True)

    if serial:
        # write combined enrichments
        if data_type == "ATAC-seq":
            meme_enr.to_csv(
                os.path.join(output_dir, output_prefix + ".meme_ame.csv"), index=False)
            homer_enr.to_csv(
                os.path.join(output_dir, output_prefix + ".homer_motifs.csv"), index=False)
            lola.to_csv(
                os.path.join(output_dir, output_prefix + ".lola.csv"), index=False)
        pathway_enr.to_csv(
            os.path.join(output_dir, output_prefix + ".enrichr.csv"), index=False)
    else:
        try:
            background = getattr(analysis, "sites").fn
        except:
            background = ""
        run_enrichment_jobs(
            analysis_name=analysis.name, results_dir=output_dir,
            genome=genome, background_bed=background)


def collect_differential_enrichment(
        differential,
        directional=True,
        data_type="ATAC-seq",
        output_dir="results/differential_analysis_{data_type}",
        output_prefix="differential_analysis",
        permissive=True):
    """
    Collect the results of enrichment analysis ran after a differential analysis.

    :param differential: Data frame with differential results as produced by
                             ``ngs_toolkit.general.differential_analysis``.
    :type differential: pandas.DataFrame
    :param directional: Whether enrichments were made in a direction-dependent way
                        (up-regulated and down-regulated features separately).
                        This implies a column named "direction" exists". Defaults to True.
    :type directional: bool, optional
    :param data_type: Data type. One of "ATAC-seq" and "RNA-seq". Defaults to "ATAC-seq".
    :type data_type: str, optional
    :param output_dir: Directory to create output files. Defaults to "results/differential_analysis_{data_type}".
    :type output_dir: str, optional
    :param output_prefix: Prefix to use when creating output files. Defaults to "differential_analysis".
    :type output_prefix: str, optional
    :param permissive: Whether to skip non-existing files, giving a warning. Defaults to True.
    :type permissive: bool, optional
    """
    import pandas as pd
    import numpy as np
    from ngs_toolkit.general import parse_ame, parse_homer
    from tqdm import tqdm

    if data_type not in ["ATAC-seq", "RNA-seq"]:
        raise AssertionError("`data_type` must match one of 'ATAC-seq' or 'RNA-seq'.")

    if "{data_type}" in output_dir:
        output_dir = output_dir.format(data_type=data_type)

    error_msg = "{} results for comparison '{}', direction '{}' were not found!"

    lola_enr = pd.DataFrame()
    meme_enr = pd.DataFrame()
    homer_enr = pd.DataFrame()
    pathway_enr = pd.DataFrame()
    # Examine each region cluster
    comps = differential['comparison_name'].drop_duplicates()
    for comp in tqdm(comps, total=len(comps)):
        if directional:
            # Separate in up/down-regulated genes
            params = list()
            if differential[(differential['comparison_name'] == comp) & (differential['log2FoldChange'] > 0)].shape[0] > 0:
                params.append('up')
            if differential[(differential['comparison_name'] == comp) & (differential['log2FoldChange'] < 0)].shape[0] > 0:
                params.append('down')
            if len(params) == 0:
                continue
        else:
            params = ["all"]

        for direction in params:
            comparison_dir = os.path.join(output_dir, "{}.{}".format(comp, direction))

            if data_type == "RNA-seq":
                # print("Collecting enrichments of comparison '{}', direction '{}'.".format(comp, direction))
                try:
                    enr = pd.read_csv(os.path.join(comparison_dir, output_prefix + ".enrichr.csv"))
                except IOError as e:
                    if permissive:
                        print(error_msg.format("Enrichr", comp, direction))
                    else:
                        raise e
                except pd.errors.EmptyDataError:
                    continue
                else:
                    enr["comparison_name"] = comp
                    enr["direction"] = direction
                    pathway_enr = pathway_enr.append(enr, ignore_index=True)
            elif data_type == "ATAC-seq":
                # print("Collecting enrichments of comparison '{}', direction '{}'.".format(comp, direction))
                # read/parse

                # MEME/AME
                try:
                    ame_motifs = parse_ame(comparison_dir).reset_index()
                    ame_motifs.columns = ["TF", "p_value"]
                except IOError as e:
                    if permissive:
                        print(error_msg.format("MEME/AME motif", comp, direction))
                    else:
                        raise e
                else:
                    ame_motifs["comparison_name"] = comp
                    ame_motifs["direction"] = direction
                    meme_enr = meme_enr.append(ame_motifs, ignore_index=True)

                # HOMER
                try:
                    homer_motifs = parse_homer(os.path.join(comparison_dir, "homerResults"))
                except IOError as e:
                    if permissive:
                        print(error_msg.format("HOMER motif", comp, direction))
                    else:
                        raise e
                else:
                    homer_motifs["comparison_name"] = comp
                    homer_motifs["direction"] = direction
                    homer_enr = homer_enr.append(homer_motifs, ignore_index=True)

                # LOLA
                try:
                    lola = pd.read_csv(os.path.join(comparison_dir, "allEnrichments.tsv"), sep="\t")
                except IOError as e:
                    if permissive:
                        print(error_msg.format("LOLA", comp, direction))
                    else:
                        raise e
                else:
                    lola["comparison_name"] = comp
                    lola["direction"] = direction
                    lola_enr = lola_enr.append(lola, ignore_index=True)

                # ENRICHR
                try:
                    enr = pd.read_csv(os.path.join(
                        comparison_dir, output_prefix + "_genes.enrichr.csv"))
                except IOError as e:
                    if permissive:
                        print(error_msg.format("Enrichr", comp, direction))
                    else:
                        raise e
                else:
                    enr["comparison_name"] = comp
                    enr["direction"] = direction
                    pathway_enr = pathway_enr.append(enr, ignore_index=True)

    # write combined enrichments
    pathway_enr.to_csv(
        os.path.join(output_dir, output_prefix + ".enrichr.csv"), index=False)
    if data_type == "ATAC-seq":
        meme_enr.to_csv(
            os.path.join(output_dir, output_prefix + ".meme_ame.csv"), index=False)
        homer_enr.to_csv(
            os.path.join(output_dir, output_prefix + ".homer_motifs.csv"), index=False)
        lola_enr.to_csv(
            os.path.join(output_dir, output_prefix + ".lola.csv"), index=False)


def plot_differential_enrichment(
        enrichment_table,
        enrichment_type,
        data_type="ATAC-seq",
        direction_dependent=True,
        output_dir="results/differential_analysis_{data_type}/enrichments",
        comp_variable="comparison_name",
        output_prefix="differential_analysis",
        to_rasterize=True,
        barplots=True,
        correlation_plots=True,
        top_n=5,
        z_score=0):
    """
    Given a table of enrichment terms across several comparisons, produced
    plots illustrating these enrichments in the various comparisons.

    # TODO: add plotting of genomic region and chromatin state enrichment.
    `comp_variable` is the column in the enrichment table that labels groups.

    :param enrichment_table: Data frame with enrichment results as produced by
                             ``ngs_toolkit.general.differential_enrichment`` or
                             ``ngs_toolkit.general.collect_differential_enrichment``.
    :type enrichment_table: pandas.DataFrame
    :param enrichment_type: One of 'lola', 'enrichr', 'motif', 'great'.
    :type enrichment_type: str
    :param data_type: Data type. One of "ATAC-seq" and "RNA-seq". Defaults to "ATAC-seq".
    :type data_type: str, optional
    :param direction_dependent: Whether enrichments were made in a direction-dependent way
                                (up-regulated and down-regulated features separately).
                                This implies a column named "direction" exists". Defaults to True.
    :type direction_dependent: bool, optional
    :param output_dir: Directory to create output files. Defaults to "results/differential_analysis_{data_type}/enrichments".
    :type output_dir: str, optional
    :param output_prefix: Prefix to use when creating output files. Defaults to "differential_analysis".
    :type output_prefix: str, optional
    :param comp_variable: Column defining which comparison enrichment terms belong to. Defaults to "comparison_name".
    :type comp_variable: str, optional
    :param barplots: Whether barplots with top enriched terms per comparison should be produced. Defaults to True.
    :type barplots: bool, optional
    :param correlation_plots: Whether correlation plots of comparisons across enriched terms should be produced. Defaults to True.
    :type correlation_plots: bool, optional
    :param top_n: Top terms to be used to make barplots. Defaults to 5
    :type top_n: number, optional
    :param z_score: Which dimention/axis to perform Z-score transformation for. Numpy/Pandas conventions are used:
                    `0` is row-wise (in this case across comparisons) and `1` is column-wise (across terms). Defaults to 0.
    :type z_score: number, optional
    :param to_rasterize: Whether or not to rasterize heatmaps for efficient plotting. Defaults to True.
    :type to_rasterize: bool, optional
    """

    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns

    if enrichment_type not in ["lola", "enrichr", "motif", 'great']:
        raise AssertionError("`enrichment_type` must be one of 'lola', 'enrichr', 'motif', 'great.")

    if "{data_type}" in output_dir:
        output_dir = output_dir.format(data_type=data_type)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if z_score == 0:
        z_score_label = "Row"
    elif z_score == 1:
        z_score_label = "Column"
    elif z_score is None:
        pass
    else:
        raise ValueError("Argument 'z_score' must be on of 0, 1 or None.")

    enrichment_table = enrichment_table.copy()
    if "direction" in enrichment_table.columns and direction_dependent:
        enrichment_table[comp_variable] = enrichment_table[comp_variable].astype(
            str) + " " + enrichment_table["direction"].astype(str)

    if enrichment_type == "lola":
        # get a unique label for each lola region set
        enrichment_table["label"] = (
            enrichment_table["description"].astype(str) + ", " +
            enrichment_table["cellType"].astype(str) + ", " +
            enrichment_table["tissue"].astype(str) + ", " +
            enrichment_table["antibody"].astype(str) + ", " +
            enrichment_table["treatment"].astype(str))
        enrichment_table["label"] = (
            enrichment_table["label"]
            .str.replace("nan", "").str.replace("None", "")
            .str.replace(", , ", "").str.replace(", $", ""))

        # Replace inf values with
        enrichment_table["pValueLog"] = enrichment_table["pValueLog"].replace(
            np.inf,
            enrichment_table.loc[enrichment_table["pValueLog"] != np.inf, "pValueLog"].max()
        )

        # Plot top_n terms of each comparison in barplots
        top_data = enrichment_table.set_index("label").groupby(
            comp_variable)["pValueLog"].nlargest(top_n).reset_index()

        if barplots:
            n = len(enrichment_table[comp_variable].drop_duplicates())
            n_side = int(np.ceil(np.sqrt(n)))

            fig, axis = plt.subplots(n_side, n_side, figsize=(
                4 * n_side, n_side * max(5, 0.12 * top_n)), sharex=False, sharey=False)
            axis = iter(axis.flatten())
            for i, comp in enumerate(top_data[comp_variable].drop_duplicates().sort_values()):
                df2 = top_data.loc[top_data[comp_variable] == comp, :]
                ax = next(axis)
                sns.barplot(
                    df2["pValueLog"], df2["label"], estimator=max,
                    orient="horizontal", ax=ax, color=sns.color_palette("colorblind")[0])
                ax.set_title(comp)
            for ax in axis:
                ax.set_visible(False)
            sns.despine(fig)
            fig.savefig(os.path.join(output_dir, output_prefix +
                                     ".lola.barplot.top_{}.svg".format(top_n)), bbox_inches="tight", dpi=300)

        # Plot heatmaps of terms for each comparison
        if len(enrichment_table[comp_variable].drop_duplicates()) < 2:
            return

        # pivot table
        lola_pivot = pd.pivot_table(enrichment_table,
                                    values="pValueLog", columns=comp_variable, index="label").fillna(0)
        lola_pivot = lola_pivot.replace(np.inf, lola_pivot[lola_pivot != np.inf].max().max())

        # plot correlation
        if correlation_plots:
            g = sns.clustermap(
                lola_pivot.corr(),
                rasterized=to_rasterize,
                xticklabels=True,
                yticklabels=True,
                cbar_kws={"label": "Correlation of enrichment\nof differential regions"})
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
            g.fig.savefig(os.path.join(output_dir, output_prefix +
                                       ".lola.correlation.svg"), bbox_inches="tight", dpi=300)

        top = enrichment_table.set_index('label').groupby(
            comp_variable)['pValueLog'].nlargest(top_n)
        top_terms = top.index.get_level_values('label').unique()

        # plot clustered heatmap
        shape = lola_pivot.loc[top_terms, :].shape
        g = sns.clustermap(
            lola_pivot.loc[top_terms, :], figsize=(
                max(6, 0.12 * shape[1]), max(6, 0.12 * shape[0])),
            xticklabels=True,
            yticklabels=True,
            rasterized=to_rasterize,
            cbar_kws={"label": "-log10(p-value) of enrichment\nof differential regions"}, metric="correlation")
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(),
                                     rotation=90, ha="right", fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(),
                                     rotation=0, fontsize="xx-small")
        g.fig.savefig(os.path.join(output_dir, output_prefix +
                                   ".lola.cluster_specific.svg"), bbox_inches="tight", dpi=300)

        # plot clustered heatmap
        if z_score is not None:
            g = sns.clustermap(
                lola_pivot.loc[top_terms, :], figsize=(
                    max(6, 0.12 * shape[1]), max(6, 0.12 * shape[0])),
                xticklabels=True,
                yticklabels=True,
                rasterized=to_rasterize,
                cmap="RdBu_r", center=0, z_score=z_score, cbar_kws={"label": "{} Z-score of enrichment\nof differential regions".format(z_score_label)}, metric="correlation")
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(),
                                         rotation=90, ha="right", fontsize="xx-small")
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(),
                                         rotation=0, fontsize="xx-small")
            g.fig.savefig(os.path.join(output_dir, output_prefix +
                                       ".lola.cluster_specific.{}_z_score.svg".format(z_score_label)), bbox_inches="tight", dpi=300)

    if enrichment_type == "motif":
        enrichment_table["log_p_value"] = (-np.log10(enrichment_table["p_value"])
                                           ).replace({np.inf: 300})
        # Plot top_n terms of each comparison in barplots
        top_data = enrichment_table.set_index("TF").groupby(
            comp_variable)["log_p_value"].nlargest(top_n).reset_index()

        if barplots:
            n = len(enrichment_table[comp_variable].drop_duplicates())
            n_side = int(np.ceil(np.sqrt(n)))

            fig, axis = plt.subplots(n_side, n_side, figsize=(
                4 * n_side, n_side * max(5, 0.12 * top_n)), sharex=False, sharey=False)
            axis = iter(axis.flatten())
            for i, comp in enumerate(top_data[comp_variable].drop_duplicates().sort_values()):
                df2 = top_data.loc[top_data[comp_variable] == comp, :]
                ax = next(axis)
                sns.barplot(
                    df2["log_p_value"], df2["TF"], estimator=max,
                    orient="horizontal", ax=ax, color=sns.color_palette("colorblind")[0])
                ax.set_title(comp)
            for ax in axis:
                ax.set_visible(False)
            sns.despine(fig)
            fig.savefig(os.path.join(output_dir, output_prefix +
                                     ".motifs.barplot.top_{}.svg".format(top_n)), bbox_inches="tight", dpi=300)

        # Plot heatmaps of terms for each comparison
        if len(enrichment_table[comp_variable].drop_duplicates()) < 2:
            return

        # pivot table
        motifs_pivot = pd.pivot_table(enrichment_table,
                                      values="log_p_value", columns="TF", index=comp_variable).fillna(0)

        # plot correlation
        if correlation_plots:
            g = sns.clustermap(motifs_pivot.T.corr(),
                               rasterized=to_rasterize,
                               xticklabels=True,
                               yticklabels=True,
                               cbar_kws={"label": "Correlation of enrichment\nof differential regions"})
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
            g.fig.savefig(os.path.join(output_dir, output_prefix +
                                       ".motifs.correlation.svg"), bbox_inches="tight", dpi=300)

        top = enrichment_table.set_index('TF').groupby(comp_variable)['log_p_value'].nlargest(top_n)
        top_terms = top.index.get_level_values('TF').unique()

        # plot clustered heatmap
        shape = motifs_pivot.loc[:, top_terms].shape
        g = sns.clustermap(
            motifs_pivot.loc[:, top_terms].T, figsize=(
                max(6, 0.12 * shape[0]), max(6, 0.12 * shape[1])),
            xticklabels=True, yticklabels=True, rasterized=to_rasterize,
            cbar_kws={"label": "-log10(p-value) of enrichment\nof differential regions"}, metric="correlation")
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(),
                                     rotation=90, ha="right", fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(),
                                     rotation=0, fontsize="xx-small")
        g.fig.savefig(os.path.join(output_dir, output_prefix +
                                   ".motifs.cluster_specific.svg"), bbox_inches="tight", dpi=300)

        # plot clustered heatmap
        shape = motifs_pivot.loc[:, top_terms].shape
        if z_score is not None:
            g = sns.clustermap(
                motifs_pivot.loc[:, top_terms].T, figsize=(
                    max(6, 0.12 * shape[0]), max(6, 0.12 * shape[1])),
                xticklabels=True, yticklabels=True, rasterized=to_rasterize,
                cmap="RdBu_r", center=0, z_score=z_score, cbar_kws={"label":  "{} Z-score of enrichment\nof differential regions".format(z_score_label)}, metric="correlation")
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(),
                                         rotation=90, ha="right", fontsize="xx-small")
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(),
                                         rotation=0, fontsize="xx-small")
            g.fig.savefig(os.path.join(output_dir, output_prefix +
                                       ".motifs.cluster_specific.{}_z_score.svg".format(z_score_label)), bbox_inches="tight", dpi=300)

    if enrichment_type == "enrichr":
        # enrichment_table["description"] = enrichment_table["description"].str.decode("utf-8")
        enrichment_table["log_p_value"] = (-np.log10(enrichment_table["p_value"])
                                           ).replace({np.inf: 300})

        for gene_set_library in enrichment_table["gene_set_library"].unique():
            print(gene_set_library)
            if gene_set_library == "Epigenomics_Roadmap_HM_ChIP-seq":
                continue

            # Plot top_n terms of each comparison in barplots
            top_data = (
                enrichment_table[enrichment_table["gene_set_library"] == gene_set_library]
                .set_index("description")
                .groupby(comp_variable)
                ["log_p_value"]
                .nlargest(top_n)
                .reset_index())

            n = len(enrichment_table[comp_variable].drop_duplicates())
            n_side = int(np.ceil(np.sqrt(n)))

            if barplots:
                fig, axis = plt.subplots(n_side, n_side, figsize=(
                    4 * n_side, n_side * max(5, 0.12 * top_n)), sharex=False, sharey=False)
                axis = iter(axis.flatten())
                for i, comp in enumerate(top_data[comp_variable].drop_duplicates().sort_values()):
                    df2 = top_data.loc[top_data[comp_variable] == comp, :]
                    ax = next(axis)
                    sns.barplot(
                        df2["log_p_value"], df2["description"], estimator=max,
                        orient="horizontal", ax=ax, color=sns.color_palette("colorblind")[0])
                    ax.set_title(comp)
                for ax in axis:
                    ax.set_visible(False)
                sns.despine(fig)
                fig.savefig(os.path.join(output_dir, output_prefix + ".enrichr.{}.barplot.top_{}.svg".format(
                    gene_set_library, top_n)), bbox_inches="tight", dpi=300)

            # Plot heatmaps of terms for each comparison
            if len(enrichment_table[comp_variable].drop_duplicates()) < 2:
                return

            # pivot table
            enrichr_pivot = pd.pivot_table(
                enrichment_table[enrichment_table["gene_set_library"] == gene_set_library],
                values="log_p_value", columns="description", index=comp_variable).fillna(0)

            # plot correlation
            if correlation_plots:
                g = sns.clustermap(enrichr_pivot.T.corr(),
                                   rasterized=to_rasterize,
                                   xticklabels=True,
                                   yticklabels=True,
                                   cbar_kws={"label": "Correlation of enrichment\nof differential genes"})
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
                g.fig.savefig(os.path.join(output_dir, output_prefix +
                                           ".enrichr.{}.correlation.svg".format(gene_set_library)), bbox_inches="tight", dpi=300)

            top = enrichment_table[enrichment_table["gene_set_library"] == gene_set_library].set_index(
                'description').groupby(comp_variable)['p_value'].nsmallest(top_n)
            top_terms = top.index.get_level_values('description').unique()
            # top_terms = top_terms[top_terms.isin(lola_pivot.columns[lola_pivot.sum() > 5])]

            # plot clustered heatmap
            shape = enrichr_pivot[list(set(top_terms))].shape
            g = sns.clustermap(
                enrichr_pivot[list(set(top_terms))].T, figsize=(
                    max(6, 0.12 * shape[0]), max(6, 0.12 * shape[1])),
                xticklabels=True, yticklabels=True, rasterized=to_rasterize,
                cbar_kws={"label": "-log10(p-value) of enrichment\nof differential genes"}, vmin=0)
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(),
                                         rotation=90, ha="right", fontsize="xx-small")
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(),
                                         rotation=0, fontsize="xx-small")
            g.fig.savefig(os.path.join(output_dir, output_prefix +
                                       ".enrichr.{}.cluster_specific.svg".format(gene_set_library)), bbox_inches="tight", dpi=300)

            # plot clustered heatmap
            shape = enrichr_pivot[list(set(top_terms))].shape
            if z_score is not None:
                g = sns.clustermap(
                    enrichr_pivot[list(set(top_terms))].T, figsize=(
                        max(6, 0.12 * shape[0]), max(6, 0.12 * shape[1])),
                    xticklabels=True, yticklabels=True, rasterized=to_rasterize,
                    cmap="RdBu_r", center=0, z_score=z_score, cbar_kws={"label":  "{} Z-score of enrichment\nof differential regions".format(z_score_label)})
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(
                ), rotation=90, ha="right", fontsize="xx-small")
                g.ax_heatmap.set_yticklabels(
                    g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
                g.fig.savefig(os.path.join(output_dir, output_prefix + ".enrichr.{}.cluster_specific.{}_z_score.svg".format(
                    gene_set_library, z_score_label)), bbox_inches="tight", dpi=300)

    if enrichment_type == "great":
        # enrichment_table["description"] = enrichment_table["description"].str.decode("utf-8")
        enrichment_table["log_q_value"] = (-np.log10(enrichment_table["HyperFdrQ"])
                                           ).replace({np.inf: 300})

        for gene_set_library in enrichment_table["Ontology"].unique():
            print(gene_set_library)

            # Plot top_n terms of each comparison in barplots
            top_data = (
                enrichment_table[enrichment_table["Ontology"] == gene_set_library]
                .set_index("Desc")
                .groupby(comp_variable)
                ["log_q_value"]
                .nlargest(top_n)
                .reset_index())

            n = len(enrichment_table[comp_variable].drop_duplicates())
            n_side = int(np.ceil(np.sqrt(n)))

            if barplots:
                fig, axis = plt.subplots(n_side, n_side, figsize=(
                    4 * n_side, n_side * max(5, 0.12 * top_n)), sharex=False, sharey=False)
                axis = iter(axis.flatten())
                for i, comp in enumerate(top_data[comp_variable].drop_duplicates().sort_values()):
                    df2 = top_data.loc[top_data[comp_variable] == comp, :]
                    ax = next(axis)
                    sns.barplot(
                        df2["log_q_value"], df2["Desc"], estimator=max,
                        orient="horizontal", ax=ax, color=sns.color_palette("colorblind")[0])
                    ax.set_title(comp)
                for ax in axis:
                    ax.set_visible(False)
                sns.despine(fig)
                fig.savefig(os.path.join(output_dir, output_prefix + ".great.{}.barplot.top_{}.svg".format(
                    gene_set_library, top_n)), bbox_inches="tight", dpi=300)

            # Plot heatmaps of terms for each comparison
            if len(enrichment_table[comp_variable].drop_duplicates()) < 2:
                return

            # pivot table
            great_pivot = pd.pivot_table(
                enrichment_table[enrichment_table["Ontology"] == gene_set_library],
                values="log_q_value", columns="Desc", index=comp_variable).fillna(0)

            # plot correlation
            if correlation_plots:
                try:
                    g = sns.clustermap(great_pivot.T.corr(),
                                       rasterized=to_rasterize,
                                       xticklabels=True,
                                       yticklabels=True,
                                       cbar_kws={"label": "Correlation of enrichment\nof differential genes"})
                    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
                    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
                    g.fig.savefig(os.path.join(output_dir, output_prefix +
                                               ".great.{}.correlation.svg".format(gene_set_library)), bbox_inches="tight", dpi=300)
                except FloatingPointError:
                    continue

            top = enrichment_table[enrichment_table["Ontology"] == gene_set_library].set_index(
                'Desc').groupby(comp_variable)['HyperP'].nsmallest(top_n)
            top_terms = top.index.get_level_values('Desc').unique()
            # top_terms = top_terms[top_terms.isin(lola_pivot.columns[lola_pivot.sum() > 5])]

            # plot clustered heatmap
            shape = great_pivot[list(set(top_terms))].shape
            g = sns.clustermap(
                great_pivot[list(set(top_terms))].T, figsize=(
                    max(6, 0.12 * shape[0]), max(6, 0.12 * shape[1])),
                xticklabels=True, yticklabels=True, rasterized=to_rasterize,
                cbar_kws={"label": "-log10(p-value) of enrichment\nof differential genes"}, metric="correlation")
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(),
                                         rotation=90, ha="right", fontsize="xx-small")
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(),
                                         rotation=0, fontsize="xx-small")
            g.fig.savefig(os.path.join(output_dir, output_prefix +
                                       ".great.{}.cluster_specific.svg".format(gene_set_library)), bbox_inches="tight", dpi=300)

            # plot clustered heatmap
            shape = great_pivot[list(set(top_terms))].shape
            if z_score is not None:
                g = sns.clustermap(
                    great_pivot[list(set(top_terms))].T, figsize=(
                        max(6, 0.12 * shape[0]), max(6, 0.12 * shape[1])),
                    xticklabels=True, yticklabels=True, rasterized=to_rasterize,
                    cmap="RdBu_r", center=0, z_score=z_score, cbar_kws={"label":  "{} Z-score of enrichment\nof differential regions".format(z_score_label)}, metric="correlation")
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(
                ), rotation=90, ha="right", fontsize="xx-small")
                g.ax_heatmap.set_yticklabels(
                    g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
                g.fig.savefig(os.path.join(output_dir, output_prefix + ".great.{}.cluster_specific.{}_z_score.svg".format(
                    gene_set_library, z_score_label)), bbox_inches="tight", dpi=300)


def chunks(l, n):
    """
    Partition iterable in `n` chunks.

    :param iterable l: Iterable (e.g. list or numpy array).
    :param int n: Number of chunks to generate.
    """
    n = max(1, n)
    return (l[i:i + n] for i in range(0, len(l), n))


def standard_score(x):
    """
    Compute a standard score, defined as (x - min(x)) / (max(x) - min(x)).

    :param numpy.array x: Numeric array.
    """
    return (x - x.min()) / (x.max() - x.min())


def z_score(x):
    """
    Compute a Z-score, defined as (x - mean(x)) / std(x).

    :param numpy.array x: Numeric array.
    """
    return (x - x.mean()) / x.std()


def signed_max(x, f=0.66):
    """
    Return maximum or minimum depending on the sign of the majority of values.
    If there isn't a clear majority (at least `f` fraction in one side), return mean of values.

    :param numpy.array x: Numeric array.
    :param float f: Threshold fraction of majority agreement.
    """
    l = float(len(x))
    neg = sum(x < 0)
    pos = sum(x > 0)
    obs_f = max(neg / l, pos / l)
    if obs_f >= f:
        if neg > pos:
            return min(x)
        else:
            return max(x)
    else:
        return np.mean(x)


# def detect_peaks(x, mph=None, mpd=1, threshold=0, edge='rising',
#                  kpsh=False, valley=False):
#     """Detect peaks in data based on their amplitude and other features.

#     Parameters
#     ----------
#     x : 1D array_like
#         data.
#     mph : {None, number}, optional (default = None)
#         detect peaks that are greater than minimum peak height.
#     mpd : positive integer, optional (default = 1)
#         detect peaks that are at least separated by minimum peak distance (in
#         number of data).
#     threshold : positive number, optional (default = 0)
#         detect peaks (valleys) that are greater (smaller) than `threshold`
#         in relation to their immediate neighbors.
#     edge : {None, 'rising', 'falling', 'both'}, optional (default = 'rising')
#         for a flat peak, keep only the rising edge ('rising'), only the
#         falling edge ('falling'), both edges ('both'), or don't detect a
#         flat peak (None).
#     kpsh : bool, optional (default = False)
#         keep peaks with same height even if they are closer than `mpd`.
#     valley : bool, optional (default = False)
#         if True (1), detect valleys (local minima) instead of peaks.
#     show : bool, optional (default = False)
#         if True (1), plot data in matplotlib figure.
#     ax : a matplotlib.axes.Axes instance, optional (default = None).

#     Returns
#     -------
#     ind : 1D array_like
#         indeces of the peaks in `x`.

#     Notes
#     -----
#     The detection of valleys instead of peaks is performed internally by simply
#     negating the data: `ind_valleys = detect_peaks(-x)`

#     The function can handle NaN's

#     See this IPython Notebook [1]_.

#     References
#     ----------
#     .. [1] http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb

#     Examples
#     --------
#     >>> from detect_peaks import detect_peaks
#     >>> x = np.random.randn(100)
#     >>> x[60:81] = np.nan
#     >>> # detect all peaks and plot data
#     >>> ind = detect_peaks(x, show=True)
#     >>> print(ind)

#     >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
#     >>> # set minimum peak height = 0 and minimum peak distance = 20
#     >>> detect_peaks(x, mph=0, mpd=20, show=True)

#     >>> x = [0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0]
#     >>> # set minimum peak distance = 2
#     >>> detect_peaks(x, mpd=2, show=True)

#     >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
#     >>> # detection of valleys instead of peaks
#     >>> detect_peaks(x, mph=0, mpd=20, valley=True, show=True)

#     >>> x = [0, 1, 1, 0, 1, 1, 0]
#     >>> # detect both edges
#     >>> detect_peaks(x, edge='both', show=True)

#     >>> x = [-2, 1, -2, 2, 1, 1, 3, 0]
#     >>> # set threshold = 2
#     >>> detect_peaks(x, threshold = 2, show=True)
#     """
#     import numpy as np

#     x = np.atleast_1d(x).astype('float64')
#     if x.size < 3:
#         return np.array([], dtype=int)
#     if valley:
#         x = -x
#     # find indices of all peaks
#     dx = x[1:] - x[:-1]
#     # handle NaN's
#     indnan = np.where(np.isnan(x))[0]
#     if indnan.size:
#         x[indnan] = np.inf
#         dx[np.where(np.isnan(dx))[0]] = np.inf
#     ine, ire, ife = np.array([[], [], []], dtype=int)
#     if not edge:
#         ine = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) > 0))[0]
#     else:
#         if edge.lower() in ['rising', 'both']:
#             ire = np.where((np.hstack((dx, 0)) <= 0) & (np.hstack((0, dx)) > 0))[0]
#         if edge.lower() in ['falling', 'both']:
#             ife = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) >= 0))[0]
#     ind = np.unique(np.hstack((ine, ire, ife)))
#     # handle NaN's
#     if ind.size and indnan.size:
#         # NaN's and values close to NaN's cannot be peaks
#         ind = ind[np.in1d(ind, np.unique(np.hstack((indnan, indnan - 1, indnan + 1))), invert=True)]
#     # first and last values of x cannot be peaks
#     if ind.size and ind[0] == 0:
#         ind = ind[1:]
#     if ind.size and ind[-1] == x.size - 1:
#         ind = ind[:-1]
#     # remove peaks < minimum peak height
#     if ind.size and mph is not None:
#         ind = ind[x[ind] >= mph]
#     # remove peaks - neighbors < threshold
#     if ind.size and threshold > 0:
#         dx = np.min(np.vstack([x[ind] - x[ind - 1], x[ind] - x[ind + 1]]), axis=0)
#         ind = np.delete(ind, np.where(dx < threshold)[0])
#     # detect small peaks closer than minimum peak distance
#     if ind.size and mpd > 1:
#         ind = ind[np.argsort(x[ind])][::-1]  # sort ind by peak height
#         idel = np.zeros(ind.size, dtype=bool)
#         for i in range(ind.size):
#             if not idel[i]:
#                 # keep peaks with the same height if kpsh is True
#                 idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
#                     & (x[ind[i]] > x[ind] if kpsh else True)
#                 idel[i] = 0  # Keep current peak
#         # remove the small peaks and sort back the indices by their occurrence
#         ind = np.sort(ind[~idel])

#     return ind


def count_jobs_running(cmd="squeue", sep="\n"):
    """
    Count running jobs on a cluster by invoquing a command that lists the jobs.
    """
    import subprocess
    return subprocess.check_output(cmd).split(sep).__len__()


def submit_job_if_possible(cmd, total_job_lim=800, refresh_time=10, in_between_time=5):
    from ngs_toolkit.general import count_jobs_running
    import time
    import os

    submit = count_jobs_running() < total_job_lim
    while not submit:
        time.sleep(refresh_time)
        submit = count_jobs_running() < total_job_lim
    os.system(cmd)
    time.sleep(in_between_time)


def sra_id2geo_id(sra_ids):
    """Query SRA ID from GEO ID"""
    import subprocess

    cmd = """esearch -db sra -query {} | efetch -format docsum | xtract -pattern DocumentSummary -element Runs |  perl -ne '@mt = ($_ =~ /SRR\d+/g); print "@mt"'"""

    geo_ids = list()
    for id in sra_ids:
        p, err = subprocess.Popen(cmd.format(id))
        geo_ids.append(p.communicate())
    return


def sra2fastq(input_sra, output_dir):
    cmd = """
\t\tfastq-dump --split-3 --outdir {} {}
    """.format(output_dir, input_sra)

    return cmd


def fastq2bam(input_fastq, output_bam, sample_name, input_fastq2=None):
    cmd = """
\t\tjava -Xmx4g -jar /cm/shared/apps/picard-tools/1.118/FastqToSam.jar"""
    cmd += " FASTQ={0}".format(input_fastq)
    cmd += " SAMPLE_NAME={0}".format(sample_name)
    if input_fastq2 is not None:
        cmd += " FASTQ2={0}".format(input_fastq2)
    cmd += """ OUTPUT={0}""".format(output_bam)

    return cmd


def download_cram(link, output_dir):
    cmd = """
    cd {}
    wget '{}'
    cd -
    """.format(output_dir, link)

    return cmd


def cram2bam(input_cram, output_bam):
    cmd = """
    samtools view -b -o {} {}
    """.format(output_bam, input_cram)

    return cmd


def download_sra(link, output_dir):
    cmd = """
    cd {}
    wget '{}'
    cd -
    """.format(output_dir, link)

    return cmd


def sra2bam_job(sra_id, base_path):
    from pypiper import NGSTk
    import textwrap
    tk = NGSTk()

    # Slurm header
    job_file = os.path.join(base_path, "%s_sra2bam.sh" % sra_id)
    log_file = os.path.join(base_path, "%s_sra2bam.log" % sra_id)

    cmd = tk.slurm_header("-".join(["sra2bam", sra_id]), log_file, cpus_per_task=2)

    # SRA to FASTQ
    cmd += sra2fastq(os.path.join(base_path, sra_id + ".sra"), base_path)

    # FASTQ to BAM
    cmd += fastq2bam(
        os.path.join(base_path, sra_id + "_1.fastq"),
        os.path.join(base_path, sra_id + ".bam"),
        sra_id,
        os.path.join(base_path, sra_id + "_2.fastq"))

    # Slurm footer
    cmd += tk.slurm_footer() + "\n"

    # Write job to file

    with open(job_file, "w") as handle:
        handle.write(textwrap.dedent(cmd))

    # Submit
    tk.slurm_submit_job(job_file)


def link2bam_job(sample_name, link, base_path):
    from pypiper import NGSTk
    import textwrap
    tk = NGSTk()

    # Slurm header
    job_file = os.path.join(base_path, "%s_link2bam.sh" % sample_name)
    log_file = os.path.join(base_path, "%s_link2bam.log" % sample_name)

    cmd = tk.slurm_header("-".join(["link2bam", sample_name]), log_file, cpus_per_task=2)

    # Download CRAM
    cmd += download_cram(
        link,
        base_path)

    # CRAM to BAM
    cmd += cram2bam(
        os.path.join(base_path, sample_name + ".cram"),
        os.path.join(base_path, sample_name + ".bam"))

    # Slurm footer
    cmd += tk.slurm_footer() + "\n"

    # Write job to file

    with open(job_file, "w") as handle:
        handle.write(textwrap.dedent(cmd))

    # Submit
    tk.slurm_submit_job(job_file)


def sralink2bam_job(sra_id, base_path):
    from pypiper import NGSTk
    import textwrap
    tk = NGSTk()

    # Slurm header
    job_file = os.path.join(base_path, "%s_sra2bam.sh" % sra_id)
    log_file = os.path.join(base_path, "%s_sra2bam.log" % sra_id)

    cmd = tk.slurm_header("-".join(["sra2bam", sra_id]), log_file, cpus_per_task=2)

    # SRA to FASTQ
    cmd += sra2fastq(sra_id, base_path)

    # FASTQ to BAM
    cmd += fastq2bam(
        os.path.join(base_path, sra_id + "_1.fastq"),
        os.path.join(base_path, sra_id + ".bam"),
        sra_id,
        os.path.join(base_path, sra_id + "_2.fastq"))

    # Slurm footer
    cmd += tk.slurm_footer() + "\n"

    # Write job to file

    with open(job_file, "w") as handle:
        handle.write(textwrap.dedent(cmd))

    # Submit
    tk.slurm_submit_job(job_file)
    print(job_file)


def series_matrix2csv(matrix_url, prefix=None):
    """
    matrix_url: gziped URL with GEO series matrix.
    """
    import gzip
    import pandas as pd

    os.system("wget {}".format(matrix_url))
    filename = matrix_url.split("/")[-1]

    with gzip.open(filename, 'rb') as f:
        file_content = f.read()

    # separate lines with only one field (project-related)
    # from lines with >2 fields (sample-related)
    prj_lines = dict()
    sample_lines = dict()

    for line in file_content.decode("utf-8").strip().split("\n"):
        line = line.strip().split("\t")
        if len(line) == 2:
            prj_lines[line[0].replace("\"", "")] = line[1].replace("\"", "")
        elif len(line) > 2:
            sample_lines[line[0].replace("\"", "")] = [x.replace("\"", "") for x in line[1:]]

    prj = pd.Series(prj_lines)
    prj.index = prj.index.str.replace("!Series_", "")

    samples = pd.DataFrame(sample_lines)
    samples.columns = samples.columns.str.replace("!Sample_", "")

    if prefix is not None:
        prj.to_csv(os.path.join(prefix + ".project_annotation.csv"), index=True)
        samples.to_csv(os.path.join(prefix + ".sample_annotation.csv"), index=False)

    return prj, samples


def project_to_geo(
        project, output_dir="geo_submission",
        samples=None, distributed=False, dry_run=False):
    """
    Prepare raw sequencing files for submission to GEO.
    Files will be copied or generated in a new directory `output_dir`.
    It will get the raw BAM file(s) of each sample, and in case of ATAC-seq/ChIP-seq samples,
    the bigWig and peak files. If multiple BAM files exist for each sample, all will be copied and sequencially named with the "file\d" suffix,
    where "\d" is the file number.
    For each copied file a md5sum will be calculated.
    A pandas DataFrame with info on the sample's files and md5sums will be returned.

    :param project: A looper Project object to process.
    :type project: looper.models.Project
    :param output_dir: Directory to create output. Will be created/overwriten if existing.
                       Defaults to "geo_submission".
    :type output_dir: str, optional
    :param samples: List of looper.models.Sample objects in project to restrict to. Defaults to all samples in project.
    :type samples: list, optional
    :param distributed: Whether processing should be distributed as jobs in a computing cluster for each sample.
                        Currently available implementation supports a SLURM cluster only. Defaults to False.
    :type distributed: bool, optional
    :param dry_run: Whether copy/execution/submisison commands should be not be run to test.
    :type dry_run: bool, optional
    :returns: pandas.DataFrame with annotation of samples and their BAM, BigWig, narrowPeak files and respective md5sums.
    """
    output_dir = os.path.abspath(output_dir)
    if samples is None:
        samples = project.samples
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    cmds = list()
    annot = pd.DataFrame()
    for sample in samples:
        various = len(sample.data_path.split(" ")) > 1
        cmd = ""
        for i, file in enumerate(sample.data_path.split(" ")):
            suffix = ".file{}".format(i) if various else ""
            # Copy raw file
            bam_file = os.path.join(output_dir, sample.name + "{}.bam".format(suffix))
            cmd += "cp {} {}; ".format(file, bam_file)
            annot.loc[sample.name, "bam_file{}".format(i)] = bam_file

            # Copy or generate md5sum
            md5_file = os.path.join(output_dir, sample.name + "{}.bam.md5".format(suffix))
            if os.path.exists(file + ".md5"):
                cmd += "cp {} {}; ".format(file + ".md5", md5_file)
            else:
                b = os.path.basename(file)
                cmd += "md5sum {} > {}; ".format(os.path.join(output_dir, b), md5_file)
            annot.loc[sample.name, "bam_file{}_md5sum".format(i)] = md5_file

        # Copy bigWig files
        if sample.library in ["ATAC-seq", "ChIP-seq"]:
            if hasattr(sample, "bigwig"):
                bigwig_file = os.path.join(output_dir, sample.name + ".bigWig")
                cmd += "cp {} {}; ".format(sample.bigwig, bigwig_file)
                annot.loc[sample.name, "bigwig_file"] = bigwig_file

                # Copy or generate md5sum
                md5_file = os.path.join(output_dir, sample.name + "{}.bigWig.md5".format(suffix))
                if os.path.exists(sample.bigwig + ".md5"):
                    cmd += "cp {} {}; ".format(sample.bigwig + ".md5", md5_file)
                else:
                    b = os.path.basename(sample.bigwig)
                    cmd += "md5sum {} > {}; ".format(os.path.join(output_dir, b), md5_file)
                annot.loc[sample.name, "bigwig_file_md5sum"] = md5_file
            else:
                print("'{}' sample '{}' does not have a 'bigwig' attribute set. Skipping bigWig file.".format(
                    sample.library, sample.name))

        # Copy peaks
        if sample.library == "ATAC-seq":
            if hasattr(sample, "peaks"):
                peaks_file = os.path.join(output_dir, sample.name + ".peaks.narrowPeak")
                cmd += "cp {} {}; ".format(sample.peaks, peaks_file)
                annot.loc[sample.name, "peaks_file"] = peaks_file

                # Copy or generate md5sum
                md5_file = os.path.join(output_dir, sample.name +
                                        "{}.peaks.narrowPeak.md5".format(suffix))
                if os.path.exists(sample.peaks + ".md5"):
                    cmd += "cp {} {}; ".format(sample.peaks + ".md5", md5_file)
                else:
                    b = os.path.basename(sample.peaks)
                    cmd += "md5sum {} > {}; ".format(os.path.join(output_dir, b), md5_file)
                annot.loc[sample.name, "peaks_file_md5sum"] = md5_file
            else:
                print("'{}' sample '{}' does not have a 'peaks' attribute set. Skipping peaks file.".format(
                    sample.library, sample.name))

        if distributed and not dry_run:
            from pypiper.ngstk import NGSTk
            import textwrap
            tk = NGSTk()

            job_name = "project_to_geo.{}".format(sample.name)
            log_file = os.path.join(output_dir, job_name + ".log")
            job_file = os.path.join(output_dir, job_name + ".sh")

            job = textwrap.dedent(tk.slurm_header(
                job_name=job_name, output=log_file,
                cpus_per_task=1, mem_per_cpu=8000))
            job += "\n" + "\n".join(cmd.split("; ")) + "\n"
            job += textwrap.dedent(tk.slurm_footer())

            with open(job_file, "w") as handle:
                handle.write(textwrap.dedent(job))
            tk.slurm_submit_job(job_file)
        else:
            cmds.append(cmd)

    if not distributed and not dry_run:
        for i, cmd in enumerate(cmds):
            print(i, cmd)
            os.system(cmd)

    return annot


def rename_sample_files(
        annotation_mapping,
        old_sample_name_column="old_sample_name",
        new_sample_name_column="new_sample_name",
        tmp_prefix="rename_sample_files",
        results_dir="results_pipeline",
        dry_run=False):
    """
    Rename existing directories with pipeline outputs for samples based on mapping of
    old/new sample names.
    All files within the directory with the old sample name will be renamed recursively.
    Old and new sample names can overlap - this procedure will handle these cases correctly
    by a 2-step process with temporary sample names with prefix `prefix`.

    NEEDS TESTING!

    :param annotation_mapping: DataFrame with mapping of
    old (column "previous_sample_name") vs new ("new_sample_name") sample names.
    :type annotation_mapping: pandas.DataFrame
    :param old_sample_name_column: Name of column with old sample names.
                        Defaults to "old_sample_name"
    :type old_sample_name_column: str, optional
    :param new_sample_name_column: Name of column with new sample names.
                        Defaults to "new_sample_name"
    :type new_sample_name_column: str, optional
    :param tmp_prefix: Prefix for temporary files to avoid overlap between old and new names.
                        Defaults to "rename_sample_files"
    :type tmp_prefix: str, optional
    :param results_dir: Pipeline output directory containing sample output directories.
                        Defaults to "results_pipeline"
    :type results_dir: str, optional
    :param dry_run: Whether to print commands instead of running them. Defaults to False
    :type dry_run: bool, optional

    :returns: None
    """
    from subprocess import call
    cmds = list()
    # 1) move to tmp name
    for i, series in annotation_mapping.iterrows():
        o = series[old_sample_name_column]
        t = "{}-{}".format(tmp_prefix, i)
        cmds.append("# Moving old sample '{}' to '{}'.".format(
            o, t))

        # directory
        cmds.append("mv {} {}".format(
            os.path.join(results_dir, o), os.path.join(results_dir, t)))
        # further directories
        cmds.append("find {} -type d -exec rename {} {} {{}} \;".format(
            os.path.join(results_dir, t), o, t))
        # files
        cmds.append("find {} -type f -exec rename {} {} {{}} \;".format(
            os.path.join(results_dir, t), o, t))

    # 2) move to new name
    for i, series in annotation_mapping.iterrows():
        t = "{}-{}".format(tmp_prefix, i)
        n = series[new_sample_name_column]
        cmds.append("# Moving old sample '{}' to '{}'.".format(
            t, n))

        # directory
        cmds.append("mv {} {}".format(
            os.path.join(results_dir, t), os.path.join(results_dir, n)))
        # further directories
        cmds.append("find {} -type d -exec rename {} {} {{}} \;".format(
            os.path.join(results_dir, n), t, n))
        # files
        cmds.append("find {} -type f -exec rename {} {} {{}} \;".format(
            os.path.join(results_dir, n), t, n))

    if dry_run:
        print("\n".join(cmds))
    else:
        for i, cmd in enumerate(cmds):
            print(i, cmd)
            if cmd.startswith("#"):
                continue
            try:
                r = call(cmd, shell=True)
            except OSError as e:
                raise e
            if r != 0:
                raise OSError("Command '{}' failed.".format(cmd))


def query_biomart(
        attributes=["ensembl_gene_id", "external_gene_name", "hgnc_id", "hgnc_symbol"],
        species="hsapiens", ensembl_version="grch37"):
    """
    Query Biomart for gene attributes.
    Returns pandas dataframe with query results.
    """
    import requests
    import pandas as pd
    import numpy as np

    # Build request XML
    url_query = "".join([
        """http://{}.ensembl.org/biomart/martservice?query=""".format(ensembl_version),
        """<?xml version="1.0" encoding="UTF-8"?>""",
        """<!DOCTYPE Query>""",
        """<Query  virtualSchemaName="default" formatter="CSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6" >""",
        """<Dataset name="{}_gene_ensembl" interface="default" >""".format(species)] +
        ["""<Attribute name="{}" />""".format(attr) for attr in attributes] +
        ["""</Dataset>""",
         """</Query>"""])
    req = requests.get(url_query, stream=True)
    content = list(req.iter_lines())
    if type(content[0]) == bytes:
        mapping = pd.DataFrame((x.decode("utf-8").strip().split(",")
                                for x in content), columns=attributes)
    else:
        mapping = pd.DataFrame((x.strip().split(",") for x in content), columns=attributes)
    return mapping.replace("", np.nan)


def subtract_principal_component(
        X, pc=1, norm=False, plot=True, plot_name="PCA_based_batch_correction.svg", pcs_to_plot=6):
    """
    Given a matrix (n_samples, n_variables), remove `pc` (1-based) from matrix.
    """
    import numpy as np
    from sklearn.decomposition import PCA

    pc -= 1

    # All regions
    if norm:
        from sklearn.preprocessing import StandardScaler
        X = StandardScaler().fit_transform(X)

    # PCA
    pca = PCA()
    X_hat = pca.fit_transform(X)

    # Remove PC
    X2 = X - np.outer(X_hat[:, pc], pca.components_[pc, :])

    # plot
    if plot:
        import matplotlib.pyplot as plt
        import seaborn as sns
        X2_hat = pca.fit_transform(X2)
        fig, axis = plt.subplots(pcs_to_plot, 2, figsize=(4 * 2, 4 * pcs_to_plot))
        for pc in range(pcs_to_plot):
            # before
            for j, sample in enumerate(X.index):
                axis[pc, 0].scatter(X_hat[j, pc], X_hat[j, pc + 1], s=50)
            axis[pc, 0].set_xlabel("PC{}".format(pc + 1))
            axis[pc, 0].set_ylabel("PC{}".format(pc + 2))
            # after
            for j, sample in enumerate(X2.index):
                axis[pc, 1].scatter(X2_hat[j, pc], X2_hat[j, pc + 1], s=35, alpha=0.8)
            axis[pc, 1].set_xlabel("PC{}".format(pc + 1))
            axis[pc, 1].set_ylabel("PC{}".format(pc + 2))
        fig.savefig(plot_name)

    return X2


def subtract_principal_component_by_attribute(df, pc=1, attributes=["CLL"]):
    """
    Given a matrix (n_samples, n_variables), remove `pc` (1-based) from matrix.
    """
    import numpy as np
    from sklearn.decomposition import PCA

    pc -= 1

    X2 = pd.DataFrame(index=df.index, columns=df.columns)
    for attr in attributes:
        print(attr)
        sel = df.index[df.index.str.contains(attr)]
        X = df.loc[sel, :]

        # PCA
        pca = PCA()
        X_hat = pca.fit_transform(X)

        # Remove PC
        X2.loc[sel, :] = X - np.outer(X_hat[:, pc], pca.components_[pc, :])
    for sample in df.index:
        if X2.loc[sample, :].isnull().all():
            X2.loc[sample, :] = df.loc[sample, :]
    return X2
