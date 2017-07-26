#!/usr/bin/env python

import os

import numpy as np
import pandas as pd


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
    Class to hold functions and data from analysis.
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
            pickle_file = os.path.join(results_dir, "analysis.{}.pickle".format(name))
        self.pickle_file = pickle_file

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
        pass

    def from_pickle(self):
        import cPickle as pickle
        return pickle.load(open(self.pickle_file, 'rb'))

    def update(self):
        self.__dict__.update(self.from_pickle().__dict__)


def count_reads_in_intervals(bam, intervals):
    """
    Counts total number of reads in a iterable holding strings
    representing genomic intervals of the type chrom:start-end.
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
    import numpy as np

    # install package
    # R
    # source('http://bioconductor.org/biocLite.R')
    # biocLite('preprocessCore')

    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()

    robjects.r('require("preprocessCore")')
    normq = robjects.r('normalize.quantiles')
    return np.array(normq(array))


def normalize_quantiles_p(df_input):
    df = df_input.copy()
    # compute rank
    dic = {}
    for col in df:
        dic.update({col : sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis = 1).tolist()
    # sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df


def deseq_analysis(
        count_matrix, experiment_matrix, comparison_table, formula,
        output_dir, output_prefix,
        overwrite=True, alpha=0.05):
    """
    Perform differential comparisons with DESeq2.
    """
    import pandas as pd
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
    experiment_matrix.to_csv(os.path.join(output_dir, output_prefix + ".experiment_matrix.tsv"), sep="\t")
    comparison_table.to_csv(os.path.join(output_dir, output_prefix + ".comparison_table.tsv"), sep="\t")

    # Run DESeq analysis
    dds = _DESeqDataSetFromMatrix(
        countData=count_matrix.astype(int),
        colData=experiment_matrix,
        design=_as_formula(formula))
    dds = _DESeq(dds, parallel=True)
    # _save(dds, file=os.path.join(output_dir, output_prefix + ".deseq_dds_object.Rdata"))

    results = pd.DataFrame()
    for comp in comparison_table["comparison_name"].drop_duplicates().sort_values():
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
        res = _as_data_frame(_results(dds, contrast=contrast, alpha=alpha, independentFiltering=False, parallel=True))

        # convert to pandas dataframe
        res2 = r2pandas_df(res)
        res2.loc[:, "comparison_name"] = comp

        # save
        res2.to_csv(out_file)
        # append
        results = results.append(res2.reset_index(), ignore_index=True)

    # save all
    results.to_csv(os.path.join(output_dir, output_prefix + ".deseq_result.all_comparisons.csv"), index=False)

    # return
    return results


def differential_analysis(
        analysis,
        comparison_table,
        data_type="ATAC-seq",
        samples=None,
        covariates=None,
        output_dir="results/differential_analysis_{data_type}",
        output_prefix="differential_analysis",
        alpha=0.05,
        overwrite=True):
    """
    Discover differential regions/genes across samples that are associated with a certain trait.

    `comparison_table` is a dataframe with 'comparison_name', 'comparison_side' and 'sample_name', 'sample_group' columns.
    """
    # Check comparisons
    # check comparison table has required columns
    req_attrs = ['comparison_name', 'comparison_side', 'sample_name', 'sample_group']
    if not all([x in comparison_table.columns for x in req_attrs]):
        raise AssertionError("Given comparison table does not have all of '{}' columns.".format("', '".join(req_attrs)))
    # check all comparisons have samples in two sides
    if not all(comparison_table.groupby("comparison_name")["comparison_side"].nunique() == 2):
        raise AssertionError("All comparisons must have samples in each side of the comparison.")
    # check if any comparison and sample group has samples disagreeing in their side
    if not all(comparison_table.groupby(['comparison_name', "sample_group"])['comparison_side'].nunique() == 1):
        raise AssertionError("Samples in same comparison and group must agree on their side of the comparison.")

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
        (s.name in count_matrix.columns)
    ]
    count_matrix = count_matrix[[s.name for s in samples]]

    # Get experiment matrix
    # by getting any other relevant covariates as required
    if covariates is not None:
        sample_table = pd.DataFrame([s.as_series() for s in samples])
        # check all covariates are in the samples and none is null
        if not all([x in sample_table.columns for x in covariates]):
            raise AssertionError("Not all of the specified covariate variables are in the selected samples.")
        if sample_table[covariates].isnull().any().any():
            raise AssertionError("None of the selected samples can have a Null value in the specified covariate variables.")
    
        # add covariates to comparison table
        comparison_table = comparison_table.set_index("sample_name").join(sample_table.set_index("sample_name")[covariates]).reset_index()

    # Make table for DESeq2
    experiment_matrix = comparison_table[
        ["sample_name", "sample_group"] + (covariates if covariates is not None else [])
    ].drop_duplicates()

    # Make formula for DESeq2
    formula = "~ {}sample_group".format(" + ".join(covariates) + " + " if covariates is not None else "")

    # Run DESeq2 analysis
    results = deseq_analysis(
        count_matrix, experiment_matrix, comparison_table,
        formula, output_dir, output_prefix, alpha=alpha, overwrite=overwrite)

    return results


def differential_overlap(
        differential,
        data_type="ATAC-seq",
        output_dir="results/differential_analysis_{data_type}",
        output_prefix="differential_analysis"
    ):
    """
    Visualize intersection of sets of differential regions/genes.
    """
    import numpy as np
    import itertools
    import matplotlib.pyplot as plt
    import matplotlib
    import seaborn as sns

    if "{data_type}" in output_dir:
        output_dir = output_dir.format(data_type=data_type)

    if data_type == "ATAC-seq":
        unit = "region"
    elif data_type == "RNA-seq":
        unit = "gene"

    if "direction" not in differential.columns:
        differential["direction"] = differential["log2FoldChange"].apply(lambda x: "up" if x > 0 else "down")

    differential.index.name = "index"
    differential["intersect"] = 1
    piv = pd.pivot_table(differential.reset_index(), index='index', columns=['comparison_name', 'direction'], values='intersect', fill_value=0)

    intersections = pd.DataFrame(columns=["group1", "group2", "dir1", "dir2", "size1", "size2", "intersection", "union"])
    for ((k1, dir1), i1), ((k2, dir2), i2) in itertools.permutations(piv.T.groupby(level=['comparison_name', 'direction']).groups.items(), 2):
        print(k1, k2)
        i1 = set(piv[i1][piv[i1] == 1].dropna().index)
        i2 = set(piv[i2][piv[i2] == 1].dropna().index)
        intersections = intersections.append(
            pd.Series(
                [k1, k2, dir1, dir2, len(i1), len(i2), len(i1.intersection(i2)), len(i1.union(i2))],
                index=["group1", "group2", "dir1", "dir2", "size1", "size2", "intersection", "union"]
            ),
            ignore_index=True
        )
    # convert to %
    intersections['intersection'] = intersections['intersection'].astype(float)
    intersections['intersection_perc'] = ((intersections['intersection'] / intersections['size2']) * 100.).astype(float)

    # save
    intersections.to_csv(os.path.join(output_dir, output_prefix + ".differential_overlap.csv"), index=False)

    # make pivot tables
    piv_up = pd.pivot_table(
        intersections[(intersections['dir1'] == "up") & (intersections['dir2'] == "up")],
        index="group1", columns="group2", values="intersection_perc").fillna(0)
    piv_down = pd.pivot_table(
        intersections[(intersections['dir1'] == "down") & (intersections['dir2'] == "down")],
        index="group1", columns="group2", values="intersection_perc").fillna(0)
    np.fill_diagonal(piv_up.values, np.nan)
    np.fill_diagonal(piv_down.values, np.nan)

    # heatmaps
    fig, axis = plt.subplots(1, 2, figsize=(8 * 2, 8))
    sns.heatmap(piv_down, square=True, cmap="summer", cbar_kws={"label": "Concordant {}s (% of group 2)".format(unit)}, ax=axis[0])
    sns.heatmap(piv_up, square=True, cmap="summer", cbar_kws={"label": "Concordant {}s (% of group 2)".format(unit)}, ax=axis[1])
    axis[0].set_title("Downregulated {}s".format(unit))
    axis[1].set_title("Upregulated {}s".format(unit))
    axis[0].set_xticklabels(axis[0].get_xticklabels(), rotation=90, ha="right")
    axis[1].set_xticklabels(axis[1].get_xticklabels(), rotation=90, ha="right")
    axis[0].set_yticklabels(axis[0].get_yticklabels(), rotation=0, ha="right")
    axis[1].set_yticklabels(axis[1].get_yticklabels(), rotation=0, ha="right")
    fig.savefig(os.path.join(output_dir, output_prefix + ".differential_overlap.up_down_split.svg"), bbox_inches="tight")

    # combined heatmap
    # with upregulated {}s in upper square matrix and downredulated in down square
    piv_combined = pd.DataFrame(np.triu(piv_up), index=piv_up.index, columns=piv_up.columns).replace(0, np.nan)
    piv_combined.update(pd.DataFrame(np.tril(piv_down), index=piv_down.index, columns=piv_down.columns).replace(0, np.nan))
    piv_combined = piv_combined.fillna(0)
    np.fill_diagonal(piv_combined.values, np.nan)

    fig, axis = plt.subplots(1, figsize=(8, 8))
    sns.heatmap(piv_combined, square=True, cmap="summer", cbar_kws={"label": "Concordant {}s (% of group 2)".format(unit)}, ax=axis)
    axis.set_xticklabels(axis.get_xticklabels(), rotation=90, ha="right")
    axis.set_yticklabels(axis.get_yticklabels(), rotation=0, ha="right")
    fig.savefig(os.path.join(output_dir, output_prefix + ".differential_overlap.up_down_together.svg"), bbox_inches="tight")

    # Observe disagreement between knockouts
    # (overlap of down-regulated with up-regulated and vice-versa)
    piv_up = pd.pivot_table(
        intersections[(intersections['dir1'] == "up") & (intersections['dir2'] == "down")],
        index="group1", columns="group2", values="intersection")
    piv_down = pd.pivot_table(
        intersections[(intersections['dir1'] == "down") & (intersections['dir2'] == "up")],
        index="group1", columns="group2", values="intersection")

    piv_disagree = pd.concat([piv_up, piv_down]).groupby(level=0).max()

    fig, axis = plt.subplots(1, 2, figsize=(16, 8))
    sns.heatmap(piv_disagree, square=True, cmap="Greens", cbar_kws={"label": "Discordant {}s".format(unit).format(unit)}, ax=axis[0])
    sns.heatmap(np.log2(1 + piv_disagree), square=True, cmap="Greens", cbar_kws={"label": "Discordant {}s (log2)".format(unit)}, ax=axis[1])

    norm = matplotlib.colors.Normalize(vmin=0, vmax=piv_disagree.max().max())
    cmap = plt.get_cmap("Greens")
    log_norm = matplotlib.colors.Normalize(vmin=0, vmax=np.log2(1 + piv_disagree).max().max())
    for i, g1 in enumerate(piv_disagree.columns):
        for j, g2 in enumerate(piv_disagree.index):
            axis[0].scatter(
                len(piv_disagree.index) - (j + 0.5), len(piv_disagree.index) - (i + 0.5),
                s=(100 ** (norm(piv_disagree.loc[g1, g2]))) - 1, color=cmap(norm(piv_disagree.loc[g1, g2])), marker="o")
            axis[1].scatter(
                len(piv_disagree.index) - (j + 0.5), len(piv_disagree.index) - (i + 0.5),
                s=(100 ** (log_norm(np.log2(1 + piv_disagree).loc[g1, g2]))) - 1, color=cmap(log_norm(np.log2(1 + piv_disagree).loc[g1, g2])), marker="o")

    for ax in axis:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="right")
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0, ha="right")
    fig.savefig(os.path.join(output_dir, output_prefix + ".differential_overlap.disagreement.svg"), bbox_inches="tight")


def plot_differential(
        analysis,
        results,
        comparison_table=None,
        samples=None,
        data_type="ATAC-seq",
        alpha=0.05,
        corrected_p_value=True,
        fold_change=None,
        output_dir="results/differential_analysis_{data_type}",
        output_prefix="differential_analysis"):
    """
    Discover differential regions across samples that are associated with a certain trait.
    The `results` matrix should be indexed by the relevant type of feature (regions/genes).
    """
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns

    req_attrs = ["log2FoldChange", "pvalue", "padj", "comparison_name"]
    if not all([x in results.columns for x in req_attrs]):
        raise AssertionError("Results dataframe must have '{}' columns.".format(", ".join(req_attrs)))

    # Extract significant based on p-value and fold-change
    if fold_change is not None:
        fc = (abs(results["log2FoldChange"]) > fold_change)
    else:
        fc = [True] * results.shape[0]
    if corrected_p_value:
        p_var = "padj"
    else:
        p_var = "pvalue"
    results.loc[(results[p_var] < alpha) & fc, "diff"] = True
    results.loc[:, "diff"] = results.loc[:, "diff"].fillna(False)

    # Annotate direction of change
    results.loc[:, "direction"] = results.loc[:, "log2FoldChange"].apply(lambda x: "up" if x >= 0 else "down")

    if results.loc[:, "diff"].sum() < 1:
        print("No significantly different regions found in any comparison.")
        return

    # Make output dir
    if "{data_type}" in output_dir:
        output_dir = output_dir.format(data_type=data_type)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get matrix and samples
    if data_type == "ATAC-seq":
        matrix = analysis.accessibility
        results.index.name = "region"
        var_name = "region"
        quantity = "Accessibility"
        unit = "RPM"
    elif data_type == "RNA-seq":
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
    if comparison_table is not None:
        samples = [s for s in samples if s.name in comparison_table['sample_name'].tolist()]
    matrix = matrix[[s.name for s in samples]]

    # PLOTS
    comparisons = sorted(results["comparison_name"].drop_duplicates())
    n_side = int(np.ceil(np.sqrt(len(comparisons))))

    # P-value distributions
    fig, axis = plt.subplots(1, 1, figsize=(4, 4))
    sns.distplot(results["pvalue"].dropna(), kde=False, ax=axis)
    axis.set_xlabel("P-value")
    axis.set_ylabel("Frequency")
    sns.despine(fig)
    fig.savefig(os.path.join(output_dir, output_prefix + ".pvalue_distribution.svg"), bbox_inches="tight")
    # per comparison
    g = sns.FacetGrid(data=results, col="comparison_name", col_wrap=n_side)
    g.map(sns.distplot, "pvalue", kde=False)
    for ax in g.axes:
        ax.set_xlabel("P-value")
        ax.set_ylabel("Frequency")
    sns.despine(g.fig)
    g.fig.savefig(os.path.join(output_dir, output_prefix + ".pvalue_distribution.per_comparison.svg"), bbox_inches="tight")

    # Number of differential vars
    n_vars = float(matrix.shape[0])
    total_diff = results.groupby(["comparison_name"])['diff'].sum().sort_values(ascending=False)
    split_diff = results.groupby(["comparison_name", "direction"])['diff'].sum().sort_values(ascending=False)
    fig, axis = plt.subplots(2, 2, figsize=(4 * 2, 4* 2))
    sns.barplot(total_diff.values, total_diff.index, orient="h", ax=axis[0][0])
    sns.barplot((total_diff.values / n_vars) * 100, total_diff.index, orient="h", ax=axis[0][1])
    sns.barplot(split_diff.values, split_diff.index, orient="h", ax=axis[1][0])
    sns.barplot((split_diff.values / n_vars) * 100, split_diff.index, orient="h", ax=axis[1][1])
    axis[0][0].set_xlabel("N. diff")
    axis[0][1].set_xlabel("N. diff (% of total)")
    axis[1][0].set_xlabel("N. diff")
    axis[1][1].set_xlabel("N. diff (% of total)")
    fig.savefig(os.path.join(output_dir, output_prefix + ".number_differential.svg"), bbox_inches="tight")

    # Pairwise scatter plots
    if comparison_table is not None:
        fig, axes = plt.subplots(n_side, n_side, figsize=(n_side * 4, n_side * 4), sharex=True, sharey=True)
        if n_side > 1 or n_side > 1:
            axes = iter(axes.flatten())
        else:
            axes = iter([axes])
        for comparison in comparisons:
            c = comparison_table[comparison_table["comparison_name"] == comparison]
            a = c.loc[c['comparison_side'] >= 1, "sample_name"]
            b = c.loc[c['comparison_side'] <= 0, "sample_name"]

            a = matrix[[s.name for s in samples if s.name in a.tolist() and s.library == data_type]].mean(axis=1)
            b = matrix[[s.name for s in samples if s.name in b.tolist() and s.library == data_type]].mean(axis=1)

            # Hexbin plot
            ax = axes.next()
            ax.hexbin(b, a, alpha=0.85, color="black", edgecolors="white", linewidths=0, bins='log', mincnt=1, rasterized=True)

            # Scatter for significant
            diff_vars = results[(results["comparison_name"] == comparison) & (results["diff"] == True)].index
            if diff_vars.shape[0] > 0:
                ax.scatter(b.loc[diff_vars], a.loc[diff_vars], alpha=0.1, color="red", s=2)
            ax.set_title(comparison)
            ax.set_xlabel("Down")
            ax.set_ylabel("Up")
            # x = y square
            lims = [np.min([ax.get_xlim(), ax.get_ylim()]), np.max([ax.get_xlim(), ax.get_ylim()])]
            ax.plot(lims, lims, linestyle='--', alpha=0.5, zorder=0, color="black")
            ax.set_aspect('equal')
            ax.set_xlim(lims)
            ax.set_ylim(lims)
        for ax in axes:
            ax.set_visible(False)
        sns.despine(fig)
        fig.savefig(os.path.join(output_dir, output_prefix + ".scatter_plots.svg"), bbox_inches="tight", dpi=300)

    # Volcano plots
    fig, axes = plt.subplots(n_side, n_side, figsize=(n_side * 4, n_side * 4), sharex=False, sharey=False)
    if n_side > 1 or n_side > 1:
        axes = iter(axes.flatten())
    else:
        axes = iter([axes])
    for comparison in comparisons:
        t = results.loc[results["comparison_name"] == comparison, :]

        # Hexbin plot
        ax = axes.next()
        ax.hexbin(
            t["log2FoldChange"], -np.log10(t["pvalue"]),
            alpha=0.85, color="black", edgecolors="white", linewidths=0, bins='log', mincnt=1, rasterized=True)

        # Scatter for significant
        diff_vars = t[t["diff"] == True].index
        if diff_vars.shape[0] > 0:
            ax.scatter(t.loc[diff_vars, "log2FoldChange"], -np.log10(t.loc[diff_vars, "pvalue"]), alpha=0.1, color="red", s=2)
        ax.set_title(comparison)
        ax.set_xlabel("log2(Fold-change)")
        ax.set_ylabel("-log10(P-value)")
        ax.axvline(0, linestyle='--', alpha=0.5, zorder=0, color="black")
        l = np.max([abs(i) for i in ax.get_xlim()])
        ax.set_xlim(-l, l)

        # Add lines of significance
        ax.axhline(-np.log10(t.loc[t["diff"] == True, p_var].max()), linestyle='--', alpha=0.5, zorder=0, color="black")
        if fold_change is not None:
            ax.axvline(-fold_change, linestyle='--', alpha=0.5, zorder=0, color="black")
            ax.axvline(fold_change, linestyle='--', alpha=0.5, zorder=0, color="black")
    for ax in axes:
        ax.set_visible(False)
    sns.despine(fig)
    fig.savefig(os.path.join(output_dir, output_prefix + ".volcano_plots.svg"), bbox_inches="tight", dpi=300)

    # MA plots
    fig, axes = plt.subplots(n_side, n_side, figsize=(n_side * 4, n_side * 4), sharex=False, sharey=False)
    if n_side > 1 or n_side > 1:
        axes = iter(axes.flatten())
    else:
        axes = iter([axes])
    for comparison in comparisons:
        t = results.loc[results["comparison_name"] == comparison, :]

        # Hexbin plot
        ax = axes.next()
        ax.hexbin(
            np.log10(t["baseMean"]), t["log2FoldChange"],
            alpha=0.85, color="black", edgecolors="white", linewidths=0, bins='log', mincnt=1, rasterized=True)

        # Scatter for significant
        diff_vars = t[t["diff"] == True].index
        if diff_vars.shape[0] > 0:
            ax.scatter(np.log10(t.loc[diff_vars, "baseMean"]), t.loc[diff_vars, "log2FoldChange"], alpha=0.1, color="red", s=2)
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
    fig.savefig(os.path.join(output_dir, output_prefix + ".ma_plots.svg"), bbox_inches="tight", dpi=300)

    # Observe values of variables across all comparisons
    all_diff = results[results["diff"] == True].index.drop_duplicates()
    if type(matrix.columns) is pd.MultiIndex:
        sample_cols = matrix.columns.get_level_values("sample_name").tolist()
    else:
        sample_cols = matrix.columns.tolist()

    if comparison_table is not None:
        if results["comparison_name"].drop_duplicates().shape[0] > 1:
            groups = pd.DataFrame()
            for sample_group in results["comparison_name"].drop_duplicates():
                c = comparison_table.loc[comparison_table["sample_group"] == sample_group, "sample_name"].drop_duplicates()
                if c.shape[0] > 0:
                    groups.loc[:, sample_group] = matrix[[d for d in c if d in sample_cols]].mean(axis=1)

            # Heatmaps
            # Comparison level
            g = sns.clustermap(
                groups.corr(),
                xticklabels=False, cbar_kws={"label": "Pearson correlation\non differential {}".format(var_name)},
                cmap="BuGn", metric="correlation", rasterized=True)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
            g.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.groups.clustermap.corr.svg".format(var_name)), bbox_inches="tight", dpi=300, metric="correlation")

            g = sns.clustermap(
                groups.loc[all_diff, :],
                yticklabels=False, cbar_kws={"label": "{} of\ndifferential {}".format(quantity, var_name)},
                metric="correlation", rasterized=True)
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
            g.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.groups.clustermap.svg".format(var_name)), bbox_inches="tight", dpi=300)

            g = sns.clustermap(
                groups.loc[all_diff, :],
                yticklabels=False, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}".format(quantity, var_name)},
                metric="correlation", rasterized=True)
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
            g.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.groups.clustermap.z0.svg".format(var_name)), bbox_inches="tight", dpi=300)

    # Fold-changes and P-values
    # pivot table of genes vs comparisons
    fold_changes = pd.pivot_table(results.loc[all_diff, :].reset_index(), index=results.index.name, columns="comparison_name", values="log2FoldChange")
    p_values = -np.log10(pd.pivot_table(results.loc[all_diff, :].reset_index(), index=results.index.name, columns="comparison_name", values="padj"))

    # fold
    if fold_changes.shape[1] > 1:
        g = sns.clustermap(fold_changes.corr(),
            xticklabels=False, cbar_kws={"label": "Pearson correlation\non fold-changes"},
            cmap="BuGn", vmin=0, vmax=1, metric="correlation", rasterized=True)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.groups.fold_changes.clustermap.corr.svg".format(var_name)), bbox_inches="tight", dpi=300, metric="correlation")

        g = sns.clustermap(fold_changes.loc[all_diff, :],
            yticklabels=False, cbar_kws={"label": "Fold-change of\ndifferential {}".format(var_name)},
            robust=True, metric="correlation", rasterized=True)
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.groups.fold_changes.clustermap.svg".format(var_name)), bbox_inches="tight", dpi=300)

    # Sample level
    if type(matrix.columns) is pd.core.indexes.multi.MultiIndex:
        matrix.columns = matrix.columns.get_level_values("sample_name")

    g = sns.clustermap(matrix.loc[all_diff, :].corr(),
        yticklabels=True, xticklabels=False,
        cbar_kws={"label": "Pearson correlation\non differential {}".format(var_name)},
        cmap="BuGn", metric="correlation", rasterized=True)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.samples.clustermap.corr.svg".format(var_name)), bbox_inches="tight", dpi=300)

    g = sns.clustermap(matrix.loc[all_diff, :],
        yticklabels=False, cbar_kws={"label": "{} of\ndifferential {}".format(quantity, var_name)},
        xticklabels=True,
        vmin=0, metric="correlation", rasterized=True)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.samples.clustermap.svg".format(var_name)), bbox_inches="tight", dpi=300)

    g = sns.clustermap(matrix.loc[all_diff, :],
        yticklabels=False, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}".format(quantity, var_name)},
        xticklabels=True,
        metric="correlation", rasterized=True)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.fig.savefig(os.path.join(output_dir, output_prefix + ".diff_{}.samples.clustermap.z0.svg".format(var_name)), bbox_inches="tight", dpi=300)


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


def enrichr(dataframe, gene_set_libraries=None, kind="genes"):
    """
    Use Enrichr on a list of genes (currently only genes supported through the API).
    """
    import json
    import requests
    import pandas as pd

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
    for gene_set_library in gene_set_libraries:
        print("Using enricher on %s gene set library." % gene_set_library)

        if kind == "genes":
            # Build payload with bed file
            attr = "\n".join(dataframe["gene_name"].dropna().tolist())
        elif kind == "regions":
            # Build payload with bed file
            attr = "\n".join(dataframe[['chrom', 'start', 'end']].apply(lambda x: "\t".join([str(i) for i in x]), axis=1).tolist())

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
            res.columns = ["rank", "description", "p_value", "z_score", "combined_score", "genes", "adjusted_p_value"]
        elif len(res.columns) == 9:
            res.columns = ["rank", "description", "p_value", "z_score", "combined_score", "genes", "adjusted_p_value", "old_p_value", "old_adjusted_p_value"]

        # Remember gene set library used
        res["gene_set_library"] = gene_set_library

        # Append to master dataframe
        results = results.append(res, ignore_index=True)

    return results


def run_enrichment_jobs(analysis_name, genome):
    """
    Submit enrichment jobs for a specifc analysis.
    """
    # LOLA
    cmds = ["""for F in `find results -name "*_regions.bed"`; do
DIR=`dirname $F`
if [ ! -f ${{DIR}}/allEnrichments.txt ]; then
echo $DIR $F
sbatch -J lola.$F -o $F.lola.log -p shortq -c 8 --mem 24000 \
--wrap "Rscript ~/jobs/run_LOLA.R $F results/{PROJECT_NAME}_peak_set.bed {GENOME}"
fi
done""".format(PROJECT_NAME=analysis_name, GENOME=genome)]

    # AME
    dbs = {
        "human": "~/resources/motifs/motif_databases/HUMAN/HOCOMOCOv10.meme",
        "mouse": "~/resources/motifs/motif_databases/MOUSE/uniprobe_mouse.meme"}
    omap = {"hg38": "human", "hg19": "human", "mm10": "mouse"}

    cmds += ["""for F in `find results -name "*_regions.fa"`; do
DIR=`dirname $F`
if [ ! -f ${{DIR}}/ame.html ]; then
echo $DIR $F
sbatch -J "meme_ame.${{F}}" -o "${{F}}.meme_ame.log" -p shortq -c 1 --mem 4000 \
--wrap "fasta-dinucleotide-shuffle -c 1 -f "${{F}}" > "${{F}}".shuffled.fa; \
ame --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 0.05 \
--control "${{F}}".shuffled.fa -o "${{DIR}}" "${{F}}" {motifs}"
fi
done""".format(motifs=dbs[omap[genome]])]

    # HOMER
    cmds += ["""for F in `find results -name "*_regions.bed"`; do
DIR=`dirname $F`
if [ ! -f ${{DIR}}/homerResults.html ]; then
echo $DIR $F
sbatch -J "homer.${{F}}" -o "${{F}}.homer.log" -p shortq -c 8 --mem 20000 \
--wrap "findMotifsGenome.pl ${{F}} {GENOME}r ${{DIR}} -size 1000 -h -p 2 -len 8,10,12,14 -noknown"
fi
done""".format(GENOME=genome)]

    # Enrichr
    cmds += ["""for F in `find results -name "*.gene_symbols.txt"`; do
if [ ! -f ${{F}}.enrichr.csv ]; then
echo $F
sbatch -J enrichr.$F -o $F.enrichr.log -p shortq -c 1 --mem 4000 \
--wrap "python ~/jobs/run_Enrichr.py --input-file "$F" --output-file "${F/gene_symbols.txt/enrichr.csv}" "
fi
done"""]
    cmds += ["""for F in `find results -name "*_genes.symbols.txt"`; do
if [ ! -f ${F/symbols.txt/enrichr.csv} ]; then
echo $F
sbatch -J enrichr.$F -o $F.enrichr.log -p shortq -c 1 --mem 4000 \
--wrap "python ~/jobs/run_Enrichr.py --input-file "$F" --output-file "${F/symbols.txt/enrichr.csv}" "
fi
done"""]

    for cmd in cmds:
        os.system(cmd)


def differential_enrichment(
        analysis,
        differential,
        data_type="ATAC-seq",
        output_dir="results/differential_analysis_{data_type}",
        output_prefix="differential_analysis",
        genome="hg19",
        max_diff=1000,
        sort_var="pvalue",
        as_jobs=True
    ):
    """
    Given a dataframe of the results of differential analysis
    (containing "comparison_name" and "log2FoldChange" columns),
    will get enrichment of gene sets, regions and TF motifs according
    to the `data_type`.

    At most will use `max_diff` regions/genes in each comparison sorted by `sort_var`.
    
    `run_mode`: one of "serial" or "job".
    """
    import pandas as pd

    serial = not as_jobs

    if data_type == "ATAC-seq":
        from ngs_toolkit.atacseq import characterize_regions_function
        matrix = analysis.coverage_annotated
        lola_enr = pd.DataFrame()
        meme_enr = pd.DataFrame()
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
    max_diff = 1000
    for comp in differential['comparison_name'].drop_duplicates():
        # Separate in up/down-regulated genes
        for f, direction, top in [(np.less, "down", "head"), (np.greater, "up", "tail")]:
            diff = differential.loc[
                (differential["comparison_name"] == comp) &
                (f(differential["log2FoldChange"], 0)), :].index

            # Handle extremes of regions
            if diff.shape[0] < 1:
                continue
            if diff.shape[0] > max_diff:
                diff = (
                    getattr(
                        differential[
                            (differential["comparison_name"] == comp) &
                            (f(differential["log2FoldChange"], 0))]
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
                .reset_index()[['gene_name']]
                .drop_duplicates()
                .sort_values()
                .to_csv(os.path.join(comparison_dir, output_prefix + ".gene_symbols.txt"), header=None, index=False))

                if serial:
                    if not os.path.exists(os.path.join(comparison_dir, output_prefix + ".enrichr.csv")):
                        enr = enrichr(comparison_df.reset_index())
                        enr.to_csv(os.path.join(comparison_dir, output_prefix + ".enrichr.csv"), index=False)
                    else:
                        enr = pd.read_csv(os.path.join(comparison_dir, output_prefix + ".enrichr.csv"))
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
                    motifs = parse_ame(comparison_dir)
                    lola = pd.read_csv(os.path.join(comparison_dir, "allEnrichments.txt"), sep="\t")
                    enr = pd.read_csv(os.path.join(comparison_dir, output_prefix + "_regions.enrichr.csv"), index=False, encoding='utf-8')
                    # label
                    for d in [lola, enr, motifs]:
                        d["comparison_name"] = comp
                        d["direction"] = direction
                        d["label"] = "{}.{}".format(comp, direction)
                    # append
                    meme_enr = meme_enr.append(motifs, ignore_index=True)
                    lola_enr = lola_enr.append(lola, ignore_index=True)
                    pathway_enr = pathway_enr.append(enr, ignore_index=True)

    if serial:
        # write combined enrichments
        if data_type == "ATAC-seq":
            meme_enr.to_csv(
                os.path.join(output_dir, output_prefix + ".meme_ame.csv"), index=False)
            lola.to_csv(
                os.path.join(output_dir, output_prefix + ".lola.csv"), index=False)
        pathway_enr.to_csv(
            os.path.join(output_dir, output_prefix + ".enrichr.csv"), index=False)
    else:
        run_enrichment_jobs(analysis_name=analysis.name, genome=genome)


def collect_differential_enrichment(
        differential,
        data_type="ATAC-seq",
        output_dir="results/differential_analysis_{data_type}",
        output_prefix="differential_analysis",
        permissive=True):
    """
    Collect the results of enrichment analysis ran after a differential analysis.
    `differential` is a dataframe of the various comparisons with columns "comparison_name" and "direction".

    If `permissive`, will skip non-existing files, giving a warning.
    """
    import pandas as pd
    import numpy as np
    from ngs_toolkit.general import parse_ame

    if data_type not in ["ATAC-seq", "RNA-seq"]:
        raise AssertionError("`data_type` must match one of 'ATAC-seq' or 'RNA-seq'.")

    if "{data_type}" in output_dir:
        output_dir = output_dir.format(data_type=data_type)

    error_msg = "{} results for comparison '{}', direction '{}' were not found!"

    lola_enr = pd.DataFrame()
    meme_enr = pd.DataFrame()
    pathway_enr = pd.DataFrame()
    # Examine each region cluster
    for comp in differential['comparison_name'].drop_duplicates():
        # Separate in up/down-regulated genes
        for f, direction, top in [(np.less, "down", "head"), (np.greater, "up", "tail")]:
            comparison_dir = os.path.join(output_dir, "{}.{}".format(comp, direction))

            if data_type == "RNA-seq":
                # print("Collecting enrichments of comparison '{}', direction '{}'.".format(comp, direction))
                try:
                    enr = pd.read_csv(os.path.join(comparison_dir, output_prefix + ".gene_symbols.txt.enrichr.csv"))
                except IOError as e:
                    if permissive:
                        print(error_msg.format("Enrichr", comp, direction))
                    else:
                        raise e
                else:
                    enr["comparison_name"] = comp
                    enr["direction"] = direction
                    pathway_enr = pathway_enr.append(enr, ignore_index=True)
            elif data_type == "ATAC-seq":
                # print("Collecting enrichments of comparison '{}', direction '{}'.".format(comp, direction))
                # read/parse
                try:
                    motifs = parse_ame(comparison_dir).reset_index()
                    motifs.columns = ["TF", "p_value"]
                except IOError as e:
                    if permissive:
                        print(error_msg.format("Motif", comp, direction))
                    else:
                        raise e
                else:
                    motifs["comparison_name"] = comp
                    motifs["direction"] = direction
                    meme_enr = meme_enr.append(motifs, ignore_index=True)

                # LOLA
                try:
                    lola = pd.read_csv(os.path.join(comparison_dir, "allEnrichments.txt"), sep="\t")
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
                    enr = pd.read_csv(os.path.join(comparison_dir, output_prefix + "_genes.enrichr.csv"))
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
        lola_enr.to_csv(
            os.path.join(output_dir, output_prefix + ".lola.csv"), index=False)


def plot_differential_enrichment(
        enrichment_table,
        enrichment_type,
        data_type="ATAC-seq",
        direction_dependent=True,
        output_dir="results/differential_analysis_{data_type}",
        comp_variable="comparison_name",
        output_prefix="differential_analysis",
        top_n=5):
    """
    Given a table of enrichment terms across several comparisons, produce
    plots illustrating these enrichments in the various comparisons.

    `enrichment_type` is one of 'lola', 'enrichr', 'motif'.
    `comp_variable` is the column in the enrichment table that labels groups.
    """
    import numpy as np
    import seaborn as sns

    if enrichment_type not in ["lola", "enrichr", "motif", 'great']:
        raise AssertionError("`enrichment_type` must be one of 'lola', 'enrichr', 'motif', 'great.")

    if "{data_type}" in output_dir:
        output_dir = output_dir.format(data_type=data_type)

    if "direction" in enrichment_table.columns and direction_dependent:
        enrichment_table[comp_variable] = enrichment_table[comp_variable].astype(str) + " " + enrichment_table["direction"].astype(str)

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
        top_data = enrichment_table.set_index("label").groupby(comp_variable)["pValueLog"].nlargest(top_n).reset_index()

        n = len(enrichment_table[comp_variable].drop_duplicates())
        n_side = int(np.ceil(np.sqrt(n)))

        fig, axis = plt.subplots(n_side, n_side, figsize=(4 * n_side, n_side * max(5, 0.12 * top_n)), sharex=False, sharey=False)
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
        fig.savefig(os.path.join(output_dir, output_prefix + ".lola.barplot.top_{}.svg".format(top_n)), bbox_inches="tight", dpi=300)

        # Plot heatmaps of terms for each comparison
        if len(enrichment_table[comp_variable].drop_duplicates()) < 2:
            return

        # pivot table
        lola_pivot = pd.pivot_table(enrichment_table,
            values="pValueLog", columns=comp_variable, index="label").fillna(0)
        lola_pivot = lola_pivot.replace(np.inf, lola_pivot[lola_pivot != np.inf].max().max())

        # plot correlation
        g = sns.clustermap(
            lola_pivot.corr(), cbar_kws={"label": "Correlation of enrichemnt\nof differential regions"})
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.fig.savefig(os.path.join(output_dir, output_prefix + ".lola.correlation.svg"), bbox_inches="tight", dpi=300)

        top = enrichment_table.set_index('label').groupby(comp_variable)['pValueLog'].nlargest(top_n)
        top_terms = top.index.get_level_values('label').unique()

        # plot clustered heatmap
        shape = lola_pivot.loc[top_terms, :].shape
        g = sns.clustermap(
            lola_pivot.loc[top_terms, :], figsize=(max(6, 0.12 * shape[1]), max(6, 0.12 * shape[0])),
            cbar_kws={"label": "-log10(p-value) of enrichment\nof differential regions"}, metric="correlation")
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.fig.savefig(os.path.join(output_dir, output_prefix + ".lola.cluster_specific.svg"), bbox_inches="tight", dpi=300)

        # plot clustered heatmap
        g = sns.clustermap(
            lola_pivot.loc[top_terms, :], figsize=(max(6, 0.12 * shape[1]), max(6, 0.12 * shape[0])),
            z_score=0, cbar_kws={"label": "Z-score of enrichment\nof differential regions"}, metric="correlation")
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.fig.savefig(os.path.join(output_dir, output_prefix + ".lola.cluster_specific.z_score.svg"), bbox_inches="tight", dpi=300)

    if enrichment_type == "motif":
        enrichment_table["log_p_value"] = (-np.log10(enrichment_table["p_value"])).replace({np.inf: 300})
        # Plot top_n terms of each comparison in barplots
        top_data = enrichment_table.set_index("TF").groupby(comp_variable)["log_p_value"].nlargest(top_n).reset_index()

        n = len(enrichment_table[comp_variable].drop_duplicates())
        n_side = int(np.ceil(np.sqrt(n)))

        fig, axis = plt.subplots(n_side, n_side, figsize=(4 * n_side, n_side * max(5, 0.12 * top_n)), sharex=False, sharey=False)
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
        fig.savefig(os.path.join(output_dir, output_prefix + ".motifs.barplot.top_{}.svg".format(top_n)), bbox_inches="tight", dpi=300)

        # Plot heatmaps of terms for each comparison
        if len(enrichment_table[comp_variable].drop_duplicates()) < 2:
            return

        # pivot table
        motifs_pivot = pd.pivot_table(enrichment_table,
            values="log_p_value", columns=comp_variable, index="TF").fillna(0)

        # plot correlation
        g = sns.clustermap(motifs_pivot.corr(), cbar_kws={"label": "Correlation of enrichemnt\nof differential regions"})
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
        g.fig.savefig(os.path.join(output_dir, output_prefix + ".motifs.correlation.svg"), bbox_inches="tight", dpi=300)

        top = enrichment_table.set_index('TF').groupby(comp_variable)['log_p_value'].nlargest(top_n)
        top_terms = top.index.get_level_values('TF').unique()

        # plot clustered heatmap
        shape = motifs_pivot.loc[top_terms, :].shape
        g = sns.clustermap(motifs_pivot.loc[top_terms, :].T, figsize=(max(6, 0.12 * shape[0]), max(6, 0.12 * shape[1])),
            cbar_kws={"label": "-log10(p-value) of enrichment\nof differential regions"}, metric="correlation")
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.fig.savefig(os.path.join(output_dir, output_prefix + ".motifs.cluster_specific.svg"), bbox_inches="tight", dpi=300)

        # plot clustered heatmap
        shape = motifs_pivot.loc[top_terms, :].shape
        g = sns.clustermap(motifs_pivot.loc[top_terms, :].T, figsize=(max(6, 0.12 * shape[0]), max(6, 0.12 * shape[1])),
            z_score=0, cbar_kws={"label": "Z-score of enrichment\nof differential regions"}, metric="correlation")
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        g.fig.savefig(os.path.join(output_dir, output_prefix + ".motifs.cluster_specific.z_score.svg"), bbox_inches="tight", dpi=300)

    if enrichment_type == "enrichr":
        # enrichment_table["description"] = enrichment_table["description"].str.decode("utf-8")
        enrichment_table["log_p_value"] = (-np.log10(enrichment_table["p_value"])).replace({np.inf: 300})

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

            fig, axis = plt.subplots(n_side, n_side, figsize=(4 * n_side, n_side * max(5, 0.12 * top_n)), sharex=False, sharey=False)
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
            fig.savefig(os.path.join(output_dir, output_prefix + ".enrichr.{}.barplot.top_{}.svg".format(gene_set_library, top_n)), bbox_inches="tight", dpi=300)

            # Plot heatmaps of terms for each comparison
            if len(enrichment_table[comp_variable].drop_duplicates()) < 2:
                return

            # pivot table
            enrichr_pivot = pd.pivot_table(
                enrichment_table[enrichment_table["gene_set_library"] == gene_set_library],
                values="log_p_value", columns="description", index=comp_variable).fillna(0)

            # plot correlation
            g = sns.clustermap(enrichr_pivot.T.corr(), cbar_kws={"label": "Correlation of enrichemnt\nof differential genes"})
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
            g.fig.savefig(os.path.join(output_dir, output_prefix + ".enrichr.{}.correlation.svg".format(gene_set_library)), bbox_inches="tight", dpi=300)

            top = enrichment_table[enrichment_table["gene_set_library"] == gene_set_library].set_index('description').groupby(comp_variable)['p_value'].nsmallest(top_n)
            top_terms = top.index.get_level_values('description').unique()
            # top_terms = top_terms[top_terms.isin(lola_pivot.columns[lola_pivot.sum() > 5])]

            # plot clustered heatmap
            shape = enrichr_pivot[list(set(top_terms))].shape
            g = sns.clustermap(enrichr_pivot[list(set(top_terms))].T, figsize=(max(6, 0.12 * shape[0]), max(6, 0.12 * shape[1])),
                cbar_kws={"label": "-log10(p-value) of enrichment\nof differential genes"}, metric="correlation")
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
            g.fig.savefig(os.path.join(output_dir, output_prefix + ".enrichr.{}.cluster_specific.svg".format(gene_set_library)), bbox_inches="tight", dpi=300)

            # plot clustered heatmap
            shape = enrichr_pivot[list(set(top_terms))].shape
            g = sns.clustermap(enrichr_pivot[list(set(top_terms))].T, figsize=(max(6, 0.12 * shape[0]), max(6, 0.12 * shape[1])),
                z_score=0, cbar_kws={"label": "Z-score of enrichment\nof differential genes"}, metric="correlation")
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
            g.fig.savefig(os.path.join(output_dir, output_prefix + ".enrichr.{}.cluster_specific.z_score.svg".format(gene_set_library)), bbox_inches="tight", dpi=300)

    if enrichment_type == "great":
        # enrichment_table["description"] = enrichment_table["description"].str.decode("utf-8")
        enrichment_table["log_q_value"] = (-np.log10(enrichment_table["HyperFdrQ"])).replace({np.inf: 300})

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

            fig, axis = plt.subplots(n_side, n_side, figsize=(4 * n_side, n_side * max(5, 0.12 * top_n)), sharex=False, sharey=False)
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
            fig.savefig(os.path.join(output_dir, output_prefix + ".great.{}.barplot.top_{}.svg".format(gene_set_library, top_n)), bbox_inches="tight", dpi=300)

            # Plot heatmaps of terms for each comparison
            if len(enrichment_table[comp_variable].drop_duplicates()) < 2:
                return

            # pivot table
            great_pivot = pd.pivot_table(
                enrichment_table[enrichment_table["Ontology"] == gene_set_library],
                values="log_q_value", columns="Desc", index=comp_variable).fillna(0)

            # plot correlation
            try:
                g = sns.clustermap(great_pivot.T.corr(), cbar_kws={"label": "Correlation of enrichemnt\nof differential genes"})
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
                g.fig.savefig(os.path.join(output_dir, output_prefix + ".great.{}.correlation.svg".format(gene_set_library)), bbox_inches="tight", dpi=300)
            except FloatingPointError:
                continue

            top = enrichment_table[enrichment_table["Ontology"] == gene_set_library].set_index('Desc').groupby(comp_variable)['HyperP'].nsmallest(top_n)
            top_terms = top.index.get_level_values('Desc').unique()
            # top_terms = top_terms[top_terms.isin(lola_pivot.columns[lola_pivot.sum() > 5])]

            # plot clustered heatmap
            shape = great_pivot[list(set(top_terms))].shape
            g = sns.clustermap(great_pivot[list(set(top_terms))].T, figsize=(max(6, 0.12 * shape[0]), max(6, 0.12 * shape[1])),
                cbar_kws={"label": "-log10(p-value) of enrichment\nof differential genes"}, metric="correlation")
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
            g.fig.savefig(os.path.join(output_dir, output_prefix + ".great.{}.cluster_specific.svg".format(gene_set_library)), bbox_inches="tight", dpi=300)

            # plot clustered heatmap
            shape = great_pivot[list(set(top_terms))].shape
            g = sns.clustermap(great_pivot[list(set(top_terms))].T, figsize=(max(6, 0.12 * shape[0]), max(6, 0.12 * shape[1])),
                z_score=0, cbar_kws={"label": "Z-score of enrichment\nof differential genes"}, metric="correlation")
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, ha="right", fontsize="xx-small")
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
            g.fig.savefig(os.path.join(output_dir, output_prefix + ".great.{}.cluster_specific.z_score.svg".format(gene_set_library)), bbox_inches="tight", dpi=300)


def chunks(l, n):
    n = max(1, n)
    return (l[i:i + n] for i in range(0, len(l), n))


def standard_scale(x):
    return (x - x.min()) / (x.max() - x.min())


def z_score(x):
    return (x - x.mean()) / x.std()


def signed_max(x, f=0.66):
    """
    Return maximum or minimum depending on the sign of the majority of values.
    If there isn't a clear majority (at least `f` fraction in one side), return mean of values.
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


def detect_peaks(x, mph=None, mpd=1, threshold=0, edge='rising',
                 kpsh=False, valley=False):
    """Detect peaks in data based on their amplitude and other features.

    Parameters
    ----------
    x : 1D array_like
        data.
    mph : {None, number}, optional (default = None)
        detect peaks that are greater than minimum peak height.
    mpd : positive integer, optional (default = 1)
        detect peaks that are at least separated by minimum peak distance (in
        number of data).
    threshold : positive number, optional (default = 0)
        detect peaks (valleys) that are greater (smaller) than `threshold`
        in relation to their immediate neighbors.
    edge : {None, 'rising', 'falling', 'both'}, optional (default = 'rising')
        for a flat peak, keep only the rising edge ('rising'), only the
        falling edge ('falling'), both edges ('both'), or don't detect a
        flat peak (None).
    kpsh : bool, optional (default = False)
        keep peaks with same height even if they are closer than `mpd`.
    valley : bool, optional (default = False)
        if True (1), detect valleys (local minima) instead of peaks.
    show : bool, optional (default = False)
        if True (1), plot data in matplotlib figure.
    ax : a matplotlib.axes.Axes instance, optional (default = None).

    Returns
    -------
    ind : 1D array_like
        indeces of the peaks in `x`.

    Notes
    -----
    The detection of valleys instead of peaks is performed internally by simply
    negating the data: `ind_valleys = detect_peaks(-x)`

    The function can handle NaN's

    See this IPython Notebook [1]_.

    References
    ----------
    .. [1] http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb

    Examples
    --------
    >>> from detect_peaks import detect_peaks
    >>> x = np.random.randn(100)
    >>> x[60:81] = np.nan
    >>> # detect all peaks and plot data
    >>> ind = detect_peaks(x, show=True)
    >>> print(ind)

    >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
    >>> # set minimum peak height = 0 and minimum peak distance = 20
    >>> detect_peaks(x, mph=0, mpd=20, show=True)

    >>> x = [0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0]
    >>> # set minimum peak distance = 2
    >>> detect_peaks(x, mpd=2, show=True)

    >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
    >>> # detection of valleys instead of peaks
    >>> detect_peaks(x, mph=0, mpd=20, valley=True, show=True)

    >>> x = [0, 1, 1, 0, 1, 1, 0]
    >>> # detect both edges
    >>> detect_peaks(x, edge='both', show=True)

    >>> x = [-2, 1, -2, 2, 1, 1, 3, 0]
    >>> # set threshold = 2
    >>> detect_peaks(x, threshold = 2, show=True)
    """
    import numpy as np

    x = np.atleast_1d(x).astype('float64')
    if x.size < 3:
        return np.array([], dtype=int)
    if valley:
        x = -x
    # find indices of all peaks
    dx = x[1:] - x[:-1]
    # handle NaN's
    indnan = np.where(np.isnan(x))[0]
    if indnan.size:
        x[indnan] = np.inf
        dx[np.where(np.isnan(dx))[0]] = np.inf
    ine, ire, ife = np.array([[], [], []], dtype=int)
    if not edge:
        ine = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) > 0))[0]
    else:
        if edge.lower() in ['rising', 'both']:
            ire = np.where((np.hstack((dx, 0)) <= 0) & (np.hstack((0, dx)) > 0))[0]
        if edge.lower() in ['falling', 'both']:
            ife = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) >= 0))[0]
    ind = np.unique(np.hstack((ine, ire, ife)))
    # handle NaN's
    if ind.size and indnan.size:
        # NaN's and values close to NaN's cannot be peaks
        ind = ind[np.in1d(ind, np.unique(np.hstack((indnan, indnan - 1, indnan + 1))), invert=True)]
    # first and last values of x cannot be peaks
    if ind.size and ind[0] == 0:
        ind = ind[1:]
    if ind.size and ind[-1] == x.size - 1:
        ind = ind[:-1]
    # remove peaks < minimum peak height
    if ind.size and mph is not None:
        ind = ind[x[ind] >= mph]
    # remove peaks - neighbors < threshold
    if ind.size and threshold > 0:
        dx = np.min(np.vstack([x[ind] - x[ind - 1], x[ind] - x[ind + 1]]), axis=0)
        ind = np.delete(ind, np.where(dx < threshold)[0])
    # detect small peaks closer than minimum peak distance
    if ind.size and mpd > 1:
        ind = ind[np.argsort(x[ind])][::-1]  # sort ind by peak height
        idel = np.zeros(ind.size, dtype=bool)
        for i in range(ind.size):
            if not idel[i]:
                # keep peaks with the same height if kpsh is True
                idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
                    & (x[ind[i]] > x[ind] if kpsh else True)
                idel[i] = 0  # Keep current peak
        # remove the small peaks and sort back the indices by their occurrence
        ind = np.sort(ind[~idel])

    return ind


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
