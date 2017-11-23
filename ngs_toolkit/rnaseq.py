
import os
from ngs_toolkit.general import Analysis, pickle_me
from collections import Counter
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybedtools
import seaborn as sns


class RNASeqAnalysis(Analysis):
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
        super(RNASeqAnalysis, self).__init__(
            name=name,
            data_dir=data_dir,
            results_dir=results_dir,
            pickle_file=pickle_file,
            samples=samples,
            prj=prj,
            from_pickle=from_pickle,
            **kwargs)

    def annotate_with_sample_metadata(
            self,
            attributes=["sample_name"]):

        samples = [s for s in self.samples if s.name in self.coverage_annotated.columns]

        attrs = list()
        for attr in attributes:
            l = list()
            for sample in samples:  # keep order of samples in matrix
                try:
                    l.append(getattr(sample, attr))
                except AttributeError:
                    l.append(np.nan)
            attrs.append(l)

        # Generate multiindex columns
        index = pd.MultiIndex.from_arrays(attrs, names=attributes)
        self.expression = self.coverage_annotated[[s.name for s in samples]]
        self.expression.columns = index

        # Save
        self.expression.to_csv(os.path.join(self.results_dir, self.name + ".expression.annotated_metadata.csv"), index=True)

    def get_level_colors(self, index=None, levels=None, pallete="Paired", cmap="RdBu_r", nan_color=(0.662745, 0.662745, 0.662745, 0.5)):
        if index is None:
            index = self.expression.columns

        if levels is not None:
            index = index.droplevel([l.name for l in index.levels if l.name not in levels])

        _cmap = plt.get_cmap(cmap)
        _pallete = plt.get_cmap(pallete)

        colors = list()
        for level in index.levels:
            # determine the type of data in each level
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
                col[np.where(index.get_level_values(level.name).to_series().isnull().tolist())] = nan_color
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

    def collect_bitseq_output(self, samples=None, permissive=True):
        """
        Collect gene expression (read counts, transcript-level) output from Bitseq into expression matrix for `samples`.
        """
        if samples is None:
            samples = self.samples

        first = True
        for i, sample in enumerate(samples):
            if first:
                try:
                    # read the "tr" file of one sample to get indexes
                    tr = pd.read_csv(
                        os.path.join(
                            sample.paths.sample_root, "bowtie1_{}".format(sample.transcriptome),
                            "bitSeq",
                            sample.name + ".tr"),
                        sep=" ", header=None, skiprows=1,
                        names=["ensembl_gene_id", "ensembl_transcript_id", "v1", "v2"])
                except IOError("Sample {} is missing.".format(sample.name)) as e:
                    if permissive:
                        print(e)
                    else:
                        raise e
                # add id index
                tr.set_index("ensembl_gene_id", append=False, inplace=True)
                tr.set_index("ensembl_transcript_id", append=True, inplace=True)
                # create dataframe
                expr = pd.DataFrame(index=tr.index)
                first = False

            # load counts
            try:
                counts = pd.read_csv(os.path.join(
                    sample.paths.sample_root, "bowtie1_{}".format(sample.transcriptome),
                    "bitSeq",
                    sample.name + ".counts"), sep=" ")
            except IOError:
                print("Sample {} is missing.".format(sample.name))
                continue
            counts.index = expr.index

            # Append
            expr[sample.name] = counts

        return expr

    def collect_esat_output(self, samples=None, permissive=True):
        """
        Collect gene expression (read counts, gene-level) output from ESAT into expression matrix for `samples`.
        """
        if samples is None:
            samples = self.samples

        first = True
        for i, sample in enumerate(samples):
            try:
                # read the "tr" file of one sample to get indexes
                c = pd.read_csv(
                    os.path.join(
                        sample.paths.sample_root, "ESAT_{}".format(sample.genome), sample.name + ".gene.txt"), sep="\t")
            except IOError("Sample {} is missing.".format(sample.name)) as e:
                if permissive:
                    print(e)
                    continue
                else:
                    raise e
            # extract only gene ID and counts
            c = c[["Symbol", "Exp1"]]
            c = c.rename(columns={"Symbol": "gene_symbol", "Exp1": sample.name}).set_index("gene_symbol")

            # Append
            if first:
                expr = c
            else:
                expr = expr.join(c)
            first = False

        return expr.sort_index()

    def get_gene_expression(self, samples=None, attributes=["sample_name"]):
        """
        Collect gene expression (read counts, transcript level) for all samples,
        annotates ensembl IDs with gene names, reduces gene expression to gene-level,
        normalizes expression matrix (quantile normalization), and
        annotates samples with given attributes (generates pandas dataframe with MultiIndex columns).
        """
        import requests

        def quantile_normalize(df_input):
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

        if samples is None:
            samples = [s for s in self.samples if s.library == "RNA-seq"]

        self.count_matrix = self.collect_bitseq_output(samples=samples)

        # Map ensembl gene IDs to gene names
        url_query = "".join([
            """http://grch37.ensembl.org/biomart/martservice?query=""",
            """<?xml version="1.0" encoding="UTF-8"?>""",
            """<!DOCTYPE Query>""",
            """<Query  virtualSchemaName = "default" formatter = "CSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >""",
            """<Dataset name = "hsapiens_gene_ensembl" interface = "default" >""",
            """<Attribute name = "ensembl_gene_id" />""",
            """<Attribute name = "external_gene_name" />""",
            """</Dataset>""",
            """</Query>"""])
        req = requests.get(url_query, stream=True)
        mapping = pd.DataFrame((x.strip().split(",") for x in list(req.iter_lines())), columns=["ensembl_gene_id", "gene_name"])
        self.count_matrix = self.count_matrix.reset_index()
        self.count_matrix['ensembl_gene_id'] = self.count_matrix['ensembl_gene_id'].str.replace("\..*", "")

        self.count_matrix = pd.merge(self.count_matrix, mapping, how="outer").dropna()
        self.count_matrix.to_csv(os.path.join(self.results_dir, self.name + ".expression_counts.csv"), index=False)

        # Quantile normalize
        self.matrix_qnorm = quantile_normalize(self.count_matrix.drop(["ensembl_transcript_id", "ensembl_gene_id", "gene_name"], axis=1))

        # Log2 TPM
        self.matrix_qnorm_log = np.log2((1 + self.matrix_qnorm) / self.matrix_qnorm.sum(axis=0) * 1e6)
        self.matrix_qnorm_log = self.count_matrix[["gene_name", "ensembl_gene_id", "ensembl_transcript_id"]].join(self.matrix_qnorm_log)
        self.matrix_qnorm_log.to_csv(os.path.join(self.results_dir, self.name + ".expression_counts.transcript_level.quantile_normalized.log2_tpm.csv"), index=False)

        # Reduce to gene-level measurements by max of transcripts
        self.expression_matrix_counts = self.count_matrix.drop(['ensembl_transcript_id'], axis=1).dropna().set_index(['ensembl_gene_id', 'gene_name']).groupby(level=["gene_name"]).max()

        # Quantile normalize
        self.expression_matrix_qnorm = quantile_normalize(self.expression_matrix_counts)

        # Log2 TPM
        self.expression = np.log2((self.expression_matrix_qnorm / self.expression_matrix_qnorm.sum(axis=0)) * 1e6)

        # Annotate with sample metadata
        _samples = [s for s in samples if s.name in self.expression.columns]
        attrs = list()
        for attr in attributes:
            l = list()
            for sample in _samples:  # keep order of samples in matrix
                try:
                    l.append(getattr(sample, attr))
                except AttributeError:
                    l.append(np.nan)
            attrs.append(l)

        # Generate multiindex columns
        index = pd.MultiIndex.from_arrays(attrs, names=attributes)
        self.expression_annotated = self.expression[[s.name for s in _samples]]
        self.expression_annotated.columns = index

        # Save
        self.expression_matrix_counts.to_csv(os.path.join(self.results_dir, self.name + ".expression_counts.gene_level.csv"), index=True)
        self.expression_matrix_qnorm.to_csv(os.path.join(self.results_dir, self.name + ".expression_counts.gene_level.quantile_normalized.csv"), index=True)
        self.expression.to_csv(os.path.join(self.results_dir, self.name + ".expression_counts.gene_level.quantile_normalized.log2_tpm.csv"), index=True)
        self.expression_annotated.to_csv(os.path.join(self.results_dir, self.name + ".expression_counts.gene_level.quantile_normalized.log2_tpm.annotated_metadata.csv"), index=True)

    def plot_expression_characteristics(self):
        fig, axis = plt.subplots(figsize=(4, 4 * np.log10(len(self.expression_matrix_counts.columns))))
        sns.barplot(
            data=self.count_matrix.drop(['gene_name', 'ensembl_gene_id', 'ensembl_transcript_id'], axis=1).sum().sort_values().reset_index(),
            y="index", x=0, orient="horiz", color=sns.color_palette("colorblind")[0], ax=axis)
        axis.set_xlabel("Transcriptome reads")
        axis.set_ylabel("Samples")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, self.name + ".expression.reads_per_sample.svg"), bbox_inches="tight")

        cov = pd.DataFrame()
        for i in [1, 2, 3, 6, 12, 24, 48, 96, 200, 300, 400, 500, 1000]:
            cov[i] = self.expression_matrix_counts.apply(lambda x: sum(x >= i))

        fig, axis = plt.subplots(1, 2, figsize=(6 * 2, 6))
        sns.heatmap(cov.drop([1, 2], axis=1), ax=axis[0], cmap="GnBu", cbar_kws={"label": "Genes covered"})
        sns.heatmap(cov.drop([1, 2], axis=1).apply(lambda x: (x - x.mean()) / x.std(), axis=0), ax=axis[1], cbar_kws={"label": "Z-score"})
        for ax in axis:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
            ax.set_ylabel("Samples")
            ax.set_xlabel("Genes with #reads")
        axis[1].set_yticklabels(ax.get_yticklabels(), visible=False)
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, self.name + ".expression.genes_with_reads.svg"), bbox_inches="tight")

        for name, matrix in [("counts", self.expression_matrix_counts), ("qnorm_TPM", self.expression)]:
            # Boxplot with values per sample
            to_plot = pd.melt(matrix, var_name="sample")

            fig, axis = plt.subplots()
            sns.boxplot(data=to_plot, y="sample", x="value", orient="horiz", ax=axis)
            axis.set_xlabel(name)
            axis.set_ylabel("Samples")
            sns.despine(fig)
            fig.savefig(os.path.join(self.results_dir, self.name + ".expression.boxplot_per_sample.{}.svg".format(name)), bbox_inches="tight")

        # # Plot gene expression along chromossomes
        # import requests
        # url_query = "".join([
        #     """http://grch37.ensembl.org/biomart/martservice?query=""",
        #     """<?xml version="1.0" encoding="UTF-8"?>""",
        #     """<!DOCTYPE Query>""",
        #     """<Query  virtualSchemaName = "default" formatter = "CSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >""",
        #     """<Dataset name = "hsapiens_gene_ensembl" interface = "default" >""",
        #     """<Attribute name = "chromosome_name" />""",
        #     """<Attribute name = "start_position" />""",
        #     """<Attribute name = "end_position" />""",
        #     """<Attribute name = "external_gene_name" />""",
        #     """</Dataset>""",
        #     """</Query>"""])
        # req = requests.get(url_query, stream=True)
        # annot = pd.DataFrame((x.strip().split(",") for x in list(req.iter_lines())), columns=["chr", "start", "end", "gene_name"])
        # gene_order = annot[annot['chr'].isin([str(x) for x in range(22)] + ['X', 'Y'])].sort_values(["chr", "start", "end"])["gene_name"]

        # # sort genes by chrom order
        # exp = self.expression.ix[gene_order].dropna()
        # # cap expression
        # exp[exp > np.percentile(exp, 75)] = np.percentile(exp, 75)
        # exp[exp < np.percentile(exp, 25)] = np.percentile(exp, 25)

        # data = ((exp.T - exp.T.mean()) / exp.T.std()).T

        # data2 = data.rolling(int(3e4), axis=0).mean().dropna()
        # plt.pcolor(data2)

    def unsupervised(
            self, quant_matrix="expression_annotated", samples=None, attributes_to_plot=["sample_name"], plot_prefix="all_genes"):
        """
        Apply unsupervised clustering (clustering of correlations) and dimentionality reduction methods (MDS, PCA) on matrix.
        Colours and labels samples by attributes in `attributes_to_plot`.
        """
        from sklearn.decomposition import PCA
        from sklearn.manifold import MDS
        from collections import OrderedDict
        import re
        import itertools
        from scipy.stats import kruskal
        from scipy.stats import pearsonr

        matrix = getattr(self, quant_matrix)

        if samples is None:
            samples = [s for s in self.samples if s.name in matrix.columns.get_level_values("sample_name")]

        # This will always be a matrix for all samples
        color_dataframe = pd.DataFrame(self.get_level_colors(index=matrix.columns, levels=attributes_to_plot), index=attributes_to_plot, columns=matrix.columns.get_level_values("sample_name"))
        # will be filtered now by the requested samples if needed
        color_dataframe = color_dataframe[[s.name for s in samples]]
        # # exclude samples if needed
        # color_dataframe = color_dataframe[[s.name for s in samples]]
        # sample_display_names = color_dataframe.columns.str.replace("ATAC-seq_", "")

        # All regions, matching samples (provided samples in matrix)
        X = matrix.loc[:, matrix.columns.get_level_values("sample_name").isin([s.name for s in samples])]

        # Pairwise correlations
        g = sns.clustermap(
            X.astype(float).corr(), xticklabels=False, annot=True,  # yticklabels=sample_display_names, 
            cmap="Spectral_r", figsize=(15, 15), cbar_kws={"label": "Pearson correlation"}, row_colors=color_dataframe.values.tolist())
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize='xx-small')
        g.ax_heatmap.set_xlabel(None, visible=False)
        g.ax_heatmap.set_ylabel(None, visible=False)
        g.fig.savefig(os.path.join(self.results_dir, "{}.{}.corr.clustermap.svg".format(self.name, plot_prefix)), bbox_inches='tight')

        # MDS
        mds = MDS(n_jobs=-1)
        x_new = mds.fit_transform(X.T)
        # transform again
        x = pd.DataFrame(x_new)
        xx = x.apply(lambda j: (j - j.mean()) / j.std(), axis=0)

        fig, axis = plt.subplots(1, len(attributes_to_plot), figsize=(4 * len(attributes_to_plot), 4 * 1))
        axis = axis.flatten()
        for i, attr in enumerate(attributes_to_plot):
            for j in range(len(xx)):
                try:
                    label = getattr(samples[j], attributes_to_plot[i])
                except AttributeError:
                    label = np.nan
                axis[i].scatter(xx.loc[j, 0], xx.loc[j, 1], s=50, color=color_dataframe.loc[attr].iloc[j], label=label)
            axis[i].set_title(attributes_to_plot[i])
            axis[i].set_xlabel("MDS 1")
            axis[i].set_ylabel("MDS 2")
            axis[i].set_xticklabels(axis[i].get_xticklabels(), visible=False)
            axis[i].set_yticklabels(axis[i].get_yticklabels(), visible=False)

            # Unique legend labels
            handles, labels = axis[i].get_legend_handles_labels()
            by_label = OrderedDict(zip(labels, handles))
            if any([type(c) in [str, unicode] for c in by_label.keys()]) and len(by_label) <= 20:
                # if not any([re.match("^\d", c) for c in by_label.keys()]):
                axis[i].legend(by_label.values(), by_label.keys())
        fig.savefig(os.path.join(self.results_dir, "{}.{}.mds.svg".format(self.name, plot_prefix)), bbox_inches="tight")

        # PCA
        pca = PCA()
        x_new = pca.fit_transform(X.T)
        # transform again
        xx = pd.DataFrame(x_new)

        # plot % explained variance per PC
        fig, axis = plt.subplots(1)
        axis.plot(
            range(1, len(pca.explained_variance_) + 1),  # all PCs
            (pca.explained_variance_ / pca.explained_variance_.sum()) * 100, 'o-')  # % of total variance
        axis.axvline(len(attributes_to_plot), linestyle='--')
        axis.set_xlabel("PC")
        axis.set_ylabel("% variance")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.{}.pca.explained_variance.svg".format(self.name, plot_prefix)), bbox_inches='tight')

        # plot
        pcs = min(xx.shape[0] - 1, 8)
        fig, axis = plt.subplots(pcs, len(attributes_to_plot), figsize=(4 * len(attributes_to_plot), 4 * pcs))
        for pc in range(pcs):
            for i, attr in enumerate(attributes_to_plot):
                for j in range(len(xx)):
                    try:
                        label = getattr(samples[j], attributes_to_plot[i])
                    except AttributeError:
                        label = np.nan
                    axis[pc, i].scatter(xx.loc[j, pc], xx.loc[j, pc + 1], s=50, color=color_dataframe.loc[attr].iloc[j], label=label)
                axis[pc, i].set_title(attributes_to_plot[i])
                axis[pc, i].set_xlabel("PC {}".format(pc + 1))
                axis[pc, i].set_ylabel("PC {}".format(pc + 2))
                axis[pc, i].set_xticklabels(axis[pc, i].get_xticklabels(), visible=False)
                axis[pc, i].set_yticklabels(axis[pc, i].get_yticklabels(), visible=False)

                # Unique legend labels
                handles, labels = axis[pc, i].get_legend_handles_labels()
                by_label = OrderedDict(zip(labels, handles))
                if any([type(c) in [str, unicode] for c in by_label.keys()]) and len(by_label) <= 20:
                    # if not any([re.match("^\d", c) for c in by_label.keys()]):
                    axis[pc, i].legend(by_label.values(), by_label.keys())
        fig.savefig(os.path.join(self.results_dir, "{}.{}.pca.svg".format(self.name, plot_prefix)), bbox_inches="tight")

        #

        # # Test association of PCs with attributes
        associations = list()
        for pc in range(pcs):
            for attr in attributes_to_plot:
                print("PC {}; Attribute {}.".format(pc + 1, attr))
                sel_samples = [s for s in samples if hasattr(s, attr)]
                sel_samples = [s for s in sel_samples if not pd.isnull(getattr(s, attr))]

                # Get all values of samples for this attr
                groups = set([getattr(s, attr) for s in sel_samples])

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
                        g1_indexes = [i for i, s in enumerate(sel_samples) if getattr(s, attr) == group1]
                        g2_indexes = [i for i, s in enumerate(sel_samples) if getattr(s, attr) == group2]

                        g1_values = xx.loc[g1_indexes, pc]
                        g2_values = xx.loc[g2_indexes, pc]

                        # Test ANOVA (or Kruskal-Wallis H-test)
                        p = kruskal(g1_values, g2_values)[1]

                        # Append
                        associations.append([pc + 1, attr, variable_type, group1, group2, p])

                elif variable_type == "numerical":
                    # It numerical, calculate pearson correlation
                    indexes = [i for i, s in enumerate(samples) if s in sel_samples]
                    pc_values = xx.loc[indexes, pc]
                    trait_values = [getattr(s, attr) for s in sel_samples]
                    p = pearsonr(pc_values, trait_values)[1]

                    associations.append([pc + 1, attr, variable_type, np.nan, np.nan, p])

        associations = pd.DataFrame(associations, columns=["pc", "attribute", "variable_type", "group_1", "group_2", "p_value"])

        # write
        associations.to_csv(os.path.join(self.results_dir, "{}.{}.pca.variable_principle_components_association.csv".format(self.name, plot_prefix)), index=False)

        # Plot
        # associations[associations['p_value'] < 0.05].drop(['group_1', 'group_2'], axis=1).drop_duplicates()
        # associations.drop(['group_1', 'group_2'], axis=1).drop_duplicates().pivot(index="pc", columns="attribute", values="p_value")
        pivot = associations.groupby(["pc", "attribute"]).min()['p_value'].reset_index().pivot(index="pc", columns="attribute", values="p_value").dropna(axis=1)

        # heatmap of -log p-values
        g = sns.clustermap(-np.log10(pivot), row_cluster=False, annot=True, cbar_kws={"label": "-log10(p_value) of association"})
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
        g.fig.savefig(os.path.join(self.results_dir, "{}.{}.pca.variable_principle_components_association.svg".format(self.name, plot_prefix)), bbox_inches="tight")

        # heatmap of masked significant
        g = sns.clustermap((pivot < 0.05).astype(int), row_cluster=False, cbar_kws={"label": "significant association"})
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
        g.fig.savefig(os.path.join(self.results_dir, "{}.{}.pca.variable_principle_components_association.masked.svg".format(self.name, plot_prefix)), bbox_inches="tight")



def knockout_plot(
    analysis=None, 
    knockout_genes=None,
    expression_matrix=None,
    output_dir=None,
    output_prefix="knockout_expression"
    ):
    """
    Plot expression of knocked-out genes in all samples.
    """

    if (analysis is None) and (expression_matrix is None):
        raise AssertionError("One of `analysis` or `expression_matrix` must be provided.")

    if (analysis is None) and (knockout_genes is None):
        raise AssertionError("If an `analysis` object is not provided, you must provide a list of `knockout_genes`.")
    elif (analysis is not None) and (knockout_genes is None):
        try:
            knockout_genes = list(set([s.knockout for s in analysis.samples]))
        except KeyError("If `knockout_genes` is not given, Samples in `analysis` must have a `knockout` attribute.") as e:
            raise e

    if expression_matrix is None:
        expression_matrix = analysis.expression

    if output_dir is None:
        if analysis is not None:
            output_dir = analysis.results_dir
        else:
            output_dir = os.path.curdir

    knockout_genes = sorted(knockout_genes)

    missing = [k for k in knockout_genes if k not in expression_matrix.index]
    if len(missing) > 0:
        print(Warning("The following `knockout_genes` were not found in the expression matrix: '{}'".format(", ".join(missing))))
    knockout_genes = [k for k in knockout_genes if k in expression_matrix.index]

    ko = expression_matrix.loc[knockout_genes, :]

    g = sns.clustermap(ko, cbar_kws={"label": "Expression"})
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + ".svg"), bbox_inches="tight")
    
    g = sns.clustermap(ko, z_score=0, cbar_kws={"label": "Expression Z-score"})
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + ".z_score.svg"), bbox_inches="tight")

    g = sns.clustermap(ko, cbar_kws={"label": "Expression"}, row_cluster=False, col_cluster=False)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + ".sorted.svg"), bbox_inches="tight")
    
    g = sns.clustermap(ko, z_score=0, cbar_kws={"label": "Expression Z-score"}, row_cluster=False, col_cluster=False)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + ".z_score.sorted.svg"), bbox_inches="tight")
