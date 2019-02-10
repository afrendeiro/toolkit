#!/usr/bin/env python


import os

import matplotlib.pyplot as plt
from ngs_toolkit import _LOGGER
from ngs_toolkit.analysis import Analysis
from ngs_toolkit.general import normalize_quantiles_p
from ngs_toolkit.general import query_biomart
import numpy as np
import pandas as pd
import scipy
from scipy.stats import zscore
import seaborn as sns
# import requests


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
            pep=False,
            **kwargs):
        super(RNASeqAnalysis, self).__init__(
            name=name,
            data_dir=data_dir,
            results_dir=results_dir,
            pickle_file=pickle_file,
            samples=samples,
            prj=prj,
            from_pickle=from_pickle,
            pep=pep,
            **kwargs)

        self.data_type = self.__data_type__ = "RNA-seq"
        self.var_names = "gene"
        self.quantity = "expression"
        self.norm_units = "RPM"
        self.raw_matrix_name = "count_matrix"
        self.norm_matrix_name = "expression"
        self.annot_matrix_name = "expression_annotated"

    def collect_bitseq_output(self, samples=None, permissive=True, expression_type="counts"):
        """
        Collect gene expression (read counts, transcript-level) output from Bitseq
        into expression matrix for `samples`.
        """
        if samples is None:
            samples = self.samples

        first = True
        for i, sample in enumerate(samples):
            if first:
                msg = "Sample {} is missing.".format(sample.name)
                try:
                    # read the "tr" file of one sample to get indexes
                    tr = pd.read_csv(
                        os.path.join(
                            sample.paths.sample_root, "bowtie1_{}".format(sample.transcriptome),
                            "bitSeq",
                            sample.name + ".tr"),
                        sep=" ", header=None, skiprows=1,
                        names=["ensembl_gene_id", "ensembl_transcript_id", "v1", "v2"])
                except IOError(msg) as e:
                    if permissive:
                        _LOGGER.warning(e)
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
                if expression_type == "counts":
                    e = pd.read_csv(os.path.join(
                            sample.paths.sample_root, "bowtie1_{}".format(sample.transcriptome),
                            "bitSeq",
                            sample.name + ".counts"), sep=" ")
                elif expression_type == "rpkm":
                    e = pd.read_csv(os.path.join(
                            sample.paths.sample_root, "bowtie1_{}".format(sample.transcriptome),
                            "bitSeq",
                            sample.name + ".mean"), sep=" ", comment="#", header=None).iloc[:, 0]
                else:
                    raise ValueError("Argument 'expression_type' must be one of {counts, rpkm}")
            except IOError:
                _LOGGER.warning("Sample {} is missing.".format(sample.name))
                continue
            e.index = expr.index

            # Append
            expr[sample.name] = e

        return expr

    def collect_esat_output(self, samples=None, permissive=True):
        """
        Collect gene expression (read counts, gene-level) output from ESAT
        into expression matrix for `samples`.
        """
        if samples is None:
            samples = self.samples

        first = True
        for i, sample in enumerate(samples):
            msg = "Sample {} is missing.".format(sample.name)
            try:
                # read the "tr" file of one sample to get indexes
                c = pd.read_csv(
                    os.path.join(
                        sample.paths.sample_root,
                        "ESAT_{}".format(sample.genome),
                        sample.name + ".gene.txt"), sep="\t")
            except IOError(msg) as e:
                if permissive:
                    _LOGGER.warning(e)
                    continue
                else:
                    raise e
            # extract only gene ID and counts
            c = c[["Symbol", "Exp1"]]
            c = c.rename(
                columns={"Symbol": "gene_symbol", "Exp1": sample.name}).set_index("gene_symbol")

            # Append
            if first:
                expr = c
            else:
                expr = expr.join(c)
            first = False

        return expr.sort_index()

    # def collect_htseqcount_output(self, samples=None, permissive=True):
    #     """
    #     Collect gene expression (read counts, gene-level) output from
    #     HTSeq-count into expression matrix for `samples`.
    #     """
    #     if samples is None:
    #         samples = self.samples

    #     first = True
    #     for i, sample in enumerate(samples):
    #         msg = "Sample {} is missing.".format(sample.name)
    #         try:
    #             # read the "tr" file of one sample to get indexes
    #             c = pd.read_csv(
    #                 os.path.join(
    #                     sample.paths.sample_root,
    #                     "tophat_{}".format(sample.genome),
    #                     sample.name + ".aln_sorted.htseq-count.tsv"), sep="\t")
    #             c.columns = ['ensembl_gene_id', 'ensembl_transcript_id', 'counts']
    #         except IOError(msg) as e:
    #             if permissive:
    #                 _LOGGER.warning(e)
    #                 continue
    #             else:
    #                 raise e
    #         # extract only gene ID and counts
    #         c = c[["Symbol", "Exp1"]]
    #         c = c.rename(
    #             columns={"Symbol": "gene_symbol", "Exp1": sample.name}).set_index("gene_symbol")

    #         # Append
    #         if first:
    #             expr = c
    #         else:
    #             expr = expr.join(c)
    #         first = False

    #     return expr.sort_index()

    def get_gene_expression(
            self, samples=None, sample_attributes=["sample_name"],
            expression_type="counts",
            genome_assembly="grch37", species="hsapiens", mul_factor=1e6):
        """
        Collect gene expression (read counts, transcript level) for all samples,
        annotates ensembl IDs with gene names, reduces gene expression to gene-level,
        normalizes expression matrix (quantile normalization), and
        annotates samples with given attributes (generates pandas dataframe with MultiIndex columns).

        Parameters
        ----------
        samples : peppy.Sample
            Samples to get expression for. Default all in analysis.

        sample_attributes : list
            Sample attributes to annotate expression matrix with.

        expression_type : str
            Type of expression quantification to get. One of "counts" or "rpkm".

        genome_assembly : str
            Genome assembly to use (e.g. "grch38") or Ensembl prefix to archive ("aug2014.archive")

        species : str
            Ensembl species name (e.g. "hsapiens", "mmusculus")

        # TODO: Rewrite to have loop perform same transformations on transcript and gene-level quantifications
        # TODO: Save all matrices of both levels with clear, consistent naming
        # TODO: Declare saved files and outputs in docstring
        """
        if samples is None:
            samples = [s for s in self.samples if s.library == "RNA-seq"]

        self.count_matrix = self.collect_bitseq_output(
            samples=samples, expression_type=expression_type)

        # Map ensembl gene IDs to gene names
        mapping = query_biomart(attributes=["ensembl_transcript_id", "external_gene_name"])

        # mapping['ensembl_transcript_id'] = mapping['ensembl_transcript_id'].str.replace("\..*", "")
        self.count_matrix = self.count_matrix.reset_index(drop=True, level="ensembl_gene_id")
        # self.count_matrix.index = self.count_matrix.index.str.replace("\..*", "")

        self.count_matrix = self.count_matrix.join(mapping.set_index("ensembl_transcript_id"))
        self.count_matrix.to_csv(
            os.path.join(self.results_dir, self.name + ".expression_{}.csv"
                         .format(expression_type)), index=True)

        # Quantile normalize
        self.matrix_qnorm = normalize_quantiles_p(
            self.count_matrix.reset_index().drop(["ensembl_transcript_id", "gene_name"], axis=1))
        self.matrix_qnorm.index = self.count_matrix.index

        # Log2 TPM
        # see if pseudocount has already been added
        if self.matrix_qnorm.min().min() <= 0:
            pseudocount = 1
        else:
            pseudocount = 0

        self.matrix_qnorm_log = np.log2(
            (pseudocount + self.matrix_qnorm) / self.matrix_qnorm.sum(axis=0)
            * mul_factor)
        self.matrix_qnorm_log.to_csv(
            os.path.join(self.results_dir, self.name +
                         ".expression_{}.transcript_level.quantile_normalized.log2_tpm.csv")
            .format(expression_type, index=False))

        # Reduce to gene-level measurements by max of transcripts
        self.expression_matrix_counts = self.count_matrix.groupby("gene_name").max()

        # Quantile normalize
        self.expression_matrix_qnorm = normalize_quantiles_p(self.expression_matrix_counts)

        # Log2 TPM
        if self.expression_matrix_qnorm.min().min() <= 0:
            pseudocount = 1
        else:
            pseudocount = 0
        self.expression = np.log2(
            ((pseudocount + self.expression_matrix_qnorm) /
                self.expression_matrix_qnorm.sum(axis=0)) * mul_factor)

        # Annotate with sample metadata
        _samples = [s for s in samples if s.name in self.expression.columns]
        attrs = list()
        for attr in sample_attributes:
            ll = list()
            for sample in _samples:  # keep order of samples in matrix
                try:
                    ll.append(getattr(sample, attr))
                except AttributeError:
                    ll.append(np.nan)
            attrs.append(ll)

        # Generate multiindex columns
        index = pd.MultiIndex.from_arrays(attrs, names=sample_attributes)
        self.expression_annotated = self.expression[[s.name for s in _samples]]
        self.expression_annotated.columns = index

        # Save
        self.expression_matrix_counts.to_csv(
            os.path.join(
                self.results_dir, self.name +
                ".expression_counts.gene_level.csv"), index=True)
        self.expression_matrix_qnorm.to_csv(
            os.path.join(
                self.results_dir, self.name +
                ".expression_counts.gene_level.quantile_normalized.csv"), index=True)
        self.expression.to_csv(
            os.path.join(
                self.results_dir, self.name +
                ".expression_counts.gene_level.quantile_normalized.log2_tpm.csv"), index=True)
        self.expression_annotated.to_csv(
            os.path.join(
                self.results_dir, self.name +
                ".expression_counts.gene_level.quantile_normalized.log2_tpm.annotated_metadata.csv"),
            index=True)

    def plot_expression_characteristics(
            self,
            output_dir="{results_dir}/expression_qc",
            output_prefix="expression_qc"):
        """
        Plot general characteristics of the gene expression distributions within and across samples.
        """
        if "{results_dir}" in output_dir:
            output_dir = output_dir.format(results_dir=self.results_dir)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        to_drop = [v for v in ['gene_name', 'ensembl_gene_id', 'ensembl_transcript_id']
                   if v in self.count_matrix.columns.tolist()]

        fig, axis = plt.subplots(figsize=(4, 4 * np.log10(len(self.expression_matrix_counts.columns))))
        sns.barplot(
            data=self.count_matrix.drop(to_drop, axis=1).sum().sort_values().reset_index(),
            y="index", x=0, orient="horiz", color=sns.color_palette("colorblind")[0], ax=axis)
        axis.set_xlabel("Transcriptome reads")
        axis.set_ylabel("Samples")
        sns.despine(fig)
        fig.savefig(
            os.path.join(output_dir, self.name + "{}.expression.reads_per_sample.svg"
                         .format(output_prefix)), bbox_inches="tight")

        cov = pd.DataFrame()
        for i in [1, 2, 3, 6, 12, 24, 48, 96, 200, 300, 400, 500, 1000]:
            cov[i] = self.expression_matrix_counts.apply(lambda x: sum(x >= i))

        fig, axis = plt.subplots(1, 2, figsize=(6 * 2, 6))
        sns.heatmap(
            cov.drop([1, 2], axis=1), ax=axis[0],
            cmap="GnBu", cbar_kws={"label": "Genes covered"})
        sns.heatmap(
            cov.drop([1, 2], axis=1).apply(lambda x: (x - x.mean()) / x.std(), axis=0),
            ax=axis[1], cbar_kws={"label": "Z-score"})
        for ax in axis:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
            ax.set_ylabel("Samples")
            ax.set_xlabel("Genes with #reads")
        axis[1].set_yticklabels(ax.get_yticklabels(), visible=False)
        sns.despine(fig)
        fig.savefig(
            os.path.join(output_dir, self.name +
                         "{}.expression.genes_with_reads.svg".format(output_prefix)),
            bbox_inches="tight")

        for name, matrix in [
                ("counts", self.expression_matrix_counts), ("qnorm_TPM", self.expression)]:
            # Boxplot with values per sample
            if name == "counts":
                matrix = np.log2(1 + matrix)
            to_plot = pd.melt(matrix, var_name="sample")

            fig, axis = plt.subplots()
            sns.boxplot(data=to_plot, y="sample", x="value", orient="horiz", ax=axis)
            axis.set_xlabel("log2({})".format(name))
            axis.set_ylabel("Samples")
            sns.despine(fig)
            fig.savefig(os.path.join(output_dir, self.name +
                        "{}.expression.boxplot_per_sample.{}.svg".format(output_prefix, name)),
                        bbox_inches="tight")

        # # Plot gene expression along chromossomes
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
            self, args, **kwargs):
        """
        RNASeqAnalysis.unsupervised is provided for backward compatibility only and will be removed
        in the future.
        Please use ngs_toolkit.general.unsupervised_analysis(RNASeqAnalysis) in the future.
        """

        _LOGGER.warning(PendingDeprecationWarning(
            "RNASeqAnalysis.unsupervised is provided for backward compatibility "
            "only and will be removed. Please use "
            "RNASeqAnalysis.unsupervised_analysis in the future."))

        self.unsupervised_analysis(args, **kwargs)


def knockout_plot(
        analysis=None,
        knockout_genes=None,
        expression_matrix=None,
        comparison_results=None,
        output_dir=None,
        output_prefix="knockout_expression",
        square=True,
        rasterized=True):
    """
    Plot expression of knocked-out genes in all samples.
    """
    if (analysis is None) and (expression_matrix is None):
        raise AssertionError("One of `analysis` or `expression_matrix` must be provided.")

    msg = "If an `analysis` object is not provided, you must provide a list of `knockout_genes`."
    if (analysis is None) and (knockout_genes is None):
        raise AssertionError(msg)
    elif (analysis is not None) and (knockout_genes is None):
        msg = "If `knockout_genes` is not given, Samples in `analysis` must have a `knockout` attribute."
        try:
            knockout_genes = list(set([s.knockout for s in analysis.samples]))
        except KeyError(msg) as e:
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
    msg = ("The following `knockout_genes` were not found in the expression matrix: '{}'"
           .format(", ".join(missing)))
    if len(missing) > 0:
        _LOGGER.warning(msg)
    knockout_genes = [k for k in knockout_genes if k in expression_matrix.index]

    ko = expression_matrix.loc[knockout_genes, :]
    msg = "None of the `knockout_genes` were found in the expression matrix.\nCannot proceed."
    if ko.empty:
        _LOGGER.warning(msg)
        return
    v = np.absolute(scipy.stats.zscore(ko, axis=1)).flatten().max()
    v += (v / 10.)

    g = sns.clustermap(
        ko, cbar_kws={"label": "Expression"}, robust=True,
        xticklabels=True, yticklabels=True, square=square, rasterized=rasterized)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + ".svg"), bbox_inches="tight")

    g = sns.clustermap(
        ko,
        z_score=0, cmap="RdBu_r", vmin=-v, vmax=v, cbar_kws={"label": "Expression Z-score"},
        rasterized=rasterized, xticklabels=True, yticklabels=True, square=square)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + ".z_score.svg"), bbox_inches="tight")

    g = sns.clustermap(
        ko, cbar_kws={"label": "Expression"}, row_cluster=False, col_cluster=False, robust=True,
        xticklabels=True, yticklabels=True, square=square, rasterized=rasterized)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + ".sorted.svg"), bbox_inches="tight")

    g = sns.clustermap(
        ko,
        z_score=0, cmap="RdBu_r", vmin=-v, vmax=v, cbar_kws={"label": "Expression Z-score"},
        row_cluster=False, col_cluster=False, rasterized=rasterized,
        xticklabels=True, yticklabels=True, square=square)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + ".z_score.sorted.svg"), bbox_inches="tight")

    # p-values and fold-changes for knockout genes
    if comparison_results is None:
        return

    # p-values
    p_table = pd.pivot_table(
        comparison_results.loc[knockout_genes, :].reset_index(),
        index="comparison_name", columns="index", values="padj")
    p_table.index.name = "Knockout gene"
    p_table.columns.name = "Gene"
    p_table = (-np.log10(p_table.loc[knockout_genes, knockout_genes].dropna()))
    p_table = p_table.replace(np.inf, p_table[p_table != np.inf].max().max())
    p_table = p_table.replace(-np.inf, 0)

    v = np.absolute(scipy.stats.zscore(p_table, axis=1)).flatten().max()
    v += (v / 10.)

    g = sns.clustermap(
        p_table, cbar_kws={"label": "-log10(FDR p-value)"},
        xticklabels=True, yticklabels=True, square=square)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + ".p_value.svg"), bbox_inches="tight")

    g = sns.clustermap(
        p_table,
        z_score=0, center=0, cmap="RdBu_r", vmax=v, cbar_kws={"label": "-log10(FDR p-value)\nZ-score"},
        xticklabels=True, yticklabels=True, square=square)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + ".p_value.z_score.svg"), bbox_inches="tight")

    g = sns.clustermap(
        p_table, cbar_kws={"label": "-log10(FDR p-value)"}, row_cluster=False, col_cluster=False,
        xticklabels=True, yticklabels=True, square=square)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + ".p_value.sorted.svg"), bbox_inches="tight")

    g = sns.clustermap(
        p_table,
        z_score=0, cmap="RdBu_r", center=0, vmax=v, cbar_kws={"label": "-log10(FDR p-value)\nZ-score"}, row_cluster=False, col_cluster=False,
        xticklabels=True, yticklabels=True, square=square)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(
        os.path.join(output_dir, output_prefix + ".p_value.z_score.sorted.svg"),
        bbox_inches="tight")

    g = sns.clustermap(
        p_table, cbar_kws={"label": "-log10(FDR p-value)"},
        row_cluster=False, col_cluster=False,
        xticklabels=True, yticklabels=True, square=square, vmax=1.3 * 5)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(
        os.path.join(output_dir, output_prefix + ".p_value.thresholded.svg"),
        bbox_inches="tight")

    # logfoldchanges
    fc_table = pd.pivot_table(
        comparison_results.loc[knockout_genes, :].reset_index(),
        index="comparison_name", columns="index", values="log2FoldChange")
    fc_table.index.name = "Knockout gene"
    fc_table.columns.name = "Gene"
    fc_table = fc_table.loc[knockout_genes, knockout_genes].dropna()

    v = np.absolute(scipy.stats.zscore(fc_table, axis=1)).flatten().max()
    v += (v / 10.)

    g = sns.clustermap(
        fc_table, center=0, cmap="RdBu_r", cbar_kws={"label": "log2(fold-change)"},
        xticklabels=True, yticklabels=True, square=square)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + ".log_fc.svg"), bbox_inches="tight")

    g = sns.clustermap(
        fc_table, center=0, cmap="RdBu_r", cbar_kws={"label": "log2(fold-change)"},
        row_cluster=False, col_cluster=False,
        xticklabels=True, yticklabels=True, square=square)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + ".log_fc.sorted.svg"), bbox_inches="tight")

    g = sns.clustermap(
        fc_table, center=0, vmin=-2, vmax=2, cmap="RdBu_r", cbar_kws={"label": "log2(fold-change)"},
        row_cluster=False, col_cluster=False,
        xticklabels=True, yticklabels=True, square=square)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(
        os.path.join(output_dir, output_prefix + ".log_fc.thresholded.sorted.svg"),
        bbox_inches="tight")


def assess_cell_cycle(analysis, quant_matrix=None, output_dir=None, output_prefix="cell_cycle_assessment"):
    """
    Predict cell cycle phase from expression data.

    :param analysis: [description]
    :type analysis: [type]
    :param quant_matrix: [description], defaults to None
    :type quant_matrix: [type], optional
    :param output_dir: [description], defaults to None
    :type output_dir: [type], optional
    :param output_prefix: [description], defaults to "cell_cycle_assessment"
    :type output_prefix: str, optional
    :returns: [description]
    :rtype: {[type]}
    """
    import anndata
    import scanpy.api as sc

    if quant_matrix is None:
        quant_matrix = analysis.expression

    # quant_matrix = pd.read_csv(os.path.join(
    #         "results",
    #         "arid1a_rnaseq.expression_counts.gene_level.quantile_normalized.log2_tpm.csv"),
    #     index_col=0)
    exp_z = pd.DataFrame(
        zscore(quant_matrix, axis=0),
        index=quant_matrix.index, columns=quant_matrix.columns)

    # Score samples for cell cycle
    cl = "https://raw.githubusercontent.com/theislab/scanpy_usage/master/"
    cl += "180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt"
    cc_prots = pd.read_csv(cl, header=None).squeeze()
    s_genes = cc_prots[:43].tolist()
    g2m_genes = cc_prots[43:].tolist()
    cc_prots = [x for x in cc_prots if x in quant_matrix.index.tolist()]

    knockout_plot(
        knockout_genes=cc_prots,
        output_prefix=output_prefix + ".gene_expression.zscore0",
        expression_matrix=exp_z)

    adata = anndata.AnnData(exp_z.T)
    sc.pp.scale(adata)
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

    preds = adata.obs

    fig, axis = plt.subplots(
        1, 2, figsize=(4 * 2, 4), sharex=True, sharey=True)
    for ax in axis:
        ax.axhline(0, alpha=0.5, linestyle="--", color="black")
        ax.axvline(0, alpha=0.5, linestyle="--", color="black")
        ax.set_xlabel("S-phase score")
        ax.set_ylabel("G2/M-phase score")
    for phase in preds['phase']:
        axis[0].scatter(
            preds.loc[preds['phase'] == phase, 'S_score'],
            preds.loc[preds['phase'] == phase, 'G2M_score'],
            label=phase, alpha=0.5, s=5)
        axis[1].scatter(
            preds.loc[(preds['phase'] == phase) & (preds.index.str.contains("HAP1")), 'S_score'],
            preds.loc[(preds['phase'] == phase) & (preds.index.str.contains("HAP1")), 'G2M_score'],
            label=phase, alpha=0.5, s=5)
    for s in preds.index:
        axis[0].text(preds.loc[s, "S_score"], preds.loc[s, "G2M_score"], s=s, fontsize=6)
    for s in preds.index[preds.index.str.contains("HAP1")]:
        axis[1].text(preds.loc[s, "S_score"], preds.loc[s, "G2M_score"], s=s, fontsize=6)
    axis[0].set_title("All samples")
    axis[1].set_title("HAP1 samples")
    fig.savefig(os.path.join(
            output_dir, output_prefix + ".cell_cycle_prediction.svg"),
        dpi=300, bbox_inches="tight")

    return preds
