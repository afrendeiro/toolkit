#!/usr/bin/env python


import os

import matplotlib.pyplot as plt
from ngs_toolkit import _LOGGER
from ngs_toolkit.analysis import Analysis
from ngs_toolkit.general import query_biomart
from ngs_toolkit.utils import normalize_quantiles_p
import numpy as np
import pandas as pd
import scipy
from scipy.stats import zscore
import seaborn as sns

# import requests


class RNASeqAnalysis(Analysis):
    """
    Class to model analysis of RNA-seq data.
    Inherits from the `ngs_toolkit.analysis.Analysis` class.

    Parameters
    ----------
    name : str, optional
        Name of the analysis.
        Defaults to ``analysis``.

    from_pep : str, optional
        PEP configuration file to initialize analysis from.
        The analysis will adopt as much attributes from the PEP as possible
        but keyword arguments passed at initialization will still have priority.
        Defaults to None (no PEP used).

    from_pickle : str, optional
        Pickle file of an existing serialized analysis object
        from which the analysis should be loaded.
        Defaults to None (will not load).

    root_dir : str, optional
        Base directory for the project.
        Defaults to current directory or to what is specified in PEP if `from_pep`.

    data_dir : str, optional
        Directory containing processed data (e.g. by looper) that will
        be input to the analysis. This is in principle not required.
        Defaults to ``data``.

    results_dir : str, optional
        Directory to contain outputs produced by the analysis.
        Defaults to ``results``.

    prj : peppy.Project, optional
        A ``peppy.Project`` object that this analysis is tied to.
        Defaults to ``None``.

    samples : list, optional
        List of ``peppy.Sample`` objects that this analysis is tied to.
        Defaults to ``None``.

    kwargs : dict, optional
        Additional keyword arguments will be passed to parent class `ngs_toolkit.analysis.Analysis`.

    """

    def __init__(
        self,
        name=None,
        from_pep=False,
        from_pickle=False,
        root_dir=None,
        data_dir="data",
        results_dir="results",
        prj=None,
        samples=None,
        **kwargs
    ):
        # The check for existance is to make sure other classes can inherit from this
        default_args = {
            "data_type": "RNA-seq",
            "__data_type__": "RNA-seq",
            "var_unit_name": "gene",
            "quantity": "expression",
            "norm_units": "RPM"}
        for k, v in default_args.items():
            if not hasattr(self, k):
                setattr(self, k, v)

        super(RNASeqAnalysis, self).__init__(
            name=name,
            from_pep=from_pep,
            from_pickle=from_pickle,
            root_dir=root_dir,
            data_dir=data_dir,
            results_dir=results_dir,
            prj=prj,
            samples=samples,
            **kwargs
        )

    def collect_bitseq_output(
        self, samples=None, permissive=True, expression_type="counts"
    ):
        """
        Collect gene expression (read counts, transcript-level) output from Bitseq
        into expression matrix for `samples`.
        """
        # TODO: drop support for legacy pipeline output and assume one input file with all required columns
        # TODO: add support for RPKM
        if samples is None:
            samples = self.samples

        if expression_type != "counts":
            raise NotImplementedError("`expression_type` must be 'counts'!")

        expr = list()
        for i, sample in enumerate(samples):
            _LOGGER.debug(
                "Reading transcriptome files for sample '{}'.".format(sample.name)
            )
            tr_file = os.path.join(
                sample.paths.sample_root,
                "bowtie1_{}".format(sample.transcriptome),
                "bitSeq",
                sample.name + ".tr",
            )
            counts_file = os.path.join(
                sample.paths.sample_root,
                "bowtie1_{}".format(sample.transcriptome),
                "bitSeq",
                sample.name + ".counts",
            )

            # read the "tr" file of one sample to get indexes
            try:
                tr = pd.read_csv(
                    tr_file,
                    sep=" ",
                    header=None,
                    skiprows=1,
                    names=["ensembl_gene_id", "ensembl_transcript_id", "v1", "v2"],
                )
            except IOError:
                msg = "Could not open file '{}'' is missing.".format(tr_file)
                if permissive:
                    _LOGGER.warning(msg)
                    continue
                else:
                    raise

            # read the "counts" file of one sample to get indexes
            try:
                e = pd.read_csv(counts_file, sep=" ")
            except IOError:
                msg = "Could not open file '{}'' is missing.".format(counts_file)
                if permissive:
                    _LOGGER.warning(msg)
                    continue
                else:
                    raise

            e = tr.drop(["v1", "v2"], axis=1).join(e)
            e.loc[:, "sample_name"] = sample.name
            expr.append(e)

        if len(expr) == 0:
            msg = "No sample had a valid expression file!"
            if permissive:
                _LOGGER.warning(msg)
                return
            else:
                _LOGGER.error(msg)
                raise IOError(msg)

        expr = (
            pd.concat(expr, axis=0, sort=False)
            .melt(id_vars=["ensembl_gene_id", "ensembl_transcript_id", "sample_name"])
            .pivot_table(
                index=["ensembl_gene_id", "ensembl_transcript_id"],
                columns="sample_name",
                values="value",
                fill_value=0,
            )
            .astype(int, downcast=True)
        )

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
                        sample.name + ".gene.txt",
                    ),
                    sep="\t",
                )
            except IOError(msg) as e:
                if permissive:
                    _LOGGER.warning(e)
                    continue
                else:
                    raise e
            # extract only gene ID and counts
            c = c[["Symbol", "Exp1"]]
            c = c.rename(
                columns={"Symbol": "gene_symbol", "Exp1": sample.name}
            ).set_index("gene_symbol")

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
        self,
        expression_type="counts",
        expression_level="gene",
        samples=None,
        sample_attributes=None,
        genome_assembly=None,
        species=None,
        pseudocount=1,
        mul_factor=1e6,
        save=True,
        assign=True,
    ):
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

        # TODO: Save all matrices of both levels with clear, consistent naming
        # TODO: Declare saved files and outputs in docstring
        # TODO: Pass correct organism, genome to query_biomart
        """
        if expression_type != "counts":
            raise NotImplementedError("`expression_type` must be 'counts'!")

        if expression_level not in ["gene", "transcript"]:
            raise NotImplementedError(
                "`expression_level` must be one of 'gene' or 'transcript'!"
            )

        if samples is None:
            samples = [s for s in self.samples if s.library == "RNA-seq"]

        if sample_attributes is None:
            if hasattr(self, "sample_attributes"):
                sample_attributes = self.sample_attributes

        transcript_counts = self.collect_bitseq_output(
            samples=samples, expression_type=expression_type
        )

        # Map ensembl gene IDs to gene names
        mapping = query_biomart(
            attributes=["ensembl_transcript_id", "external_gene_name"]
        )
        mapping.columns = ["ensembl_transcript_id", "gene_name"]

        # Join gene names to existing Ensembl
        transcript_counts = transcript_counts.reset_index(
            drop=True, level="ensembl_gene_id"
        )
        transcript_counts = (
            transcript_counts.join(mapping.set_index("ensembl_transcript_id"))
            .set_index(["gene_name"], append=True)
            .sort_index(axis=0)
        )

        gene_counts = transcript_counts.groupby("gene_name").max()

        if expression_level == "gene":
            matrix = gene_counts
        elif expression_level == "transcript":
            matrix = transcript_counts

        # Quantile normalize
        matrix_qnorm = normalize_quantiles_p(matrix)
        # Log2 TPM
        # # make matrix non-negative
        if matrix_qnorm.min().min() <= 0:
            matrix_qnorm += np.absolute(matrix_qnorm.min().min())
        # # for some reason bitSeq kind of already adds a pseudocount or something :/
        # # so, if minimum is 1, we'll skip that
        if matrix_qnorm.min().min() == 1:
            pseudocount = 0

        # # Log2 transform
        matrix_qnorm_log = np.log2(pseudocount + matrix_qnorm)

        # Annotate with sample metadata
        _samples = [s for s in samples if s.name in matrix_qnorm_log.columns]
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
        expression_annotated = matrix_qnorm_log[[s.name for s in _samples]]
        expression_annotated.columns = index

        # Save
        if assign:
            self.transcript_counts = transcript_counts
            self.gene_counts = gene_counts
            self.expression = matrix_qnorm_log
            self.expression_annotated = expression_annotated
        if save:
            transcript_counts.to_csv(
                os.path.join(
                    self.results_dir,
                    self.name
                    + ".expression_{}.transcript_level.csv".format(expression_type),
                ),
                index=True,
            )
            gene_counts.to_csv(
                os.path.join(
                    self.results_dir,
                    self.name + ".expression_{}.gene_level.csv".format(expression_type),
                ),
                index=True,
            )
            matrix_qnorm_log.to_csv(
                os.path.join(
                    self.results_dir,
                    self.name
                    + ".expression_{}.{}_level.quantile_normalized.log2_tpm.csv".format(
                        expression_type, expression_level
                    ),
                ),
                index=True,
            )
            expression_annotated.to_csv(
                os.path.join(
                    self.results_dir,
                    self.name
                    + ".expression_{}.{}_level.quantile_normalized.log2_tpm.annotated_metadata.csv".format(
                        expression_type, expression_level
                    ),
                ),
                index=True,
            )

    def plot_expression_characteristics(
        self, output_dir="{results_dir}/expression_qc", output_prefix="expression_qc"
    ):
        """
        Plot general characteristics of the gene expression distributions within and across samples.
        """
        if "{results_dir}" in output_dir:
            output_dir = output_dir.format(results_dir=self.results_dir)

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        to_drop = [
            v
            for v in ["gene_name", "ensembl_gene_id", "ensembl_transcript_id"]
            if v in self.count_matrix.columns.tolist()
        ]

        fig, axis = plt.subplots(
            figsize=(4, 4 * np.log10(len(self.expression_matrix_counts.columns)))
        )
        sns.barplot(
            data=self.count_matrix.drop(to_drop, axis=1)
            .sum()
            .sort_values()
            .reset_index(),
            y="index",
            x=0,
            orient="horiz",
            color=sns.color_palette("colorblind")[0],
            ax=axis,
        )
        axis.set_xlabel("Transcriptome reads")
        axis.set_ylabel("Samples")
        sns.despine(fig)
        fig.savefig(
            os.path.join(
                output_dir,
                self.name + "{}.expression.reads_per_sample.svg".format(output_prefix),
            ),
            bbox_inches="tight",
        )

        cov = pd.DataFrame()
        for i in [1, 2, 3, 6, 12, 24, 48, 96, 200, 300, 400, 500, 1000]:
            cov[i] = self.expression_matrix_counts.apply(lambda x: sum(x >= i))

        fig, axis = plt.subplots(1, 2, figsize=(6 * 2, 6))
        sns.heatmap(
            cov.drop([1, 2], axis=1),
            ax=axis[0],
            cmap="GnBu",
            cbar_kws={"label": "Genes covered"},
        )
        sns.heatmap(
            cov.drop([1, 2], axis=1).apply(lambda x: (x - x.mean()) / x.std(), axis=0),
            ax=axis[1],
            cbar_kws={"label": "Z-score"},
        )
        for ax in axis:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
            ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
            ax.set_ylabel("Samples")
            ax.set_xlabel("Genes with #reads")
        axis[1].set_yticklabels(ax.get_yticklabels(), visible=False)
        sns.despine(fig)
        fig.savefig(
            os.path.join(
                output_dir,
                self.name + "{}.expression.genes_with_reads.svg".format(output_prefix),
            ),
            bbox_inches="tight",
        )

        for name, matrix in [
            ("counts", self.expression_matrix_counts),
            ("qnorm_TPM", self.expression),
        ]:
            # Boxplot with values per sample
            if name == "counts":
                matrix = np.log2(1 + matrix)
            to_plot = pd.melt(matrix, var_name="sample")

            fig, axis = plt.subplots()
            sns.boxplot(data=to_plot, y="sample", x="value", orient="horiz", ax=axis)
            axis.set_xlabel("log2({})".format(name))
            axis.set_ylabel("Samples")
            sns.despine(fig)
            fig.savefig(
                os.path.join(
                    output_dir,
                    self.name
                    + "{}.expression.boxplot_per_sample.{}.svg".format(
                        output_prefix, name
                    ),
                ),
                bbox_inches="tight",
            )

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

    def unsupervised(self, args, **kwargs):
        """
        RNASeqAnalysis.unsupervised is provided for backward compatibility only and will be removed
        in the future.
        Please use ngs_toolkit.general.unsupervised_analysis(RNASeqAnalysis) in the future.
        """

        _LOGGER.warning(
            PendingDeprecationWarning(
                "RNASeqAnalysis.unsupervised is provided for backward compatibility "
                "only and will be removed. Please use "
                "RNASeqAnalysis.unsupervised_analysis in the future."
            )
        )

        self.unsupervised_analysis(args, **kwargs)


def knockout_plot(
    analysis=None,
    knockout_genes=None,
    expression_matrix=None,
    comparison_results=None,
    output_dir=None,
    output_prefix="knockout_expression",
    square=True,
    rasterized=True,
):
    """
    Plot expression of knocked-out genes in all samples.
    """
    if (analysis is None) and (expression_matrix is None):
        raise AssertionError(
            "One of `analysis` or `expression_matrix` must be provided."
        )

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
    msg = "The following `knockout_genes` were not found in the expression matrix: '{}'".format(
        ", ".join(missing)
    )
    if len(missing) > 0:
        _LOGGER.warning(msg)
    knockout_genes = [k for k in knockout_genes if k in expression_matrix.index]

    ko = expression_matrix.loc[knockout_genes, :]
    msg = "None of the `knockout_genes` were found in the expression matrix.\nCannot proceed."
    if ko.empty:
        _LOGGER.warning(msg)
        return
    v = np.absolute(scipy.stats.zscore(ko, axis=1)).flatten().max()
    v += v / 10.0

    g = sns.clustermap(
        ko,
        cbar_kws={"label": "Expression"},
        robust=True,
        xticklabels=True,
        yticklabels=True,
        square=square,
        rasterized=rasterized,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join(output_dir, output_prefix + ".svg"), bbox_inches="tight")

    g = sns.clustermap(
        ko,
        z_score=0,
        cmap="RdBu_r",
        vmin=-v,
        vmax=v,
        cbar_kws={"label": "Expression Z-score"},
        rasterized=rasterized,
        xticklabels=True,
        yticklabels=True,
        square=square,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(
        os.path.join(output_dir, output_prefix + ".z_score.svg"), bbox_inches="tight"
    )

    g = sns.clustermap(
        ko,
        cbar_kws={"label": "Expression"},
        row_cluster=False,
        col_cluster=False,
        robust=True,
        xticklabels=True,
        yticklabels=True,
        square=square,
        rasterized=rasterized,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(
        os.path.join(output_dir, output_prefix + ".sorted.svg"), bbox_inches="tight"
    )

    g = sns.clustermap(
        ko,
        z_score=0,
        cmap="RdBu_r",
        vmin=-v,
        vmax=v,
        cbar_kws={"label": "Expression Z-score"},
        row_cluster=False,
        col_cluster=False,
        rasterized=rasterized,
        xticklabels=True,
        yticklabels=True,
        square=square,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(
        os.path.join(output_dir, output_prefix + ".z_score.sorted.svg"),
        bbox_inches="tight",
    )

    # p-values and fold-changes for knockout genes
    if comparison_results is None:
        return

    # p-values
    p_table = pd.pivot_table(
        comparison_results.loc[knockout_genes, :].reset_index(),
        index="comparison_name",
        columns="index",
        values="padj",
    )
    p_table.index.name = "Knockout gene"
    p_table.columns.name = "Gene"
    p_table = -np.log10(p_table.loc[knockout_genes, knockout_genes].dropna())
    p_table = p_table.replace(np.inf, p_table[p_table != np.inf].max().max())
    p_table = p_table.replace(-np.inf, 0)

    v = np.absolute(scipy.stats.zscore(p_table, axis=1)).flatten().max()
    v += v / 10.0

    g = sns.clustermap(
        p_table,
        cbar_kws={"label": "-log10(FDR p-value)"},
        xticklabels=True,
        yticklabels=True,
        square=square,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(
        os.path.join(output_dir, output_prefix + ".p_value.svg"), bbox_inches="tight"
    )

    g = sns.clustermap(
        p_table,
        z_score=0,
        center=0,
        cmap="RdBu_r",
        vmax=v,
        cbar_kws={"label": "-log10(FDR p-value)\nZ-score"},
        xticklabels=True,
        yticklabels=True,
        square=square,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(
        os.path.join(output_dir, output_prefix + ".p_value.z_score.svg"),
        bbox_inches="tight",
    )

    g = sns.clustermap(
        p_table,
        cbar_kws={"label": "-log10(FDR p-value)"},
        row_cluster=False,
        col_cluster=False,
        xticklabels=True,
        yticklabels=True,
        square=square,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(
        os.path.join(output_dir, output_prefix + ".p_value.sorted.svg"),
        bbox_inches="tight",
    )

    g = sns.clustermap(
        p_table,
        z_score=0,
        cmap="RdBu_r",
        center=0,
        vmax=v,
        cbar_kws={"label": "-log10(FDR p-value)\nZ-score"},
        row_cluster=False,
        col_cluster=False,
        xticklabels=True,
        yticklabels=True,
        square=square,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(
        os.path.join(output_dir, output_prefix + ".p_value.z_score.sorted.svg"),
        bbox_inches="tight",
    )

    g = sns.clustermap(
        p_table,
        cbar_kws={"label": "-log10(FDR p-value)"},
        row_cluster=False,
        col_cluster=False,
        xticklabels=True,
        yticklabels=True,
        square=square,
        vmax=1.3 * 5,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(
        os.path.join(output_dir, output_prefix + ".p_value.thresholded.svg"),
        bbox_inches="tight",
    )

    # logfoldchanges
    fc_table = pd.pivot_table(
        comparison_results.loc[knockout_genes, :].reset_index(),
        index="comparison_name",
        columns="index",
        values="log2FoldChange",
    )
    fc_table.index.name = "Knockout gene"
    fc_table.columns.name = "Gene"
    fc_table = fc_table.loc[knockout_genes, knockout_genes].dropna()

    v = np.absolute(scipy.stats.zscore(fc_table, axis=1)).flatten().max()
    v += v / 10.0

    g = sns.clustermap(
        fc_table,
        center=0,
        cmap="RdBu_r",
        cbar_kws={"label": "log2(fold-change)"},
        xticklabels=True,
        yticklabels=True,
        square=square,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(
        os.path.join(output_dir, output_prefix + ".log_fc.svg"), bbox_inches="tight"
    )

    g = sns.clustermap(
        fc_table,
        center=0,
        cmap="RdBu_r",
        cbar_kws={"label": "log2(fold-change)"},
        row_cluster=False,
        col_cluster=False,
        xticklabels=True,
        yticklabels=True,
        square=square,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(
        os.path.join(output_dir, output_prefix + ".log_fc.sorted.svg"),
        bbox_inches="tight",
    )

    g = sns.clustermap(
        fc_table,
        center=0,
        vmin=-2,
        vmax=2,
        cmap="RdBu_r",
        cbar_kws={"label": "log2(fold-change)"},
        row_cluster=False,
        col_cluster=False,
        xticklabels=True,
        yticklabels=True,
        square=square,
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(
        os.path.join(output_dir, output_prefix + ".log_fc.thresholded.sorted.svg"),
        bbox_inches="tight",
    )


def assess_cell_cycle(
    analysis, matrix=None, output_dir=None, output_prefix="cell_cycle_assessment"
):
    """
    Predict cell cycle phase from expression data.

    :param analysis: [description]
    :type analysis: [type]
    :param matrix: [description], defaults to None
    :type matrix: [type], optional
    :param output_dir: [description], defaults to None
    :type output_dir: [type], optional
    :param output_prefix: [description], defaults to "cell_cycle_assessment"
    :type output_prefix: str, optional
    :returns: [description]
    :rtype: {[type]}
    """
    import anndata
    import scanpy.api as sc

    if matrix is None:
        matrix = analysis.expression

    # matrix = pd.read_csv(os.path.join(
    #         "results",
    #         "arid1a_rnaseq.expression_counts.gene_level.quantile_normalized.log2_tpm.csv"),
    #     index_col=0)
    exp_z = pd.DataFrame(
        zscore(matrix, axis=0), index=matrix.index, columns=matrix.columns
    )

    # Score samples for cell cycle
    cl = "https://raw.githubusercontent.com/theislab/scanpy_usage/master/"
    cl += "180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt"
    cc_prots = pd.read_csv(cl, header=None).squeeze()
    s_genes = cc_prots[:43].tolist()
    g2m_genes = cc_prots[43:].tolist()
    cc_prots = [x for x in cc_prots if x in matrix.index.tolist()]

    knockout_plot(
        knockout_genes=cc_prots,
        output_prefix=output_prefix + ".gene_expression.zscore0",
        expression_matrix=exp_z,
    )

    adata = anndata.AnnData(exp_z.T)
    sc.pp.scale(adata)
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

    preds = adata.obs

    fig, axis = plt.subplots(1, 2, figsize=(4 * 2, 4), sharex=True, sharey=True)
    for ax in axis:
        ax.axhline(0, alpha=0.5, linestyle="--", color="black")
        ax.axvline(0, alpha=0.5, linestyle="--", color="black")
        ax.set_xlabel("S-phase score")
        ax.set_ylabel("G2/M-phase score")
    for phase in preds["phase"]:
        axis[0].scatter(
            preds.loc[preds["phase"] == phase, "S_score"],
            preds.loc[preds["phase"] == phase, "G2M_score"],
            label=phase,
            alpha=0.5,
            s=5,
        )
        axis[1].scatter(
            preds.loc[
                (preds["phase"] == phase) & (preds.index.str.contains("HAP1")),
                "S_score",
            ],
            preds.loc[
                (preds["phase"] == phase) & (preds.index.str.contains("HAP1")),
                "G2M_score",
            ],
            label=phase,
            alpha=0.5,
            s=5,
        )
    for s in preds.index:
        axis[0].text(
            preds.loc[s, "S_score"], preds.loc[s, "G2M_score"], s=s, fontsize=6
        )
    for s in preds.index[preds.index.str.contains("HAP1")]:
        axis[1].text(
            preds.loc[s, "S_score"], preds.loc[s, "G2M_score"], s=s, fontsize=6
        )
    axis[0].set_title("All samples")
    axis[1].set_title("HAP1 samples")
    fig.savefig(
        os.path.join(output_dir, output_prefix + ".cell_cycle_prediction.svg"),
        dpi=300,
        bbox_inches="tight",
    )

    return preds
