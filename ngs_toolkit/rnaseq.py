#!/usr/bin/env python


import os

import numpy as np
import pandas as pd
from ngs_toolkit import _LOGGER
from ngs_toolkit.analysis import Analysis


class RNASeqAnalysis(Analysis):
    """
    Class to model analysis of RNA-seq data.
    Inherits from the :class:`~ngs_toolkit.analysis.Analysis` class.

    Parameters
    ----------
    name : :obj:`str`, optional
        Name of the analysis.

        Defaults to "analysis".
    from_pep : :obj:`str`, optional
        PEP configuration file to initialize analysis from.
        The analysis will adopt as much attributes from the PEP as possible
        but keyword arguments passed at initialization will still have priority.

        Defaults to :obj:`None` (no PEP used).
    from_pickle : :obj:`str`, optional
        Pickle file of an existing serialized analysis object
        from which the analysis should be loaded.

        Defaults to :obj:`None` (will not load from pickle).
    root_dir : :obj:`str`, optional
        Base directory for the project.

        Defaults to current directory or to what is specified in PEP if
        :attr:`~ngs_toolkit.analysis.Analysis.from_pep`.
    data_dir : :obj:`str`, optional
        Directory containing processed data (e.g. by looper) that will
        be input to the analysis. This is in principle not required.

        Defaults to "data".
    results_dir : :obj:`str`, optional
        Directory to contain outputs produced by the analysis.

        Defaults to "results".
    prj : :class:`peppy.Project`, optional
        A :class:`peppy.Project` object that this analysis is tied to.

        Defaults to :obj:`None`.
    samples : :obj:`list`, optional
        List of :class:`peppy.Sample` objects that this analysis is tied to.

        Defaults to :obj:`None`.
    kwargs : :obj:`dict`, optional
        Additional keyword arguments will be passed to parent class
        :class:`~ngs_toolkit.analysis.Analysis`.
    """

    _data_type = "RNA-seq"

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
            "norm_units": "RPM",
        }
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

    def collect_bitseq_output(self, samples=None, permissive=True, expression_type="counts"):
        """
        Collect gene expression (read counts, transcript-level) output from Bitseq
        into expression matrix for `samples`.
        """
        # TODO: drop support for legacy pipeline output and assume one input file with all required columns
        # TODO: add support for RPKM
        if samples is None:
            samples = self.samples

        if expression_type != "counts":
            msg = "`expression_type` must be 'counts'!"
            _LOGGER.error(msg)
            raise NotImplementedError(msg)

        expr = list()
        for i, sample in enumerate(samples):
            _LOGGER.debug("Reading transcriptome files for sample '{}'.".format(sample.name))
            tr_file = os.path.join(
                sample.sample_root,
                "bowtie1_{}".format(sample.transcriptome),
                "bitSeq",
                sample.name + ".tr",
            )
            counts_file = os.path.join(
                sample.sample_root,
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
        for sample in samples:
            try:
                c = pd.read_csv(
                    os.path.join(
                        sample.sample_root,
                        "ESAT_{}".format(sample.genome),
                        sample.name + ".gene.txt",
                    ),
                    sep="\t",
                )
            except IOError:
                if permissive:
                    _LOGGER.warning("Sample '%s' is missing file: %s", sample.name, sample.counts)
                    continue
                else:
                    raise
            # extract only gene ID and counts
            c = c[["Symbol", "Exp1"]]
            c = c.rename(columns={"Symbol": "gene_symbol", "Exp1": sample.name}).set_index(
                "gene_symbol"
            )

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
    #                     sample.sample_root,
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
        reduction_func=max,
        quantification_prog="bitseq",
        samples=None,
        save=True,
        assign=True,
        output_file=None,
        permissive=False,
        species=None,
        ensembl_version=None,
    ):
        """
        Collect gene expression (read counts per transcript or gene) for all samples.

        If `expression_level` is "gene", then, transcripts will be reduced per gene ID
        using `reduction_func` (defaults to `max`) and features will be named with gene symbols.

        Parameters
        ----------
        expression_type : :obj:`str`, optional
            Type of expression quantification to get. One of "counts" or "rpkm".

            Defaults to "counts".
        expression_level : :obj:`str`, optional
            Type of expression quantification to get. One of "transcript" or "gene".

            Defaults to "gene".
        reduction_func : func, optional
            Function to reduce gene expression between transcript and gene if `expression_level` is "gene".

            Defaults to `max`.
        quantification_prog : :obj:`str`, optional
            Name of program used to produce the quantification of gene expression. One of "bitseq", "htseq" or "esat".

            Defaults to "bitseq".
        samples : list[peppy.Sample], optional
            Subset of samples to get expression for.

            Defaults to all in analysis.
        save: :obj:`bool`, optional
            Whether to save output as CSV.

            Default is :obj:`None`.
        assign: :obj:`bool`, optional
            Whether to assign output to `matrix_raw`.

            Default is :obj:`None`.
        output_file : :obj:`str`, optional
            Path of resulting file if `save` is `True`.

            Defaults to "{results_dir}/{name}.matrix_raw.csv".
        permissive: :obj:`bool`, optional
            Whether to skip samples with non-existing gene expression quantification.

            Default is `False`.
        species : :obj:`str`, optional
            Ensembl species name (e.g. "hsapiens", "mmusculus")

            Defaults to analysis' organism.
        ensembl_version : :obj:`str`, optional
            Ensembl version of annotation to use (e.g. "grch38", "grcm38")

            Defaults to analysis' genome.

        Attributes
        ----------
        matrix_raw : :obj:`pandas.DataFrame`
            DataFrame with gene expression.
        """
        from ngs_toolkit.general import query_biomart
        from ngs_toolkit import constants

        if expression_type != "counts":
            msg = "`expression_type` must be 'counts'!"
            _LOGGER.error(msg)
            raise NotImplementedError(msg)

        if expression_level not in ["gene", "transcript"]:
            raise NotImplementedError("`expression_level` must be one of 'gene' or 'transcript'!")

        if samples is None:
            samples = self.samples

        # Check which samples to use (dependent on permissive)
        samples = self._get_samples_with_input_file(
            "counts", permissive=permissive, samples=samples
        )

        if quantification_prog == "bitseq":
            transcript_counts = self.collect_bitseq_output(
                samples=samples, expression_type=expression_type
            )
        else:
            msg = "Only implemented for `quantification_prog`='bitseq'"
            _LOGGER.error(msg)
            raise NotImplementedError(msg)

        # Map ensembl gene IDs to gene names
        mapping = query_biomart(
            attributes=["ensembl_transcript_id", "external_gene_name"],
            species=species or constants.organism_to_species_mapping[self.organism],
            ensembl_version=ensembl_version or constants.genome_to_ensembl_mapping[self.genome],
        )
        mapping.columns = ["ensembl_transcript_id", "gene_name"]

        # Join gene names to existing Ensembl
        transcript_counts = transcript_counts.reset_index(drop=True, level="ensembl_gene_id")
        transcript_counts = (
            transcript_counts.join(mapping.set_index("ensembl_transcript_id"))
            .set_index(["gene_name"], append=True)
            .sort_index(axis=0)
        )
        if expression_level == "transcript":
            matrix_raw = transcript_counts
        elif expression_level == "gene":
            if reduction_func is max:
                matrix_raw = transcript_counts.groupby("gene_name").max()
            else:
                matrix_raw = transcript_counts.groupby("gene_name").apply(reduction_func)

        if save:
            if output_file is not None:
                matrix_raw.to_csv(output_file)
            else:
                matrix_raw.to_csv(
                    os.path.join(self.results_dir, self.name + ".matrix_raw.csv"), index=True,
                )
        if assign:
            self.matrix_raw = matrix_raw
        return matrix_raw

    def plot_expression_characteristics(
        self,
        matrix_raw=None,
        matrix_norm=None,
        samples=None,
        output_dir="{results_dir}/quality_control",
        output_prefix="quality_control",
    ):
        """
        Plot general characteristics of the gene expression
        distributions within and across samples.

        matrix_raw : {str, pandas.DataFrame}, optional
            Name of analysis attribute with raw expression
            values or pandas dataframe.

            Defaults to analysis' `matrix_raw`.
        matrix_norm : {str, pandas.DataFrame}, optional
            Name of analysis attribute with normalized expression
            values or pandas dataframe.

            Defaults to analysis' `matrix_norm`.
        samples : :obj:`list`, optional
            List of samples to include.

            Defaults to all samples in analysis
        output_dir : :obj:`str`, optional
            Directory for output files.

            Defaults to "{results_dir}/quality_control"
        output_prefix : :obj:`str`, optional
            Prefix for output files.

            Defaults to "quality_control"
        """
        # TODO: Plot gene expression along chromossomes
        import matplotlib.pyplot as plt
        import seaborn as sns

        output_dir = self._format_string_with_attributes(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        if matrix_raw is None:
            matrix_raw = self.get_matrix("matrix_raw", samples=samples)
        else:
            if samples is not None:
                matrix_raw = matrix_raw.loc[:, [s.name for s in samples]]

        if matrix_norm is None:
            matrix_norm = self.get_matrix("matrix_raw", samples=samples)
        else:
            if samples is not None:
                matrix_norm = matrix_norm.loc[:, [s.name for s in samples]]

        fig, axis = plt.subplots(figsize=(4, 4 * np.log10(len(matrix_raw.columns))))
        sns.barplot(
            data=matrix_raw.sum().sort_values().reset_index(),
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
                output_dir, self.name + "{}.expression.reads_per_sample.svg".format(output_prefix),
            ),
            bbox_inches="tight",
        )

        cov = pd.DataFrame()
        for i in [1, 2, 3, 6, 12, 24, 48, 96, 200, 300, 400, 500, 1000]:
            cov[i] = matrix_raw.apply(lambda x: sum(x >= i))

        fig, axis = plt.subplots(1, 2, figsize=(6 * 2, 6))
        sns.heatmap(
            cov.drop([1, 2], axis=1), ax=axis[0], cmap="GnBu", cbar_kws={"label": "Genes covered"},
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
                output_dir, self.name + "{}.expression.genes_with_reads.svg".format(output_prefix),
            ),
            bbox_inches="tight",
        )

        for name, matrix in [
            ("counts", matrix_raw),
            ("normalized", self.matrix_norm),
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
                    + "{}.expression.boxplot_per_sample.{}.svg".format(output_prefix, name),
                ),
                bbox_inches="tight",
            )


def plot_features(
    analysis=None,
    knockout_genes=None,
    matrix="matrix_norm",
    samples=None,
    differential_results=None,
    output_dir=None,
    output_prefix="knockout_expression",
):
    """
    Plot expression of genes in samples or sample groups.

    Parameters
    ----------

    analysis : :class:`ngs_toolkit.RNASeqAnalysis`, optional
        Analysis object.

        Not required if `matrix` is given.
    knockout_genes : :obj:`list`, optional
        List of perturbed genes to plot.

        Defaults to the set of `knockout` attributes in the analysis'
        samples if `analysis` is given. Otherwise must be given.
    matrix : str, optional
        Matrix with expression values to use.

        Defaults to "matrix_norm"
    samples : [type], optional
        [description]

        Defaults to :obj:`None`.
    differential_results : [type], optional
        [description]

        Defaults to :obj:`None`.
    output_dir : [type], optional
        [description]

        Defaults to :obj:`None`.
    output_prefix : str, optional
        Prefix for output files.

        Defaults to "knockout_expression"
    """
    from ngs_toolkit.graphics import clustermap_varieties

    if (analysis is None) and (matrix is None):
        raise AssertionError("One of `analysis` or `matrix` must be provided.")

    msg = "If an `analysis` object is not provided, you must provide a list of `knockout_genes`."
    if (analysis is None) and (knockout_genes is None):
        raise AssertionError(msg)
    elif (analysis is not None) and (knockout_genes is None):
        msg = "If `knockout_genes` is not given, Samples in `analysis` must have a `knockout` attribute."
        try:
            knockout_genes = list(set([s.knockout for s in analysis.samples]))
        except KeyError(msg) as e:
            raise e

    matrix = analysis.get_matrix(matrix=matrix, samples=samples)

    if output_dir is None:
        if analysis is not None:
            output_dir = analysis.results_dir
        else:
            output_dir = os.path.curdir

    knockout_genes = sorted(knockout_genes)

    missing = [k for k in knockout_genes if k not in matrix.index]
    msg = "Some `knockout_genes` were not found in the expression matrix: '%s'"
    if len(missing) > 0:
        _LOGGER.warning(msg % ", ".join(missing))
    knockout_genes = [k for k in knockout_genes if k in matrix.index]

    ko = matrix.loc[knockout_genes, :]
    msg = "None of the `knockout_genes` were found in the expression matrix.\nCannot proceed."
    if ko.empty:
        _LOGGER.warning(msg)
        return

    # expression values
    clustermap_varieties(ko, output_dir=output_dir, output_prefix=output_prefix)

    # p-values and fold-changes for knockout genes
    if differential_results is None:
        differential_results = getattr(analysis, "differential_results", None)
    if differential_results is None:
        return

    if len(differential_results["comparison_name"].unique()) <= 1:
        msg = "Could not plot values per comparison as only one found!"
        _LOGGER.warning(msg)
        return

    # p-values
    p_table = pd.pivot_table(
        differential_results.loc[knockout_genes, :].reset_index(),
        index="comparison_name",
        columns="index",
        values="padj",
    )
    p_table.index.name = "Knockout gene"
    p_table.columns.name = "Gene"
    p_table = -np.log10(p_table.loc[:, knockout_genes].dropna())
    p_table = p_table.replace(np.inf, p_table[p_table != np.inf].max().max())
    p_table = p_table.replace(-np.inf, 0)

    clustermap_varieties(
        p_table,
        output_dir=output_dir,
        output_prefix=output_prefix + ".p_value",
        quantity="-log10(FDR p-value)",
    )
    clustermap_varieties(
        p_table,
        output_dir=output_dir,
        output_prefix=output_prefix + ".p_value.thresholded",
        steps=["base", "sorted"],
        quantity="-log10(FDR p-value)",
        vmax=1.3 * 5,
    )

    # logfoldchanges
    fc_table = pd.pivot_table(
        differential_results.loc[knockout_genes, :].reset_index(),
        index="comparison_name",
        columns="index",
        values="log2FoldChange",
    )
    fc_table.index.name = "Knockout gene"
    fc_table.columns.name = "Gene"
    fc_table = fc_table.loc[:, knockout_genes].dropna()

    clustermap_varieties(
        fc_table,
        output_dir=output_dir,
        output_prefix=output_prefix + "log_fc",
        steps=["base", "sorted"],
        quantity="log2(fold-change)",
    )
    clustermap_varieties(
        fc_table,
        output_dir=output_dir,
        output_prefix=output_prefix + "log_fc.thresholded",
        steps=["base", "sorted"],
        quantity="log2(fold-change)",
        vmin=-2,
        vmax=2,
    )


def assess_cell_cycle(
    analysis, matrix=None, output_dir=None, output_prefix="cell_cycle_assessment"
):
    """
    Predict cell cycle phase from expression data.
    """
    import anndata
    import scanpy.api as sc
    import matplotlib.pyplot as plt
    from scipy.stats import zscore

    if matrix is None:
        matrix = analysis.expression

    # matrix = pd.read_csv(os.path.join(
    #         "results",
    #         "arid1a_rnaseq.expression_counts.gene_level.quantile_normalized.log2_tpm.csv"),
    #     index_col=0)
    exp_z = pd.DataFrame(zscore(matrix, axis=0), index=matrix.index, columns=matrix.columns)

    # Score samples for cell cycle
    cl = "https://raw.githubusercontent.com/theislab/scanpy_usage/master/"
    cl += "180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt"
    cc_prots = pd.read_csv(cl, header=None).squeeze()
    s_genes = cc_prots[:43].tolist()
    g2m_genes = cc_prots[43:].tolist()
    cc_prots = [x for x in cc_prots if x in matrix.index.tolist()]

    plot_features(
        knockout_genes=cc_prots,
        output_prefix=output_prefix + ".gene_expression.zscore0",
        matrix=exp_z,
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
            preds.loc[(preds["phase"] == phase) & (preds.index.str.contains("HAP1")), "S_score",],
            preds.loc[(preds["phase"] == phase) & (preds.index.str.contains("HAP1")), "G2M_score",],
            label=phase,
            alpha=0.5,
            s=5,
        )
    for s in preds.index:
        axis[0].text(preds.loc[s, "S_score"], preds.loc[s, "G2M_score"], s=s, fontsize=6)
    for s in preds.index[preds.index.str.contains("HAP1")]:
        axis[1].text(preds.loc[s, "S_score"], preds.loc[s, "G2M_score"], s=s, fontsize=6)
    axis[0].set_title("All samples")
    axis[1].set_title("HAP1 samples")
    fig.savefig(
        os.path.join(output_dir, output_prefix + ".cell_cycle_prediction.svg"),
        dpi=300,
        bbox_inches="tight",
    )

    return preds
