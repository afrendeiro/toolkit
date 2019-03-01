#!/usr/bin/env python

import os
import string
import warnings

import matplotlib.pyplot as plt
from ngs_toolkit import _LOGGER
from ngs_toolkit.analysis import Analysis
from ngs_toolkit.general import query_biomart
from ngs_toolkit.general import subtract_principal_component
import numpy as np
import pandas as pd
import pybedtools
import seaborn as sns
from tqdm import tqdm


class CNVAnalysis(Analysis):
    """
    Class to model analysis of CNV data.
    Inherits from the `ngs_toolkit.general.Analysis` class.

    Parameters
    ----------
    name : str, optional
        Name to give analysis object.
        Default is ``analysis`` or if ``from_pep`` the name of the PEP.

    samples : list, optional
        Iterable of peppy.Sample objects use in analysis.
        Default is samples from PEP if ``from_pep``, otherwise an empty list.

    prj : peppy.Project, optional
        Project to tie analysis to.
        Default is the PEP project if ``from_pep``, otherwise None.

    data_dir : str, optional
        Directory containing relevant data for analysis.
        Default is `data`.

    results_dir : str, optional
        Directory to output relevant analysis results.
        Default is `results`.

    pickle_file : str, optional
        File path to use to save serialized object in `pickle` format.
        Default is "{results_dir}/{name}.pickle"

    from_pickle : bool, optional
        If the analysis should be loaded from an existing pickle object.
        Default is `False.

    from_pep : str, optional
        PEP configuration file to initialize analysis from.
        Defaults to None.

    kwargs : dict, optional
        Additional keyword arguments will be passed to parent class `ngs_toolkit.analysis.Analysis`.

    :Example:

    .. code-block:: python
        from ngs_toolkit.cnv import CNVAnalysis

        pep = "metadata/project_config.yaml"
        a = CNVAnalysis(from_pep=pep)

        # Get consensus peak set from all samples
        a.get_cnv_data()

        # Normalize
        a.normalize(method="median")

        # Segmentation
        a.segment_genome()

        # General plots
        a.plot_all_data()
        a.plot_segmentation_stats()

        # Unsupervised analysis
        a.unsupervised_analysis()

        # Save object
        a.to_pickle()
    """
    def __init__(
            self,
            name=None,
            samples=None,
            prj=None,
            data_dir="data",
            results_dir="results",
            pickle_file=None,
            from_pickle=False,
            pep=False,
            **kwargs):

        self.data_type = self.__data_type__ = "CNV"
        self.var_unit_name = "bin"
        self.quantity = "copy_number"
        self.norm_units = "log2(ratio)"

        super(CNVAnalysis, self).__init__(
            name=name,
            data_dir=data_dir,
            results_dir=results_dir,
            pickle_file=pickle_file,
            samples=samples,
            prj=prj,
            from_pickle=from_pickle,
            pep=pep,
            **kwargs)

        if not hasattr(self, "resolutions"):
            self.resolutions = ['1000kb', '100kb', '20kb', '10kb'],

    def load_data(
            self,
            output_map=None,
            only_these_keys=None,
            prefix="{results_dir}/{name}",
            permissive=True):
        """
        Load the output files of the major functions of the Analysis.

        Parameters
        ----------
        output_map : dict
            Dictionary with {attribute_name: (file_path, kwargs)} to load the files.
            The kwargs in the tuple will be passed to pandas.read_csv.
            The default is the required to read the keys in `only_these_keys`.

        only_these_keys : list, optional
            Iterable of analysis attributes to load up.
            Possible attributes:
                "matrix_raw"
                "matrix_norm"
                "matrix_features"
                "differential_results"

        prefix : str, optional
            String prefix of files to load.
            Variables in curly braces will be formated with attributes of analysis.
            Defaults to "{results_dir}/{name}".

        bool : permissive, optional
            Whether an error should be ignored if reading a file causes IOError.
            Default is True.

        Attributes
        ----------
        pandas.DataFrame
            Dataframes holding the respective data, available as attributes described
            in the `only_these_keys` parameter.

        Raises
        ----------
        IOError
            If not permissive and a file is not found
        """
        from ngs_toolkit.utils import fix_dataframe_header

        prefix = self._format_string_with_attributes(prefix)

        if output_map is None:
            kwargs = {"index_col": 0}
            output_map = {
                "matrix_raw":
                    (os.path.join(self.results_dir, self.name + ".{}.matrix_raw.csv"), kwargs),
                "matrix_norm":
                    (os.path.join(self.results_dir, self.name + ".{}.matrix_norm.csv"), kwargs),
                "segmentation":
                    (os.path.join(self.results_dir, self.name + ".{}.segmentation.csv"), {}),
                "segmentation_annot":
                    (os.path.join(self.results_dir, self.name + ".{}.segmentation.annotated.csv"), {})}
        if only_these_keys is None:
            only_these_keys = list(output_map.keys())

        output_map = {k: v for k, v in output_map.items() if k in only_these_keys}

        for resolution in self.resolutions:
            for name, (file, kwargs) in output_map.items():
                file = file.format(resolution)
                _LOGGER.info("Loading '{}' analysis attribute.".format(name))
                if not hasattr(self, name):
                    setattr(self, name, {})
                try:
                    getattr(self, name)[resolution] = pd.read_csv(file, **kwargs)
                    # Fix possible multiindex for matrix_norm
                    if name == "matrix_norm":
                        getattr(self, name)[resolution] = fix_dataframe_header(getattr(self, name)[resolution])
                except IOError as e:
                    if not permissive:
                        raise e
                    else:
                        _LOGGER.warning(e)

    def get_cnv_data(
            self, resolutions=None, samples=None, save=True, assign=True, permissive=False):
        """
        Load CNV data from ATAC-seq CNV pipeline and create CNV matrix at various resolutions.

        Parameters
        ----------
        resolutions : list, optional
            Resolutions of analysis.
            Defaults to resolutions in Analysis object.

        samples : list, optional
            Samples to restrict analysis to.
            Defaults to samples in Analysis object.

        save : bool, optional
            Whether results should be saved to disc.
            Defaults to True

        assign : bool, optional
            Whether results should be assigned to an attribute in the Analsyis object.
            Defaults to True

        permissive : bool, optional
            Whether missing files should be allowed.
            Defaults to False

        Returns
        -------
        dict
            Dictionary with CNV matrices for each resolution.

        Raises
        -------
        IOError
            If not permissive and input files can't be read.

        Attributes
        ----------
        matrix : dict
            Sets a `matrix` dictionary with CNV matrices for each resolution.
        """
        if resolutions is None:
            resolutions = self.resolutions

        if samples is None:
            samples = self.samples

        matrix = dict()

        for resolution in tqdm(resolutions, total=len(resolutions), desc="Resolution"):
            matrix[resolution] = pd.DataFrame()

            for sample in tqdm(samples, total=len(samples), desc="Sample"):
                # Read log2 file
                sample.copywriter_output_file = os.path.join(
                    self.data_dir, sample.name + "_" + resolution,
                    "CNAprofiles", "log2_read_counts.igv")
                try:
                    cov = pd.read_csv(
                        sample.copywriter_output_file,
                        sep="\t", comment="#").set_index("Feature")
                except IOError as e:
                    e = "Sample {} does not have a 'log2_read_counts.igv' file: '{}'.".format(
                        sample.name, sample.copywriter_output_file)
                    if permissive:
                        print(e)
                        continue
                    else:
                        raise IOError(e)

                cov.columns = (
                    cov.columns.str.replace("log2.", "")
                    .str.replace(".trimmed.bowtie2.filtered.bam", "")
                    .str.replace(".merged.sorted.subsample.bam", ""))

                # normalize signal to control
                matrix[resolution][sample.name] = (
                    np.log2(((0.1 + (2 ** cov.loc[:, sample.name])) /
                             (0.1 + (2 ** cov.iloc[:, -1]))))
                )

            matrix[resolution].index = (
                cov['Chromosome'] + ":" +
                cov['Start'].astype(int).astype(str) + "-" +
                cov['End'].astype(int).astype(str))
            matrix[resolution].index.name = "index"

            if save:
                matrix[resolution].to_csv(
                    os.path.join(self.results_dir, self.name + "{}.matrix_raw.csv".format(
                        resolution)),
                    index=True)

        if assign:
            self.matrix = matrix

        return matrix

    def normalize(
            self,
            method="median",
            matrix="matrix_raw", samples=None,
            save=True, assign=True,
            **kwargs):
        """
        Normalization of dictionary of matrices with (n_features, n_samples).

        Parameters
        ----------
        resolutions : list, optional
            Resolutions of analysis.
            Defaults to resolutions in Analysis object.

        method : str
            Normalization method to apply. One of:
                a) `median` (subtract median value across all samples;
                b) `pca` (subtraction of principal component number `pc`).

        pc : int
            Principal Component to remove. 1-based.
            Must be specified if `method=="pca"`.

        matrix : str, optional
            Attribute name of dictionary of matrices to normalize.
            Defaults to `matrix_raw`.

        samples : list
            Iterable of peppy.Sample objects to restrict matrix to.
            Default is all in analysis.

        save : bool, optional
            Whether results should be saved to disc.
            Defaults to True

        assign : bool, optional
            Whether results should be assigned to an attribute in the Analsyis object.
            Defaults to True

        kwargs : dict, optional
            Additional kwargs are passed to the respective normalization method.

        Returns
        -------
        dict
            Dictionary with normalized CNV matrices for each resolution.

        Attributes
        ----------
        matrix_norm : dict
            Sets a `matrix_norm` dictionary with CNV matrices for each resolution.
        """
        matrix_norm = dict()
        if "resolutions" not in kwargs:
            resolutions = self.resolutions

        if matrix is None:
            matrix = self.matrix_raw

        for resolution in resolutions:
            if method == "median":
                matrix_norm[resolution] = self.normalize_median(
                    matrix=matrix[resolution],
                    samples=samples,
                    save=False, assign=False)
            elif method == "pca":
                if "pc" not in kwargs:
                    raise ValueError("If method is 'pca', the value of 'pc' must be given.")
                matrix_norm[resolution] = self.normalize_pca(
                    matrix=matrix[resolution],
                    samples=samples, pc=kwargs["pc"],
                    save=False, assign=False)
            else:
                raise ValueError("Requested method '{}' is not known.".format(method))

            if save:
                matrix_norm[resolution].to_csv(
                    os.path.join(
                        self.results_dir,
                        self.name + "{}.matrix_norm.csv".format),
                    index=True).to_csv(os.path.join(self.results_dir))
        if assign:
            self.matrix_norm = matrix_norm
        return matrix_norm

    def plot_all_data(
            self, matrix=None, resolutions=None, samples=None,
            output_dir=None,
            output_prefix="{analysis_name}.all_data",
            robust=True, vmin=None, vmax=None, rasterized=True, dpi=300,
            sample_labels=True):
        """
        Visualize CNV data genome-wide using heatmaps.
        Will be done independently for each specified resolution.

        Parameters
        ----------
        matrix : str, optional
            Attribute name of dictionary of matrices to normalize.
            Defaults to `matrix_norm`.

        resolutions : list, optional
            Resolutions of analysis.
            Defaults to resolutions in Analysis object.

        samples : list, optional
            Samples to restrict analysis to.
            Defaults to samples in Analysis object.

        output_dir : str, optional
            Output directory.
            Defaults to Analysis results directory.

        output_prefix : str, optional
            Prefix to add to plots.
            Defaults to "{analysis_name}.all_data"

        robust : bool, optional
            Whether to scale the color scale robustly (to quantiles rather than extremes).
            Defaults to True

        vmin : float, optional
            Minimum value of color scale.

        vmax : float, optional
            Maximum value of color scale.
            Defaults to None

        rasterized : bool, optional
            Whether to rasterize main heatmap.
            Defaults to True

        dpi : int, optional
            DPI resolution of rasterized image.
            Defaults to 300

        sample_labels : bool, optional
            Whether to label samples with their name.
            Defaults to True
        """
        # TODO: add support for group colours
        if matrix is None:
            matrix = self.matrix_norm
        if resolutions is None:
            resolutions = self.resolutions
        if samples is None:
            samples = self.samples
        if output_dir is None:
            output_dir = self.results_dir
        if "{analysis_name}" in output_prefix:
            output_prefix = output_prefix.format(analysis_name=self.name)
        names = [s.name for s in samples]

        # Plot mean and variationper chromosome
        for resolution in tqdm(resolutions, desc="Resolution"):
            # Plot all data
            fig, axis = plt.subplots(1, 1, figsize=(4 * 2, 4 * 1))
            sns.heatmap(
                matrix[resolution].T,
                cmap="RdBu_r", center=0, robust=robust, ax=axis, vmin=vmin, vmax=vmax,
                rasterized=rasterized, xticklabels=False,
                yticklabels=names if sample_labels else False,
                cbar_kws={"label": "log2(change)"})
            axis.set_xlabel("Chromosome position")
            axis.set_ylabel("Sample")
            fig.savefig(
                os.path.join(
                    self.results_dir,
                    self.name + "{}.{}.full_data.heatmap.svg".format(resolution, output_prefix)),
                bbox_inches="tight", dpi=dpi)

            grid = sns.clustermap(
                matrix[resolution].fillna(0).T, metric="correlation",
                cmap="RdBu_r", center=0, col_cluster=False,
                robust=robust, vmin=vmin, vmax=vmax,
                rasterized=rasterized, xticklabels=False,
                yticklabels=names if sample_labels else False,
                cbar_kws={"label": "log2(change)"}, figsize=(4 * 2, 4 * 1))
            grid.ax_heatmap.set_xlabel("Chromosome position")
            grid.ax_heatmap.set_ylabel("Sample")
            grid.savefig(
                os.path.join(
                    self.results_dir,
                    self.name + "{}.{}.full_data.fillna.clustermap.svg"
                                .format(resolution, output_prefix)),
                bbox_inches="tight", dpi=dpi)

    def plot_stats_per_chromosome(
            self, matrix="matrix_norm", resolutions=None, samples=None,
            output_dir=None,
            output_prefix="{analysis_name}.all_data",
            robust=True, rasterized=True, dpi=300,
            sample_labels=True):
        """
        Visualize mean and variation of CNV data for each chromosome using heatmaps.
        Will be done independently for each specified resolution.
        Will also be done once for all chromosomes and another time without sex chromosomes.

        Parameters
        ----------
        matrix : str, optional
            Attribute name of dictionary of matrices to normalize.
            Defaults to `matrix_norm`.

        resolutions : list, optional
            Resolutions of analysis.
            Defaults to resolutions in Analysis object.

        samples : list, optional
            Samples to restrict analysis to.
            Defaults to samples in Analysis object.

        output_dir : str, optional
            Output directory.
            Defaults to Analysis results directory.

        output_prefix : str, optional
            Prefix to add to plots.
            Defaults to "{analysis_name}.all_data"

        robust : bool, optional
            Whether to scale the color scale robustly (to quantiles rather than extremes).
            Defaults to True

        rasterized : bool, optional
            Whether to rasterize main heatmap.
            Defaults to True

        dpi : int, optional
            DPI resolution of rasterized image.
            Defaults to 300

        sample_labels : bool, optional
            Whether to label samples with their name.
            Defaults to True
        """
        matrix = self.get_matrix(matrix, samples=samples)
        if resolutions is None:
            resolutions = self.resolutions
        if samples is None:
            samples = self.samples
        if output_dir is None:
            output_dir = self.results_dir
        if "{analysis_name}" in output_prefix:
            output_prefix = output_prefix.format(analysis_name=self.name)

        # Plot mean and variationper chromosome
        for resolution in tqdm(resolutions, desc="Resolution"):
            for label, function in tqdm(
                    [("variation", np.std), ("mean", np.mean)], desc="metric"):
                to_plot = matrix[resolution].copy()
                names = [s.name for s in samples if s.name in to_plot.columns]
                to_plot = to_plot.loc[:, names]
                to_plot['chr'] = map(lambda x: x[0], to_plot.index.str.split(":"))

                p = to_plot[names + ["chr"]].groupby("chr").apply(function)
                grid = sns.clustermap(
                    p + p.min().min(),
                    cbar_kws={"label": label}, cmap="Greens",
                    xticklabels=sample_labels, yticklabels=True, rasterized=rasterized)
                grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0)
                grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), rotation=90)
                grid.savefig(
                    os.path.join(
                        output_dir,
                        self.name + ".{}.".format(resolution) + output_prefix +
                        ".{}_per_chrom.svg".format(label)),
                    bbox_inches="tight", dpi=dpi)

                p = to_plot.loc[
                    ~to_plot['chr'].str.contains("X|Y"),
                    names + ["chr"]].groupby("chr").apply(function)
                grid = sns.clustermap(
                    p + abs(p.min().min()),
                    cbar_kws={"label": label}, cmap="Greens",
                    xticklabels=sample_labels,
                    yticklabels=True, rasterized=rasterized, robust=robust)
                grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0)
                grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), rotation=90)
                grid.savefig(
                    os.path.join(
                        output_dir,
                        self.name + ".{}.".format(resolution) +
                        output_prefix + ".{}_per_chrom.no_sex_chroms.svg".format(label)),
                    bbox_inches="tight", dpi=dpi)

                grid = sns.clustermap(
                    p,
                    cbar_kws={"label": label + " (Z-score)"}, cmap="RdBu_r", center=0, z_score=1,
                    xticklabels=sample_labels, yticklabels=True,
                    rasterized=rasterized, robust=robust)
                grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0)
                grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), rotation=90)
                grid.savefig(
                    os.path.join(
                        output_dir,
                        self.name + ".{}.".format(resolution) +
                        output_prefix + ".{}_per_chrom.no_sex_chroms.zscore.svg".format(label)),
                    bbox_inches="tight", dpi=dpi)

    def segment_genome(
            self, matrix="matrix_norm", resolutions=None, samples=None, save=True, assign=True):
        """
        Segment CNV data to create calls of significant deviations.
        Will be done independently for each specified resolution.

        Requires the R package "DNAcopy" to be installed:
            >>> source('http://bioconductor.org/biocLite.R')
            >>> biocLite('DNAcopy')

        Parameters
        ----------
        matrix : str, optional
            Attribute name of dictionary of matrices to segment.
            Defaults to `matrix_norm`.

        resolutions : list, optional
            Resolutions of analysis.
            Defaults to resolutions in Analysis object.

        samples : list, optional
            Samples to restrict analysis to.
            Defaults to samples in Analysis object.

        save : bool, optional
            Whether results should be saved to disc.
            Defaults to True

        assign : bool, optional
            Whether results should be assigned to an attribute in the Analsyis object.
            Defaults to True

        Returns
        -------
        dict
            Dictionary with segmentation for each resolution.

        Attributes
        ----------
        segmentation : dict
            Dictionary with CNV matrices for each resolution.
        """
        # TODO: implement as_job
        from rpy2.rinterface import RRuntimeWarning
        from rpy2.robjects import numpy2ri, pandas2ri
        import rpy2.robjects as robjects
        warnings.filterwarnings("ignore", category=RRuntimeWarning)
        numpy2ri.activate()
        pandas2ri.activate()

        robjects.r('require("DNAcopy")')
        _CNA = robjects.r('CNA')
        _segment = robjects.r('segment')
        _segments_p = robjects.r('segments.p')
        _segments_summary = robjects.r('segments.summary')
        # _plotSample = robjects.r('plotSample')

        if matrix is None:
            matrix = self.matrix_norm
        if resolutions is None:
            resolutions = self.resolutions
        if samples is None:
            samples = self.samples

        segmentation = dict()
        for resolution in tqdm(resolutions, desc="Resolution"):
            chrom = np.array(map(lambda x: x[0], matrix[resolution].index.str.split(":")))
            start = np.array(map(lambda x: int(x[1].split("-")[0]), matrix[resolution].index.str.split(":")))

            names = [s.name for s in samples]
            cna = _CNA(matrix[resolution].loc[:, names], chrom=chrom, maploc=start)
            cna = _segment(cna)
            seg = _segments_p(cna)
            seg = pd.DataFrame(np.asarray(seg), index=seg.names).T
            summary = _segments_summary(cna)
            summary_seg = pd.DataFrame(np.asarray(summary), index=summary.names).T

            seg = pd.merge(seg, summary_seg)
            for col in seg.columns[1:5]:
                seg[col] = seg[col].replace("NA", np.nan).astype(float).astype(int)
                # columns LCL and UCL are also INT but may contain NA
            for col in seg.columns[5:]:
                seg[col] = seg[col].replace("NA", np.nan).astype(float)

            sample_labels = pd.Series(dict(zip(
                ["Sample.{}".format(i) for i in range(1, len(names) + 1)],
                names)), name="sample_name").to_frame()
            sample_labels.index.name = "ID"

            segmentation[resolution] = pd.merge(sample_labels.reset_index(), seg).drop("ID", axis=1)

            segmentation[resolution].columns = (
                ['sample_name', 'chrom', 'start', 'end', 'bin_size'] +
                ['segment_mean', 'B_stat', 'p_value', "lcl", "ucl", 'segment_std', "segment_median", "segment_mad"])
            segmentation[resolution]["segment_length"] = (
                segmentation[resolution]['end'] -
                segmentation[resolution]['start']) + int(resolution.replace("kb", "")) * 1000

            if save:
                segmentation[resolution].to_csv(
                    os.path.join(self.results_dir, self.name + ".{}.segmentation.csv".format(resolution)),
                    index=False)
        if assign:
            self.segmentation = segmentation

        return segmentation

    def annotate_with_chrom_bands(
            self, segmentation=None, resolutions=None, save=True, assign=True):
        """
        Annotate segmentation with chromosome bands and overlaping genes.
        Will be done independently for each specified resolution.

        Parameters
        ----------
        segmentation : str, optional
            Attribute name of dictionary of segmentation results.
            Defaults to `segmentation`.

        resolutions : list, optional
            Resolutions of analysis.
            Defaults to resolutions in Analysis object.

        samples : list, optional
            Samples to restrict analysis to.
            Defaults to samples in Analysis object.

        save : bool, optional
            Whether results should be saved to disc.
            Defaults to True

        assign : bool, optional
            Whether results should be assigned to an attribute in the Analsyis object.
            Defaults to True

        Returns
        -------
        dict
            Dictionary with annotated segmentation for each resolution.

        Attributes
        ----------
        segmentation_annot : dict
            Dictionary with CNV matrices for each resolution.
        """
        if resolutions is None:
            resolutions = self.resolutions

        if segmentation is None:
            segmentation = self.segmentation

        annotation = query_biomart(
            attributes=["chromosome_name", "start_position", "end_position", "external_gene_name", "band"])

        annotation["chromosome_name"] = "chr" + annotation["chromosome_name"]
        annotation = annotation[~annotation['chromosome_name'].str.contains("_")]
        annotation['band'] = annotation['chromosome_name'] + "_" + annotation['band']
        annotation_bed = pybedtools.BedTool.from_dataframe(annotation.fillna("NA"))

        segmentation_annot = dict()
        for resolution in resolutions:
            seg = segmentation[resolution]
            seg['chrom'] = "chr" + (seg['chrom'].astype(int).astype(str)
                                    .replace("23", "X").replace("24", "Y"))
            seg['start'] = seg['start'].astype(int)
            seg['end'] = seg['end'].astype(int)
            seg_bed = pybedtools.BedTool.from_dataframe(seg[['chrom', 'start', 'end']])

            inter = seg_bed.intersect(annotation_bed, wa=True, wb=True).to_dataframe()
            annot = inter.groupby(['chrom', 'start', 'end'])['thickStart'].apply(
                lambda x: string.join(set(x)))
            annot_band = inter.groupby(['chrom', 'start', 'end'])['thickEnd'].apply(
                lambda x: string.join(set(x)))
            annot_band.name = "chromosome_band"
            annot = annot.to_frame(name="gene").join(annot_band)

            seg_annot = seg.set_index(['chrom', 'start', 'end']).join(annot).reset_index()
            segmentation_annot[resolution] = seg_annot.reindex(
                seg.columns.tolist() + ['chromosome_band', 'gene'], axis=1)

            if save:
                segmentation_annot[resolution].to_csv(
                    os.path.join(
                        self.results_dir,
                        self.name + ".{}.segmentation.annotated.csv".format(resolution)),
                    index=False)

        if assign:
            self.segmentation_annot = segmentation_annot

    def plot_segmentation_stats(
            self, segmentation=None, resolutions=None, per_sample=False,
            output_dir=None, output_prefix="{resolution}.segmentation_metrics"):
        """
        Visualize distribution of statistics of CNV data segmentation.
        Will be done independently for each specified resolution.

        Parameters
        ----------
        segmentation : str, optional
            Dictionary of segmentation results.
            Defaults to `segmentation`.

        resolutions : list, optional
            Resolutions of analysis. Defaults to resolutions in Analysis object.

        per_sample : bool, optional
            Whether plots should be made for each sample too.
            Defaults to False

        output_dir : str, optional
            Output directory.
            Defaults to Analysis results directory.

        output_prefix : str, optional
            Prefix to add to plots.
            Defaults to "{resolution}.segmentation_metrics"
        """
        # TODO: plot stats per sample too
        if segmentation is None:
            segmentation = self.segmentation

        if resolutions is None:
            resolutions = self.resolutions

        if output_dir is None:
            output_dir = self.results_dir

        metric_vars = [
            'bin_size', 'segment_mean', 'B_stat', 'p_value', 'lcl', 'ucl',
            'segment_std', 'segment_median', 'segment_mad', 'segment_length']

        for resolution in resolutions:
            segments = segmentation[resolution]
            segments['log_p_value'] = -np.log10(segments['p_value'].astype(float))

            grid = sns.pairplot(
                segments[metric_vars + ['log_p_value']].dropna(),
                vars=metric_vars + ['log_p_value'],
                plot_kws=dict(s=3, alpha=0.2, linewidth=0, rasterized=True),
                diag_kws=dict(rasterized=True))
            if "{resolution}" in output_prefix:
                op = output_prefix.format(resolution=resolution)
            grid.savefig(
                os.path.join(output_dir, op + ".all_samples.svg"),
                dpi=300, bbox_inches="tight")

        # if per_sample:


def all_to_igv(matrix, output_prefix, **kwargs):
    """
    Convert dictionary of DataFrame with CNV data in several resolutions to IGV format.

    Parameters
    ----------
    matrix : pandas.DataFrame
        DataFrame with CNV data to convert.

    optional output_prefix : str,
        Prefix to add to plots.

    optional kwargs : dict,
        Additional parameters will be passed to ngs_toolkit.cnv.to_igv

    Returns
    -------
    dict
        Dictionary of CNV data in IGV format for each resolution.
    """
    igvs = dict()
    resolutions = matrix.keys()
    for resolution in tqdm(resolutions, total=len(resolutions), desc="Resolution"):
        print("Making IGV visualization for resolution '{}'.".format(resolution))
        igvs[resolution] = to_igv(
            matrix[resolution],
            output_file="{}_{}".format(output_prefix, resolution),
            **kwargs)

    return igvs


def to_igv(matrix, output_file=None, save=True, view_limits=(-2, 2)):
    """
    Convert DataFrame with CNV data to IGV format.

    Parameters
    ----------
    matrix : pandas.DataFrame
        DataFrame with CNV data to convert.

    optional output_file : str,
        Output file.
        Required is `save` is True.

    optional save : bool,
        Whether results should be saved to disc.
        Defaults to True

    optional view_limits : tuple,
        Extreme values (min, max) of color scale used to visualize in IGV.
        Defaults to (-2, 2).

    Returns
    -------
    pandas.DataFrame
        CNV data in IGV format.

    Raises
    -------
    ValueError
    """
    print("Making IGV visualization")

    # as IGV file
    igv = pd.DataFrame(index=matrix.index)
    igv.loc[:, 'Chromosome'] = map(lambda x: x[0], matrix.index.str.split(":"))
    igv.loc[:, 'Start'] = map(lambda x: int(x[1].split("-")[0]), matrix.index.str.split(":"))
    igv.loc[:, 'End'] = map(lambda x: int(x[1].split("-")[1]), matrix.index.str.split(":"))
    igv.loc[:, 'Name'] = igv.index

    igv = igv.join(matrix).reset_index(drop=True)

    if save:
        if output_file is None:
            raise ValueError("If the 'save' option is specified, 'output_file' must also be!")
        open(output_file, 'w')
        output_handle = open(output_file, 'a')
        output_handle.write('#track viewLimits=-{}:{} graphType=heatmap color=255,0,0\n'
                            .format(*view_limits))
        igv.to_csv(output_handle, sep="\t", index=False)
        output_handle.close()

    return igv
