#!/usr/bin/env python


import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm

from ngs_toolkit.general import Analysis


class CNVAnalysis(Analysis):
    """
    Class to model analysis of CNV data.
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
            resolutions=['1000kb', '100kb', '20kb', '10kb'],
            **kwargs):
        super(CNVAnalysis, self).__init__(
            name=name,
            data_dir=data_dir,
            results_dir=results_dir,
            pickle_file=pickle_file,
            samples=samples,
            prj=prj,
            from_pickle=from_pickle,
            **kwargs)

        self.data_type = self.__data_type__ = "CNV"
        self._var_names = "region"
        self._quantity = "copy_number"
        self._norm_units = "log2(ratio)"
        self._raw_matrix_name = "coverage"
        self._norm_matrix_name = "coverage_norm"
        self._annot_matrix_name = "cnv"

        self.resolutions = resolutions

    def load_data(self, resolutions=None, permissive=True, norm_method="median"):
        """
        Load data from a previous analysis.

        :param resolutions: Resolutions of analysis, defaults to resolutions in
                            Analysis object.
        :type resolutions: list, optional
        :param permissive: Whether missing files should be allowed. Defaults to True
        :type permissive: bool, optional
        :param norm_method: Normalization method used in the previous files.
                            Defaults to "median".
        :type norm_method: str, optional
        :raises: ValueError, e
        """
        if resolutions is None:
            resolutions = self.resolutions

        if norm_method == "median":
            suf = "coverage_median"
        elif norm_method == "pca":
            suf = "coverage_pcanorm"
        else:
            raise ValueError("Normalization method unknown!")

        self.coverage = dict()
        self.coverage_norm = dict()
        self.segmentation = dict()
        self.segmentation_annot = dict()
        to_load = [
            ("coverage",
                os.path.join(self.results_dir, self.name + ".{}.raw_coverage.csv"),
                "Raw coverage", {"index_col": 0}),
            ("coverage_norm",
                os.path.join(self.results_dir, self.name + ".{}" + ".{}.csv".format(suf)),
                "Normalized coverage", {"index_col": 0}),
            ("segmentation",
                os.path.join(self.results_dir, self.name + ".{}.segmentation.csv"),
                "Segmentation", {}),
            ("segmentation_annot",
                os.path.join(self.results_dir, self.name + ".{}.segmentation.annotated.csv"),
                "Annotated segmentation", {})]
        for resolution in self.resolutions:
            for attr, f, desc, params in to_load:
                f = f.format(resolution)
                try:
                    getattr(self, attr)[resolution] = pd.read_csv(f, **params)
                except IOError as e:
                    if permissive:
                        print("{} file '{}' could not be read.".format(desc, f))
                    else:
                        raise e

    def get_cnv_data(
            self, resolutions=None, samples=None, save=True, assign=True, permissive=False):
        """
        Load CNV data from ATAC-seq CNV pipeline and create CNV matrix at various resolutions.

        :param resolutions: Resolutions of analysis. Defaults to resolutions in
                            Analysis object.
        :type resolutions: list, optional
        :param samples: Samples to restrict analysis to. Defaults to samples in
                        Analysis object.
        :type samples: list, optional
        :param save: Whether results should be saved to disc. Defaults to True
        :type save: bool, optional
        :param assign: Whether results should be assigned to an attribute in the
                       Analsyis object. Defaults to True
        :type assign: bool, optional
        :param permissive: :param permissive: Whether missing files should be allowed.
                                              Defaults to False
        :type permissive: bool, optional
        :returns: Dictionary with CNV matrices for each resolution.
        :rtype: dict
        :raises: IOError

        :var dict coverage: Sets a `coverage` dictionary with CNV matrices for
                            each resolution.
        """
        if resolutions is None:
            resolutions = self.resolutions

        if samples is None:
            samples = self.samples

        coverage = dict()

        for resolution in tqdm(resolutions, total=len(resolutions), desc="Resolution"):
            coverage[resolution] = pd.DataFrame()

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
                coverage[resolution][sample.name] = (
                    np.log2(((0.1 + (2 ** cov.loc[:, sample.name])) /
                             (0.1 + (2 ** cov.iloc[:, -1]))))
                )

            coverage[resolution].index = (
                cov['Chromosome'] + ":" +
                cov['Start'].astype(int).astype(str) + "-" +
                cov['End'].astype(int).astype(str))
            coverage[resolution].index.name = "index"

            if save:
                coverage[resolution].to_csv(
                    os.path.join(self.results_dir, self.name + "{}.raw_coverage.csv".format(
                        resolution)),
                    index=True)

        if assign:
            self.coverage = coverage

        return coverage

    def normalize_coverage_median(
            self, matrix=None, resolutions=None,
            samples=None, function=np.nanmedian, fillna=True,
            save=True, assign=True):
        """
        Normalization of matrices of (n_features, n_samples) by subtracting the
        median from each sample/feature.

        :param matrix: Attribute name of dictionary of matrices to normalize.
                       Defaults to `coverage`.
        :type matrix: str, optional
        :param resolutions: Resolutions of analysis. Defaults to resolutions in
                            Analysis object.
        :type resolutions: list, optional
        :param samples: Samples to restrict analysis to. Defaults to samples in
                        Analysis object.
        :type samples: list, optional
        :param function: An alternative function to calculate across samples.
                         Data will be subtracted by this. Defaults to numpy.nanmedian
        :type function: function, optional
        :param fillna: Whether to fill NaN with zero. Defaults to True
        :type fillna: bool, optional
        :param save: Whether results should be saved to disc. Defaults to True
        :type save: bool, optional
        :param assign: Whether results should be assigned to an attribute in the
                       Analsyis object. Defaults to True
        :type assign: bool, optional

        :returns: Dictionary with normalized CNV matrices for each resolution.
        :rtype: dict
        :var dict coverage_norm: Sets a `coverage_norm` dictionary with CNV matrices
                                 for each resolution.
        """
        if matrix is None:
            matrix = self.coverage
        if resolutions is None:
            resolutions = self.resolutions
        if samples is None:
            samples = self.samples

        coverage_norm = dict()
        for resolution in tqdm(resolutions, desc="Resolution"):
            to_norm = matrix[resolution].loc[
                :, [s.name for s in samples if s.name in matrix[resolution].columns]]
            coverage_norm[resolution] = (to_norm.T - function(to_norm, axis=1)).T
            if fillna:
                coverage_norm[resolution] = coverage_norm[resolution].fillna(0)
            if save:
                coverage_norm[resolution].to_csv(
                    os.path.join(
                        self.results_dir,
                        self.name + "{}.coverage_median.csv".format(resolution)),
                    index=True)
        if assign:
            self.coverage_norm = coverage_norm

        return coverage_norm

    def normalize_coverage_pca(
            self, matrix=None, resolutions=None,
            samples=None, pc=None,
            save=True, assign=True):
        """
        Normalization of matrices of (n_features, n_samples) by subtracting the
        contribution of a Principal Component from each sample/feature.

        :param matrix: Attribute name of dictionary of matrices to normalize.
                       Defaults to `coverage`.
        :type matrix: str, optional
        :param resolutions: Resolutions of analysis. Defaults to resolutions in
                            Analysis object.
        :type resolutions: list, optional
        :param samples: Samples to restrict analysis to. Defaults to samples in
                        Analysis object.
        :type samples: list, optional
        :param pc: Principal Component to remove. 1-based. Must be specified.
        :type pc: int
        :param save: Whether results should be saved to disc. Defaults to True
        :type save: bool, optional
        :param assign: Whether results should be assigned to an attribute in the
                       Analsyis object. Defaults to True
        :type assign: bool, optional

        :returns: Dictionary with normalized CNV matrices for each resolution.
        :rtype: dict
        :var dict coverage_norm: Sets a `coverage_norm` dictionary with CNV matrices
                                 for each resolution.
        """
        from ngs_toolkit.general import subtract_principal_component

        if pc is None:
            raise ValueError("Principal Component to remove must be specified!")

        if matrix is None:
            matrix = self.coverage
        if resolutions is None:
            resolutions = self.resolutions
        if samples is None:
            samples = self.samples

        # first make sure data is centered
        to_norm = self.normalize_coverage_median(
            matrix, resolutions=resolutions, samples=samples, save=False, assign=False)

        coverage_norm = dict()
        for resolution in tqdm(resolutions, desc="Resolution"):
            # then remove the PC
            coverage_norm[resolution] = subtract_principal_component(
                to_norm[resolution].T.fillna(0), pc=pc).T

            if save:
                coverage_norm[resolution].to_csv(
                    os.path.join(
                        self.results_dir,
                        self.name + "{}.coverage_pcanorm.csv".format(resolution)),
                    index=True)
        if assign:
            self.coverage_norm = coverage_norm

        return coverage_norm

    def normalize(
            self, resolutions=None, method="median",
            pc=None, matrix=None, samples=None, save=True, assign=True):
        """
        Normalization of dictionary of matrices with (n_features, n_samples).

        :param resolutions: Resolutions of analysis. Defaults to resolutions in
                            Analysis object.
        :type resolutions: list, optional
        :param method: Normalization method to apply. One of: a) `median`
                       (subtract median value across all samples; b) `pca`
                       (subtraction of principal component number `pc`).
        :type method: str
        :param pc: Principal Component to remove. 1-based. Must be specified if
                   `method=="pca"`.
        :type pc: int
        :param matrix: Attribute name of dictionary of matrices to normalize.
                       Defaults to `coverage`.
        :type matrix: str, optional
        :param samples: Iterable of peppy.Sample objects to restrict matrix to.
                        If not provided (`None` is passed) the matrix will not be subsetted.
        :type samples: list
        :param save: Whether results should be saved to disc. Defaults to True
        :type save: bool, optional
        :param assign: Whether results should be assigned to an attribute in the
                       Analsyis object. Defaults to True
        :type assign: bool, optional

        :returns: Dictionary with normalized CNV matrices for each resolution.
        :rtype: dict
        :var dict coverage_norm: Sets a `coverage_norm` dictionary with CNV matrices
                                 for each resolution.
        """
        if method == "median":
            return self.normalize_coverage_median(
                matrix=matrix, resolutions=resolutions, samples=samples,
                save=save, assign=assign)
        elif method == "pca":
            if pc is None:
                raise ValueError("If method is 'pca', the value of 'pc' must be given.")
            return self.normalize_coverage_pca(
                matrix=matrix, resolutions=resolutions, samples=samples, pc=pc,
                save=save, assign=assign)
        else:
            raise ValueError("Requested method '{}' is not known.".format(method))

    def plot_all_data(
            self, matrix=None, resolutions=None, samples=None,
            output_dir=None,
            output_prefix="{analysis_name}.all_data",
            robust=True, vmin=None, vmax=None, rasterized=True, dpi=300,
            sample_labels=True):
        """
        Visualize CNV data genome-wide using heatmaps.
        Will be done independently for each specified resolution.

        :param matrix: Attribute name of dictionary of matrices to normalize.
                       Defaults to `coverage_norm`.
        :type matrix: str, optional
        :param resolutions: Resolutions of analysis. Defaults to resolutions
                            in Analysis object.
        :type resolutions: list, optional
        :param samples: Samples to restrict analysis to. Defaults to samples
                            in Analysis object.
        :type samples: list, optional
        :param output_dir: Output directory. Defaults to Analysis results directory.
        :type output_dir: str, optional
        :param output_prefix: Prefix to add to plots. Defaults to "{analysis_name}.all_data"
        :type output_prefix: str, optional
        :param robust: Whether to scale the color scale robustly (to quantiles rather
                       than extremes). Defaults to True
        :type robust: bool, optional
        :param vmin: Minimum value of color scale.
        :type vmin: float, optional
        :param vmax: Maximum value of color scale. Defaults to None
        :type vmax: float, optional
        :param rasterized: Whether to rasterize main heatmap. Defaults to True
        :type rasterized: bool, optional
        :param dpi: DPI resolution of rasterized image. Defaults to 300
        :type dpi: int, optional
        :param sample_labels: Whether to label samples with their name. Defaults to True
        :type sample_labels: bool, optional
        """

        # TODO:
        # add support for group colours

        if matrix is None:
            matrix = self.coverage_norm
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
            self, matrix=None, resolutions=None, samples=None,
            output_dir=None,
            output_prefix="{analysis_name}.all_data",
            robust=True, rasterized=True, dpi=300,
            sample_labels=True):
        """
        Visualize mean and variation of CNV data for each chromosome using heatmaps.
        Will be done independently for each specified resolution.
        Will also be done once for all chromosomes and another time without sex chromosomes.

        :param matrix: Attribute name of dictionary of matrices to normalize.
                       Defaults to `coverage_norm`.
        :type matrix: str, optional
        :param resolutions: Resolutions of analysis. Defaults to resolutions
                            in Analysis object.
        :type resolutions: list, optional
        :param samples: Samples to restrict analysis to. Defaults to samples
                            in Analysis object.
        :type samples: list, optional
        :param output_dir: Output directory. Defaults to Analysis results directory.
        :type output_dir: str, optional
        :param output_prefix: Prefix to add to plots. Defaults to "{analysis_name}.all_data"
        :type output_prefix: str, optional
        :param robust: Whether to scale the color scale robustly (to quantiles
                       rather than extremes). Defaults to True
        :type robust: bool, optional
        :param rasterized: Whether to rasterize main heatmap. Defaults to True
        :type rasterized: bool, optional
        :param dpi: DPI resolution of rasterized image. Defaults to 300
        :type dpi: int, optional
        :param sample_labels: Whether to label samples with their name. Defaults to True
        :type sample_labels: bool, optional
        """

        if matrix is None:
            matrix = self.coverage_norm
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
            self, matrix=None, resolutions=None, samples=None, save=True, assign=True):
        """
        Segment CNV data to create calls of significant deviations.
        Will be done independently for each specified resolution.

        Requires the R package "DNAcopy" to be installed:
            >>> source('http://bioconductor.org/biocLite.R')
            >>> biocLite('DNAcopy')

        :param matrix: Attribute name of dictionary of matrices to segment.
                       Defaults to `coverage_norm`.
        :type matrix: str, optional
        :param resolutions: Resolutions of analysis. Defaults to resolutions in
                            Analysis object.
        :type resolutions: list, optional
        :param samples: Samples to restrict analysis to. Defaults to samples in
                        Analysis object.
        :type samples: list, optional
        :param save: Whether results should be saved to disc. Defaults to True
        :type save: bool, optional
        :param assign: Whether results should be assigned to an attribute in the
                       Analsyis object. Defaults to True
        :type assign: bool, optional

        :returns: Dictionary with segmentation for each resolution.
        :rtype: dict
        :var dict segmentation: Sets a `segmentation` dictionary with CNV matrices
                                for each resolution.
        """

        # TODO: implement as_job

        import rpy2
        from rpy2.robjects import numpy2ri, pandas2ri
        import rpy2.robjects as robjects
        from rpy2.rinterface import RRuntimeError
        numpy2ri.activate()
        pandas2ri.activate()

        robjects.r('require("DNAcopy")')
        _CNA = robjects.r('CNA')
        _segment = robjects.r('segment')
        _segments_p = robjects.r('segments.p')
        _segments_summary = robjects.r('segments.summary')
        # _plotSample = robjects.r('plotSample')

        if matrix is None:
            matrix = self.coverage_norm
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

        :param segmentation: Attribute name of dictionary of segmentation results.
                             Defaults to `segmentation`.
        :type segmentation: str, optional
        :param resolutions: Resolutions of analysis. Defaults to resolutions in
                            Analysis object.
        :type resolutions: list, optional
        :param samples: Samples to restrict analysis to. Defaults to samples in
                        Analysis object.
        :type samples: list, optional
        :param save: Whether results should be saved to disc. Defaults to True
        :type save: bool, optional
        :param assign: Whether results should be assigned to an attribute in the
                       Analsyis object. Defaults to True
        :type assign: bool, optional

        :returns: Dictionary with annotated segmentation for each resolution.
        :rtype: dict
        :var dict segmentation_annot: Sets a `segmentation_annot` dictionary with
                                      CNV matrices for each resolution.
        """
        from ngs_toolkit.general import query_biomart
        import pybedtools
        import string

        if resolutions is None:
            resolutions = self.resolutions

        if segmentation is None:
            segmentation = self.segmentation

        # "http://grch37.ensembl.org/biomart/martview/3ee51bd6fe6cad2220e6b26d21547a1e?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.chromosome_name|hsapiens_gene_ensembl.default.feature_page.start_position|hsapiens_gene_ensembl.default.feature_page.end_position|hsapiens_gene_ensembl.default.feature_page.external_gene_name|hsapiens_gene_ensembl.default.feature_page.band&FILTERS=&VISIBLEPANEL=resultspanel"
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

        :param segmentation: Dictionary of segmentation results.
                             Defaults to `segmentation`.
        :type segmentation: str, optional
        :param resolutions: Resolutions of analysis. Defaults to resolutions in Analysis object.
        :type resolutions: list, optional
        :param per_sample: Whether plots should be made for each sample too. Defaults to False
        :type per_sample: bool, optional
        :param output_dir: Output directory. Defaults to Analysis results directory.
        :type output_dir: str, optional
        :param output_prefix: Prefix to add to plots. Defaults to
                              "{resolution}.segmentation_metrics"
        :type output_prefix: str, optional
        """

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
            fig, axis = plt.subplots(1, 1, figsize=(4 * 1, 4 * 1))

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

        # TODO: plot stats per sample too
        # if per_sample:


def all_to_igv(coverage, output_prefix, **kwargs):
    """
    Convert dictionary of DataFrame with CNV data in several resolutions to IGV format.

    :param coverage: DataFrame with CNV data to convert.
    :type coverage: pandas.DataFrame
    :param output_prefix: Prefix to add to plots.
    :type output_prefix: str, optional
    :param kwargs: Additional parameters will be passed to ngs_toolkit.cnv.to_igv
    :type kwargs: dict, optional
    :returns: Dictionary of CNV data in IGV format for each resolution.
    :rtype: dict
    """
    from ngs_toolkit.cnv import to_igv

    igvs = dict()
    resolutions = coverage.keys()
    for resolution in tqdm(resolutions, total=len(resolutions), desc="Resolution"):
        print("Making IGV visualization for resolution '{}'.".format(resolution))
        igvs[resolution] = to_igv(
            coverage[resolution],
            output_file="{}_{}".format(output_prefix, resolution),
            **kwargs)

    return igvs


def to_igv(matrix, output_file=None, save=True, viewLimits=(-2, 2)):
    """
    Convert DataFrame with CNV data to IGV format.

    :param matrix: DataFrame with CNV data to convert.
    :type matrix: pandas.DataFrame
    :param output_file: Output file. Required is `save` is True.
    :type output_file: str, optional
    :param save: Whether results should be saved to disc. Defaults to True
    :type save: bool, optional
    :param viewLimits: Extreme values (min, max) of color scale used to visualize in IGV.
                       Defaults to (-2, 2).
    :type viewLimits: tuple, optional
    :returns: CNV data in IGV format.
    :rtype: pandas.DataFrame
    :raises: ValueError
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
                            .format(*viewLimits))
        igv.to_csv(output_handle, sep="\t", index=False)
        output_handle.close()

    return igv
