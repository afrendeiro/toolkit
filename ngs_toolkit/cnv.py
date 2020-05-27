#!/usr/bin/env python


import os

import numpy as np
import pandas as pd

from ngs_toolkit import _LOGGER
from ngs_toolkit.analysis import Analysis
from ngs_toolkit.decorators import check_has_attributes
from ngs_toolkit.utils import warn_or_raise

from ngs_toolkit.demo.data_generator import DEFAULT_CNV_RESOLUTIONS


class CNVAnalysis(Analysis):
    """
    Class to model analysis of CNV data.
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

        Defaults to current directory or to what is specified in PEP if :attr:`~ngs_toolkit.analysis.Analysis.from_pep`.
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
        Additional keyword arguments will be passed to parent class :class:`~ngs_toolkit.analysis.Analysis`.

    Examples
    --------
    >>> from ngs_toolkit.cnv import CNVAnalysis

    This is an example of a CNV analysis:

    >>> pep = "metadata/project_config.yaml"
    >>> a = CNVAnalysis(from_pep=pep)
    >>> # Get consensus peak set from all samples
    >>> a.get_cnv_data()
    >>> # Normalize
    >>> a.normalize(method="median")
    >>> # Segmentation
    >>> a.segment_genome()
    >>> # General plots
    >>> a.plot_all_data()
    >>> a.plot_segmentation_stats()
    >>> # Unsupervised analysis
    >>> a.unsupervised_analysis()
    >>> # Save object
    >>> a.to_pickle()
    """

    _data_type = "CNV"

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
            "data_type": "CNV",
            "__data_type__": "CNV",
            "var_unit_name": "bin",
            "quantity": "copy_number",
            "norm_units": "log2(ratio)",
        }
        for k, v in default_args.items():
            if not hasattr(self, k):
                setattr(self, k, v)

        super(CNVAnalysis, self).__init__(
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

        if not hasattr(self, "resolutions"):
            self.resolutions = DEFAULT_CNV_RESOLUTIONS

    def load_data(
        self,
        output_map=None,
        only_these_keys=None,
        resolutions=None,
        prefix="{results_dir}/{name}",
        permissive=True,
    ):
        """
        Load the output files of the major functions of the Analysis.

        Parameters
        ----------
        output_map : :obj:`dict`
            Dictionary with {attribute_name: (file_path, kwargs)} to load the files.
            The kwargs in the tuple will be passed to :meth:`pandas.read_csv`.

            Defaults to the required to read the keys in ``only_these_keys``.
        only_these_keys : :obj:`list`, optional
            Iterable of analysis attributes to load up.
            Possible attributes:

                * "matrix_raw"
                * "matrix_norm"
                * "matrix_features"
                * "differential_results"

            Defaults to all of the above.

        resolutions: :obj:`list`
            List of resolution strings to get data for.

            Defaults to value of ``resolutions`` attribute of Analysis.
        prefix : :obj:`str`, optional
            String prefix of files to load.
            Variables in curly braces will be formated with attributes of analysis.

            Defaults to "{results_dir}/{name}".
        permissive : :obj:`bool`, optional
            Whether an error should be ignored if reading a file causes IOError.

            Default is :obj:`True`.

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

        if resolutions is None:
            resolutions = self.resolutions

        if output_map is None:
            kwargs = {"index_col": 0}
            output_map = {
                "matrix_raw": {
                    r: (prefix + ".{}.matrix_raw.csv".format(r), kwargs) for r in resolutions
                },
                "matrix_norm": {
                    r: (prefix + ".{}.matrix_norm.csv".format(r), kwargs) for r in resolutions
                },
                "segmentation": {
                    r: (prefix + ".{}.segmentation.csv".format(r), {}) for r in resolutions
                },
                "segmentation_annot": {
                    r: (prefix + ".{}.segmentation.annotated.csv".format(r), {})
                    for r in resolutions
                },
            }
        if only_these_keys is None:
            only_these_keys = list(output_map.keys())

        output_map = {k: v for k, v in output_map.items() if k in only_these_keys}

        for name, f in output_map.items():
            for resolution, (file, kwargs) in f.items():
                file = file.format(resolution)
                _LOGGER.info(
                    "Loading '{}' analysis attribute for resolution '{}'.".format(name, resolution)
                )
                if not hasattr(self, name):
                    setattr(self, name, {resolution: None})
                try:
                    getattr(self, name)[resolution] = pd.read_csv(file, **kwargs)
                    # Fix possible multiindex for matrix_norm
                    if name == "matrix_norm":
                        getattr(self, name)[resolution] = fix_dataframe_header(
                            getattr(self, name)[resolution]
                        )
                except IOError as e:
                    if not permissive:
                        raise e
                    else:
                        _LOGGER.warning(e)

    def _copy_cnv_profile_plots(
        self,
        output_dir="{results_dir}/cnv_profiles",
        output_prefix="log2_profile",
        resolutions=None,
        samples=None,
        permissive=True,
    ):
        """
        Convenience to copy output plots from runnning several samples independently
        to a given directory.

        Parameters
        ----------
        output_dir : :obj:`str`, optional
            Directory to copy to.

            Defaults to "{results_dir}/cnv_profiles".
        output_prefix : :obj:`str`, optional
            Prefix for copied files.

            Defaults to "log2_profile".
        resolutions : :obj:`list`, optional
            Resolutions of analysis.

            Defaults to resolutions in Analysis object.
        samples : :obj:`list`, optional
            Samples to restrict analysis to.

            Defaults to samples in Analysis object.
        permissive: :obj:`bool`, optional
            Whether missing files should raise an error.

            Defaults to :obj:`True.`
        """
        from tqdm import tqdm
        from glob import glob
        from shutil import copyfile

        if resolutions is None:
            resolutions = self.resolutions

        if samples is None:
            samples = self.samples

        output_dir = self._format_string_with_attributes(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        for resolution in tqdm(resolutions, total=len(resolutions), desc="Resolution"):

            for sample in tqdm(samples, total=len(samples), desc="Sample"):
                # Read log2 file
                if not hasattr(sample, "log2_read_counts"):
                    sample.log2_read_counts = os.path.join(
                        self.data_dir,
                        sample.sample_root,
                        sample.name + "_{resolution}",
                        "CNAprofiles",
                        "log2_read_counts.igv",
                    )
                if "{resolution}" in sample.log2_read_counts:
                    input_file = sample.log2_read_counts.format(resolution=resolution)
                f = glob(input_file)
                if len(f) == 1:
                    f = f[0]
                else:
                    msg = "Sample '{}' does not have a PDF file!".format(sample.name)
                    if permissive:
                        _LOGGER.warning(msg)
                        continue
                    else:
                        raise OSError(msg)

                d = os.path.join(
                    output_dir, sample.name + "." + resolution + "." + output_prefix + ".pdf"
                )
                try:
                    copyfile(f, d)
                except OSError:
                    msg = "Could not copy file '{}' to '{}'!".format(f, d)
                    if permissive:
                        _LOGGER.warning(msg)
                    else:
                        raise OSError(msg)

    def get_cnv_data(
        self, resolutions=None, samples=None, save=True, assign=True, permissive=False
    ):
        """
        Load CNV data from ATAC-seq CNV pipeline and create CNV matrix at various resolutions.

        Parameters
        ----------
        resolutions : :obj:`list`, optional
            Resolutions of analysis.

            Defaults to resolutions in Analysis object.
        samples : :obj:`list`, optional
            Samples to restrict analysis to.

            Defaults to samples in Analysis object.
        save: :obj:`bool`, optional
            Whether results should be saved to disc.

            Defaults to :obj:`True`
        assign: :obj:`bool`, optional
            Whether results should be assigned to an attribute in the Analsyis object.

            Defaults to :obj:`True`
        permissive: :obj:`bool`, optional
            Whether missing files should be allowed.

            Defaults to :obj:`False`

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
        matrix : :obj:`dict`
            Sets a `matrix` dictionary with CNV matrices for each resolution.
        """
        # TODO: figure out a way of having the input file specified before hand
        from tqdm import tqdm
        from ngs_toolkit.utils import bed_to_index

        if resolutions is None:
            resolutions = self.resolutions

        if samples is None:
            samples = self.samples

        matrix_raw = dict()

        for resolution in tqdm(resolutions, total=len(resolutions), desc="Resolution"):
            matrix_raw[resolution] = pd.DataFrame()

            for sample in tqdm(samples, total=len(samples), desc="Sample"):
                # Read log2 file
                if not hasattr(sample, "log2_read_counts"):
                    msg = "Sample does not have a 'log2_read_counts' attribute."
                    warn_or_raise(AttributeError(msg), permissive)

                input_file = sample.log2_read_counts[resolution].format(resolution=resolution)
                try:
                    cov = pd.read_csv(input_file, sep="\t", comment="#").set_index("Feature")
                except IOError as e:
                    e = IOError(
                        "Sample %s does not have a 'log2_read_counts' file: '%s'."
                        % (sample.name, input_file)
                    )
                    warn_or_raise(e, permissive)

                # TODO: this is specific to CopyWriter, should be removed later
                # and probably replaced with the column position
                cov.columns = (
                    cov.columns.str.replace("log2.", "")
                    .str.replace(".trimmed.bowtie2.filtered.bam", "")
                    .str.replace(".merged.sorted.subsample.bam", "")
                )

                # normalize signal to control
                # # TODO: check whether there was a reason I was previously
                # # undoing and redoing the log
                # matrix_raw[resolution][sample.name] = np.log2(
                #     (
                #         (0.1 + (2 ** cov.loc[:, sample.name]))
                #         / (0.1 + (2 ** cov.iloc[:, -1]))
                #     )
                # )
                matrix_raw[resolution][sample.name] = cov.loc[:, sample.name] - cov.iloc[:, -1]
            if "cov" not in locals():
                msg = "None of the samples had a valid 'log2_read_counts' file."
                _LOGGER.error(msg)
                raise ValueError(msg)

            c = cov.columns.tolist()
            c[:3] = ["chrom", "start", "end"]
            cov.columns = c
            matrix_raw[resolution].index = bed_to_index(cov)

            if save:
                matrix_raw[resolution].to_csv(
                    os.path.join(
                        self.results_dir, self.name + ".{}.matrix_raw.csv".format(resolution),
                    ),
                    index=True,
                )

        if assign:
            self.matrix_raw = matrix_raw

        return matrix_raw

    def normalize(
        self, method="median", matrix="matrix_raw", samples=None, save=True, assign=True, **kwargs
    ):
        """
        Normalization of dictionary of matrices with (n_features, n_samples).

        Parameters
        ----------
        resolutions : :obj:`list`, optional
            Resolutions of analysis.

            Defaults to resolutions in Analysis object.
        method : :obj:`str`
            Normalization method to apply.

            Defaults to "median".
        matrix : :obj:`str`, optional
            Attribute name of dictionary of matrices to normalize.

            Defaults to `matrix_raw`.
        samples : :obj:`list`
            Iterable of peppy.Sample objects to restrict matrix to.

            Defaults to all in analysis.
        save: :obj:`bool`, optional
            Whether results should be saved to disc.

            Defaults to :obj:`True`
        assign: :obj:`bool`, optional
            Whether results should be assigned to an attribute in the Analsyis object.

            Defaults to :obj:`True`
        kwargs : :obj:`dict`, optional
            Additional kwargs are passed to the respective normalization method.

        Returns
        -------
        dict
            Dictionary with normalized CNV matrices for each resolution.

        Attributes
        ----------
        matrix_norm : :obj:`dict`
            Sets a `matrix_norm` dictionary with CNV matrices for each resolution.
        """
        matrix_norm = dict()
        if "resolutions" not in kwargs:
            resolutions = self.resolutions
        else:
            resolutions = kwargs["resolutions"]

        if matrix is None:
            matrix = self.matrix_raw
        elif isinstance(matrix, str):
            matrix = getattr(self, matrix)

        for resolution in resolutions:
            if method == "median":
                matrix_norm[resolution] = self.normalize_median(
                    matrix=matrix[resolution], samples=samples, save=False, assign=False
                )
            elif method == "pca":
                if "pc" not in kwargs:
                    raise ValueError("If method is 'pca', the value of 'pc' must be given.")
                matrix_norm[resolution] = self.normalize_pca(
                    matrix=matrix[resolution],
                    samples=samples,
                    pc=kwargs["pc"],
                    save=False,
                    assign=False,
                )
            else:
                raise ValueError("Requested method '{}' is not known.".format(method))

            if save:
                matrix_norm[resolution].to_csv(
                    os.path.join(
                        self.results_dir, self.name + ".{}.matrix_norm.csv".format(resolution)
                    ),
                    index=True,
                )
        if assign:
            self.matrix_norm = matrix_norm
        return matrix_norm

    def plot_all_data(
        self,
        matrix="matrix_norm",
        resolutions=None,
        samples=None,
        output_dir=None,
        output_prefix="all_data",
        sample_labels=True,
        **kwargs
    ):
        """
        Visualize CNV data genome-wide using heatmaps.
        Will be done independently for each specified resolution.

        Parameters
        ----------
        matrix : :obj:`str`, optional
            Attribute name of dictionary of matrices to normalize.

            Defaults to `matrix_norm`.
        resolutions : :obj:`list`, optional
            Resolutions of analysis.

            Defaults to resolutions in Analysis object.
        samples : :obj:`list`, optional
            Samples to restrict analysis to.

            Defaults to samples in Analysis object.
        output_dir : :obj:`str`, optional
            Output directory.

            Defaults to Analysis results directory.
        output_prefix : :obj:`str`, optional
            Prefix to add to plots.

            Defaults to "{analysis_name}.all_data"
        sample_labels : :obj:`bool`, optional
            Whether to label samples with their name.

            Defaults to :obj:`True`
        **kwargs : :obj:`dict`
            Additional kwargs will be passed to `seaborn.clustermap`.
        """
        # TODO: add support for group colours
        from tqdm import tqdm
        import matplotlib.pyplot as plt
        import seaborn as sns
        from ngs_toolkit.graphics import savefig

        matrix = self.get_matrix(matrix)
        if resolutions is None:
            resolutions = self.resolutions
        if samples is None:
            samples = self.samples
        if output_dir is None:
            output_dir = self.results_dir
        names = [s.name for s in samples]

        # Plot mean and variationper chromosome
        for resolution in tqdm(resolutions, desc="Resolution"):
            r_names = [n for n in names if n in matrix[resolution].columns]
            p = os.path.join(
                self.results_dir, self.name + ".{}.{}.full_data.".format(resolution, output_prefix)
            )
            # Plot all data
            kws = dict(
                cmap="RdBu_r",
                robust=True,
                rasterized=True,
                xticklabels=False,
                yticklabels=r_names if sample_labels else False,
                cbar_kws={"label": "log2(change)"},
            )
            kws.update(kwargs)

            fig, axis = plt.subplots(1, 1, figsize=(4 * 2, 4 * 1))
            sns.heatmap(matrix[resolution].T, center=0, ax=axis, **kws)
            axis.set_xlabel("Chromosome position")
            axis.set_ylabel("Sample")
            savefig(fig, p + "heatmap.svg")

            grid = sns.clustermap(
                matrix[resolution].fillna(0).T,
                metric="correlation",
                center=0,
                col_cluster=False,
                figsize=(4 * 2, 4 * 1),
                **kws
            )
            grid.ax_heatmap.set_xlabel("Chromosome position")
            grid.ax_heatmap.set_ylabel("Sample")
            savefig(grid, p + "fillna.clustermap.svg")

    def plot_stats_per_chromosome(
        self,
        matrix="matrix_norm",
        resolutions=None,
        samples=None,
        output_dir="{results_dir}",
        output_prefix="all_data",
        sample_labels=True,
        **kwargs
    ):
        """
        Visualize mean and variation of CNV data for each chromosome using heatmaps.
        Will be done independently for each specified resolution.
        Will also be done once for all chromosomes and another time without sex chromosomes.

        Parameters
        ----------
        matrix : :obj:`str`, optional
            Attribute name of dictionary of matrices to normalize.

            Defaults to `matrix_norm`.
        resolutions : :obj:`list`, optional
            Resolutions of analysis.

            Defaults to resolutions in Analysis object.
        samples : :obj:`list`, optional
            Samples to restrict analysis to.

            Defaults to samples in Analysis object.
        output_dir : :obj:`str`, optional
            Output directory.

            Defaults to Analysis results directory.
        output_prefix : :obj:`str`, optional
            Prefix to add to plots.

            Defaults to "{analysis_name}.all_data"
        sample_labels: :obj:`bool`, optional
            Whether to label samples with their name.

            Defaults to :obj:`True`
        **kwargs : :obj:`dict`
            Additional kwargs will be passed to `seaborn.clustermap`.
        """
        import seaborn as sns
        from tqdm import tqdm
        from ngs_toolkit.graphics import clustermap_fix_label_orientation, savefig

        matrix = self.get_matrix(matrix)
        if resolutions is None:
            resolutions = self.resolutions
        if samples is None:
            samples = self.samples
        output_dir = self._format_string_with_attributes(output_dir)
        if "{analysis_name}" in output_prefix:
            output_prefix = output_prefix.format(analysis_name=self.name)

        # Plot mean and variation per chromosome
        for resolution in tqdm(resolutions, desc="Resolution"):
            names = [s.name for s in samples if s.name in matrix[resolution].columns]
            to_plot = matrix[resolution].loc[:, names]

            to_plot["chr"] = list(map(lambda x: x[0], to_plot.index.str.split(":")))

            for label, function in tqdm([("variation", np.std), ("mean", np.mean)], desc="metric"):
                prefix = os.path.join(
                    output_dir,
                    self.name + ".{}.".format(resolution) + output_prefix + ".{}".format(label),
                )
                kws = {"yticklabels": False, "xticklabels": sample_labels}
                kws.update(kwargs)

                p = to_plot[names + ["chr"]].groupby("chr").apply(function)
                grid = sns.clustermap(
                    p + p.min().min(), cbar_kws={"label": label}, cmap="Greens", **kws
                )
                clustermap_fix_label_orientation(grid)
                savefig(grid.fig, prefix + "_per_chrom.svg")

                p = (
                    to_plot.loc[~to_plot["chr"].str.contains("X|Y"), names + ["chr"]]
                    .groupby("chr")
                    .apply(function)
                )
                grid = sns.clustermap(
                    p + abs(p.min().min()), cbar_kws={"label": label}, cmap="Greens", **kws
                )
                clustermap_fix_label_orientation(grid)
                savefig(grid.fig, prefix + "_per_chrom.no_sex_chroms.svg")

                grid = sns.clustermap(
                    p,
                    cbar_kws={"label": label + " (Z-score)"},
                    cmap="RdBu_r",
                    center=0,
                    z_score=1,
                    **kws
                )
                clustermap_fix_label_orientation(grid)
                savefig(grid.fig, prefix + "_per_chrom.no_sex_chroms.zscore.svg")

    def segment_genome(
        self, matrix="matrix_norm", resolutions=None, samples=None, save=True, assign=True,
    ):
        """
        Segment CNV data to create calls of significant deviations.
        Will be done independently for each specified resolution.

        Requires the R package "DNAcopy" to be installed:
            >>> source('http://bioconductor.org/biocLite.R')
            >>> biocLite('DNAcopy')

        Parameters
        ----------
        matrix : :obj:`str`, optional
            Attribute name of dictionary of matrices to segment.

            Defaults to `matrix_norm`.
        resolutions : :obj:`list`, optional
            Resolutions of analysis.

            Defaults to resolutions in Analysis object.
        samples : :obj:`list`, optional
            Samples to restrict analysis to.

            Defaults to samples in Analysis object.
        save: :obj:`bool`, optional
            Whether results should be saved to disc.

            Defaults to :obj:`True`
        assign: :obj:`bool`, optional
            Whether results should be assigned to an attribute in the Analsyis object.

            Defaults to :obj:`True`

        Returns
        -------
        dict
            Dictionary with segmentation for each resolution.

        Attributes
        ----------
        segmentation : :obj:`dict`
            Dictionary with CNV matrices for each resolution.
        """
        # TODO: implement distributed mode
        import warnings

        from tqdm import tqdm
        import rpy2
        from rpy2.rinterface import RRuntimeWarning
        from rpy2.robjects import numpy2ri, pandas2ri
        from rpy2.robjects.packages import STAP

        warnings.filterwarnings("ignore", category=RRuntimeWarning)
        numpy2ri.activate()
        pandas2ri.activate()

        DNAcopy = rpy2.robjects.packages.importr("DNAcopy")
        if matrix is None:
            matrix = self.matrix_norm
        elif isinstance(matrix, str):
            matrix = getattr(self, matrix)
        if resolutions is None:
            resolutions = self.resolutions
        if samples is None:
            samples = self.samples

        segmentation = dict()
        for resolution in tqdm(resolutions, desc="Resolution"):
            chrom = np.array(list(map(lambda x: x[0], matrix[resolution].index.str.split(":"))))
            start = np.array(
                list(
                    map(lambda x: int(x[1].split("-")[0]), matrix[resolution].index.str.split(":"),)
                )
            )
            names = [s.name for s in samples if s.name in matrix[resolution].columns]
            df = matrix[resolution].reindex(names, axis=1)

            code = """
            run = function(df, chrom, start, names){
                return(
                    DNAcopy::segment(
                        DNAcopy::CNA(
                            df, chrom=chrom, maploc=start, sampleid=names)))
            }
            """
            run = STAP(code, name="run").run
            res = run(df, chrom, start, names)
            seg = DNAcopy.segments_p(res)
            summary = DNAcopy.segments_summary(res)
            segmentation[resolution] = seg.merge(summary)
            segmentation[resolution].columns = [
                "sample_name",
                "chrom",
                "start",
                "end",
                "bin_size",
                "segment_mean",
                "B_stat",
                "p_value",
                "lcl",
                "ucl",
                "segment_std",
                "segment_median",
                "segment_mad",
            ]
            segmentation[resolution]["segment_length"] = (
                segmentation[resolution]["end"] - segmentation[resolution]["start"]
            ) + int(resolution.replace("kb", "")) * 1000

            if save:
                segmentation[resolution].to_csv(
                    os.path.join(
                        self.results_dir, self.name + ".{}.segmentation.csv".format(resolution),
                    ),
                    index=False,
                )
        if assign:
            self.segmentation = segmentation

        return segmentation

    @check_has_attributes(["organism", "genome"])
    def annotate_with_chrom_bands(
        self, segmentation=None, resolutions=None, save=True, assign=True
    ):
        """
        Annotate segmentation with chromosome bands and overlaping genes.
        Will be done independently for each specified resolution.

        Parameters
        ----------
        segmentation : :obj:`str`, optional
            Attribute name of dictionary of segmentation results.

            Defaults to `segmentation`.
        resolutions : :obj:`list`, optional
            Resolutions of analysis.

            Defaults to resolutions in Analysis object.
        samples : :obj:`list`, optional
            Samples to restrict analysis to.

            Defaults to samples in Analysis object.
        save: :obj:`bool`, optional
            Whether results should be saved to disc.

            Defaults to :obj:`True`
        assign: :obj:`bool`, optional
            Whether results should be assigned to an attribute in the Analsyis object.

            Defaults to :obj:`True`

        Returns
        -------
        dict
            Dictionary with annotated segmentation for each resolution.

        Attributes
        ----------
        segmentation_annot : :obj:`dict`
            Dictionary with CNV matrices for each resolution.
        """
        from ngs_toolkit.general import query_biomart
        import pybedtools

        if resolutions is None:
            resolutions = self.resolutions

        if segmentation is None:
            segmentation = self.segmentation

        organisms = {
            "human": {"species": "hsapiens", "ensembl_version": "grch37"},
            "mouse": {"species": "mmusculus", "ensembl_version": "grcm38"},
            "yeast": {"species": "scerevisiae", "ensembl_version": "R64"},
        }

        annotation = query_biomart(
            attributes=[
                "chromosome_name",
                "start_position",
                "end_position",
                "external_gene_name",
                "band",
            ],
            species=organisms[self.organism]["species"],
            ensembl_version=organisms[self.organism]["ensembl_version"],
        )

        annotation["chromosome_name"] = "chr" + annotation["chromosome_name"]
        annotation = annotation[~annotation["chromosome_name"].str.contains("_")]
        annotation["band"] = annotation["chromosome_name"] + "_" + annotation["band"]
        annotation_bed = pybedtools.BedTool.from_dataframe(annotation.fillna("NA"))

        segmentation_annot = dict()
        for resolution in resolutions:
            seg = segmentation[resolution]
            seg_bed = pybedtools.BedTool.from_dataframe(seg[["chrom", "start", "end"]])

            inter = seg_bed.intersect(annotation_bed, wa=True, wb=True).to_dataframe()
            annot = inter.groupby(["chrom", "start", "end"])["thickStart"].apply(
                lambda x: str.join(" ", set(x))
            )
            annot_band = inter.groupby(["chrom", "start", "end"])["thickEnd"].apply(
                lambda x: str.join(" ", set(x))
            )
            annot_band.name = "chromosome_band"
            annot = annot.to_frame(name="gene").join(annot_band)

            seg_annot = seg.set_index(["chrom", "start", "end"]).join(annot).reset_index()
            segmentation_annot[resolution] = seg_annot.reindex(
                seg.columns.tolist() + ["chromosome_band", "gene"], axis=1
            )

            if save:
                segmentation_annot[resolution].to_csv(
                    os.path.join(
                        self.results_dir,
                        self.name + ".{}.segmentation.annotated.csv".format(resolution),
                    ),
                    index=False,
                )

        if assign:
            self.segmentation_annot = segmentation_annot

    def plot_segmentation_stats(
        self,
        segmentation=None,
        resolutions=None,
        per_sample=False,
        output_dir="{results_dir}/segmentation",
        output_prefix="{resolution}.segmentation_metrics",
    ):
        """
        Visualize distribution of statistics of CNV data segmentation.
        Will be done independently for each specified resolution.

        Parameters
        ----------
        segmentation : :obj:`str`, optional
            Dictionary of segmentation results.

            Defaults to "segmentation".
        resolutions : :obj:`list`, optional
            Resolutions of analysis.

            Defaults to resolutions in Analysis object.
        per_sample: :obj:`bool`, optional
            Whether plots should be made for each sample too.

            Defaults to :obj:`False`
        output_dir : :obj:`str`, optional
            Output directory.

        output_prefix : :obj:`str`, optional
            Prefix to add to plots.

            Defaults to "{resolution}.segmentation_metrics"
        """
        # TODO: plot stats per sample too
        from ngs_toolkit.utils import log_p_value
        import seaborn as sns
        from ngs_toolkit.graphics import savefig

        if segmentation is None:
            segmentation = self.segmentation

        if resolutions is None:
            resolutions = self.resolutions

        output_dir = self._format_string_with_attributes(output_dir)

        metric_vars = [
            "bin_size",
            "segment_mean",
            "B_stat",
            "p_value",
            "lcl",
            "ucl",
            "segment_std",
            "segment_median",
            "segment_mad",
            "segment_length",
        ]

        for resolution in resolutions:
            segments = segmentation[resolution]
            segments["log_p_value"] = log_p_value(segments["p_value"])

            grid = sns.pairplot(
                segments[metric_vars + ["log_p_value"]].dropna(),
                vars=metric_vars + ["log_p_value"],
                plot_kws=dict(s=3, alpha=0.2, linewidth=0, rasterized=True),
                diag_kws=dict(rasterized=True),
            )
            if "{resolution}" in output_prefix:
                op = output_prefix.format(resolution=resolution)
            savefig(grid.fig, os.path.join(output_dir, op + ".all_samples.svg"))

        # if per_sample:


def all_to_igv(matrix, output_prefix, **kwargs):
    """
    Convert dictionary of DataFrame with CNV data in several resolutions to IGV format.

    Parameters
    ----------
    matrix : :obj:`pandas.DataFrame`
        DataFrame with CNV data to convert.

    output_prefix : str
        Prefix to add to plots.

    **kwargs : :obj:`dict`, optional
        Additional parameters will be passed to ngs_toolkit.cnv.to_igv

    Returns
    -------
    dict
        Dictionary of CNV data in IGV format for each resolution.
    """
    from tqdm import tqdm

    igvs = dict()
    resolutions = matrix.keys()
    for resolution in tqdm(resolutions, total=len(resolutions), desc="Resolution"):
        _LOGGER.info("Making IGV visualization for resolution '{}'.".format(resolution))
        igvs[resolution] = to_igv(
            matrix[resolution], output_file="{}_{}.igv".format(output_prefix, resolution), **kwargs
        )

    return igvs


def to_igv(matrix, output_file=None, save=True, view_limits=(-2, 2)):
    """
    Convert DataFrame with CNV data to IGV format.

    Parameters
    ----------
    matrix : :obj:`pandas.DataFrame`
        DataFrame with CNV data to convert.
    output_file : str, optional
        Output file.

        Required if `save` is True.
    save: :obj:`bool`, optional
        Whether results should be saved to disc.

        Defaults to :obj:`True`.
    view_limits : tuple, optional
        Extreme values (min, max) of color scale used to visualize in IGV.

        Defaults to (-2, 2).

    Returns
    -------
    pandas.DataFrame
        CNV data in IGV format.

    Raises
    -------
    ValueError:
        If `save` is True but `output_file` is None.
    """
    _LOGGER.info("Making IGV visualization")

    # as IGV file
    igv = pd.DataFrame(index=matrix.index)
    igv.loc[:, "Chromosome"] = list(map(lambda x: x[0], matrix.index.str.split(":")))
    igv.loc[:, "Start"] = list(map(lambda x: int(x[1].split("-")[0]), matrix.index.str.split(":")))
    igv.loc[:, "End"] = list(map(lambda x: int(x[1].split("-")[1]), matrix.index.str.split(":")))
    igv.loc[:, "Name"] = igv.index

    igv = igv.join(matrix).reset_index(drop=True)

    if save:
        if output_file is None:
            raise ValueError("If the 'save' option is specified, 'output_file' must also be!")
        open(output_file, "w")
        output_handle = open(output_file, "a")
        output_handle.write(
            "#track viewLimits=-{}:{} graphType=heatmap color=255,0,0\n".format(*view_limits)
        )
        igv.to_csv(output_handle, sep="\t", index=False)
        output_handle.close()

    return igv
