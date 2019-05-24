#!/usr/bin/env python


import os

import numpy as np
import pandas as pd

from ngs_toolkit import _LOGGER
from ngs_toolkit.atacseq import ATACSeqAnalysis


class ChIPSeqAnalysis(ATACSeqAnalysis):
    """
    Class to model analysis of ChIP-seq data.
    Inherits from the `ngs_toolkit.atacseq.ATACSeqAnalysis` class.

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

        Defaults to :obj:`None` (will not load).
    root_dir : :obj:`str`, optional
        Base directory for the project.

        Defaults to current directory or to what is specified in PEP if `from_pep`.
    data_dir : :obj:`str`, optional
        Directory containing processed data (e.g. by looper) that will
        be input to the analysis. This is in principle not required.

        Defaults to "data".
    results_dir : :obj:`str`, optional
        Directory to contain outputs produced by the analysis.

        Defaults to "results".
    prj : :obj:`peppy.Project`, optional
        A ``peppy.Project`` object that this analysis is tied to.

        Defaults to :obj:`None`.
    samples : :obj:`list`, optional
        List of ``peppy.Sample`` objects that this analysis is tied to.

        Defaults to :obj:`None`.
    kwargs : :obj:`dict`, optional
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
            "data_type": "ChIP-seq",
            "__data_type__": "ChIP-seq",
            "var_unit_name": "region",
            "quantity": "binding",
            "norm_units": "RPM"}
        for k, v in default_args.items():
            if not hasattr(self, k):
                setattr(self, k, v)

        super(ChIPSeqAnalysis, self).__init__(
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

    def call_peaks_from_comparisons(
        self,
        comparison_table=None,
        output_dir="{results_dir}/chipseq_peaks",
        permissive=True,
        overwrite=True,
        distributed=True,
    ):
        """
        Call peaks for ChIP-seq samples using an annotation of which samples
        belong in each comparison and which samples represent signal or background.

        Parameters
        ----------
        comparison_table : :obj:`pandas.DataFrame`
            Comparison table with the following required columns:
            "comparison_name", "sample_name", "comparison_side", "sample_group".

            Defaults to analysis' own `comparison_table`.
        output_dir : :obj:`str`
            Parent directory where peaks will be created.

            Will be created if does not exist.
        permissive: :obj:`bool`
            If incomplete/incoherent comparisons should be skipped or an error should be thrown.

            Default is :obj:`True`.
        overwrite: :obj:`bool`
            If incomplete/incoherent comparisons should be skipped or an error should be thrown.

            Default is :obj:`True`.
        distributed: :obj:`bool`
            Whether peak calling should be run in serial or in distributed mode as jobs.

            Default is :obj:`True`.

        Raises
        ----------
        ValueError
            If not `permissive` and incomplete/incoherent comparisons are detected.
        """
        import subprocess

        from ngs_toolkit.utils import (
            macs2_call_chipseq_peak,
            homer_call_chipseq_peak_job,
        )
        from tqdm import tqdm

        if comparison_table is None:
            comparison_table = self.comparison_table
        req_columns = [
            "comparison_name",
            "sample_name",
            "comparison_side",
            "sample_group",
        ]
        msg = "Comparison table is missing some of the following columns: '{}'.".format(
            ",".join(req_columns)
        )
        if not all([col in comparison_table.columns for col in req_columns]):
            _LOGGER.error(msg)
            raise AssertionError(msg)

        # Complement default `output_dir`
        if "{results_dir}" in output_dir:
            output_dir = os.path.abspath(
                output_dir.format(results_dir=self.results_dir)
            )
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # For each comparison
        comps = comparison_table["comparison_name"].drop_duplicates().sort_values()
        for comparison in tqdm(comps, total=len(comps), desc="Comparison"):
            # If there aren't two sides to each comparison, skip it or throw error
            if (
                len(
                    set(
                        comparison_table[
                            (comparison_table["comparison_name"] == comparison)
                        ]["comparison_side"].tolist()
                    )
                )
                != 2
            ):
                error = "Comparison '{}' does not contain two sides.".format(comparison)
                if permissive:
                    _LOGGER.warning(error)
                    continue
                else:
                    _LOGGER.error(error)
                    raise ValueError(error)

            # Get the sample names of samples in each side
            pos_names = comparison_table[
                (comparison_table["comparison_name"] == comparison)
                & (comparison_table["comparison_side"] == 1)
            ]["sample_name"].tolist()
            neg_names = comparison_table[
                (comparison_table["comparison_name"] == comparison)
                & (comparison_table["comparison_side"] < 1)
            ]["sample_name"].tolist()

            # Now get the actual samples
            signal_samples = [s for s in self.samples if s.name in pos_names]
            control_samples = [s for s in self.samples if s.name in neg_names]

            if (len(signal_samples) == 0) or (len(control_samples) == 0):
                error = "Comparison side for '{}' comparison does not contain samples.".format(
                    comparison
                )
                if permissive:
                    print(error)
                    continue
                else:
                    raise ValueError(error)

            print(
                "Doing comparison '{}' with positive samples '{}' and background samples '{}'".format(
                    comparison,
                    [s.name for s in signal_samples],
                    [s.name for s in control_samples],
                )
            )
            # Call peaks
            cmds = list()
            kwargs = {
                "signal_samples": signal_samples, "control_samples": control_samples,
                "output_dir": output_dir, "name": comparison, "distributed": distributed}
            if overwrite:
                cmds += [
                    macs2_call_chipseq_peak(**kwargs), homer_call_chipseq_peak_job(**kwargs)
                ]
            else:
                if not os.path.exists(
                    os.path.join(
                        output_dir, comparison, comparison + "_peaks.narrowPeak"
                    )
                ):
                    cmds += [
                        macs2_call_chipseq_peak(**kwargs)
                    ]
                if not os.path.exists(
                    os.path.join(
                        output_dir, comparison, comparison + "_homer_peaks.narrowPeak"
                    )
                ):
                    cmds += [
                        homer_call_chipseq_peak_job(**kwargs)
                    ]
                else:
                    _LOGGER.warning(
                        "Peak files for comparison '{}' already exist. Skipping.".format(
                            comparison
                        )
                    )

            if not distributed:
                for cmd in cmds:
                    _LOGGER.info(
                        "Calling peaks for comparison '{}' with command: '{}'.\n".format(
                            comparison, cmd
                        )
                    )
                    subprocess.call(cmd.split(" "))

    def filter_peaks(
        self,
        comparison_table=None,
        filter_bed=None,
        peaks_dir="{results_dir}/chipseq_peaks",
    ):
        """
        Filter peak calls for various comparisons for entries that do not overlap another BED file.

        Parameters
        ----------
        comparison_table : :obj:`pandas.DataFrame`, optional
            Comparison table with the following required columns:
            "comparison_name", "sample_name", "comparison_side", "sample_group".

            Defaults to analysis' own `comparison_table`.
        filter_bed : :obj:`str`
            BED file with entries to filter out from the BED files of each comparison.

            Defaults to the set of Blacklisted regions from the analysis' genome.
            In that case it will be fetched if not present.
        peaks_dir : :obj:`str`
            Parent directory where peak calls for each comparison exist.
            Will be created if does not exist.

            Defaults to "{results_dir}/chipseq_peaks".

        Raises
        ----------
        AttributeError
            If `filter_bed` is not given and failes to be retrieved.
        """
        import pybedtools
        from ngs_toolkit.utils import homer_peaks_to_bed

        def _filter(input_bed, filter_bed, output_bed):
            """
            Filter BED file for entries that overlap another BED file.

            Parameters
            ----------
            input_bed : :obj:`str`
                BED file to filter.

            filter_bed : :obj:`str`
                BED file with entries to filter from input_bed.

            output_bed : :obj:`str`
                Output BED file.
            """
            (
                pybedtools.BedTool(input_bed)
                .intersect(pybedtools.BedTool(filter_bed), v=True)
                .saveas(output_bed)
            )

        if comparison_table is None:
            comparison_table = self.comparison_table

        if filter_bed is None:
            _LOGGER.info("Blacklist file not provided. Downloading...")
            try:
                filter_bed = self.get_resources(steps=["blacklist"])[
                    "blacklist_file"
                ]
            except AttributeError:
                msg = "Blacklist file was not provided and cannot be"
                msg += " get one without analysis having `organism` and `genome` set."
                _LOGGER.error(msg)
                raise AttributeError(msg)

        peaks_dir = self._format_string_with_attributes(peaks_dir)
        if not os.path.exists(peaks_dir):
            os.makedirs(peaks_dir)

        for comparison in comparison_table["comparison_name"].unique():
            # MACS2
            prefix = os.path.join(peaks_dir, comparison, comparison + "_peaks")
            _filter(prefix + ".narrowPeak", filter_bed, prefix + ".filtered.bed")
            # HOMER
            prefix = os.path.join(peaks_dir, comparison, comparison + "_homer_peaks")
            homer_peaks_to_bed(prefix + ".factor.narrowPeak", prefix + ".factor.bed")
            _filter(prefix + ".factor.bed", filter_bed, prefix + ".factor.filtered.bed")
            homer_peaks_to_bed(prefix + ".histone.narrowPeak", prefix + ".histone.bed")
            _filter(
                prefix + ".histone.bed", filter_bed, prefix + ".histone.filtered.bed"
            )

    def summarize_peaks_from_comparisons(
        self,
        comparison_table=None,
        output_dir="{results_dir}/chipseq_peaks",
        filtered=True,
        permissive=True,
    ):
        """
        Call peaks for ChIP-seq samples using an annotation of which samples belong in each comparison and which
        samples represent signal or background.

        Parameters
        ----------
        comparison_table : :obj:`pandas.DataFrame`, optional
            Comparison table with the following required columns:
            "comparison_name", "sample_name", "comparison_side", "sample_group".

            Defaults to analysis' own `comparison_table`.
        output_dir : :obj:`str`
            Parent directory where peaks will be created. Will be created if does not exist.
        permissive: :obj:`bool`
            If incomplete/incoherent comparisons should be skipped or an error should be thrown.

        Raises
        ----------
        ValueError
            Will be raised if not `permissive` and incomplete/incoherent comparisons are detected.
        """
        from ngs_toolkit.utils import homer_peaks_to_bed

        if comparison_table is None:
            comparison_table = self.comparison_table

        req_columns = [
            "comparison_name",
            "sample_name",
            "comparison_side",
            "sample_group",
        ]
        msg = "Comparison table is missing some of the following columns: '{}'.".format(
            ",".join(req_columns)
        )
        if not all([col in comparison_table.columns for col in req_columns]):
            _LOGGER.error(msg)
            raise AssertionError(msg)

        # Complement default `output_dir`
        if "{results_dir}" in output_dir:
            output_dir = os.path.abspath(
                output_dir.format(results_dir=self.results_dir)
            )

        # For each comparison, count called peaks
        peak_counts = pd.DataFrame()
        for comparison in (
            comparison_table["comparison_name"].drop_duplicates().sort_values()
        ):
            _LOGGER.info(comparison)
            ending = ".filtered.bed" if filtered else ".narrowPeak"
            for peak_type, file in [
                (
                    "macs",
                    os.path.join(
                        output_dir, comparison, comparison + "_peaks" + ending
                    ),
                ),
                (
                    "homer_factor",
                    os.path.join(
                        output_dir,
                        comparison,
                        comparison + "_homer_peaks.factor" + ending,
                    ),
                ),
                (
                    "homer_histone",
                    os.path.join(
                        output_dir,
                        comparison,
                        comparison + "_homer_peaks.histone" + ending,
                    ),
                ),
            ]:
                error = "Peak files for comparison '{}' with '{}' parameters don't exist.".format(
                    comparison, peak_type
                )

                if "homer" in peak_type and not filtered:
                    try:
                        homer_peaks_to_bed(file, file.replace("narrowPeak", "bed"))
                    except IOError:
                        if permissive:
                            _LOGGER.warning(error)
                            peak_counts = peak_counts.append(
                                pd.Series([comparison, peak_type, np.nan]),
                                ignore_index=True,
                            )
                            continue
                        else:
                            raise
                    except pd.errors.EmptyDataError:
                        peak_counts = peak_counts.append(
                            pd.Series([comparison, peak_type, 0.0]), ignore_index=True
                        )

                    file = file.replace("narrowPeak", "bed")
                try:
                    df = pd.read_csv(file, sep="\t")
                except IOError:
                    if permissive:
                        _LOGGER.warning(error)
                        peak_counts = peak_counts.append(
                            pd.Series([comparison, peak_type, np.nan]),
                            ignore_index=True,
                        )
                        continue
                    else:
                        raise
                except pd.errors.EmptyDataError:
                    peak_counts = peak_counts.append(
                        pd.Series([comparison, peak_type, 0.0]), ignore_index=True
                    )

                peak_counts = peak_counts.append(
                    pd.Series([comparison, peak_type, df.shape[0]]), ignore_index=True
                )
        peak_counts.columns = ["comparison_name", "peak_type", "peak_counts"]

        return peak_counts  # .fillna(0)

    def get_consensus_sites(self, **kwargs):
        """
        Get consensus (union) of enriched sites (peaks) across all comparisons.
        If `region_type` is "summits, regions used will be peak summits which will be extended by `extension`
        before union. Otherwise sample peaks will be used with no modification.

        Parameters
        ----------
        comparison_table : :obj:`pandas.DataFrame`, optional
            DataFrame with signal/background combinations used to call peaks

            Defaults to analysis' own `comparison_table`.
        peak_dir : :obj:`str`, optional
            Path to peaks output directory.

            Defaults to `{analysis.results_dir}/chipseq_peaks`.
        region_type : :obj:`str`, optional
            Type of region to use.
            One of "summits" or "peaks".

            Default is "peaks".
        extension : :obj:`int`, optional
            Width to extend peak summits is `region_type` is "summits".

            Default is 250.
        blacklist_bed : :obj:`str`, optional
            A BED file with genomic positions to exclude from consensus peak set.

            Will try to download if not given.

        Attributes
        ----------
        sites : :class:`pybedtools.BedTool`
            Bedtool with consensus sites.
        """
        import re
        from ngs_toolkit.general import get_blacklist_annotations
        import pybedtools
        from tqdm import tqdm

        if "comparison_table" not in kwargs:
            comparison_table = self.comparison_table
        else:
            comparison_table = kwargs["comparison_table"]

        if "peak_dir" not in kwargs:
            peak_dir = os.path.join(self.results_dir, "chipseq_peaks")

        if "region_type" not in kwargs:
            region_type = "peaks"
        if region_type not in ["summits", "peaks"]:
            msg = "`region_type` attribute must be one of 'summits' or 'peaks'!"
            _LOGGER.error(msg)
            raise ValueError(msg)

        if "extension" not in kwargs:
            extension = 250

        if "permissive" not in kwargs:
            permissive = False

        if "blacklist_bed" not in kwargs:
            _LOGGER.info("Blacklist file not provided. Downloading...")
            try:
                blacklist_bed = get_blacklist_annotations(self.organism, self.genome)
            except AttributeError:
                msg = "Blacklist file was not provided and cannot"
                msg += " get one without analysis having `organism` and `genome` set."
                _LOGGER.error(msg)
                raise AttributeError(msg)

        first = True
        comps = comparison_table["comparison_name"].drop_duplicates()
        for comparison in tqdm(comps, total=len(comps), desc="Comparison"):
            peak_files = [
                os.path.join(peak_dir, comparison, comparison + "_peaks.narrowPeak"),
                os.path.join(
                    peak_dir, comparison, comparison + "_homer_peaks.factor.bed"
                ),
                os.path.join(
                    peak_dir, comparison, comparison + "_homer_peaks.histone.bed"
                ),
            ]
            for peak_file in peak_files:
                genome = (
                    comparison_table.loc[
                        comparison_table["comparison_name"] == comparison,
                        "comparison_genome",
                    ]
                    .drop_duplicates()
                    .squeeze()
                )

                msg = "Could not determine genome of comparison '{}'.".format(
                    comparison
                )
                if not isinstance(genome, str):
                    _LOGGER.error(msg)
                    raise AssertionError(msg)

                # Get peaks
                if region_type == "summits":
                    try:
                        f = re.sub("_peaks.narrowPeak", "_summits.bed", peak_file)
                        peaks = pybedtools.BedTool(f).slop(b=extension, genome=genome)
                    except ValueError:
                        _LOGGER.warning(
                            "Summits for comparison {} ({}) not found!".format(
                                comparison, f
                            )
                        )
                        if not permissive:
                            raise
                else:
                    try:
                        peaks = pybedtools.BedTool(peak_file)
                    except ValueError:
                        _LOGGER.warning(
                            "Peaks for comparison {} ({}) not found!".format(
                                comparison, peak_file
                            )
                        )
                        if not permissive:
                            raise
                # Merge overlaping peaks within a comparison
                peaks = peaks.merge()
                if first:
                    sites = peaks
                    first = False
                else:
                    # Concatenate all peaks
                    sites = sites.cat(peaks)

        # Merge overlaping peaks across comparisons
        sites = sites.merge()

        # Filter
        # remove blacklist regions
        blacklist = pybedtools.BedTool(
            os.path.join(self.data_dir, "external", blacklist_bed)
        )
        # remove chrM peaks and save
        sites.intersect(v=True, b=blacklist).filter(lambda x: x.chrom != "chrM").saveas(
            os.path.join(self.results_dir, self.name + ".peak_set.bed")
        )

        # Read up again
        self.sites = pybedtools.BedTool(
            os.path.join(self.results_dir, self.name + ".peak_set.bed")
        )

    def calculate_peak_support(self, **kwargs):
        """
        Calculate a measure of support for each region in peak set
        (i.e. ratio of samples containing a peak overlapping region in union set of peaks).

        Parameters
        ----------
        comparison_table : :obj:`pandas.DataFrame`, optional
            DataFrame with signal/background combinations used to call peaks

            Defaults to analysis' own `comparison_table`.
        peak_dir : :obj:`str`, optional
            Path to peaks output directory.
            Defaults to {analysis.results_dir}/chipseq_peaks

        Attributes
        ----------
        support : :obj:`pandas.DataFrame`
            DataFrame with signal/background combinations used to call peaks
        """
        import pybedtools
        from tqdm import tqdm

        if "comparison_table" not in kwargs:
            comparison_table = self.comparison_table
        else:
            comparison_table = kwargs["comparison_table"]
        if "peak_dir" not in kwargs:
            peak_dir = "{results_dir}/chipseq_peaks"

        if "{results_dir}" in peak_dir:
            peak_dir = os.path.abspath(peak_dir.format(results_dir=self.results_dir))

        # get index
        index = self.sites.to_dataframe()
        index = (
            index["chrom"]
            + ":"
            + index["start"].astype(str)
            + "-"
            + index["end"].astype(str)
        )

        # calculate support (number of samples overlaping each merged peak)
        comps = comparison_table["comparison_name"].drop_duplicates()
        support = pd.DataFrame(index=index)
        for comparison in tqdm(comps, total=comps.shape[0], desc="Comparison"):
            peak_files = [
                (
                    "MACS",
                    os.path.join(
                        peak_dir, comparison, comparison + "_peaks.narrowPeak"
                    ),
                ),
                (
                    "HOMER_factor",
                    os.path.join(
                        peak_dir, comparison, comparison + "_homer_peaks.factor.bed"
                    ),
                ),
                (
                    "HOMER_histone",
                    os.path.join(
                        peak_dir, comparison, comparison + "_homer_peaks.histone.bed"
                    ),
                ),
            ]
            for peak_type, peak_file in peak_files:
                try:
                    sample_support = self.sites.intersect(
                        peak_file, wa=True, c=True
                    ).to_dataframe()
                except (
                    ValueError,
                    pybedtools.MalformedBedLineError,
                    pybedtools.helpers.BEDToolsError,
                ):
                    _LOGGER.warning(
                        "Peaks for comparison {} ({}) not found!".format(
                            comparison, peak_file
                        )
                    )
                    continue
                sample_support.index = index
                support[(comparison, peak_type)] = sample_support.iloc[:, 3]

        # Make multiindex labeling comparisons and peak type
        support.columns = pd.MultiIndex.from_tuples(
            support.columns, names=["comparison", "peak_type"]
        )
        support.to_csv(
            os.path.join(
                self.results_dir, self.name + "_peaks.binary_overlap_support.csv"
            ),
            index=True,
        )

        # divide sum (of unique overlaps) by total to get support value between 0 and 1
        support["support"] = support.astype(bool).astype(int).sum(axis=1) / float(
            support.shape[1]
        )

        # save
        support.to_csv(
            os.path.join(self.results_dir, self.name + "_peaks.support.csv"), index=True
        )

        self.support = support

    def get_supported_peaks(self, **kwargs):
        """
        Get mask of sites with 0 support in the given samples.
        Requires support matrix produced by `ngs_toolkit.atacseq.ATACSeqAnalysis.calculate_peak_support`.

        Parameters
        ----------
        comparisons : :obj:`list`
            Iterable of comparison names to restrict to.
            Must match name of comparisons used in comparison_table.

        Returns
        -------
        pd.Series
            Boolean Pandas Series with sites with at least one of the given samples having a peak called.
        """
        if "comparisons" not in kwargs:
            msg = "Requires a keyword argument `comparisons`."
            _LOGGER.error(msg)
            raise ValueError(msg)
        return self.support[[c for c in kwargs["comparisons"]]].sum(1) != 0

    def normalize_by_background(
            self,
            comparison_table=None,
            reduction_func=np.mean,
            comparison_func=np.subtract,
            by_group=False,
            matrix="matrix_norm",
            samples=None):
        """
        Normalize values in matrix by background samples in a comparison-specific way as specified in `comparison_table`.

        The background samples will be pooled by the `reduction_func` and their values wil be removed
        from the signal samples using the `comparison_func`.

        Parameters
        ----------
        comparison_table : :class:`pandas.DataFrame`
            Table with comparisons from which peaks were called.

            Defaults to analysis' `comparison_table`.
        reduction_func : func
            Function to reduce the region to gene values to.

            Defaults to :obj:`numpy.mean`.
        comparison_func : func
            Function to use for normalization of signal samples against background samples.
            You can also try for example :obj:`numpy.divide`.

            Defaults to :obj:`numpy.subtract`.
        by_group : :obj:`bool`
            Whether output should be by group (:obj:`True`) or for each sample (:obj:`False`).

            Default is :obj:`False`.
        matrix : {:class:`pandas.DataFrame`, :obj:`str`, optional}
            Name of attribute or pandas DataFrame to use.

            Defaults to "matrix_norm".
        samples : :obj:`list`, optional
            Subset of samples to consider.

            Defaults to all samples in analysis.

        Returns
        -------
        :class:`pandas.DataFrame`
            Dataframe with values normalized by background samples.
        """
        if comparison_table is None:
            comparison_table = self.comparison_table
        if samples is None:
            samples = self.samples
        matrix = self.get_matrix(matrix, samples=samples)

        comparisons = set(comparison_table['comparison_name'])
        s_names = [s.name for s in samples]
        if not by_group:
            res = pd.DataFrame(index=matrix.index)
            for comparison in comparisons:
                comp = comparison_table.query("comparison_name == '{}'".format(comparison))
                signal_samples = [s for s in comp.query("comparison_side == 1")['sample_name'] if s in s_names]
                background_samples = [s for s in comp.query("comparison_side == 0")['sample_name'] if s in s_names]

                for s in signal_samples:
                    res.loc[:, s] = comparison_func(
                        matrix.loc[:, s],
                        matrix.loc[:, background_samples].apply(reduction_func, axis=1)
                        if reduction_func != np.mean
                        else matrix.loc[:, background_samples].mean(axis=1))
        if by_group:
            comparisons = set(comparison_table['comparison_name'])
            s_names = [s.name for s in samples]
            res = pd.DataFrame(index=matrix.index, columns=comparisons)
            for comparison in comparisons:
                comp = comparison_table.query("comparison_name == '{}'".format(comparison))
                signal_samples = [s for s in comp.query("comparison_side == 1")['sample_name'] if s in s_names]
                background_samples = [s for s in comp.query("comparison_side == 0")['sample_name'] if s in s_names]

                res.loc[:, comparison] = comparison_func(
                    matrix.loc[:, signal_samples].apply(reduction_func, axis=1)
                    if reduction_func != np.mean
                    else matrix.loc[:, signal_samples].mean(axis=1),
                    matrix.loc[:, background_samples].apply(reduction_func, axis=1)
                    if reduction_func != np.mean
                    else matrix.loc[:, background_samples].mean(axis=1))
        return res
