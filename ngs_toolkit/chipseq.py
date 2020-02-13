#!/usr/bin/env python


import os

import numpy as np
import pandas as pd

from ngs_toolkit import _LOGGER
from ngs_toolkit.atacseq import ATACSeqAnalysis
from ngs_toolkit.decorators import check_has_attributes


class ChIPSeqAnalysis(ATACSeqAnalysis):
    """
    Class to model analysis of ChIP-seq data.
    Inherits from the :class:`~ngs_toolkit.atacseq.ATACSeqAnalysis` class.

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
        Additional keyword arguments will be passed to parent class :class:`~ngs_toolkit.atacseq.ATACSeqAnalysis`.
    """
    _data_type = "ChIP-seq"

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

        if hasattr(self, "comparison_table"):
            self.set_comparisons()
        else:
            msg = "No comparison table was given. Will not prefill peak calling comparisons."
            _LOGGER.warning(msg)

    def set_comparisons(self, comparison_table=None, peak_dir="{results_dir}/chipseq_peaks"):
        """
        Set up an attribute containing information about the
        sample comparisons necessary for peak calling.

        Structure:

            * comparison_name:
                * signal_samples
                * background_samples
                * directory
                * prefix
                * resulting_files
                    * macs
                    * homer_histone
                    * homer_factor

        Parameters
        ----------
        comparison_table : :obj:`str`, optional
            Comparison table wit peak comparisons.

            Defaults to one from PEP project if available.
        peak_dir : :obj:`str`, optional
            Directory with peak calls.

            Defaults to "{results_dir}/chipseq_peaks".

        Returns
        -------
        :obj:`dict`
            The dictionary with the attributes.

        Attributes
        ----------
        :obj:`dict`
            The dictionary with the attributes.

        Raises
        ------
        ValueError
            If comparisons are not correctly specified.
        """
        if comparison_table is None:
            comparison_table = self.comparison_table

        comparison_names = (
            comparison_table.loc[
                comparison_table['comparison_type'] == 'peaks',
                "comparison_name"]
            .drop_duplicates().sort_values()).tolist()
        if not comparison_names:
            _LOGGER.warning("Could not find any comparisons of type 'peak'.")

        peak_dir = os.path.abspath(self._format_string_with_attributes(peak_dir))

        self.comparisons = dict()
        for name in comparison_names:
            _LOGGER.info("Setting comparison '%s' up", name)

            # If there aren't two sides to each comparison, skip it or throw error
            if len(set(comparison_table.query("comparison_name == '{}'".format(name))["comparison_side"])) != 2:
                error = "Comparison '{}' does not contain two sides.".format(name)
                _LOGGER.error(error)
                raise ValueError(error)

            # Get the sample names of samples in each side
            pos_names = comparison_table.loc[
                (comparison_table["comparison_name"] == name) & (comparison_table["comparison_side"] == 1),
                "sample_name"].tolist()
            neg_names = comparison_table.loc[
                (comparison_table["comparison_name"] == name) & (comparison_table["comparison_side"] < 1),
                "sample_name"].tolist()

            signal_samples = [s for s in self.samples if s.name in pos_names]
            control_samples = [s for s in self.samples if s.name in neg_names]

            co = dict()
            co['signal_samples'] = signal_samples
            co['control_samples'] = control_samples

            # Additional info
            co['output_dir'] = os.path.join(peak_dir, name)
            co['prefix'] = os.path.join(co['output_dir'], name)
            g = comparison_table.query("comparison_name == '{}'".format(name))['comparison_genome'].drop_duplicates().squeeze()
            if not isinstance(g, str):
                msg = "Could not determine genome of comparison '%s'." % g
                _LOGGER.error(msg)
                raise AssertionError(msg)
            co['genome'] = g

            # resulting files files
            res = dict()
            res['macs'] = co['prefix'] + "_peaks.narrowPeak"
            res["homer_factor"] = co['prefix'] + "_homer_peaks.factor.narrowPeak"
            res["homer_histone"] = co['prefix'] + "_homer_peaks.histone.narrowPeak"
            co['peak_calls'] = dict()
            co['peak_calls']["original"] = res
            co['peak_calls']["filtered"] = {
                k: v.replace(".narrowPeak", ".filtered.bed")
                for k, v in res.items()}

            self.comparisons[name] = co
        return self.comparisons

    @check_has_attributes(["comparisons"], [dict])
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
            filter_kwargs_by_callable
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
        for name, comp in tqdm(self.comparisons.items(), total=len(self.comparisons), desc="Comparison"):

            _LOGGER.info(
                "Doing comparison '{}' with positive samples '{}' and background samples '{}'".format(
                    name,
                    [s.name for s in comp['signal_samples']],
                    [s.name for s in comp['control_samples']],
                )
            )
            # Call peaks
            cmds = list()
            bkws = filter_kwargs_by_callable(comp, macs2_call_chipseq_peak)
            kwargs = {
                "name": name, "distributed": distributed, **bkws}
            if overwrite:
                cmds += [macs2_call_chipseq_peak(**kwargs), homer_call_chipseq_peak_job(**kwargs)]
            else:
                if not os.path.exists(comp['peak_calls']['original']['macs']):
                    cmds += [macs2_call_chipseq_peak(**kwargs)]
                if not os.path.exists(comp['peak_calls']['original']['homer_factor']):
                    cmds += [homer_call_chipseq_peak_job(**kwargs)]
                else:
                    _LOGGER.warning("Peak files for comparison '%s' already exist. Skipping.", name)
            if not distributed:
                for cmd in cmds:
                    _LOGGER.info("Calling peaks for comparison '%s' with command: '%s'.\n", (name, cmd))
                    subprocess.call(cmd.split(" "))

    @check_has_attributes(["comparisons"], [dict])
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
        from ngs_toolkit.utils import homer_peaks_to_bed, filter_bed_file

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

        for name, comp in self.comparisons.items():
            # MACS2
            filter_bed_file(comp['peak_calls']['original']['macs'], filter_bed, comp['peak_calls']['filtered']['macs'])
            # HOMER
            tmp_bed = comp['peak_calls']['original']['homer_factor'].replace(".narrowPeak", ".bed")
            homer_peaks_to_bed(comp['peak_calls']['original']['homer_factor'], tmp_bed)
            filter_bed_file(tmp_bed, filter_bed, comp['peak_calls']['filtered']['homer_factor'])
            os.remove(tmp_bed)

            tmp_bed = comp['peak_calls']['original']['homer_histone'].replace(".narrowPeak", ".bed")
            homer_peaks_to_bed(comp['peak_calls']['original']['homer_histone'], tmp_bed)
            filter_bed_file(tmp_bed, filter_bed, comp['peak_calls']['filtered']['homer_histone'])
            os.remove(tmp_bed)

    @check_has_attributes(["comparisons"], [dict])
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
        peak_type = "filtered" if filtered else "original"
        peak_counts = list()
        for name, comp in self.comparisons.items():
            _LOGGER.info(name)

            for peak_caller, file in comp['peak_calls'][peak_type].items():
                error = "Peak files for comparison '%s' with '%s' parameters don't exist."

                if "homer" in peak_caller and not filtered:
                    try:
                        homer_peaks_to_bed(file, file.replace("narrowPeak", "bed"))
                    except IOError:
                        if permissive:
                            _LOGGER.warning(error, (name, peak_caller))
                            peak_counts.append([name, peak_caller, np.nan])
                            continue
                        else:
                            raise
                    except pd.errors.EmptyDataError:
                        peak_counts.append([name, peak_caller, 0.0])
                    file = file.replace("narrowPeak", "bed")
                try:
                    df = pd.read_csv(file, sep="\t")
                except IOError:
                    if permissive:
                        _LOGGER.warning(error, (name, peak_caller))
                        peak_counts.append([name, peak_caller, np.nan])
                        continue
                    else:
                        raise
                except pd.errors.EmptyDataError:
                    peak_counts.append([name, peak_caller, 0.0])
                peak_counts.append([name, peak_caller, df.shape[0]])
        peak_counts = pd.DataFrame(peak_counts, columns=["comparison_name", "peak_caller", "peak_counts"])

        return peak_counts  # .fillna(0)

    def get_consensus_sites(
            self,
            samples=None,
            region_type="summits",
            peak_type="filtered",
            extension=250,
            blacklist_bed=None,
            filter_chroms=True,
            permissive=False,
            save=True,
            assign=True,
            **kwargs):
        """
        Get consensus (union) of enriched sites (peaks) across all comparisons.
        There are two modes possible, defined by the value of ``region_type``:

         * peaks: simple union of all sites;
         * summits: peak summits are extended by ``extension`` and a union is made.

        For ChIP-seq, the ``comparison_table`` keyword argument or a
        ``comparison_table`` attribute set is required. Peaks/summits will be
        aggregated for the peaks called in each sample comparison.

        Parameters
        ----------
        samples : :obj:`list`
            Iterable of :class:`peppy.Sample` objects to restrict to.
            Must have a ``peaks`` attribute set.

            Defaults to all samples in the analysis (``samples`` attribute).
        region_type : :obj:`str`
            The type of region to use to create the consensus region set
            - one of "summits" or "peaks".
            If "summits", peak summits will be extended by ``extension``
            before union.
            If "peaks", sample peaks will be used with no modification prior to
            union.

            Default is "summits".
        extension : :obj:`int`
            Amount to extend peaks summits by in both directions.

            Default is 250.
        blacklist_bed : {:obj:`False`, :obj:`str`}
            Either :obj:`False` or a path to a BED file with genomic positions
            to exclude from consensus peak set.

            Default is to use a blacklist file for the analysis ``genome``.
        filter_chroms : {:obj:`list`, :obj:`str`}
            A list of chromosomes to filter out or
            a string with a pattern to match to exclude chromosomes.
            Uses Pandas string methods :class:`pandas.Series.str.match`.
            Pass for example `'.*_.*|chrM'` to filter out chromosomes with a "_"
            character and a "chrM" chromosome.

            Default is not to filter anything.
        permissive : :obj:`bool`
            Whether Samples that which ``region_type`` attribute file
            does not exist should be simply skipped or an error thrown.

        comparison_table : :obj:`pandas.DataFrame`, optional
            DataFrame with signal/background combinations used to call peaks.
            Part of kwargs.

            Defaults to analysis own ``comparison_table``.
        peak_dir : :obj:`str`, optional
            Path to peaks output directory. Part of kwargs.

            Defaults to "{analysis.results_dir}/chipseq_peaks".

        Attributes
        ----------
        sites : :class:`pybedtools.BedTool`
            Bedtool with consensus sites.
        """
        import re
        from ngs_toolkit.general import get_blacklist_annotations
        import pybedtools
        from tqdm import tqdm
        import tempfile

        if "comparison_table" not in kwargs:
            # TODO: allow not requiring peak_dir to be passed if specifying a new table
            self.set_comparisons(kwargs["comparison_table"], peak_dir=kwargs["peak_dir"])

        if region_type not in ["summits", "peaks"]:
            msg = "`region_type` attribute must be one of 'summits' or 'peaks'!"
            _LOGGER.error(msg)
            raise ValueError(msg)

        if blacklist_bed is None:
            _LOGGER.info("Blacklist file not provided. Downloading...")
            try:
                blacklist_bed = get_blacklist_annotations(self.organism, self.genome)
            except AttributeError:
                msg = "Blacklist file was not provided and cannot"
                msg += " get one without analysis having `organism` and `genome` set."
                _LOGGER.error(msg)
                raise AttributeError(msg)

        # Simply concatenate all peaks in one file
        f = tempfile.NamedTemporaryFile()
        with open(f.name, "a") as handle:
            for name, comp in tqdm(self.comparisons.items(), total=len(self.comparisons), desc="Comparison"):
                for peak_caller, peak_file in comp['peak_calls'][peak_type].items():
                    try:
                        # TODO: check if homer has summits and they match this pattern
                        summit = re.sub("_peaks.narrowPeak", "_summits.bed", peak_file)
                        file = (
                            pybedtools.BedTool(summit).slop(b=extension, genome=comp['genome']).fn
                            if region_type == "summits"
                            else peak_file)
                    except (ValueError, FileNotFoundError):
                        _LOGGER.warning("Input file for comparison {} ({}) not found!", (name, f))
                        if not permissive:
                            raise

                    for line in open(file, 'r'):
                        handle.write(line)

        # Merge overlaping peaks across comparisons
        sites = pybedtools.BedTool(f.name).sort().merge()

        # Filter
        # # remove blacklist regions
        if blacklist_bed is not False:
            if not isinstance(blacklist_bed, pybedtools.BedTool):
                blacklist = pybedtools.BedTool(blacklist_bed)
            sites = sites.intersect(v=True, b=blacklist)

        # # filter requested chromosomes
        if filter_chroms is not None:
            if isinstance(filter_chroms, list):
                sites = sites.filter(lambda x: x.chrom not in filter_chroms).saveas()
            elif isinstance(filter_chroms, str):
                s = sites.to_dataframe()
                sites = pybedtools.BedTool.from_dataframe(s.loc[~s['chrom'].str.match(filter_chroms)])

        # Save and assign
        if save:
            output_file = os.path.join(self.results_dir, self.name + ".peak_set.bed")
            sites.saveas(output_file)
            sites = pybedtools.BedTool(output_file)
        if assign:
            self.sites = sites
        return sites

    def calculate_peak_support(
            self, samples=None, region_type="summits", peak_type="filtered", permissive=True,
            comparison_table=None, peak_dir="{results_dir}/chipseq_peaks"):
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
        samples: :obj:`list`
            Not used. Provided for compatibility with ATACSeqAnalysis class.
        region_type: :obj:`str`
            Not used. Provided for compatibility with ATACSeqAnalysis class.
        permissive: :obj:`bool`
            Not used. Provided for compatibility with ATACSeqAnalysis class.

        Attributes
        ----------
        support : :obj:`pandas.DataFrame`
            DataFrame with signal/background combinations used to call peaks
        """
        import pybedtools
        from tqdm import tqdm
        from ngs_toolkit.utils import bed_to_index

        if comparison_table is None:
            comparison_table = self.comparison_table

        peak_dir = os.path.abspath(self._format_string_with_attributes(peak_dir))

        # get index
        index = bed_to_index(self.sites.to_dataframe())

        # calculate support (number of samples overlaping each merged peak)
        support = pd.DataFrame(index=index)
        for name, comp in tqdm(self.comparisons.items(), total=len(self.comparisons), desc="Comparison"):
            for peak_caller, peak_file in comp['peak_calls'][peak_type].items():
                try:
                    sample_support = self.sites.intersect(peak_file, wa=True, c=True).to_dataframe()
                except (
                    ValueError,
                    pybedtools.MalformedBedLineError,
                    pybedtools.helpers.BEDToolsError,
                ):
                    _LOGGER.warning(
                        "Peaks for comparison %s (%s) not found!", (name, peak_file))
                    if permissive:
                        continue
                    else:
                        raise
                sample_support.index = index
                support[(name, peak_caller)] = sample_support.iloc[:, 3]

        # Make multiindex labeling comparisons and peak type
        support.columns = pd.MultiIndex.from_tuples(
            support.columns, names=["comparison", "peak_caller"]
        )
        support.to_csv(
            os.path.join(
                self.results_dir, self.name + "_peaks.binary_overlap_support.csv"
            ),
            index=True,
        )

        # divide sum (of unique overlaps) by total to get support value between 0 and 1
        support["support"] = support.astype(bool).sum(axis=1) / float(support.shape[1])

        # save
        support.to_csv(
            os.path.join(self.results_dir, self.name + "_peaks.support.csv"), index=True
        )

        self.support = support

    def get_supported_peaks(self, samples=None, **kwargs):
        """
        Get mask of sites with 0 support in the given samples.
        Requires support matrix produced by `ngs_toolkit.atacseq.ATACSeqAnalysis.calculate_peak_support`.

        Parameters
        ----------
        samples : :obj:`list`
            Not used. Provided for compatibility with ATACSeqAnalysis class.

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

        # TODO: I believe the indexing of a dataframe with the sample objects themselves will not work
        # one would need to get their names
        if not by_group:
            res = pd.DataFrame(index=matrix.index)
            for name, comp in self.comparisons.items():
                for s in comp['signal_samples']:
                    res.loc[:, s] = comparison_func(
                        matrix.loc[:, s],
                        matrix.loc[:, comp['background_samples']].apply(reduction_func, axis=1)
                        if reduction_func != np.mean
                        else matrix.loc[:, comp['background_samples']].mean(axis=1))
        if by_group:
            comparisons = set(comparison_table['comparison_name'])
            res = pd.DataFrame(index=matrix.index, columns=comparisons)
            for name, comp in self.comparisons.items():
                res.loc[:, name] = comparison_func(
                    matrix.loc[:, comp['signal_samples']].apply(reduction_func, axis=1)
                    if reduction_func != np.mean
                    else matrix.loc[:, comp['signal_samples']].mean(axis=1),
                    matrix.loc[:, comp['background_samples']].apply(reduction_func, axis=1)
                    if reduction_func != np.mean
                    else matrix.loc[:, comp['background_samples']].mean(axis=1))
        return res
