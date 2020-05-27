#!/usr/bin/env python


import os

import numpy as np
import pandas as pd

from ngs_toolkit import _LOGGER
from ngs_toolkit.analysis import Analysis
from ngs_toolkit.decorators import check_has_attributes


class ATACSeqAnalysis(Analysis):
    """
    Class to model analysis of ATAC-seq data.
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
    >>> from ngs_toolkit.atacseq import ATACSeqAnalysis

    This is an example of the beginning of an ATAC-seq analysis:

    >>> pep = "metadata/project_config.yaml"
    >>> a = ATACSeqAnalysis(from_pep=pep)
    >>> # Get consensus peak set from all samples
    >>> a.get_consensus_sites(a.samples)
    >>> # Annotate regions
    >>> a.get_peak_gene_annotation()
    >>> a.get_peak_genomic_location()
    >>> # Get coverage values for each peak in each sample of ATAC-seq
    >>> a.measure_coverage()
    >>> # Normalize jointly (quantile normalization + GC correction)
    >>> a.normalize(method="gc_content")
    >>> # Annotate quantified peaks with previously calculated metrics and features
    >>> a.annotate_features()
    >>> # Annotate with sample metadata
    >>> a.annotate_samples()
    >>> # Save object
    >>> a.to_pickle()
    """

    _data_type = "ATAC-seq"

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
        **kwargs,
    ):
        # The check for existance is to make sure other classes can inherit from this
        default_args = {
            "data_type": "ATAC-seq",
            "__data_type__": "ATAC-seq",
            "var_unit_name": "region",
            "quantity": "accessibility",
            "norm_units": "RPM",
        }
        for k, v in default_args.items():
            if not hasattr(self, k):
                setattr(self, k, v)

        super(ATACSeqAnalysis, self).__init__(
            name=name,
            from_pep=from_pep,
            from_pickle=from_pickle,
            root_dir=root_dir,
            data_dir=data_dir,
            results_dir=results_dir,
            prj=prj,
            samples=samples,
            **kwargs,
        )

    def load_data(
        self, output_map=None, only_these_keys=None, prefix="{results_dir}/{name}", permissive=True,
    ):
        """
        Load the output files of the major functions of the Analysis.

        Parameters
        ----------
        output_map : :obj:`dict`
            Dictionary with "attribute name": "path prefix" to load the files.

        only_these_keys : :obj:`list`, optional
            Iterable of analysis attributes to load up.
            Possible attributes:

                * "matrix_raw"
                * "matrix_norm"
                * "matrix_features"
                * "sites"
                * "support"
                * "nuc"
                * "coverage_gc_corrected"
                * "gene_annotation"
                * "region_annotation"
                * "region_annotation_b"
                * "chrom_state_annotation"
                * "chrom_state_annotation_b"
                * "stats"
                * "differential_results"

            Default is all of the above.
        prefix : :obj:`str`, optional
            String prefix of files to load.
            Variables in curly braces will be formated with attributes of analysis.

            Defaults to "{results_dir}/{name}".
        bool : permissive, optional
            Whether an error should be ignored if reading a file causes IOError.

            Default is :obj:`True`.

        Attributes
        ----------
        sites : :class:`pybedtools.bedtool.BedTool`
            Sets a `sites` variable.

        pandas.DataFrame
            Dataframes holding the respective data, available as attributes described
            in the `only_these_keys` parameter.

        Raises
        ----------
        IOError
            If not permissive and a file is not found
        """
        import pybedtools

        prefix = self._format_string_with_attributes(prefix)

        if output_map is None:
            kwargs = {"index_col": 0}
            output_map = {
                "sites": (prefix + ".peak_set.bed", {}),
                "matrix_raw": (prefix + ".matrix_raw.csv", kwargs),
                "matrix_norm": (prefix + ".matrix_norm.csv", kwargs),
                "matrix_features": (prefix + ".matrix_features.csv", kwargs),
                "support": (prefix + ".support.csv", kwargs),
                "nuc": (prefix + ".gccontent_length.csv", kwargs),
                "gene_annotation": (prefix + ".gene_annotation.csv", kwargs),
                "closest_tss_distances": (prefix + ".closest_tss_distances.csv", kwargs,),
                "region_annotation": (prefix + ".region_annotation.csv", kwargs),
                "region_annotation_b": (prefix + ".region_annotation_background.csv", kwargs,),
                "region_annotation_mapping": (prefix + ".region_annotation_mapping.csv", kwargs,),
                "region_annotation_b_mapping": (
                    prefix + ".region_annotation_background_mapping.csv",
                    kwargs,
                ),
                "chrom_state_annotation": (prefix + ".chrom_state_annotation.csv", kwargs,),
                "chrom_state_annotation_b": (
                    prefix + ".chrom_state_annotation_background.csv",
                    kwargs,
                ),
                "chrom_state_annotation_mapping": (
                    prefix + ".chrom_state_annotation_mapping.csv",
                    kwargs,
                ),
                "chrom_state_annotation_b_mapping": (
                    prefix + ".chrom_state_annotation_background_mapping.csv",
                    kwargs,
                ),
                "stats": (prefix + ".stats_per_feature.csv", kwargs),
                "differential_results": (
                    os.path.join(
                        self.results_dir,
                        "differential_analysis_{}".format(self.data_type),
                        "differential_analysis.deseq_result.all_comparisons.csv",
                    ),
                    kwargs,
                ),
            }

        if only_these_keys is None:
            only_these_keys = list(output_map.keys())

        # Use the parent method just with an updated output_map dictionary
        Analysis.load_data(
            self,
            output_map={k: v for k, v in output_map.items() if k != "sites"},
            only_these_keys=only_these_keys,
            prefix=prefix,
            permissive=permissive,
        )

        # Special case
        if "sites" in only_these_keys:
            _LOGGER.info("Loading 'sites' analysis attribute.")
            file = os.path.join(self.results_dir, self.name + ".peak_set.bed")
            try:
                setattr(self, "sites", pybedtools.BedTool(file))
            except IOError as e:
                if not permissive:
                    raise e
                else:
                    _LOGGER.warning(e)

    @staticmethod
    def check_region_index(matrix):
        return matrix.index.str.contains(":").all() and matrix.index.str.contains("-").all()

    @staticmethod
    def set_region_index(matrix, force=False):
        from ngs_toolkit.utils import bed_to_index

        if (not force) and ATACSeqAnalysis.check_region_index(matrix):
            _LOGGER.warning("Matrix already has well-formatted index.")
            return matrix
        else:
            req = ["chrom", "start", "end"]
            if all([x in matrix.columns for x in req]):
                matrix.index = bed_to_index(matrix)
            else:
                raise ValueError(
                    "Could not format matrix index. Missing '{}' columns.".format(",".join(req))
                )

    def get_consensus_sites(
        self,
        samples=None,
        region_type="summits",
        extension=250,
        blacklist_bed=None,
        filter_chroms=None,
        permissive=False,
        save=True,
        assign=True,
        **kwargs,
    ):
        """
        Get consensus (union) of enriched sites (peaks) across samples.
        There are two modes possible, defined by the value of ``region_type``:

         * peaks: simple union of all sites;
         * summits: peak summits are extended by ``extension`` and a union is made.

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
            Pass for example `".*_.*|chrM"` to filter out chromosomes with a "_"
            character and a "chrM" chromosome.

            Default is not to filter anything.
        permissive : :obj:`bool`
            Whether Samples that which ``region_type`` attribute file does not exist
            should be simply skipped or an error thrown.

            Default is :obj:`True`.
        **kwargs
            Not used. Provided for compatibility with :class:`ngs_toolkit.ChIPSeqAnalysis` class.

        Raises
        ----------
        :obj:`ValueError`
            If not ``permissive`` and either the ``peaks`` or ``summits`` file
            of a sample is not readable, or if ``permissive``
            but none of the samples has an existing file.
        :obj:`AttributeError`
            If analysis does not have ``organism`` and ``genome`` attributes.

        Attributes
        ----------
        sites : :class:`pybedtools.bedtool.BedTool`
            Sets a ``sites`` variable with the consensus peak set.

        Returns
        ----------
        sites : :class:`pybedtools.bedtool.BedTool`
            The consensus peak set.
        """
        from tqdm import tqdm
        import pybedtools
        import tempfile

        if region_type not in ["summits", "peaks"]:
            msg = "`region_type` attribute must be one of 'summits' or 'peaks'!"
            _LOGGER.error(msg)
            raise ValueError(msg)

        if samples is None:
            samples = self.samples

        # Check which samples to run (dependent on permissive)
        samples = self._get_samples_with_input_file(
            region_type, permissive=permissive, samples=samples
        )

        if blacklist_bed is not False and blacklist_bed is None:
            _LOGGER.info("Blacklist file not provided. Downloading...")
            try:
                blacklist_bed = self.get_resources(steps=["blacklist"])["blacklist_file"]
            except AttributeError:
                msg = "Blacklist file was not provided and cannot be"
                msg += " get one without `organism` and `genome` set."
                _LOGGER.error(msg)
                raise AttributeError(msg)

        f = tempfile.NamedTemporaryFile()
        with open(f.name, "a") as handle:
            for sample in tqdm(samples, total=len(samples), desc="Sample"):
                try:
                    file = (
                        pybedtools.BedTool(sample.summits)
                        .slop(b=extension, genome=sample.genome)
                        .fn
                        if region_type == "summits"
                        else sample.peaks
                    )
                except (ValueError, FileNotFoundError):
                    if not permissive:
                        raise
                    else:
                        _LOGGER.warning(
                            "Peaks for sample {} ({}) not found!".format(sample, sample.peaks)
                        )
                        continue
                for line in open(file, "r"):
                    handle.write(line)

        # NCBI genome FASTA files are sorted naturally while Ensembl are not
        # depending on that you might want to sort the resulting BED file
        # accordingly with the following:
        #     sites = sort_bed(f.name).merge()
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
                sites = pybedtools.BedTool.from_dataframe(
                    s.loc[~s["chrom"].str.match(filter_chroms)]
                )

        # Save and assign
        if save:
            output_file = os.path.join(self.results_dir, self.name + ".peak_set.bed")
            sites.saveas(output_file)
            self.record_output_file(output_file, "consensus_sites")
            sites = pybedtools.BedTool(output_file)
        if assign:
            self.sites = sites
        return sites

    def set_consensus_sites(self, bed_file, overwrite=True):
        """
        Set consensus (union) sites across samples given a BED file.

        Parameters
        ----------
        bed_file : :obj:`str`
            BED file to use as consensus sites.

        overwrite : :obj:`bool`
            Whether a possibly existing file with a consensus peak set
            for this analysis should be overwritten in disk.

        Attributes
        ----------
        sites : :class:`~pybedtools.BedTool`
            Sets a `sites` variable with consensus peak set.
        """
        import pybedtools

        self.sites = pybedtools.BedTool(bed_file)
        if overwrite:
            default_sites = os.path.join(self.results_dir, self.name + ".peak_set.bed")
            # pybedtools will pipe to the input file!
            if bed_file == default_sites:
                self.sites.saveas(default_sites + ".new")
                os.rename(default_sites + ".new", default_sites)
            else:
                self.sites.saveas(default_sites)
        # TODO: warn if not overwrite and file exists already

    @check_has_attributes(["sites"])
    def calculate_peak_support(
        self,
        samples=None,
        region_type="summits",
        permissive=False,
        comparison_table=None,
        peak_dir=None,
    ):
        """
        Count number of called peaks per sample in the consensus region set.
        In addition calculate a measure of peak support (or ubiquitouness) by
        observing the ratio of samples containing a peak overlapping each region.

        Parameters
        ----------
        samples : :obj:`list`
            Iterable of :class:`peppy.Sample` objects to restrict to.
            Must have a ``peaks`` attribute set.

            Defaults to all samples in the analysis (``samples`` attribute).
        region_type : :obj:`str`
            The type of region to use to create the consensus region set.
            One of "summits" or "peaks".
            If `summits`, peak summits will be extended by ``extension`` before union.
            Otherwise sample peaks will be used with no modification.

            Default is "summits".
        permissive: :obj:`bool`
            Whether Samples that which `region_type` attribute file does
            not exist should be simply skipped or an error thrown.

        comparison_table: :obj:`pandas.DataFrame`
            Not used. Provided for compatibility with ChIPSeqAnalysis class.
        peak_dir: :obj:`str`
            Not used. Provided for compatibility with ChIPSeqAnalysis class.

        Raises
        ----------
        IOError
            If not ``permissive`` and either the ``peaks`` or ``summits`` file of a sample is not readable.
            Or if ``permissive`` but none of the samples has an existing file.

        Attributes
        ----------
        support : :obj:`pandas.DataFrame`
            A dataframe with counts of peaks overlapping each feature of consensus set.
        """
        # TODO: Implement distributed
        from tqdm import tqdm
        import pybedtools
        import tempfile
        from ngs_toolkit.utils import bed_to_index

        if samples is None:
            samples = self.samples

        # Check which samples to run (dependent on permissive)
        samples = self._get_samples_with_input_file(
            region_type, permissive=permissive, samples=samples
        )

        # calculate support (number of samples overlaping each merged peak)
        for i, sample in tqdm(enumerate(samples), total=len(samples), desc="Sample"):
            if region_type == "summits":
                peaks = sample.summits
            else:
                peaks = sample.peaks
            # print(sample, peaks)
            if i == 0:
                support = self.sites.intersect(peaks, wa=True, c=True)
            else:
                support = support.intersect(peaks, wa=True, c=True)

        try:
            support = support.to_dataframe()
        except (
            ValueError,
            pybedtools.MalformedBedLineError,
            pybedtools.helpers.BEDToolsError,
            OverflowError,
        ):
            _LOGGER.debug(
                "Could not convert support intersection directly to dataframe,"
                "saving/reading temporary file."
            )
            t = tempfile.NamedTemporaryFile()
            support.saveas(t.name)
            support = pd.read_csv(t.name, sep="\t", header=None)

        support.columns = ["chrom", "start", "end"] + [sample.name for sample in samples]
        support.index = bed_to_index(support)
        support.to_csv(
            os.path.join(self.results_dir, self.name + ".binary_overlap_support.csv"), index=True,
        )

        # divide sum (of unique overlaps) by total to get support value between 0 and 1
        support.loc[:, "support"] = support[[sample.name for sample in samples]].apply(
            lambda x: sum([i if i <= 1 else 1 for i in x]) / float(len(samples)), axis=1
        )
        # save
        support.to_csv(os.path.join(self.results_dir, self.name + ".support.csv"), index=True)

        setattr(self, "support", support)
        return self.support

    def get_supported_peaks(self, samples=None, **kwargs):
        """
        Get mask of sites with 0 support in the given samples.
        Requires support matrix produced by `ngs_toolkit.atacseq.ATACSeqAnalysis.calculate_peak_support`.

        Parameters
        ----------
        samples : :obj:`list`
            Iterable of :class:`peppy.Sample` objects to restrict to.

        **kwargs
            Not used. Provided for compatibility with ChIPSeqAnalysis class.

        Returns
        -------
        pd.Series
            Boolean Pandas Series with sites with at least one of the given samples having a peak called.
        """
        if samples is None:
            samples = self.samples
        return self.support.loc[:, [s.name for s in samples]].sum(1) != 0

    def measure_coverage(
        self,
        samples=None,
        sites=None,
        save=True,
        assign=True,
        peak_set_name="peak_set",
        output_file="{results_dir}/{name}.matrix_raw.csv",
        permissive=False,
        distributed=False,
        overwrite=True,
        **kwargs,
    ):
        """
        Measure read coverage (counts) of each sample in each region
        in consensus sites.
        Uses parallel computing using the :class:`parmap` library.
        However, for many samples (hundreds), parallelization in a
        computing cluster is possible with the `distributed` option.

        Parameters
        ----------
        samples : :obj:`list`
            Iterable of :class:`peppy.Sample` objects to restrict to.
            Must have a ``aligned_filtered_bam`` attribute set.

            Defaults to all samples in the analysis (``samples`` attribute).
        sites : {:class:`pybedtools.bedtool.BedTool`, :class:`pandas.DataFrame`, :obj:`str`}
            Sites in the genome to quantify, usually a :class:`pybedtools.bedtool.BedTool`
            from :class:`ngs_toolkit.atacseq.ATACSeqAnalysis.get_consensus_sites`.
            If a DataFrame, will try to convert to BED format assuming first
            three columns are chr,start,end.
            If a string assumes a path to a BED file.

            Defaults to ``sites`` attribute of analysis object.
        save : :obj:`bool`
            Whether to save to disk the coverage matrix with filename ``output_file``.

            Default is :obj:`True`.
        assign : :obj:`bool`
            Whether to assign the matrix to an attribute named ``coverage``.

            Default is :obj:`True`.
        peak_set_name : :obj:`bool`
            Suffix to files containing coverage of ``distributed`` is True.

            Defaults to "peak_set".
        output_file : :obj:`str`
            A path to a CSV file with coverage output.

            Default is "{results_dir}/{name}.raw_coverage.csv".
        permissive : :obj:`bool`
            Whether Samples for which ``region_type`` attribute file does not exist
            should be simply skipped or an error thrown.

            Default is :obj:`False`.
        distributed : :obj:`bool`
            Whether it should be run as jobs for each sample
            separately in parallel.

            Default is :obj:`False`.
        overwrite : :obj:`bool`
            Whether to overwrite existing files if ``distributed`` is True.

            Default is :obj:`True`.
        **kwargs : :obj:`dict`
            Additional keyword arguments will be passed to
            `ngs_toolkit.utils.submit_job` if `distributed` is True,
            and on to a divvy submission template.
            Pass for example: computing_configuration="slurm", jobname="job",
            cores=2, mem=8000, partition="longq".

        Raises
        ----------
        IOError
            If not ``permissive`` and the 'aligned_filtered_bam'
            file attribute of a sample is not readable.
            Or if ``permissive`` but none of the samples has an existing file.

        Attributes
        ----------
        matrix_raw : :class:`pandas.DataFrame`
            The dataframe of raw coverage values (counts) of
            shape (n_features, m_samples).

        Returns
        -------
        :class:`pandas.DataFrame`
            Pandas DataFrame with read counts of shape (n_sites, m_samples).
        """
        import sys
        import multiprocessing
        import parmap

        from ngs_toolkit.utils import count_reads_in_intervals, submit_job, to_bed_index

        if samples is None:
            samples = self.samples

        output_file = self._format_string_with_attributes(output_file)

        # Check which samples to run (dependent on permissive)
        samples = self._get_samples_with_input_file(
            "aligned_filtered_bam", samples=samples, permissive=permissive
        )

        if sites is None:
            sites = self.sites

        if not distributed:
            # Count reads with pysam
            # make strings with intervals
            sites_str = to_bed_index(sites)
            # count, create dataframe
            matrix_raw = pd.DataFrame(
                map(
                    pd.Series,
                    parmap.map(
                        count_reads_in_intervals,
                        [sample.aligned_filtered_bam for sample in samples],
                        sites_str,
                        pm_parallel=True,
                    ),
                ),
                index=[sample.name for sample in samples],
            ).T

            if assign:
                self.matrix_raw = matrix_raw
            if save:
                matrix_raw.to_csv(output_file, index=True)
            return matrix_raw

        else:
            for s in samples:
                output_dir = os.path.join(s.sample_root, "coverage")
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)

                job_name = "{}.{}_coverage".format(peak_set_name, s.name)
                prefix = os.path.join(output_dir, s.name + ".{}_coverage".format(peak_set_name))
                output_file = prefix + ".bed"
                log_file = prefix + ".log"
                job_file = prefix + ".sh"
                if not overwrite:
                    if os.path.exists(output_file):
                        continue

                cmd = "\\\n".join(
                    [
                        "{executable} -m ngs_toolkit.recipes.coverage ",
                        "--no-overwrite" if not overwrite else "",
                        "{input_bed} {input_bam} {output_bed}",
                    ]
                ).format(
                    executable=sys.executable,
                    input_bed=sites.fn,
                    input_bam=s.aligned_filtered_bam,
                    output_bed=output_file,
                )
                for k, v in [("cores", 1), ("mem", 8000), ("time", "04:00:00")]:
                    if k not in kwargs:
                        kwargs[k] = v

                submit_job(cmd, job_file, jobname=job_name, logfile=log_file, **kwargs)
            if getattr(kwargs, "computing_configuration", None) in ["localhost", "default"]:
                _LOGGER.info("Collecting job results.")
                return self.collect_coverage(
                    samples=samples,
                    save=save,
                    assign=assign,
                    output_file=output_file,
                    permissive=permissive,
                )

    def collect_coverage(
        self,
        samples=None,
        save=True,
        assign=True,
        output_file=None,
        permissive=False,
        peak_set_name="peak_set",
        fast_and_unsafe=False,
    ):
        """
        Collect read coverage (counts) of each sample in each region in consensus sites from existing files.
        Useful after runnning analysis.measure_coverage() in distributed mode.

        Parameters
        ----------
        samples : :obj:`list`
            Iterable of :class:`peppy.Sample` objects to restrict to.
            If not provided (`None` is passed) if will default to all samples in the analysis (``samples`` attribute).

        save: :obj:`bool`
            Whether to save to disk the coverage matrix with filename ``output_file``.

        assign: :obj:`bool`
            Whether to assign the matrix to an attribute of self named ``coverage``.

        output_file : :obj:`str`
            A path to a CSV file with coverage output.

            Default is "{results_dir}/{name}.raw_coverage.csv".
        permissive: :obj:`bool`
            Whether Samples without an existing coverage file does not exist
            should be simply skipped or an error thrown.

        peak_set_name: :obj:`bool`
            Suffix to files containing coverage.
            Defaults to "peak_set".

        fast_and_unsafe: :obj:`bool`
            Whether to use a faster but unsafer method to concatenate the data.
            If the order of all rows in all samples is the same then the result should be the same.
            The default, slower method assures that all rows are matched and is therefore slower.

            Defaults to :obj:`False`.

        Raises
        ----------
        IOError
            If not ``permissive`` and the coverage file of a sample is not readable or is empty.
            Or if ``permissive`` but none of the samples has an existing file or are empty.

        Attributes
        ----------
        matrix_raw : :class:`pandas.DataFrame`
            The dataframe of raw coverage values (counts) of shape (n_features, m_samples).

        Returns
        -------
        :class:`pandas.DataFrame`
            Pandas DataFrame with read counts of shape (n_sites, m_samples).
        """
        from tqdm import tqdm

        from ngs_toolkit.utils import bed_to_index

        if samples is None:
            samples = self.samples

        for sample in samples:
            setattr(
                sample,
                "_coverage",
                os.path.join(
                    sample.sample_root,
                    "coverage",
                    sample.name + ".{}_coverage.bed".format(peak_set_name),
                ),
            )

        # Check which samples to run (dependent on permissive)
        samples = self._get_samples_with_input_file(
            "_coverage", permissive=permissive, samples=samples
        )

        # Read in counts
        matrix_raw = list()
        for sample in tqdm(samples, total=len(samples)):
            cov = pd.read_csv(
                sample._coverage,
                sep="\t",
                header=None,
                names=["chrom", "start", "end", sample.name],
            )

            if cov.empty:
                msg = "Coverage file for sample '{}' is empty!".format(sample.name)
                if not permissive:
                    _LOGGER.error(msg)
                    raise IOError(msg)
                else:
                    _LOGGER.warning(msg)
                    continue
            matrix_raw.append(cov)

        if len(matrix_raw) == 0:
            msg = "No sample had a valid coverage file!"
            if permissive:
                _LOGGER.warning(msg)
                return
            else:
                _LOGGER.error(msg)
                raise IOError(msg)

        if fast_and_unsafe:
            _LOGGER.warning("Using a concatenation method that is not 100% safe.")
            matrix_raw = pd.concat(matrix_raw, axis=1, sort=False).dropna().sort_index()
            matrix_raw = (
                matrix_raw.loc[:, ~matrix_raw.columns.duplicated()]
                .set_index(["chrom", "start", "end"])
                .astype(int)
            )
        else:
            matrix_raw = (
                pd.concat(matrix_raw, axis=0, sort=False)
                .melt(id_vars=["chrom", "start", "end"])
                .pivot_table(
                    index=["chrom", "start", "end"],
                    columns="variable",
                    values="value",
                    fill_value=0,
                )
                .astype(int)
            )
        matrix_raw.index = bed_to_index(matrix_raw.index.to_frame())

        if assign:
            self.matrix_raw = matrix_raw
        if save:
            if output_file is not None:
                matrix_raw.to_csv(output_file, index=True)
            else:
                self.matrix_raw.to_csv(
                    os.path.join(self.results_dir, self.name + ".matrix_raw.csv"), index=True,
                )
        return matrix_raw

    @check_has_attributes(["organism", "genome"])
    def get_peak_gccontent_length(self, bed_file=None, fasta_file=None):
        """
        Get length and GC content of features in region set.

        bed_file : :obj:`str`
            A BED file with regions to calculate GC content on. Must be a 3-column BED!
            If not provided the calculation will be for the analysis `sites` attribute.

        genome : :obj:`str`
            Genome assembly.

        fasta_file : :obj:`str`
            Fasta file of `genome`. Preferably indexed. If not given, will try to download.

        :var nuc:
            DataFrame with nucleotide content and length of each region.

        Returns
        -------
        pandas.DataFrame
            DataFrame with nucleotide content and length of each region.

        Attributes
        ----------
        nuc : :class:`pandas.DataFrame`
            Dataframe with length and GC-content of each feature.
        """
        import pybedtools
        from ngs_toolkit.utils import bed_to_index

        if bed_file is None:
            sites = self.sites
        else:
            sites = pybedtools.BedTool(bed_file)

        if fasta_file is None:
            _LOGGER.info("Reference genome FASTA file was not given, will try to get it.")
            _LOGGER.info(
                "Getting genome FASTA file for organism '{}', genome '{}'. ".format(
                    self.organism, self.genome
                )
            )
            fasta_file = self.get_resources(steps=["genome"])["genome_file"]["fasta"]

        nuc = sites.nucleotide_content(fi=fasta_file).to_dataframe(comment="#")[
            ["score", "blockStarts"]
        ]
        nuc.columns = ["gc_content", "length"]
        nuc.index = bed_to_index(sites)

        # get only the sites matching the coverage (not overlapping blacklist)
        self.nuc = nuc.loc[self.matrix_raw.index]

        self.nuc.to_csv(
            os.path.join(self.results_dir, self.name + ".gccontent_length.csv"), index=True,
        )

        return self.nuc

    def normalize_cqn(self, matrix="matrix_raw", samples=None, save=True, assign=True):
        """
        Conditional quantile normalization (CQN) of a matrix.
        It uses GC content and length of regulatory elements as covariates.

        Requires the R package "cqn" to be installed:

        .. highlight:: R
        .. code-block:: R

            if (!requireNamespace("BiocManager", quietly = TRUE))
                install.packages("BiocManager")
            BiocManager::install("cqn")

        Parameters
        ----------
        matrix : :obj:`str`
            Attribute name of matrix to normalize.

            Defaults to "matrix_raw".
        samples : :obj:`list`
            Iterable of :class:`peppy.Sample` objects to restrict matrix to.

            Defaults to all samples in analysis.
        save: :obj:`bool`
            Whether to write normalized DataFrame to disk.

            Default is :obj:`True`.
        assign: :obj:`bool`
            Whether to assign the normalized DataFrame to an attribute `matrix_norm`.

            Default is :obj:`True`.

        Attributes
        ----------
        matrix_norm : :class:`pandas.DataFrame`
            If `assign`, the dataframe with normalized values.
        norm_method : :obj:`str`
            If ``assign``, it is the name of method used to normalize: "cqn".
        """
        from ngs_toolkit.utils import cqn

        # Perform quantile normalization first
        if not hasattr(self, "nuc"):
            self.normalize_quantiles(matrix=matrix, samples=samples)

        # Get GC content and length of each feature
        if not hasattr(self, "nuc"):
            self.get_peak_gccontent_length()

        matrix_norm = cqn(
            matrix=self.get_matrix(matrix=matrix, samples=samples),
            gc_content=self.nuc["gc_content"],
            lengths=self.nuc["length"],
        )

        if save:
            matrix_norm.to_csv(
                os.path.join(self.results_dir, self.name + ".matrix_norm.csv"), index=True,
            )
        if assign:
            self.matrix_norm = matrix_norm
            self.norm_method = "cqn"

        return matrix_norm

    @check_has_attributes(["sites"])
    def get_peak_gene_annotation(
        self, tss_file=None, max_dist=100000, save=True, output_prefix="", assign=True
    ):
        """
        Annotates peaks with closest gene.
        The annotation reference can either be given in the `tss_file` parameter
        but if ommited, it will be fetched if analysis has `genome` and `organism`
        attributes.
        A dataframe with each feature's distance to the nearest gene is also saved.

        Parameters
        ----------
        tss_file : :obj:`str`, optional
            A valid BED file where the name field (4th column) identifies the gene
            and the strand column (6th column). Other fields will not be used.

            Default is to get gene position annotations.
        max_dist : :obj:`int`, optional
            Maximum absolute distance allowed to perform associations.
            Regions with no genes within the range will have NaN values.

            Default is 100000.
        save: :obj:`bool`, optional
            Whether to write the annotated DataFrame to disk.

            Default is :obj:`True`.
        output_prefix: :obj:`str`, optional
            Prefix to add to output file when save is True.

            Default is "" (empty string).
        assign: :obj:`bool`, optional
            Whether to assign the DataFrames to `Attributes`.

            Default is :obj:`True`.

        Attributes
        ----------
        gene_annotation : :class:`pandas.DataFrame`
            A pandas DataFrame containing the genome annotations of the region features.
            If a feature overlaps more than one gene, the two gene values will be concatenated with a comma.

        closest_tss_distances : :class:`pandas.DataFrame`
            A pandas DataFrame containing unique region->gene associations.
            In contrast to gene_annotation dataframe, this contains one row per region->gene assignment.

        Returns
        -------
        pandas.DataFrame
            A dataframe with genes annotated for the peak set.
        """
        import pybedtools
        from ngs_toolkit.utils import bed_to_index, get_this_file_or_timestamped

        cols = [6, 8, -1]  # gene_name, strand, distance

        if tss_file is None:
            _LOGGER.info("Reference TSS file was not given, will try to get TSS annotations.")
            _LOGGER.info(
                "Getting TSS annotations for organism '{}', genome '{}'.".format(
                    self.organism, self.genome
                )
            )
            tss_file = self.get_resources(steps=["tss"])["tss_file"]

        # extract only relevant columns
        tss_file = get_this_file_or_timestamped(tss_file)
        if isinstance(self.sites, str):
            self.sites = pybedtools.BedTool(self.sites)

        # get closest TSS of each region
        tss = pybedtools.BedTool(tss_file)
        columns = ["chrom", "start", "end", "gene_name", "strand", "distance"]
        closest_tss_distances = self.sites.closest(tss, D="b").to_dataframe()

        closest_tss_distances = closest_tss_distances.iloc[:, [0, 1, 2] + cols]
        closest_tss_distances.columns = columns

        # set NaN to distance without assignment (rather than the default '-1' from bedtools)
        closest_tss_distances.loc[closest_tss_distances["gene_name"] == ".", "distance"] = np.nan

        # set NaN to assignments out of range
        closest_tss_distances.loc[
            closest_tss_distances["distance"].abs() > max_dist, ["gene_name", "strand", "distance"],
        ] = np.nan

        # aggregate annotation per peak, concatenate various genes (comma-separated)
        gene_annotation = (
            closest_tss_distances.groupby(["chrom", "start", "end"])
            .aggregate(lambda x: ",".join(set([str(i) for i in x if i != "."])))
            .reset_index()
        )
        closest_tss_distances.index = bed_to_index(closest_tss_distances)
        gene_annotation.index = bed_to_index(gene_annotation)

        # save to disk
        if save:
            if output_prefix != "":
                output_prefix += "."
            closest_tss_distances.to_csv(
                os.path.join(
                    self.results_dir,
                    self.name + ".closest_tss_distances.{}csv".format(output_prefix),
                ),
                index=True,
            )
            gene_annotation.to_csv(
                os.path.join(
                    self.results_dir, self.name + ".gene_annotation.{}csv".format(output_prefix)
                ),
                index=True,
            )
        if assign:
            self.closest_tss_distances = closest_tss_distances
            self.gene_annotation = gene_annotation
        return gene_annotation

    @check_has_attributes(["organism", "genome", "sites"])
    def get_peak_genomic_location(
        self, genomic_context_file=None, save=True, output_prefix="", assign=True
    ):
        """
        Annotates a consensus peak set (``sites`` attribute of analysis) with their genomic context.
        The genomic context is mostly gene-centric, which includes overlap with
        gene promoters, UTRs, exons, introns and remaining intergenic space.

        If no reference genomic annotation file is given (genomic_context_file kwarg),
        it will use the ngs_toolkit.general.get_genomic_context function to get
        such data. For more customization of the annotations, use that function directly and pass
        the output file to this function.

        Parameters
        ----------
        genomic_context_file : :obj:`str`
            A 4 column BED file (chrom, start, end, feature), where feature is a string with the type of region.
            If not provided will be get with the get_genomic_context function.
        save: :obj:`bool`, optional
            Whether to write the annotated DataFrame to disk.

            Default is :obj:`True`.
        output_prefix: :obj:`str`, optional
            Prefix to add to output file when save is True.

            Default is "" (empty string).
        assign: :obj:`bool`, optional
            Whether to assign the DataFrames to `Attributes`.

            Default is :obj:`True`.

        Attributes
        ----------
        region_annotation, region_annotation_b : :class:`pandas.DataFrame`
            A DataFrame with the genome annotations of the region features or genome background.

        region_annotation_mapping, region_annotation_b_mapping : :class:`pandas.DataFrame`
            A DataFrame with one row for each chromatin state-region mapping or genome background.

        Returns
        -------
        pandas.DataFrame
            The genomic context annotation for the peak set.
        """
        import pybedtools

        from ngs_toolkit.utils import bed_to_index, get_this_file_or_timestamped

        if genomic_context_file is None:
            _LOGGER.info("Reference genomic context file was not given, will try to get it.")
            _LOGGER.info(
                "Getting genomic context annotations for organism '{}', genome '{}'. ".format(
                    self.organism, self.genome
                )
            )
            genomic_context_file = self.get_resources(steps=["genomic_context"])[
                "genomic_context_file"
            ]

        context = pybedtools.BedTool(get_this_file_or_timestamped(genomic_context_file))

        if isinstance(self.sites, str):
            self.sites = pybedtools.BedTool(self.sites)

        # create background
        # shuffle regions in genome to create background (keep them in the same chromossome)
        background = self.sites.shuffle(genome=self.genome, chrom=True)

        cols = [0, 1, 2, -1]
        for label, attr, bed in [
            ("background", "region_annotation_b", background),
            ("real", "region_annotation", self.sites),
        ]:
            annot = (
                bed.intersect(context, wa=True, wb=True, f=0.2).sort().to_dataframe().iloc[:, cols]
            )
            annot.index = bed_to_index(annot)
            annot.columns = ["chrom", "start", "end", "genomic_region"]

            # remove duplicates (there shouldn't be anyway)
            annot = annot.drop_duplicates()
            # join various annotations per peak
            annot_comp = (
                annot.groupby(["chrom", "start", "end"])
                .aggregate(lambda x: ",".join(set([str(i) for i in x])))
                .reset_index()
            )
            annot_comp.index = bed_to_index(annot_comp)
            annot_comp.columns = ["chrom", "start", "end", "genomic_region"]
            # save to disk
            if save:
                a = "" if (label == "real") else ("_" + label)
                if output_prefix != "":
                    output_prefix += "."
                annot.to_csv(
                    os.path.join(
                        self.results_dir,
                        self.name + ".region_annotation{}_mapping.{}csv".format(a, output_prefix),
                    ),
                    index=True,
                )
                annot_comp.to_csv(
                    os.path.join(
                        self.results_dir,
                        self.name + ".region_annotation{}.{}csv".format(a, output_prefix),
                    ),
                    index=True,
                )

            if assign:
                setattr(self, attr, annot_comp)
                setattr(self, attr + "_mapping", annot)
        return self.region_annotation

    @check_has_attributes(["organism", "genome", "sites"])
    def get_peak_chromatin_state(
        self, chrom_state_file, frac=0.2, save=True, output_prefix="", assign=True
    ):
        """
        Annotates a consensus peak set (``sites`` attribute of analysis)
        with their chromatin state context.
        This would be given, for example by a chromatin state segmentation
        file from projects such as Roadmap Epigenomics.

        See examples of such files for various cell types/assemblies
        `in this website <https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/>`_
        (the "dense.bed.gz" files are optimal).

        Parameters
        ----------
        chrom_state_file : :obj:`str`
            A 4 column BED file (chrom, start, end, feature), where feature is a string with the type of region.
            Additional columns are ignored.

        frac : float
            Minimal fraction of region to overlap with a feature.

            Defaults to 0.2.
        save: :obj:`bool`, optional
            Whether to write the annotated DataFrame to disk.

            Default is :obj:`True`.
        output_prefix: :obj:`str`, optional
            Prefix to add to output file when save is True.

            Default is "" (empty string).
        assign: :obj:`bool`, optional
            Whether to assign the DataFrames to `Attributes`.

            Default is :obj:`True`.

        Returns
        ----------
        pandas.DataFrame
            The chromatin state annotation for the peak set.

        Attributes
        ----------
        chrom_state_annotation, chrom_state_annotation_b : :class:`pandas.DataFrame`
            A DataFrame with the chromatin state annotations of
            the region features or of the genome background.

        chrom_state_annotation_mapping, chrom_state_annotation_b_mapping : :class:`pandas.DataFrame`
            A DataFrame with one row for each chromatin state-region mapping
            or for the genome background.
        """
        import pybedtools
        from ngs_toolkit.utils import bed_to_index

        states = pybedtools.BedTool(chrom_state_file)

        if isinstance(self.sites, str):
            self.sites = pybedtools.BedTool(self.sites)

        # create background
        # shuffle regions in genome to create background (keep them in the same chromossome)
        background = self.sites.shuffle(genome=self.genome, chrom=True)

        for label, attr, bed in [
            ("real", "chrom_state_annotation", self.sites),
            ("background", "chrom_state_annotation_b", background),
        ]:
            _LOGGER.debug("Overlapping chromatin state annotation with {} regions.".format(label))
            annot = bed.intersect(states, wa=True, wb=True, f=frac, loj=True)
            try:
                annot = annot.to_dataframe(usecols=[0, 1, 2, 6])
            except pd.errors.EmptyDataError:
                _LOGGER.error("Could not annotate region set.")
                raise

            annot.iloc[:, 3] = annot.iloc[:, 3].astype(str)

            # remove duplicates (there shouldn't be anyway)
            annot = annot.drop_duplicates()
            annot.index = bed_to_index(annot)
            annot.columns = ["chrom", "start", "end", "chromatin_state"]
            # join various annotations per peak
            annot_comp = (
                annot.groupby(["chrom", "start", "end"])
                .aggregate(lambda x: ",".join(set([str(i) for i in x])))
                .reset_index()
            )
            annot_comp.index = bed_to_index(annot_comp)
            annot_comp.columns = ["chrom", "start", "end", "chromatin_state"]
            # save to disk
            if save:
                a = "" if (label == "real") else ("_" + label)
                if output_prefix != "":
                    output_prefix += "."
                annot.to_csv(
                    os.path.join(
                        self.results_dir,
                        self.name
                        + ".chrom_state_annotation{}_mapping.{}csv".format(a, output_prefix),
                    ),
                    index=True,
                )
                annot_comp.to_csv(
                    os.path.join(
                        self.results_dir,
                        self.name + ".chrom_state_annotation{}.{}csv".format(a, output_prefix),
                    ),
                    index=True,
                )

            if assign:
                setattr(self, attr, annot_comp)
                setattr(self, attr + "_mapping", annot)
        return self.chrom_state_annotation

    def get_sex_chrom_ratio(
        self,
        matrix="matrix_norm",
        sex_chroms=["chrX", "chrY"],
        output_dir="{results_dir}",
        output_prefix="sex_chrom_ratio",
        plot=True,
    ):
        """
        Get ratio of signal between sex chromosomes.
        Useful to quickly assign sex to samples.

        Parameters
        ----------
        matrix : :obj:`pandas.DataFrame`, optional
            Matrix to use.
            Defaults to `matrix_norm`.

        sex_chroms : :obj:`list`, optional
            Names of the two sex chromosomes to use.

        output_dir : :obj:`str`, optional
            Directory to write output to.

        output_prefix : :obj:`str`, optional
            String to prefix output with.

        plot: :obj:`bool`, optional
            Whether to produce illustrative plots.

        Returns
        ----------
        pd.Series
            Ratio of sex chromosomes defined as
            `sex_chroms[1] - sex_chroms[0]`.
        """
        import matplotlib.pyplot as plt
        import seaborn as sns
        from natsort import natsorted
        from ngs_toolkit.graphics import savefig

        matrix = self.get_matrix(matrix)

        output_dir = self._format_string_with_attributes(output_dir)

        # remove per sample mean
        matrix -= matrix.mean()

        chroms = matrix.index.str.extract(r"^(.*):\d+-\d+$").iloc[:, 0].values
        if not all([x in chroms for x in sex_chroms]):
            msg = f"Requested sex chromosomes {', '.join(sex_chroms)} not found in matrix."
            _LOGGER.error(msg)
            raise ValueError
        matrix = matrix.assign(chrom=chroms)
        m = matrix.groupby("chrom").mean()
        # remove per chromosome mean
        m = (m.T - m.mean(1)).T
        # order chromosomes
        m = m.reindex(natsorted(m.index))

        # calculate ratio
        ratio = m.loc[sex_chroms[1]] - m.loc[sex_chroms[0]]
        ratio.name = "{}_to_{}_ratio".format(sex_chroms[1], sex_chroms[0])
        ratio.to_csv(
            os.path.join(output_dir, self.name + "." + output_prefix + ".csv"), header=True
        )

        if plot:
            ratio.sort_values(inplace=True)
            m = m.reindex(ratio.index, axis=1)

            # Clustermap
            if isinstance(ratio.index, pd.MultiIndex):
                cols = m.columns.get_level_values("sample_name")
            else:
                cols = m.columns
            grid = sns.clustermap(
                m.T,
                z_score=1,
                center=0,
                cmap="RdBu_r",
                figsize=(m.shape[0] * 0.3, m.shape[1] * 0.3),
                row_cluster=False,
                col_cluster=False,
                cbar_kws={"label": "Deviation from mean\nchromosome accessibility"},
                yticklabels=cols,
            )
            grid.ax_heatmap.set_xlabel("Chromosomes")
            grid.ax_heatmap.set_ylabel("Samples")
            savefig(
                grid, os.path.join(output_dir, self.name + "." + output_prefix + ".clustermap.svg")
            )

            # Value vs Rank
            fig, axis = plt.subplots(1, figsize=(3, ratio.shape[0] * 0.3))
            axis.scatter(ratio, ratio.rank(), linestyle="-")
            axis.axvline(0, linestyle="--", color="grey")
            axis.set_yticks(range(1, ratio.shape[0] + 1))
            if isinstance(ratio.index, pd.MultiIndex):
                axis.set_yticklabels(ratio.index.get_level_values("sample_name"))
            else:
                axis.set_yticklabels(ratio.index)
            axis.set_ylabel("Samples")
            axis.set_xlabel(ratio.name.replace("_", " "))
            v = ratio.abs().max()
            v += v * 0.1
            axis.set_xlim((-v, v))
            sns.despine(fig)
            savefig(
                fig,
                os.path.join(output_dir, self.name + "." + output_prefix + ".rank_vs_ratio.svg"),
            )

        return ratio.sort_index()

    def get_gene_level_matrix(
        self,
        matrix="matrix_norm",
        reduce_func=np.mean,
        assign=True,
        save=True,
        output_file="{results_dir}/{name}.gene_coverage.csv",
    ):
        """
        Get gene-level measurements of coverage.

        Requires a 'gene_annotation' or 'closest_tss_distances' attribute to be set
        containing a mapping between the index of `matrix` and genes
        (produced from `get_peak_gene_annotation`).

        Parameters
        ----------
        matrix : :obj:`str`, optional
            Quantification matrix to use (e.g. 'matrix_raw' or 'matrix_norm')

            Default is "matrix_norm".
        reduce_func : func
            Function to apply to reduce values.

            Default is mean.
        assign: :obj:`bool`
            Whether to assign the matrix to an attribute of self named `matrix_gene`.

            Default is :obj:`True`.
        save: :obj:`bool`
            Whether to save to disk the coverage matrix with filename `output_file`.

            Default is :obj:`True`.
        output_file : :obj:`str`
            Path to save a CSV file with coverage output if `save` is `True`.

            Default is `self.results_dir/self.name + ".raw_coverage.csv"`.

        Returns
        ---------
        pandas.DataFrame
            Coverage values reduced per gene.

        Attributes
        ---------
        matrix_gene : :obj:`pandas.DataFrame`
            Coverage values reduced per gene.
        """
        msg = "Analysis object lacks a 'gene_annotation' or 'closest_tss_distances' dataframe."
        hint = " Call 'analysis.get_peak_gene_annotation' to have region-gene associations."
        if not (hasattr(self, "gene_annotation") or hasattr(self, "closest_tss_distances")):
            _LOGGER.error(msg + hint)
            raise AssertionError(msg)

        matrix = self.get_matrix(matrix).copy()

        if hasattr(self, "closest_tss_distances"):
            matrix2 = matrix.join(self.closest_tss_distances[["gene_name"]])
            matrix2 = matrix2.set_index("gene_name", append=True)
        else:
            g = self.gene_annotation["gene_name"].str.split(",").apply(pd.Series).stack()
            g.index = g.index.droplevel(1)
            g.name = "gene_name"
            matrix2 = matrix.join(g).drop("gene_name", axis=1)
            matrix2.index = matrix.join(g).reset_index().set_index(["index", "gene_name"]).index

        matrix2.columns = matrix.columns
        matrix3 = matrix2.groupby(level="gene_name").apply(reduce_func)
        matrix3 = matrix3.loc[:, ~matrix3.isnull().all()]
        if assign:
            self.matrix_gene = matrix3
        if save:
            matrix3.to_csv(self._format_string_with_attributes(output_file))
        return matrix3

    def get_gene_level_changes(self, differential_results=None, reduce_func=np.mean):
        """
        Redcuce changes in regulatory elements to gene-level by
        aggregating across regulatory elements.
        Requires a 'gene_annotation' attribute to be set containing
        a mapping between the index of `matrix` and genes (produced
        from `get_peak_gene_annotation`).

        Parameters
        ----------
        differential_results : :obj:`pandas.DataFrame`
            Matrix with differential results to use.
            Default is a 'differential_results' attribute of self.

        reduce_func : func
            Function to apply to reduce values. Default is mean

        Returns
        ---------
        pandas.DataFrame
            Changes in chromatin accessibility (log2FoldChanges)
            reduced per gene.
        """
        msg = "Analysis object lacks a 'gene_annotation' or 'closest_tss_distances' dataframe."
        hint = " Call 'analysis.get_peak_gene_annotation' to have region-gene associations."
        if not (hasattr(self, "gene_annotation") or hasattr(self, "closest_tss_distances")):
            _LOGGER.error(msg + hint)
            raise AssertionError(msg)

        if differential_results is None:
            differential_results = self.differential_results

        if hasattr(self, "closest_tss_distances"):
            dr2 = differential_results.join(self.closest_tss_distances[["gene_name"]])
            dr2 = dr2.set_index("gene_name", append=True)
        else:
            g = self.gene_annotation["gene_name"].str.split(",").apply(pd.Series).stack()
            g.index = g.index.droplevel(1)
            g.name = "gene_name"

            dr2 = differential_results.join(g).drop("gene_name", axis=1)
            dr2.index = (
                differential_results.join(g).reset_index().set_index(["index", "gene_name"]).index
            )

        dr2.columns = differential_results.columns
        dr3 = dr2.reset_index().groupby(["gene_name", "comparison_name"]).apply(reduce_func)

        return dr3.loc[:, ~dr3.isnull().all()]

    def plot_peak_characteristics(
        self,
        samples=None,
        by_attribute=None,
        genome_space=3e9,
        output_dir="{results_dir}/peak_characteristics",
        output_prefix="{name}",
    ):
        """
        Several diagnostic plots on the analysis' consensus peak set
        and the sample's signal on them.

        Provides plots with samples grouped `by_attribute` if given (a string or a list of strings).

        Parameters
        ----------
        samples : :obj:`list`, optional
            List of samples to restrict analysis to.

        by_attribute : {str, list}, optional
            Attribute or list of sample attributes to groupby samples by when plotting.
            This is done in addition to the plots with individual values per sample.

        genome_space : :obj:`int`
            Length of genome.

            Defaults to 3e9 basepairs (human genome).

        output_dir : :obj:`str`
            Directory to output files. Will be formated with variables from Analysis.

            Defaults to "peak_characteristics" under the Analysis "results_dir".
        output_prefix : :obj:`str`
            Prefix to add to output files.

            Defaults to the Analysis' name.
        """
        # TODO: abstract genome space kwarg
        import multiprocessing
        import parmap
        import matplotlib.pyplot as plt
        import seaborn as sns

        from ngs_toolkit.utils import (
            count_bam_file_length,
            count_lines,
            get_total_region_area,
            get_region_lengths,
            get_regions_per_chromosomes,
        )
        from ngs_toolkit.graphics import savefig

        if samples is None:
            samples = self.samples

        output_dir = self._format_string_with_attributes(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_prefix = self._format_string_with_attributes(output_prefix)

        reads = parmap.map(count_bam_file_length, [s.aligned_filtered_bam for s in samples])
        peaks = list(map(count_lines, [s.peaks for s in samples]))
        open_chrom = list(map(get_total_region_area, [s.peaks for s in samples]))

        stats = pd.DataFrame(
            [reads, peaks, open_chrom],
            index=["reads_used", "peak_number", "open_chromatin"],
            columns=[s.name for s in samples],
        ).T

        stats["peaks_norm"] = (stats["peak_number"] / stats["reads_used"]) * 1e3
        stats["open_chromatin_norm"] = stats["open_chromatin"] / stats["reads_used"]
        stats.to_csv(
            os.path.join(output_dir, "{}.open_chromatin_space.csv".format(output_prefix)),
            index=True,
        )
        # stats = pd.read_csv(os.path.join(
        #         output_dir,
        #         "{}.open_chromatin_space.csv"
        #         .format(output_prefix)),
        #     index_col=0)

        # median lengths per sample (split-apply-combine)
        if by_attribute is not None:
            stats = stats.join(self.get_sample_annotation())
            stats = pd.merge(
                stats,
                stats.groupby(by_attribute)["open_chromatin"]
                .median()
                .to_frame(name="group_open_chromatin")
                .reset_index(),
            ).set_index(stats.index)
            stats = pd.merge(
                stats,
                stats.groupby(by_attribute)["open_chromatin_norm"]
                .median()
                .to_frame(name="group_open_chromatin_norm")
                .reset_index(),
            ).set_index(stats.index)

        # plot
        stats = stats.sort_values("open_chromatin_norm")
        fig, axis = plt.subplots(2, 1, figsize=(1 * 3, 2 * 3))
        for ax, var in zip(axis, ["open_chromatin", "open_chromatin_norm"]):
            sns.barplot(
                y="index", x=var, orient="horiz", data=stats.reset_index(), palette="summer", ax=ax,
            )
            ax.set_ylabel("Sample name")
        axis[0].set_xlabel("Total open chromatin space (bp)")
        axis[1].set_xlabel("Total open chromatin space (normalized)")
        sns.despine(fig)
        savefig(
            fig,
            os.path.join(
                output_dir, "{}.total_open_chromatin_space.per_sample.svg".format(output_prefix),
            ),
        )

        if by_attribute is not None:
            fig, axis = plt.subplots(2, 1, figsize=(4 * 2, 6 * 1))
            stats = stats.sort_values("group_open_chromatin")
            sns.barplot(
                y=by_attribute,
                x="open_chromatin",
                orient="horiz",
                data=stats.reset_index(),
                palette="summer",
                ax=axis[0],
            )
            sns.stripplot(
                y=by_attribute,
                x="open_chromatin",
                orient="horiz",
                data=stats.reset_index(),
                palette="summer",
                ax=axis[0],
            )
            stats = stats.sort_values("group_open_chromatin_norm")
            sns.barplot(
                y=by_attribute,
                x="open_chromatin_norm",
                orient="horiz",
                data=stats.reset_index(),
                palette="summer",
                ax=axis[1],
            )
            sns.stripplot(
                y=by_attribute,
                x="open_chromatin_norm",
                orient="horiz",
                data=stats.reset_index(),
                palette="summer",
                ax=axis[1],
            )
            axis[0].set_xlabel("Total open chromatin space (bp)")
            axis[1].set_xlabel("Total open chromatin space (normalized)")
            sns.despine(fig)
            savefig(
                fig,
                os.path.join(
                    output_dir,
                    "{}.total_open_chromatin_space.per_{}.svg".format(self.name, by_attribute),
                ),
            )

        # plot distribution of peak lengths
        sample_peak_lengths = map(get_region_lengths, [s.peaks for s in samples])
        lengths = pd.melt(
            pd.DataFrame(sample_peak_lengths, index=[s.name for s in samples]).T,
            value_name="peak_length",
            var_name="sample_name",
        ).dropna()

        # median lengths per sample (split-apply-combine)
        lengths = pd.merge(
            lengths,
            lengths.groupby("sample_name")["peak_length"]
            .median()
            .to_frame(name="mean_peak_length")
            .reset_index()
            .sort_values("mean_peak_length"),
        )

        fig, axis = plt.subplots(2, 1, figsize=(3 * 1, 3 * 2), sharex=False)
        for ax in axis:
            sns.boxplot(
                y="sample_name",
                x="peak_length",
                orient="horiz",
                data=lengths,
                palette="summer",
                ax=ax,
                showfliers=False,
            )
            ax.set_ylabel("Sample name")
            ax.set_xlabel("Peak length (bp)")
        axis[1].set_xscale("log")
        sns.despine(fig)
        savefig(
            fig, os.path.join(output_dir, "{}.peak_lengths.per_sample.svg".format(output_prefix)),
        )

        if by_attribute is not None:
            lengths = lengths.merge(self.get_sample_annotation())
            lengths = pd.merge(
                lengths,
                lengths.groupby(by_attribute)["peak_length"]
                .median()
                .to_frame(name="group_mean_peak_length")
                .reset_index()
                .sort_values("group_mean_peak_length"),
            ).set_index(lengths.index)

            fig, axis = plt.subplots(2, 1, figsize=(3 * 1, 3 * 2), sharex=False)
            for ax in axis:
                sns.boxplot(
                    y=by_attribute,
                    x="peak_length",
                    orient="horiz",
                    data=lengths,
                    palette="summer",
                    ax=ax,
                    showfliers=False,
                )
                ax.set_xlabel("Peak length (bp)")
            axis[1].set_xscale("log")
            sns.despine(fig)
            savefig(
                fig,
                os.path.join(
                    output_dir, "{}.peak_lengths.per_{}.svg".format(output_prefix, by_attribute),
                ),
            )

        # peaks per chromosome per sample
        chroms = (
            pd.DataFrame(
                map(get_regions_per_chromosomes, [s.peaks for s in samples]),
                index=[s.name for s in samples],
            )
            .fillna(0)
            .T
        )
        chroms_norm = (chroms / chroms.sum(axis=0)) * 100

        if hasattr(self, "organism"):
            chroms_norm = chroms_norm.loc[~chroms_norm.index.str.contains("_"), :]

            fig, axis = plt.subplots(1, 1, figsize=(8 * 1, 8 * 1))
            sns.heatmap(
                chroms_norm,
                square=True,
                cmap="summer",
                xticklabels=True,
                yticklabels=True,
                cbar_kws={"label": "Normalized accessibility"},
                ax=axis,
            )
            axis.set_xlabel("Sample name")
            axis.set_ylabel("Chromosome")
            axis.set_xticklabels(axis.get_xticklabels(), rotation=90, ha="right")
            axis.set_yticklabels(axis.get_yticklabels(), rotation=0, ha="right")
            savefig(
                fig,
                os.path.join(output_dir, "{}.peak_location.per_sample.svg".format(output_prefix)),
            )

        # Peak set across samples:
        # interval lengths
        fig, axis = plt.subplots(1, 1, figsize=(3, 3))
        sns.distplot(
            [interval.length for interval in self.sites if interval.length < 2000],
            hist=False,
            kde=True,
            ax=axis,
        )
        axis.set_xlabel("Peak width (bp)")
        axis.set_ylabel("Density")
        sns.despine(fig)
        savefig(fig, os.path.join(output_dir, "{}.lengths.svg".format(output_prefix)))

        # plot support
        if hasattr(self, "support"):
            fig, axis = plt.subplots()
            sns.distplot(self.support["support"], bins=40, ax=axis)
            axis.set_ylabel("frequency")
            sns.despine(fig)
            savefig(fig, os.path.join(output_dir, "{}.support.svg".format(output_prefix)))

        # Plot distance to nearest TSS
        if hasattr(self, "closest_tss_distances"):
            fig, axis = plt.subplots(2, 1, figsize=(3 * 1, 3 * 2), sharex=False, sharey=False)
            for i, ax in enumerate(axis):
                sns.distplot(
                    self.closest_tss_distances["distance"],
                    bins=1000,
                    kde=True,
                    hist=False if (i % 2 == 0) else True,
                    ax=ax,
                )
                ax.set_xlabel("Distance to nearest TSS (bp)")
                ax.set_ylabel("Density")
            axis[1].set_yscale("log")
            sns.despine(fig)
            savefig(
                fig, os.path.join(output_dir, "{}.tss_distance.svg".format(output_prefix)),
            )

        # Plot genomic regions
        datas = list()
        for name, attr, attr_b in [
            ("genomic_region", "region_annotation_mapping", "region_annotation_b_mapping",),
            (
                "chromatin_state",
                "chrom_state_annotation_mapping",
                "chrom_state_annotation_b_mapping",
            ),
        ]:
            if hasattr(self, attr):
                f = getattr(self, attr)
                b = getattr(self, attr_b)
                # count region frequency
                data = f[name].value_counts().sort_values(ascending=False)
                background = b[name].value_counts().sort_values(ascending=False)
                data = data.to_frame(name="foreground").join(background.to_frame(name="background"))
                data["fold_change"] = np.log2(data["foreground"] / data["background"])
                data.index.name = "region"

                # plot also % of genome space "used"
                f["length"] = f["end"] - f["start"]
                b["length"] = b["end"] - b["start"]

                s = f.groupby(name)["length"].sum()
                data.loc[:, "percent_genome_space"] = (s / genome_space).dropna() * 100
                datas.append(data)

                # plot together
                g = sns.FacetGrid(
                    data=pd.melt(data.reset_index(), id_vars="region"),
                    col="variable",
                    col_wrap=2,
                    sharex=False,
                    sharey=True,
                    size=2,
                    aspect=1.2,
                )
                g.map(sns.barplot, "value", "region", orient="horiz")
                sns.despine(fig)
                savefig(
                    g, os.path.join(output_dir, "{}.{}s.svg".format(output_prefix, name)),
                )

        # plot together
        if len(datas) > 1:
            data = pd.concat(datas)
            g = sns.FacetGrid(
                data=pd.melt(data.reset_index(), id_vars="region").sort_values("value"),
                col="variable",
                col_wrap=2,
                sharex=False,
                sharey=True,
                size=3,
                aspect=1.2,
            )
            g.map(sns.barplot, "value", "region", orient="horiz")
            sns.despine(fig)
            savefig(
                g,
                os.path.join(
                    output_dir, "{}.genomic_region_and_chromatin_states.svg".format(output_prefix),
                ),
            )

        # distribution of count statistics
        if hasattr(self, "stats"):
            for attr in self.stats.columns:
                fig, axis = plt.subplots(1, 1, figsize=(3, 3))
                sns.distplot(self.stats.loc[:, attr], hist=False, kde=True, ax=axis)
                sns.despine(fig)
                savefig(
                    fig, os.path.join(output_dir, "{}.{}.distplot.svg".format(output_prefix, attr)),
                )
        if hasattr(self, "support"):
            attr = "support"
            fig, axis = plt.subplots(1, 1, figsize=(3, 3))
            sns.distplot(self.support.loc[:, attr], hist=False, kde=True, ax=axis)
            sns.despine(fig)
            savefig(
                fig, os.path.join(output_dir, "{}.{}.distplot.svg".format(output_prefix, attr)),
            )

        # Pairwise against mean
        if hasattr(self, "stats"):
            stats = self.stats.copy()
            if hasattr(self, "support"):
                stats = stats.join(self.support.loc[:, "support"])
            for attr in stats.columns:
                if attr == "mean":
                    continue
                p = stats[(stats["mean"] > 0) & (stats[attr] < np.percentile(stats[attr], 99) * 3)]
                g = sns.jointplot(p["mean"], p[attr], s=1, alpha=0.1, rasterized=True, height=3)
                savefig(
                    g.fig,
                    os.path.join(output_dir, "{}.mean_vs_{}.svg".format(output_prefix, attr)),
                )

    def plot_raw_coverage(self, samples=None, by_attribute=None):
        """
        Diagnostic plots on the Sample's signal.
        Provides plots with Samples grouped `by_attribute` if given (a string or a list of strings).

        Parameters
        ----------
        samples : :obj:`list`
            List of peppy.Samples objects to use for plotting.

        by_attribute : :obj:`str`, optional
            Attribute of samples to group by.
            Values will be aggregated across samples by that attribute.
        """
        # TODO: get matrix as input, move to graphics
        import matplotlib.pyplot as plt
        import seaborn as sns
        from ngs_toolkit.graphics import savefig

        if samples is None:
            samples = self.samples

        if by_attribute is None:
            cov = pd.melt(
                np.log2(1 + self.matrix_raw[[s.name for s in samples]]),
                var_name="Sample name",
                value_name="Raw counts (log2)",
            )
            fig, axis = plt.subplots(1, 1, figsize=(6, 1 * 4))
            sns.violinplot(
                "Raw counts (log2)",
                "Sample name",
                orient="horizontal",
                palette="tab20",
                data=cov,
                ax=axis,
            )
            sns.despine(fig)
            savefig(fig, os.path.join(self.results_dir, self.name + ".raw_counts.violinplot.svg"))
        else:
            attrs = set([getattr(s, by_attribute) for s in samples])
            fig, axis = plt.subplots(len(attrs), 1, figsize=(8, len(attrs) * 6))
            for i, attr in enumerate(attrs):
                _LOGGER.info(attr)
                cov = pd.melt(
                    np.log2(
                        1
                        + self.matrix_raw[
                            [s.name for s in samples if getattr(s, by_attribute) == attr]
                        ]
                    ),
                    var_name="Sample name",
                    value_name="Raw counts (log2)",
                )
                sns.violinplot(
                    "Raw counts (log2)",
                    "Sample name",
                    orient="horizontal",
                    palette="tab20",
                    data=cov,
                    ax=axis[i],
                )
                axis[i].set_title(attr)
                axis[i].set_xticklabels(axis[i].get_xticklabels(), rotation=90)
            sns.despine(fig)
            savefig(
                fig,
                os.path.join(
                    self.results_dir,
                    self.name + ".raw_counts.violinplot.by_{}.svg".format(by_attribute),
                ),
            )

    def plot_coverage(self):
        import matplotlib.pyplot as plt
        import seaborn as sns

        # TODO: add plots for overal genome
        # TODO: add raw counts too

        data = self.matrix_norm.copy()
        # (rewrite to avoid putting them there in the first place)
        variables = ["gene_name", "genomic_region", "chromatin_state"]

        for variable in variables:
            d = (
                data[variable].str.split(",").apply(pd.Series).stack()
            )  # separate comma-delimited fields
            d.index = d.index.droplevel(
                1
            )  # returned a multiindex Series, so get rid of second index level (first is from original row)
            data = data.drop([variable], axis=1)  # drop original column so there are no conflicts
            d.name = variable
            data = data.join(d)  # joins on index

        variables = [
            "chrom",
            "start",
            "end",
            "ensembl_transcript_id",
            "distance",
            "ensembl_gene_id",
            "support",
            "mean",
            "variance",
            "std_deviation",
            "dispersion",
            "qv2",
            "amplitude",
            "gene_name",
            "genomic_region",
            "chromatin_state",
        ]
        # Plot
        data_melted = pd.melt(data, id_vars=variables, var_name="sample", value_name="norm_counts")

        # transform dispersion
        data_melted["dispersion"] = np.log2(1 + data_melted["dispersion"])

        # Together in same violin plot
        fig, axis = plt.subplots(1)
        sns.violinplot("genomic_region", "norm_counts", data=data_melted, ax=axis)
        fig.savefig(
            os.path.join(
                self.results_dir, self.name + ".norm_counts.per_genomic_region.violinplot.svg",
            ),
            bbox_inches="tight",
        )

        # dispersion
        fig, axis = plt.subplots(1)
        sns.violinplot("genomic_region", "dispersion", data=data_melted, ax=axis)
        fig.savefig(
            os.path.join(
                self.results_dir,
                self.name + ".norm_counts.dispersion.per_genomic_region.violinplot.svg",
            ),
            bbox_inches="tight",
        )

        # dispersion
        fig, axis = plt.subplots(1)
        sns.violinplot("genomic_region", "qv2", data=data_melted, ax=axis)
        fig.savefig(
            os.path.join(
                self.results_dir, self.name + ".norm_counts.qv2.per_genomic_region.violinplot.svg",
            ),
            bbox_inches="tight",
        )

        fig, axis = plt.subplots(1)
        sns.violinplot("chromatin_state", "norm_counts", data=data_melted, ax=axis)
        fig.savefig(
            os.path.join(
                self.results_dir, self.name + ".norm_counts.chromatin_state.violinplot.svg",
            ),
            bbox_inches="tight",
        )

        fig, axis = plt.subplots(1)
        sns.violinplot("chromatin_state", "dispersion", data=data_melted, ax=axis)
        fig.savefig(
            os.path.join(
                self.results_dir,
                self.name + ".norm_counts.dispersion.chromatin_state.violinplot.svg",
            ),
            bbox_inches="tight",
        )

        fig, axis = plt.subplots(1)
        sns.violinplot("chromatin_state", "qv2", data=data_melted, ax=axis)
        fig.savefig(
            os.path.join(
                self.results_dir, self.name + ".norm_counts.qv2.chromatin_state.violinplot.svg",
            ),
            bbox_inches="tight",
        )

        # separated by variable in one grid
        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "mean", hist=False, rug=False)
        g.fig.savefig(
            os.path.join(
                self.results_dir, self.name + ".norm_counts.mean.per_genomic_region.distplot.svg",
            ),
            bbox_inches="tight",
        )

        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "dispersion", hist=False, rug=False)
        g.fig.savefig(
            os.path.join(
                self.results_dir,
                self.name + ".norm_counts.dispersion.per_genomic_region.distplot.svg",
            ),
            bbox_inches="tight",
        )

        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "qv2", hist=False, rug=False)
        g.fig.savefig(
            os.path.join(
                self.results_dir, self.name + ".norm_counts.qv2.per_genomic_region.distplot.svg",
            ),
            bbox_inches="tight",
        )

        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "support", hist=False, rug=False)
        g.fig.savefig(
            os.path.join(
                self.results_dir,
                self.name + ".norm_counts.support.per_genomic_region.distplot.svg",
            ),
            bbox_inches="tight",
        )

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "mean", hist=False, rug=False)
        g.fig.savefig(
            os.path.join(
                self.results_dir, self.name + ".norm_counts.mean.chromatin_state.distplot.svg",
            ),
            bbox_inches="tight",
        )

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "dispersion", hist=False, rug=False)
        g.fig.savefig(
            os.path.join(
                self.results_dir,
                self.name + ".norm_counts.dispersion.chromatin_state.distplot.svg",
            ),
            bbox_inches="tight",
        )

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "qv2", hist=False, rug=False)
        g.fig.savefig(
            os.path.join(
                self.results_dir, self.name + ".norm_counts.qv2.chromatin_state.distplot.svg",
            ),
            bbox_inches="tight",
        )

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "support", hist=False, rug=False)
        g.fig.savefig(
            os.path.join(
                self.results_dir, self.name + ".norm_counts.support.chromatin_state.distplot.svg",
            ),
            bbox_inches="tight",
        )
        plt.close("all")

    def region_context_enrichment(
        self,
        regions,
        steps=["genomic_region", "chromatin_state"],
        background="region_set",
        prefix="region_type_enrichment",
        output_dir="{results_dir}",
    ):
        """
        Characterize a subset of the regions (e.g. differential regions)
        in terms of their genomic context.

        Parameters
        ----------
        regions : {:obj:`list`, :obj:`pandas.DataFrame`, :obj:`pandas.Index`}
            Subset of regions of interest to analysis.
            Must be a subset of the universe (i.e. ``sites`` attribute).

        steps : :obj:`list`, optional
            Steps of enrichment to perform.
            Defaults to all available: 'genomic_region' and 'chromatin_state'.

        background : :obj:`str`, optional
            Which set to consider as backgroud.
            Options are:

                * "region_set": the consensus region_set of the analysis
                * "genome": a randomized set of size as region_set across the genome

        prefix : :obj:`str`, optional
            Prefix for saved files.

            Default is "region_type_enrichment".
        output_dir : :obj:`str`, optional
            Directory to write results to.
            Can be formatted with Analysis attributes.

            Default is "{results_dir}".

        Returns
        ---------
        pandas.DataFrame
            Enrichment results
        """
        from ngs_toolkit.utils import log_pvalues
        from scipy.stats import fisher_exact

        output_dir = self._format_string_with_attributes(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        if isinstance(regions, pd.DataFrame):
            regions = regions.index.tolist()
        elif isinstance(regions, pd.Index):
            regions = regions.tolist()
        elif isinstance(regions, list):
            pass
        else:
            msg = "Input not understood."
            _LOGGER.error(msg)
            raise ValueError(msg)

        options = ["region_set", "genome"]
        if background not in options:
            msg = "Option `background` must be one of '{}'.".format("', '".join(options))
            raise ValueError(msg)

        # compare genomic regions and chromatin_states
        enr = list()
        msg = "'{}' step selected, but analysis does not have '{}'."
        msg2 = "'genome' selected, but analysis does not have '{}'."
        for step, matrix, matrix_b in [
            ("genomic_region", "region_annotation_mapping", "region_annotation_b_mapping",),
            (
                "chromatin_state",
                "chrom_state_annotation_mapping",
                "chrom_state_annotation_b_mapping",
            ),
        ]:
            if step not in steps:
                continue
            if not hasattr(self, matrix):
                _LOGGER.warning(msg.format(step, matrix))
                continue

            _LOGGER.debug("Getting enrichment of regions in '{}'.".format(step))

            # Count foreground occurences
            annot = getattr(self, matrix)
            res = annot.loc[set(regions), step].value_counts().to_frame("foreground")

            # Count background occurences
            # # in case of background == 'genome', we simply replace the dataframes
            if background == "genome":
                try:
                    annot = getattr(self, matrix_b)
                except AttributeError:
                    _LOGGER.warning(msg2.format(matrix_b))
                    continue
            # # account for foreground regions not annotated
            res_b = annot.loc[:, step].value_counts().to_frame(name="universe")
            if not all([x in res_b.index for x in res.index]):
                m = res.index[~res.index.isin(res_b.index)]
                msg3 = "Foreground regions contains type of {} not in background: {}".format(
                    step, "', '".join(m)
                )
                msg3 += " Continuing without those."
                _LOGGER.warning(msg3)
            res = res.reindex(res_b.index)

            # # join
            res = res.join(res_b, how="outer").fillna(0).astype(int)

            # Calculate log fold enrichment:
            # # normalize to total:
            res.loc[:, "foreground_fraction"] = res["foreground"] / res["foreground"].sum()
            res.loc[:, "universe_fraction"] = res["universe"] / res["universe"].sum()
            res.loc[:, "log2_fold_change"] = np.log2(
                res["foreground_fraction"] / res["universe_fraction"]
            )
            # Calculate overlap p-value:
            for feature in res["foreground"].index:
                a = res.loc[feature, "foreground"]
                b = res.loc[:, "foreground"].drop(feature).sum()
                c = annot.loc[(~annot.index.isin(regions)), step].value_counts()[feature]
                d = annot.loc[(~annot.index.isin(regions)), step].value_counts().drop(feature).sum()
                res.loc[feature, "odds_ratio"], res.loc[feature, "p_value"] = fisher_exact(
                    [[a, c], [b, d]], alternative="two-sided"
                )
            res.loc[:, "log2_odds_ratio"] = np.log2(res["odds_ratio"])
            res.loc[:, "-log10(p-value)"] = log_pvalues(res["p_value"])
            res.loc[:, "region_type"] = step
            # Append
            enr.append(res)

        # save
        enr = pd.concat(enr)
        enr.index.name = "region"
        enr.to_csv(os.path.join(output_dir, prefix + ".csv"), index=True)
        return enr

    def characterize_regions_function(
        self,
        differential,
        output_dir,
        prefix,
        universe_file=None,
        run=True,
        genome=None,
        steps=["region", "lola", "meme", "homer", "enrichr"],
    ):
        """
        Performs a range of functional enrichments of a set of regions
        given in ``differential`` (a dataframe which is typically a subset
        of an annotated coverage dataframe). Will extract regions,
        their underlying sequence, associated genes, perform enrichment of
        genomic regions, chromatin states against a background, motif enrichment,
        location overlap analysis (LOLA), and gene set enrichment (using the Enrichr API).

        This requires several programs and R libraries:
            - MEME suite (AME)
            - HOMER suite (findMotifsGenome.pl)
            - LOLA (R library)

        Additionally, some genome-specific databases are needed
        to run these programs.

        Parameters
        ----------
        differential : :obj:`pandas.DataFrame`
            Results of differential analysis for a given comparison of interest.

        output_dir : :obj:`str`
            Directory to output results to.

        prefix : :obj:`str`
            Prefix to use for output files.

        universe_file : :obj:`str`, optional
            Path to BED file with set of universe regions where
            differential were selected from.

            Default is ``sites`` attribute of Analysis.
        run: :obj:`bool`, optional
            Whether to run enrichment commands now or to simply
            prepare the input files for it.
            Default is :obj:`True`.

        genome : :obj:`str`, optional
            Genome assembly of analysis.
            Default is genome genome assembly of analysis (``genome`` attribute).

        steps : :obj:`list`, optional
            Which steps of the analysis to perform.
            Default is all: ['region', 'lola', 'meme', 'homer', 'enrichr'].
        """
        from ngs_toolkit.general import meme_ame, homer_motifs, lola, enrichr
        from ngs_toolkit.utils import (
            bed_to_fasta,
            standard_score,
            location_index_to_bed,
            get_this_file_or_timestamped,
        )

        # use all sites as universe
        if universe_file is None:
            try:
                universe_file = getattr(self, "sites").fn
                _LOGGER.info(
                    "Using default background region set from 'analysis.sites': {}'.".format(
                        universe_file
                    )
                )
            except AttributeError as e:
                _LOGGER.error("Background region set 'analysis.sites' is not set! Cannot run LOLA!")
                raise e

        if genome is None:
            genome = self.genome

        # make output dirs
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        bed = location_index_to_bed(differential.index)
        # save to bed
        bed_file = os.path.join(output_dir, "{}_regions.bed".format(prefix))
        bed.to_csv(bed_file, sep="\t", header=False, index=False)
        # save as tsv
        tsv_file = os.path.join(output_dir, "{}_regions.tsv".format(prefix))
        bed.reset_index().to_csv(tsv_file, sep="\t", header=False, index=False)

        # export gene names
        clean_gene = (
            differential["gene_name"].str.split(",").apply(pd.Series, 1).stack().drop_duplicates()
        )
        clean_gene = clean_gene[~clean_gene.isin([".", "nan", ""])]
        clean_gene.to_csv(
            os.path.join(output_dir, "{}.gene_symbols.txt".format(prefix)),
            header=False,
            index=False,
        )
        if "ensembl_gene_id" in differential.columns:
            # export ensembl gene names
            clean = (
                differential["ensembl_gene_id"]
                .str.split(",")
                .apply(pd.Series, 1)
                .stack()
                .drop_duplicates()
            )
            clean.to_csv(
                os.path.join(output_dir, "{}_genes.ensembl.txt".format(prefix)),
                header=False,
                index=False,
            )

        # export gene symbols with scaled absolute fold change
        if "log2FoldChange" in differential.columns:
            differential["score"] = standard_score(abs(differential["log2FoldChange"]))
            differential["abs_fc"] = abs(differential["log2FoldChange"])

            d = differential[["gene_name", "score"]].sort_values("score", ascending=False)

            # split gene names from score if a reg.element was assigned to more than one gene
            a = d["gene_name"].str.split(",").apply(pd.Series, 1).stack()
            a.index = a.index.droplevel(1)
            a.name = "gene_name"
            d = d[["score"]].join(a)
            # reduce various ranks to mean per gene
            d = d.groupby("gene_name").mean().reset_index()
            d.to_csv(
                os.path.join(output_dir, "{}.gene_symbols.score.csv".format(prefix)), index=False,
            )

        # get fasta file with sequence underlying region
        if ("meme" in steps) or ("homer" in steps):
            hint = " Will not do motif enrichment analysis."
            fasta_file = os.path.join(output_dir, "{}_regions.fa".format(prefix))

            resources = self.get_resources(steps=["genome"])["genome_file"]
            if "fasta" not in resources:
                reason = "Could not get genome sequence file in either FASTA or 2bit format."
                _LOGGER.warning(reason + hint)
            else:
                try:
                    bed_to_fasta(
                        input_bed=bed_file, output_fasta=fasta_file, genome_file=resources["fasta"],
                    )
                except EnvironmentError:
                    reason = "Could not get FASTA sequence for regions."
                    _LOGGER.warning(reason + hint)

        if not run:
            return

        if "region" in steps:
            self.region_context_enrichment(differential, output_dir=output_dir)

        # MEME
        if "meme" in steps:
            _LOGGER.info("Running MEME-AME for '{}'".format(prefix))
            omap = {"hg38": "human", "hg19": "human", "mm10": "mouse"}
            meme_ame(fasta_file, output_dir, organism=omap[genome])

        # HOMER
        if "homer" in steps:
            _LOGGER.info("Running HOMER for '{}'".format(prefix))
            homer_motifs(bed_file, output_dir, genome_assembly=genome)

        # LOLA
        if "lola" in steps:
            _LOGGER.info("Running LOLA for '{}'".format(prefix))
            try:
                lola(
                    get_this_file_or_timestamped(bed_file),
                    get_this_file_or_timestamped(universe_file),
                    output_dir,
                    genome=genome,
                )
            except:
                _LOGGER.error("LOLA analysis for '{}' failed!".format(prefix))

        # Enrichr
        if "enrichr" in steps:
            _LOGGER.info("Running Enrichr for '{}'".format(prefix))
            results = enrichr(clean_gene.to_frame(name="gene_name"))
            results.to_csv(
                os.path.join(output_dir, "{}.enrichr.csv".format(prefix)),
                index=False,
                encoding="utf-8",
            )
