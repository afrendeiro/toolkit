#!/usr/bin/env python


import multiprocessing
import os

import matplotlib.pyplot as plt
from ngs_toolkit import _LOGGER
from ngs_toolkit.analysis import Analysis
from ngs_toolkit.decorators import (
    check_organism_genome, check_has_sites)
from ngs_toolkit.general import (
    get_blacklist_annotations,
    meme_ame, homer_motifs, lola, enrichr)
from ngs_toolkit.utils import (
    bed_to_index, bed_to_fasta,
    log_pvalues, standard_score,
    count_reads_in_intervals,
    normalize_quantiles_r,
    normalize_quantiles_p)
import numpy as np
import pandas as pd
import parmap
import pybedtools
import pysam
from rpy2 import robjects
from rpy2.rinterface import RRuntimeWarning
from scipy.stats import fisher_exact
import seaborn as sns
from tqdm import tqdm


class ATACSeqAnalysis(Analysis):
    """
    Class to model analysis of ATAC-seq data.
    Inherits from the `ngs_toolkit.general.Analysis` class.

    Parameters
    ----------
    name : str
        Name to give analysis object.
        Default is ``atacseq_analysis``.

    samples : list
        Iterable of peppy.Sample objects use in analysis.
        If not provided (`None` is passed) and `prj` is.
        Defaults to all samples in the `prj` object (`samples` attribute).

    prj : peppy.Project
        Project to tie analysis to.

    data_dir : str
        Directory containing relevant data for analysis.
        Default is `data`.

    results_dir : str
        Directory to output relevant analysis results.
        Default is `results`.

    pickle_file : str
        File path to use to save serialized object in `pickle` format.

    from_pickle : bool
        If the analysis should be loaded from an existing pickle object.
        Default is `False.

    Remaining keyword arguments will be passed to parent class `ngs_toolkit.general.Analysis`.

    :Example:

    .. code-block:: python
        from peppy import Project
        import os
        from ngs_toolkit.atacseq import ATACSeqAnalysis

        prj = Project(os.path.join("metadata", "project_config.yaml"))
        atac_analysis = ATACSeqAnalysis(
            name=prj.project_name, prj=prj,
            samples=[s for s in prj.samples if s.protocol == "ATAC-seq"])

        # Get consensus peak set from all samples
        atac_analysis.get_consensus_sites(atac_analysis.samples)

        # Annotate regions
        atac_analysis.calculate_peak_support(atac_analysis.samples)
        atac_analysis.get_peak_gene_annotation()
        atac_analysis.get_peak_genomic_location()
        atac_analysis.get_peak_chromatin_state(
            os.path.join(atac_analysis.data_dir, "external", "E032_15_coreMarks_mnemonics.bed"))

        # Get coverage values for each peak in each sample of ATAC-seq
        atac_analysis.measure_coverage(atac_analysis.samples)

        # Normalize jointly (quantile normalization + GC correction)
        atac_analysis.normalize(method="gc_content")

        # Annotate quantified peaks with previously calculated metrics and features
        atac_analysis.annotate()

        # Annotate with sample metadata
        atac_analysis.accessibility = atac_analysis.annotate_with_sample_metadata(
            quant_matrix="coverage_annotated",
            attributes=atac_analysis.sample_variables)

        # Save object
        atac_analysis.to_pickle()

    """
    def __init__(
            self,
            name="atacseq_analysis",
            samples=None,
            prj=None,
            data_dir="data",
            results_dir="results",
            pickle_file=None,
            from_pickle=False,
            pep=False,
            **kwargs):
        super(ATACSeqAnalysis, self).__init__(
            name=name,
            data_dir=data_dir,
            results_dir=results_dir,
            pickle_file=pickle_file,
            samples=samples,
            prj=prj,
            from_pickle=from_pickle,
            pep=pep,
            **kwargs)

        self.data_type = self.__data_type__ = "ATAC-seq"
        self.var_names = "region"
        self.quantity = "accessibility"
        self.norm_units = "RPM"
        self.raw_matrix_name = "coverage"
        self.norm_matrix_name = "coverage_qnorm"  # or coverage_gc_corrected
        # TODO: have normalization method reset the above value, with info to user
        # or just overwrite
        self.annot_matrix_name = "accessibility"

    def load_data(self, output_mapping=None, only_these_keys=None, permissive=True, n_header_vars=None):
        """
        Load the output files of the major functions of the Analysis.

        Parameters
        ----------
        output_mapping : dict
            Dictionary with "attribute name": "path prefix" to load the files.

        only_these_keys : list, optional
            Iterable of analysis attributes to load up.
            Possible attributes:
                "sites", "support", "coverage", "coverage_qnorm",
                "nuc", "coverage_gc_corrected", "gene_annotation",
                "region_annotation", "region_annotation_b",
                "chrom_state_annotation", "chrom_state_annotation_b",
                "coverage_annotated", "accessibility", "differential_results".

        bool : permissive
            Whether an error should be thrown if reading a file causes IOError.

        Attributes
        ----------
        sites : pybedtools.BedTool
            Sets a `sites` variable.

        pandas.DataFrame
            Dataframes holding the respective data, available as attributes described
            in the `only_these_keys` parameter.

        Raises
        ----------
        IOError
            If not permissive and a file is not found
        """
        if only_these_keys is None:
            only_these_keys = [
                "sites", "support", "coverage", "coverage_qnorm",
                "nuc", "coverage_gc_corrected", "closest_tss_distances",
                "gene_annotation",
                "region_annotation", "region_annotation_b",
                "region_annotation_mapping", "region_annotation_b_mapping",
                "chrom_state_annotation", "chrom_state_annotation_b",
                "chrom_state_annotation_mapping", "chrom_state_annotation_b_mapping",
                "stats",
                "coverage_annotated", "accessibility",
                "differential_results"]

        # Figure out how many levels of MultiIndex does the 'accessibility' dataframe has
        if n_header_vars is None:
            if hasattr(self, "sample_attributes"):
                n_header_vars = len(self.sample_attributes)
            else:
                msg = "`n_header_vars` was not given and analysis does not have `sample_attributes`."
                msg += " Cannot load `accessibility` attribute."
                hint = " Other attributes loaded though."
                if permissive:
                    _LOGGER.warn(msg + hint)
                    try:
                        only_these_keys.pop(only_these_keys.index("accessibility"))
                    except ValueError:
                        pass
                else:
                    _LOGGER.error(msg + hint)
                    raise ValueError(msg)

        prefix = os.path.join(self.results_dir, self.name)
        kwargs = {"index_col": 0}
        if output_mapping is None:
            # TODO: get default mapping by having functions declare what they output
            # perhaps also with a dict of kwargs to pass to pandas.read_csv
            output_mapping = {
                "support":
                    (prefix + "_peaks.support.csv", kwargs),
                "coverage":
                    (prefix + "_peaks.raw_coverage.csv", kwargs),
                "coverage_qnorm":
                    (prefix + "_peaks.coverage_qnorm.csv", kwargs),
                "nuc":
                    (prefix + "_peaks.gccontent_length.csv", kwargs),
                "coverage_gc_corrected":
                    (prefix + "_peaks.coverage_gc_corrected.csv", kwargs),
                "gene_annotation":
                    (prefix + "_peaks.gene_annotation.csv", kwargs),
                "closest_tss_distances":
                    (prefix + "_peaks.closest_tss_distances.csv", kwargs),
                "region_annotation":
                    (prefix + "_peaks.region_annotation.csv", kwargs),
                "region_annotation_b":
                    (prefix + "_peaks.region_annotation_background.csv", kwargs),
                "region_annotation_mapping":
                    (prefix + "_peaks.region_annotation_mapping.csv", kwargs),
                "region_annotation_b_mapping":
                    (prefix + "_peaks.region_annotation_background_mapping.csv", kwargs),
                "chrom_state_annotation":
                    (prefix + "_peaks.chrom_state_annotation.csv", kwargs),
                "chrom_state_annotation_b":
                    (prefix + "_peaks.chrom_state_annotation_background.csv", kwargs),
                "chrom_state_annotation_mapping":
                    (prefix + "_peaks.chrom_state_annotation_mapping.csv", kwargs),
                "chrom_state_annotation_b_mapping":
                    (prefix + "_peaks.chrom_state_annotation_background_mapping.csv", kwargs),
                "stats":
                    (prefix + "_peaks.stats_per_region.csv", kwargs),
                "coverage_annotated":
                    (prefix + "_peaks.coverage_qnorm.annotated.csv", kwargs),
                "accessibility":
                    (prefix + ".accessibility.annotated_metadata.csv",
                        {"index_col": 0, "header": list(range(n_header_vars))}),
                "differential_results":
                    (os.path.join(self.results_dir, "differential_analysis_{}".format(self.data_type),
                                  "differential_analysis.deseq_result.all_comparisons.csv"), kwargs)}

        output_mapping = {k: v for k, v in output_mapping.items() if k in only_these_keys}

        for name, (file, kwargs) in output_mapping.items():
            _LOGGER.info("Loading '{}' analysis attribute.".format(name))
            try:
                setattr(self, name, pd.read_csv(file, **kwargs))
            except IOError as e:
                if not permissive:
                    raise e
                else:
                    _LOGGER.warning(e)

        # Special cases
        if "sites" in only_these_keys:
            file = os.path.join(self.results_dir, self.name + "_peak_set.bed")
            try:
                setattr(self, "sites", pybedtools.BedTool(file))
            except IOError as e:
                if not permissive:
                    raise e
                else:
                    _LOGGER.warning(e)

    @staticmethod
    def check_region_index(matrix):
        return (
            matrix.index.str.contains(":").all()
            and
            matrix.index.str.contains("-").all())

    @staticmethod
    def set_region_index(matrix):
        if ATACSeqAnalysis.check_region_index(matrix):
            _LOGGER.warning("Matrix already has well-formatted index.")
            return matrix
        else:
            req = ['chrom', 'start', 'end']
            if all([x in matrix.columns for x in req]):
                matrix.index = (
                    matrix['chrom'].astype(str) +
                    ":" +
                    matrix['start'].astype(str) +
                    "-" +
                    matrix['end'].astype(str))
            else:
                raise ValueError()

    def get_consensus_sites(
            self, samples=None, region_type="summits", extension=250,
            blacklist_bed=None,
            filter_mito_chr=True,
            permissive=False,
            **kwargs):
        """
        Get consensus (union) of enriched sites (peaks) across samples.
        There are two modes possible, defined by the value of ``region_type``:
         - peaks: simple union of all sites
         - summits: peak summits are extended by ``extension`` and a union is made,

        Parameters
        ----------
        samples : list
            Iterable of peppy.Sample objects to restrict to.
            Must have a `peaks` attribute set.
            Defaults to all samples in the analysis (`samples` attribute).

        region_type : str
            The type of region to use to create the consensus region set - one of `summits` or `peaks`.
            If `summits`, peak summits will be extended by `extension` before union.
            Otherwise sample peaks will be used with no modification.

        extension : int
            Amount to extend peaks summits by in both directions.

        blacklist_bed : str
            A (3 column) BED file with genomic positions to exclude from consensus peak set.

        filter_mito_chr : bool
            Whether to exclude 'chrM' from peak set.

        permissive : bool
            Whether Samples that which `region_type` attribute file does not exist
            should be simply skipped or an error thrown.

        **kwargs
            Not used. Provided for compatibility with ChIPSeqAnalysis class.

        Raises
        ----------
        IOError
            If not `permissive` and either the `peaks` or `summits` file
            of a sample is not readable, or if `permissive`
            but none of the samples has an existing file.

        Attributes
        ----------
        sites : pybedtools.BedTool
            Sets a `sites` variable with consensus peak set.
        """
        if region_type not in ['summits', 'peaks']:
            msg = "`region_type` attribute must be one of 'summits' or 'peaks'!"
            _LOGGER.error(msg)
            raise ValueError(msg)

        if samples is None:
            samples = self.samples

        # Check which samples to run (dependent on permissive)
        samples = self._get_samples_with_input_file(region_type, permissive=permissive, samples=samples)

        if blacklist_bed is None:
            _LOGGER.info("Blacklist file not provided. Downloading...")
            try:
                blacklist_bed = get_blacklist_annotations(self.organism, self.genome)
            except AttributeError:
                msg = "Blacklist file was not provided and cannot"
                msg += " get one without analysis having `organism` and `genome` set."
                _LOGGER.error(msg)
                raise AttributeError(msg)

        for i, sample in tqdm(enumerate(samples), total=len(samples), desc="Sample"):
            # print(sample.name)
            if region_type == "summits":
                try:
                    peaks = pybedtools.BedTool(sample.summits).slop(b=extension, genome=sample.genome)
                except ValueError as e:
                    if permissive:
                        _LOGGER.warning(
                            "Summits for sample {} ({}) not found!".format(sample, sample.summits))
                        continue
                    else:
                        raise e
            else:
                try:
                    peaks = pybedtools.BedTool(sample.peaks)
                except ValueError as e:
                    if permissive:
                        _LOGGER.warning(
                            "Peaks for sample {} ({}) not found!".format(sample, sample.peaks))
                        continue
                    else:
                        raise e
            # Merge overlaping peaks within a sample
            peaks = peaks.merge()
            if i == 0:
                sites = peaks
            else:
                # Concatenate all peaks
                sites = sites.cat(peaks)

        if "sites" not in locals():
            msg = "Couldn't read peak file for any sample."
            _LOGGER.error(msg)
            raise ValueError(msg)
        # Merge overlaping peaks across samples
        sites = sites.merge()

        # Filter
        # # remove blacklist regions
        if blacklist_bed is not False:
            if not isinstance(blacklist_bed, pybedtools.BedTool):
                blacklist = pybedtools.BedTool(blacklist_bed)
                sites = sites.intersect(v=True, b=blacklist)
        # # remove chrM peaks
        if filter_mito_chr:
            sites = sites.filter(lambda x: x.chrom != 'chrM')

        # Save
        sites.saveas(os.path.join(self.results_dir, self.name + "_peak_set.bed"))

        # Read up again
        self.sites = pybedtools.BedTool(os.path.join(self.results_dir, self.name + "_peak_set.bed"))

    def set_consensus_sites(self, bed_file, overwrite=True):
        """
        Set consensus (union) sites across samples given a BED file.

        Parameters
        ----------
        bed_file : str
            BED file to use as consensus sites.

        overwrite : book
            Whether a possibly existing file with a consensus peak set
            for this analysis should be overwritten in disk.

        Attributes
        ----------
        sites : pybedtools.BedTool
            Sets a `sites` variable with consensus peak set.
        """
        self.sites = pybedtools.BedTool(bed_file)
        if overwrite:
            default_sites = os.path.join(self.results_dir, self.name + "_peak_set.bed")
            # pybedtools will pipe to the input file!
            if bed_file == default_sites:
                self.sites.saveas(default_sites + ".new")
                os.rename(default_sites + ".new", default_sites)
            else:
                self.sites.saveas(default_sites)
        # TODO: warn if not overwrite and file exists already

    @check_has_sites
    def calculate_peak_support(
            self,
            samples=None,
            region_type="summits",
            permissive=False,
            **kwargs):
        """
        Count number of called peaks per sample in the consensus region set.
        In addition calculate a measure of peak support (or ubiquitouness) by
        observing the ratio of samples containing a peak overlapping each region.

        Parameters
        ----------
        samples : list
            Iterable of peppy.Sample objects to restrict to.
            Must have a `peaks` attribute set.
            If not provided (`None` is passed) if will default to all samples in the analysis (`samples` attribute).

        region_type : str
            The type of region to use to create the consensus region set.
            One of `summits` or `peaks`.
            If `summits`, peak summits will be extended by `extension` before union.
            Otherwise sample peaks will be used with no modification.

        permissive : bool
            Whether Samples that which `region_type` attribute file does
            not exist should be simply skipped or an error thrown.

        **kwargs
            Not used. Provided for compatibility with ChIPSeqAnalysis class.

        Raises
        ----------
        IOError
            If not `permissive` and either the `peaks` or `summits` file of a sample is not readable.
            Or if `permissive` but none of the samples has an existing file.

        Attributes
        ----------
        support : pandas.DataFrame
            Sets a `support` variable with peak set overlap.
        """
        if samples is None:
            samples = self.samples

        # Check which samples to run (dependent on permissive)
        samples = self._get_samples_with_input_file(region_type, permissive=permissive, samples=samples)

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
        except (ValueError, pybedtools.MalformedBedLineError, pybedtools.helpers.BEDToolsError, OverflowError):
            _LOGGER.debug("Could not convert support intersection directly to dataframe, saving temporarly file.")
            support.saveas("_tmp.peaks.bed")
            support = pd.read_csv("_tmp.peaks.bed", sep="\t", header=None)
            os.remove("_tmp.peaks.bed")

        support.columns = ["chrom", "start", "end"] + [sample.name for sample in samples]
        support.index = pd.Index(
            support['chrom'] + ":" + support['start'].astype(str) + "-" + support['end'].astype(str), name="index")
        support.to_csv(os.path.join(self.results_dir, self.name + "_peaks.binary_overlap_support.csv"), index=True)

        # divide sum (of unique overlaps) by total to get support value between 0 and 1
        support["support"] = (
            support[[sample.name for sample in samples]]
            .apply(lambda x: sum([i if i <= 1 else 1 for i in x]) / float(len(samples)), axis=1))
        # save
        support.to_csv(os.path.join(self.results_dir, self.name + "_peaks.support.csv"), index=True)

        setattr(self, 'support', support)
        return self.support

    def get_supported_peaks(
            self, samples=None,
            **kwargs):
        """
        Get mask of sites with 0 support in the given samples.
        Requires support matrix produced by `ngs_toolkit.atacseq.ATACSeqAnalysis.calculate_peak_support`.

        Parameters
        ----------
        samples : list
            Iterable of peppy.Sample objects to restrict to.

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
            self, samples=None, sites=None, assign=True,
            save=True, output_file=None, permissive=False, distributed=False):
        """
        Measure read coverage (counts) of each sample in each region in consensus sites.
        Uses parallel computing using the `parmap` library.
        However, for many samples (hundreds), parallelization in a computing cluster is possible
        with the `distributed` option. Only supports SLURM clusters fow now though.

        Parameters
        ----------
        samples : list
            Iterable of peppy.Sample objects to restrict to. Must have a `filtered` attribute set.
            If not provided (`None` is passed) if will default to all samples in the analysis (`samples` attribute).

        sites : pybedtools.BedTool,pd.DataFrame,str
            Sites in the genome to quantify, usually a pybedtools.BedTool from analysis.get_consensus_sites()
            If a DataFrame, will try to convert to BED format assuming first three columns are chr,start,end.
            If a string assumes a path to a BED file.
            If `None` the object's `sites` attribute will be used.

        assign : bool
            Whether to assign the matrix to an attribute of self named `coverage`.

        save : bool
            Whether to save to disk the coverage matrix with filename `output_file`.

        output_file : str
            A path to a CSV file with coverage output.
            Default is `self.results_dir/self.name + "_peaks.raw_coverage.csv"`.

        permissive : bool
            Whether Samples that which `region_type` attribute file does not exist
            should be simply skipped or an error thrown.

        distributed : bool
            Whether it should be run as jobs for each sample separately in parallel.
            Currently only implemented for a SLURM cluster.
            Default False.

        Raises
        ----------
        IOError
            If not `permissive` and the 'aligned_filtered_bam' file attribute of a sample is not readable.
            Or if `permissive` but none of the samples has an existing file.

        Attributes
        ----------
        coverage : pd.DataFrame
            Sets a `coverage` variable with DataFrame with read counts of shape (n_sites, m_samples).

        Returns
        -------
        pd.DataFrame
            Pandas DataFrame with read counts of shape (n_sites, m_samples).
        """
        from pypiper.ngstk import NGSTk
        if samples is None:
            samples = self.samples

        # Check which samples to run (dependent on permissive)
        samples = self._get_samples_with_input_file("aligned_filtered_bam", permissive=permissive, samples=samples)

        if sites is None:
            sites = self.sites

        if not distributed:
            # Count reads with pysam
            # make strings with intervals
            if isinstance(sites, pybedtools.BedTool):
                sites_str = [
                    str(i.chrom) + ":" +
                    str(i.start) + "-" +
                    str(i.stop) for i in self.sites]
            elif isinstance(sites, pd.core.frame.DataFrame):
                sites_str = (
                    sites.iloc[:, 0] + ":" +
                    sites.iloc[:, 1].astype(str) + "-" +
                    sites.iloc[:, 2].astype(str)).astype(str).tolist()
            elif isinstance(sites, str):
                sites_str = [
                    str(i.chrom) + ":" +
                    str(i.start) + "-" +
                    str(i.stop) for i in pybedtools.BedTool(sites)]

            # count, create dataframe
            coverage = pd.DataFrame(
                map(
                    lambda x:
                        pd.Series(x),
                        parmap.map(
                            count_reads_in_intervals,
                            [sample.aligned_filtered_bam for sample in samples],
                            sites_str,
                            parallel=True
                        )
                ),
                index=[sample.name for sample in samples]
            ).T

            if assign:
                self.coverage = coverage
            if save:
                if output_file is not None:
                    coverage.to_csv(output_file, index=True)
                else:
                    self.coverage.to_csv(os.path.join(
                        self.results_dir, self.name + "_peaks.raw_coverage.csv"), index=True)
            return coverage

        else:
            import textwrap
            tk = NGSTk()
            for s in samples:
                output_dir = os.path.join(s.paths.sample_root, "coverage")
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)

                job_name = "peak_set_coverage.{}".format(s.name)
                output_file = os.path.join(output_dir, s.name + ".peak_set_coverage.bed")
                log_file = os.path.join(output_dir, s.name + ".peak_set_coverage.log")
                job_file = os.path.join(output_dir, s.name + ".peak_set_coverage.sh")

                cmd = tk.slurm_header(
                    job_name=job_name, output=log_file,
                    cpus_per_task=1, mem_per_cpu=8000)
                cmd += "bedtools coverage -counts"
                cmd += " -abam {}".format(s.aligned_filtered_bam)
                cmd += " -b {}".format(sites.fn)
                cmd += " > {}".format(output_file)
                cmd += " \n"
                cmd += tk.slurm_footer()

                with open(job_file, "w") as handle:
                    handle.write(textwrap.dedent(cmd).replace("\n ", "\n").replace("  ", ""))

                tk.slurm_submit_job(job_file)

    def collect_coverage(
            self, samples=None, sites=None, assign=True,
            save=True, output_file=None, permissive=False):
        """
        Collect read coverage (counts) of each sample in each region in consensus sites from existing files.
        Useful after runnning analysis.measure_coverage() in distributed mode.

        Importantly, it assumes all coverage files are ordered equally and that this order matches the order in `sites`!

        Parameters
        ----------
        samples : list
            Iterable of peppy.Sample objects to restrict to.
            If not provided (`None` is passed) if will default to all samples in the analysis (`samples` attribute).

        sites : pybedtools.BedTool,pd.DataFrame,str
            Sites in the genome to quantify, usually a pybedtools.BedTool from analysis.get_consensus_sites()
            If a DataFrame, will try to convert to BED format assuming first three columns are chr,start,end.
            If a string assumes a path to a BED file.
            If `None` the object's `sites` attribute will be used.

        assign : bool
            Whether to assign the matrix to an attribute of self named `coverage`.

        save : bool
            Whether to save to disk the coverage matrix with filename `output_file`.

        output_file : str
            A path to a CSV file with coverage output.
            Default is `self.results_dir/self.name + "_peaks.raw_coverage.csv"`.

        permissive : bool
            Whether Samples without an existing coverage file does not exist
            should be simply skipped or an error thrown.

        Raises
        ----------
        IOError
            If not `permissive` and the coverage file of a sample is not readable.
            Or if `permissive` but none of the samples has an existing file.

        Attributes
        ----------
        coverage : pd.DataFrame
            Sets a `coverage` variable with DataFrame with read counts of shape (n_sites, m_samples).

        Returns
        -------
        pd.DataFrame
            Pandas DataFrame with read counts of shape (n_sites, m_samples).
        """
        if samples is None:
            samples = self.samples

        for sample in samples:
            setattr(
                sample, "_coverage",
                os.path.join(sample.paths.sample_root, "coverage", sample.name + ".peak_set_coverage.bed"))

        # Check which samples to run (dependent on permissive)
        samples = self._get_samples_with_input_file("_coverage", permissive=permissive, samples=samples)

        # Read in counts
        for i, sample in tqdm(enumerate(samples), total=len(samples)):
            if i == 0:
                cov = pd.read_csv(
                    sample._coverage, sep="\t", header=None,
                    names=['chrom', 'start', 'end', sample.name])
                coverage = pd.DataFrame(
                    index=bed_to_index(cov),
                    columns=[s.name for s in samples])
                coverage.loc[:, sample.name] = cov[sample.name].values
            else:
                coverage.loc[:, sample.name] = pd.read_csv(
                    sample._coverage, sep="\t", header=None,
                    names=['chrom', 'start', 'end', sample.name]).loc[:, sample.name].values

        if assign:
            self.coverage = coverage
        if save:
            if output_file is not None:
                coverage.to_csv(output_file, index=True)
            else:
                self.coverage.to_csv(os.path.join(
                    self.results_dir, self.name + "_peaks.raw_coverage.csv"), index=True)
        return coverage

    def normalize_coverage_rpm(
            self, matrix=None, samples=None,
            mult_factor=1e6, log_transform=True, pseudocount=1,
            save=True, assign=True):
        """
        Normalization of matrix of (n_features, n_samples) by total in each sample.

        Parameters
        ----------
        matrix : str
            Attribute name of matrix to normalize. Defaults to 'coverage'.

        samples : list
            Iterable of peppy.Sample objects to restrict matrix to.
            If not provided (`None` is passed) the matrix will not be subsetted.

        mult_factor : float
            A constant to multiply values for.

        log_transform : bool
            Whether to log transform values or not.

        pseudocount : int|float
            A constant to add to values.

        save : bool
            Whether to write normalized DataFrame to disk.

        assign : bool
            Whether to assign the normalized DataFrame to an attribute ``.
        """
        to_norm = self.get_matrix(matrix=matrix, samples=samples, matrix_name="coverage")
        # apply normalization over total
        coverage_rpm = ((pseudocount + to_norm) / (pseudocount + to_norm).sum()) * mult_factor

        # Log2 transform
        if log_transform:
            coverage_rpm = np.log2(coverage_rpm)

        # coverage_rpm = coverage_rpm.join(self.coverage[['chrom', 'start', 'end']])
        if save:
            coverage_rpm.to_csv(os.path.join(self.results_dir, self.name + "_peaks.coverage_rpm.csv"), index=True)
        if assign:
            self.coverage_rpm = coverage_rpm

        return coverage_rpm

    def normalize_coverage_quantiles(
            self, matrix=None, samples=None, implementation="Python",
            log_transform=True, pseudocount=1, save=True, assign=True):
        """
        Quantile normalization of matrix of (n_features, n_samples).

        Parameters
        ----------
        matrix : str
            Attribute name of matrix to normalize. Defaults to 'coverage'.

        samples : list
            Iterable of peppy.Sample objects to restrict matrix to.
            If not provided (`None` is passed) the matrix will not be subsetted.

        implementation : str
            One of `"R"` or `"Python"`. Dictates which implementation is to be used.
            The R implementation comes from the `preprocessCore` package,
            and the Python one is from https://github.com/ShawnLYU/Quantile_Normalize.

        log_transform : bool
            Whether to log transform values or not.

        pseudocount : float
            A constant to add before log transformation.

        save : bool
            Whether to write normalized DataFrame to disk.

        assign : bool
            Whether to assign the normalized DataFrame to an attribute `coverage_qnorm`.
        """
        if matrix is None:
            to_norm = self.get_matrix(matrix=matrix, samples=samples, matrix_name="coverage")
        else:
            to_norm = matrix

        if implementation == "R":
            coverage_qnorm = pd.DataFrame(
                normalize_quantiles_r(to_norm.values),
                index=to_norm.index,
                columns=to_norm.columns)
        elif implementation == "Python":
            coverage_qnorm = normalize_quantiles_p(to_norm)
        else:
            msg = "Implementation of quantile normalization must be one of 'R' of 'Python'"
            _LOGGER.error(msg)
            raise ValueError(msg)

        # Log2 transform
        if log_transform:
            coverage_qnorm = np.log2(pseudocount + coverage_qnorm)

        if coverage_qnorm.min().min() <= 0:
            coverage_qnorm += coverage_qnorm.abs().min().min()

        # # Add back postition columns
        # coverage_qnorm = coverage_qnorm.join(self.coverage[['chrom', 'start', 'end']])
        if save:
            coverage_qnorm.to_csv(os.path.join(
                self.results_dir, self.name + "_peaks.coverage_qnorm.csv"), index=True)
        if assign:
            self.coverage_qnorm = coverage_qnorm

        return coverage_qnorm

    @check_organism_genome
    def get_peak_gccontent_length(
            self, bed_file=None, fasta_file=None):
        """
        Get length and GC content of features in region set.

        bed_file : str
            A BED file with regions to calculate GC content on. Must be a 3-column BED!
            If not provided the calculation will be for the analysis `sites` attribute.

        genome : str
            Genome assembly.

        fasta_file : str
            Fasta file of `genome`. Preferably indexed. If not given, will try to download.

        :var nuc:
            DataFrame with nucleotide content and length of each region.

        Returns
        -------
        pandas.DataFrame
            DataFrame with nucleotide content and length of each region.
        """
        if bed_file is None:
            sites = self.sites
        else:
            sites = pybedtools.BedTool(bed_file)

        if fasta_file is None:
            _LOGGER.info("Reference genome FASTA file was not given, will try to get it.")
            _LOGGER.info("Getting genome FASTA file for organism '{}', genome '{}'. "
                         .format(self.organism, self.genome))
            fasta_file = self.get_annotations(steps=['genome'])['genome_file']['fasta']

        nuc = sites.nucleotide_content(fi=fasta_file).to_dataframe(comment="#")[["score", "blockStarts"]]
        nuc.columns = ["gc_content", "length"]
        nuc.index = [str(i.chrom) + ":" + str(i.start) + "-" + str(i.stop) for i in sites]

        # get only the sites matching the coverage (not overlapping blacklist)
        self.nuc = nuc.ix[self.coverage.index]

        self.nuc.to_csv(os.path.join(
            self.results_dir, self.name + "_peaks.gccontent_length.csv"), index=True)

        return self.nuc

    def normalize_gc_content(self, matrix=None, samples=None, save=True, assign=True):
        """
        Quantile normalization of matrix of (n_features, n_samples) followed by GC content
        correction by regression.

        Requires the R package "cqn" to be installed:
            >>> source('http://bioconductor.org/biocLite.R')
            >>> biocLite('cqn')

        Parameters
        ----------
        matrix : str
            Attribute name of matrix to normalize. Defaults to 'coverage'.

        samples : list
            Iterable of peppy.Sample objects to restrict matrix to.
            If not provided (`None` is passed) the matrix will not be subsetted.

        save : bool
            Whether to write normalized DataFrame to disk.

        assign : bool
            Whether to assign the normalized DataFrame to an attribute ``.
        """
        def cqn(cov, gc_content, lengths):
            # install R package
            # source('http://bioconductor.org/biocLite.R')
            # biocLite('cqn')
            import rpy2.robjects.pandas2ri
            import warnings
            warnings.filterwarnings("ignore", category=RRuntimeWarning)
            robjects.numpy2ri.deactivate()
            rpy2.robjects.pandas2ri.activate()

            robjects.r('require("cqn")')
            cqn = robjects.r('cqn')

            cqn_out = cqn(cov, x=gc_content, lengths=lengths)

            y_r = cqn_out[list(cqn_out.names).index('y')]
            y = pd.DataFrame(
                np.array(y_r),
                index=cov.index,
                columns=cov.columns)
            offset_r = cqn_out[list(cqn_out.names).index('offset')]
            offset = pd.DataFrame(
                np.array(offset_r),
                index=cov.index,
                columns=cov.columns)

            return y + offset

        # Perform quantile normalization first
        if not hasattr(self, "nuc"):
            self.normalize_coverage_quantiles(samples)

        # Get GC content and length of each feature
        if not hasattr(self, "nuc"):
            self.get_peak_gccontent_length()

        to_norm = self.get_matrix(matrix=matrix, samples=samples, matrix_name="coverage")
        coverage_gc_corrected = (
            cqn(cov=to_norm, gc_content=self.nuc["gc_content"], lengths=self.nuc["length"])
            # .join(self.coverage[['chrom', 'start', 'end']])
        )

        if save:
            coverage_gc_corrected.to_csv(os.path.join(
                self.results_dir, self.name + "_peaks.coverage_gc_corrected.csv"), index=True)
        if assign:
            self.coverage_gc_corrected = coverage_gc_corrected

        return coverage_gc_corrected

    def normalize(self, method="quantile", matrix=None, samples=None, save=True, assign=True):
        """
        Normalization of matrix of (n_features, n_samples).

        Parameters
        ----------
        method : str
            Normalization method to apply. One of:
             - `total`: Reads per total normalization (RPM).
             - `quantile`: Quantile normalization and log2 transformation.
             - `gc_content`: Quantile normalization followed by GC content
                             correction by regression (cqn R package)
                             and log2 transformation.

        matrix : str
            Attribute name of matrix to normalize.

        samples : list
            Iterable of peppy.Sample objects to restrict matrix to.
            If not provided (`None` is passed) the matrix will not be subsetted.

        save : bool
            Whether to write normalized DataFrame to disk.

        assign : bool
            Whether to assign the normalized DataFrame to an attribute
            (see variables below for each respective normalization type).

        Attributes
        ----------
        {coverage_rpm, coverage_qnorm, coverage_gc_corrected} : pd.DataFrame
            If `assign` is True, a pandas DataFrame normalized with respective method.

        Returns
        -------
        pd.DataFrame
            Normalized pandas DataFrame.
        """
        if method == "total":
            return self.normalize_coverage_rpm(
                matrix=matrix, samples=samples, save=save, assign=assign)
        elif method == "quantile":
            return self.normalize_coverage_quantiles(
                matrix=matrix, samples=samples, save=save, assign=assign)
        elif method == "gc_content":
            return self.normalize_gc_content(
                matrix=matrix, samples=samples, save=save, assign=assign)
        else:
            msg = "Requested normalization method is not available!"
            _LOGGER.error(msg)
            raise ValueError(msg)

    @check_has_sites
    def get_peak_gene_annotation(self, tss_file=None, max_dist=100000):
        """
        Annotates peaks with closest gene.
        The annotation reference can either be given in the `tss_file` parameter
        but if ommited, it will be fetched if analysis has `genome` and `organism`
        attributes.
        A dataframe with each feature's distance to the nearest gene is also saved.

        Parameters
        ----------
        tss_file : str
            A valid BED file where the name field (4th column) identifies the gene
            and the strand column (6th column). Other fields will not be used.

        max_dist : int
            Maximum absolute distance allowed to perform associations.
            Regions with no genes within the range will have NaN values.

        Attributes
        ----------
        gene_annotation : pd.DataFrame
            A pandas DataFrame containing the genome annotations of the region features.
            If a feature overlaps more than one gene, the two gene values will be concatenated with a comma.

        closest_tss_distances : pd.DataFrame
            A pandas DataFrame containing unique region->gene associations.
            In contrast to gene_annotation dataframe, this contains one row per region->gene assignment.

        Returns
        -------
        pandas.DataFrame
            A dataframe with genes annotated for the peak set.
        """
        cols = [6, 8, 9]

        if tss_file is None:
            _LOGGER.info("Reference TSS file was not given, will try to get TSS annotations.")
            _LOGGER.info("Getting TSS annotations for organism '{}', genome '{}'. "
                         .format(self.organism, self.genome))
            tss_file = self.get_annotations(steps=['tss'])['tss_file']

        # extract only relevant columns
        tss = pd.read_csv(tss_file, header=None, sep="\t")
        tss = tss.iloc[:, list(range(6))]
        tss = pybedtools.BedTool.from_dataframe(tss)

        if isinstance(self.sites, str):
            self.sites = pybedtools.BedTool(self.sites)

        # get closest TSS of each region
        self.closest_tss_distances = self.sites.closest(tss, D="b").to_dataframe()

        # rename
        self.closest_tss_distances = (self.closest_tss_distances[
            ['chrom', 'start', 'end'] +
            self.closest_tss_distances.columns[cols].tolist()])
        self.closest_tss_distances.columns = [
            'chrom', 'start', 'end', 'gene_name', "strand", 'distance']

        # set NaN to distance without assignment (rather than the default '-1' from bedtools)
        self.closest_tss_distances.loc[
            self.closest_tss_distances['gene_name'] == '.', 'distance'] = np.nan

        # set NaN to assignments out of range
        self.closest_tss_distances.loc[
            self.closest_tss_distances['distance'].abs() > max_dist,
            ['gene_name', 'strand', 'distance']] = np.nan

        # aggregate annotation per peak, concatenate various genes (comma-separated)
        self.gene_annotation = (
            self.closest_tss_distances
            .groupby(['chrom', 'start', 'end'])
            .aggregate(
                lambda x: ",".join(set([str(i) for i in x if i != '.'])))
            .reset_index())
        self.closest_tss_distances.index = bed_to_index(self.closest_tss_distances)
        self.gene_annotation.index = bed_to_index(self.gene_annotation)

        # save to disk
        self.closest_tss_distances.to_csv(
            os.path.join(self.results_dir,
                         self.name + "_peaks.closest_tss_distances.csv"),
            index=True)
        self.gene_annotation.to_csv(
            os.path.join(self.results_dir,
                         self.name + "_peaks.gene_annotation.csv"),
            index=True)
        return self.gene_annotation

    @check_organism_genome
    @check_has_sites
    def get_peak_genomic_location(
            self, genomic_context_file=None):
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
        genomic_context_file : str
            A 4 column BED file (chrom, start, end, feature), where feature is a string with the type of region.
            If not provided will be get with the get_genomic_context function.

        Attributes
        ----------
        region_annotation, region_annotation_b : pd.DataFrame
            A DataFrame with the genome annotations of the region features or genome background.

        region_annotation_mapping, region_annotation_b_mapping : pd.DataFrame
            A DataFrame with one row for each chromatin state-region mapping or genome background.

        Returns
        -------
        pandas.DataFrame
            The genomic context annotation for the peak set.
        """

        if genomic_context_file is None:
            _LOGGER.info("Reference genomic context file was not given, will try to get it.")
            _LOGGER.info("Getting genomic context annotations for organism '{}', genome '{}'. "
                         .format(self.organism, self.genome))
            genomic_context_file = self.get_annotations(steps=['genomic_context'])['genomic_context_file']

        context = pybedtools.BedTool(genomic_context_file)

        if isinstance(self.sites, str):
            self.sites = pybedtools.BedTool(self.sites)

        # create background
        # shuffle regions in genome to create background (keep them in the same chromossome)
        background = self.sites.shuffle(genome=self.genome, chrom=True)

        cols = [0, 1, 2, 6]
        for label, attr, bed in [
                ("background", "region_annotation_b", background),
                ("real", "region_annotation", self.sites)]:
            annot = (
                bed.intersect(context, wa=True, wb=True, f=0.2)
                .sort().to_dataframe(usecols=cols))
            annot.index = bed_to_index(annot)
            annot.columns = ['chrom', 'start', 'end', 'genomic_region']

            # remove duplicates (there shouldn't be anyway)
            annot = annot.drop_duplicates()
            # join various annotations per peak
            annot_comp = (
                annot.groupby(['chrom', 'start', 'end'])
                .aggregate(lambda x: ",".join(set([str(i) for i in x]))).reset_index())
            annot_comp.index = bed_to_index(annot_comp)
            annot_comp.columns = ['chrom', 'start', 'end', 'genomic_region']
            # save to disk
            a = "" if (label == "real") else ("_" + label)
            annot.to_csv(os.path.join(
                self.results_dir, self.name + "_peaks.region_annotation{}_mapping.csv".format(a)),
                index=True)
            annot_comp.to_csv(os.path.join(
                self.results_dir, self.name + "_peaks.region_annotation{}.csv".format(a)),
                index=True)

        setattr(self, attr, annot_comp)
        setattr(self, attr + "_mapping", annot)
        return self.region_annotation

    @check_organism_genome
    @check_has_sites
    def get_peak_chromatin_state(
            self, chrom_state_file, frac=0.2):
        """
        Annotates a consensus peak set (``sites`` attribute of analysis) with their chromatin
        state context. This would be given, for example by a chromatin state segmentation
        file from projects such as Roadmap Epigenomics.

        See examples of such files for various cell types/assemblies here:
        https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/
        (the *_dense.bed.gz files are optimal).

        Parameters
        ----------
        chrom_state_file : str
            A 4 column BED file (chrom, start, end, feature), where feature is a string with the type of region.

        frac : float
            Minimal fraction of region to overlap with a feature.

        Returns
        ----------
        pandas.DataFrame
            The chromatin state annotation for the peak set.

        Attributes
        ----------
        chrom_state_annotation, chrom_state_annotation_b : pd.DataFrame
            A DataFrame with the chromatin state annotations of
            the region features or of the genome background.

        chrom_state_annotation_mapping, chrom_state_annotation_b_mapping : pd.DataFrame
            A DataFrame with one row for each chromatin state-region mapping
            or for the genome background.
        """
        states = pybedtools.BedTool(chrom_state_file)

        if isinstance(self.sites, str):
            self.sites = pybedtools.BedTool(self.sites)

        # create background
        # shuffle regions in genome to create background (keep them in the same chromossome)
        background = self.sites.shuffle(genome=self.genome, chrom=True)

        for label, attr, bed in [
                ("real", "chrom_state_annotation", self.sites),
                ("background", "chrom_state_annotation_b", background)]:
            _LOGGER.debug("Getting chromatinn state annotation for {} regions."
                          .format(label))
            annot = (
                bed.intersect(states, wa=True, wb=True, f=frac, loj=True)
                .to_dataframe(usecols=[0, 1, 2, 6]))
            annot.iloc[:, 3] = annot.iloc[:, 3].astype(str)

            # remove duplicates (there shouldn't be anyway)
            annot = annot.drop_duplicates()
            annot.index = bed_to_index(annot)
            annot.columns = ['chrom', 'start', 'end', 'chromatin_state']
            # join various annotations per peak
            annot_comp = (
                annot.groupby(['chrom', 'start', 'end'])
                .aggregate(lambda x: ",".join(set([str(i) for i in x]))).reset_index())
            annot_comp.index = bed_to_index(annot_comp)
            annot_comp.columns = ['chrom', 'start', 'end', 'chromatin_state']
            # save to disk
            a = "" if (label == "real") else ("_" + label)
            annot.to_csv(os.path.join(
                self.results_dir, self.name + "_peaks.chrom_state_annotation{}_mapping.csv".format(a)),
                index=True)
            annot_comp.to_csv(os.path.join(
                self.results_dir, self.name + "_peaks.chrom_state_annotation{}.csv".format(a)),
                index=True)

            setattr(self, attr, annot_comp)
            setattr(self, attr + "_mapping", annot)
        return self.chrom_state_annotation

    def get_matrix_stats(self, quant_matrix=None, samples=None):
        """
        Gets a matrix of feature-wise (i.e. for every reg. element) statistics such
        across samples such as mean, variance, deviation, dispersion and amplitude.

        Parameters
        ----------
        matrix : str
            Attribute name of matrix to normalize. Defaults to 'coverage'.

        Returns
        -------
        pandas.DataFrame
            Statistics for each feature.

        Attributes
        ----------
        stats : pd.DataFrame
            A DataFrame with statistics for each feature.
        """
        if samples is None:
            samples = self.samples
        if quant_matrix is None:
            quant_matrix = "coverage_gc_corrected"
        quant_matrix = getattr(self, quant_matrix)

        quant_matrix = quant_matrix.loc[:, [s.name for s in samples]]

        matrix = pd.DataFrame(index=pd.Index(quant_matrix.index, name="region"))
        # calculate mean coverage
        matrix.loc[:, 'mean'] = quant_matrix.mean(axis=1)
        # calculate coverage variance
        matrix.loc[:, 'variance'] = quant_matrix.var(axis=1)
        # calculate std deviation (sqrt(variance))
        matrix.loc[:, 'std_deviation'] = np.sqrt(matrix.loc[:, 'variance'])
        # calculate dispersion (variance / mean)
        matrix.loc[:, 'dispersion'] = matrix.loc[:, 'variance'] / matrix.loc[:, 'mean']
        # calculate qv2 (std / mean) ** 2
        matrix.loc[:, 'qv2'] = (matrix.loc[:, 'std_deviation'] / matrix.loc[:, 'mean']) ** 2
        # calculate "amplitude" (max - min)
        matrix.loc[:, 'amplitude'] = (quant_matrix.max(axis=1) - quant_matrix.min(axis=1))
        # calculate interquantile range
        matrix.loc[:, 'iqr'] = (quant_matrix.quantile(0.75, axis=1) - quant_matrix.quantile(0.25, axis=1))
        matrix.index.name = "index"
        matrix.to_csv(os.path.join(
            self.results_dir, self.name + "_peaks.stats_per_region.csv"),
            index=True)

        self.stats = matrix
        return self.stats

    def annotate(self, samples=None, quant_matrix=None, permissive=True):
        """
        Annotates analysis regions by aggregating region-wise annotations
        (region, chromatin state, gene annotations and statistics - if present).

        The numeric matrix to be used is specified in `quant_matrix`.
        If two annotation dataframes have equally named columns (e.g. chrom, start, end),
        the value of the first is kept.

        Parameters
        ----------
        samples : list
            Iterable of peppy.Sample objects to restrict matrix to.
            If not provided (`None` is passed) the matrix will not be subsetted.
            Calculated metrics will be restricted to these samples.

        quant_matrix : str
            Attribute name of matrix to annotate.

        permissive : bool
            Whether DataFrames that do not exist should be simply skipped or an error will be thrown.

        Raises
        ----------
        AttributeError
            If not `permissive` a required DataFrame does not exist as an object attribute.

        Attributes
        ----------
        coverage_annotated : pd.DataFrame
            A pandas DataFrame containing annotations of the region features.
        """
        if samples is None:
            samples = self.samples
        # TODO: come up with a resonable iterative way to figure out which matrix to use by default
        if quant_matrix is None:
            quant_matrix = "coverage_gc_corrected"
        quant_matrix = getattr(self, quant_matrix)

        next_matrix = quant_matrix
        # add closest gene
        msg = "`{}` attribute does not exist."

        # TODO: decide if this function should call the others to produce the annotations first

        for matrix_name in ['gene_annotation', 'region_annotation', 'chrom_state_annotation', 'support', 'stats']:
            if hasattr(self, matrix_name):
                matrix = getattr(self, matrix_name)
                self.coverage_annotated = pd.merge(
                    next_matrix,
                    matrix[matrix.columns.difference(next_matrix.columns)],
                    left_index=True, right_index=True, how="left")
                next_matrix = self.coverage_annotated
            else:
                if not permissive:
                    _LOGGER.error(msg.format(matrix_name))
                    raise AttributeError(msg.format(matrix_name))
                else:
                    _LOGGER.debug(msg.format(matrix_name) + " Proceeding anyway.")

        if not hasattr(self, "coverage_annotated"):
            self.coverage_annotated = next_matrix

        # Pair indexes
        msg = "Annotated matrix does not have same feature length as coverage matrix."
        if not self.coverage.shape[0] == self.coverage_annotated.shape[0]:
            _LOGGER.error(msg)
            raise AssertionError(msg)
        self.coverage_annotated.index = self.coverage.index
        self.coverage_annotated.index.name = "index"

        # TODO: call annotate with samples_metadata maybe?

        # Save
        self.coverage_annotated.to_csv(os.path.join(self.results_dir, self.name + "_peaks.coverage_qnorm.annotated.csv"), index=True)

    def get_gene_level_accessibility(self, matrix="accessibility", reduce_func=np.mean):
        """
        Get gene-level measurements of coverage.
        Requires a 'gene_annotation' attribute to be set containing a mapping between the
        index of `matrix` and genes (produced from `get_peak_gene_annotation`).

        Parameters
        ----------
        matrix : str
            Quantification matrix to be used (e.g. 'coverage' or 'accessibility')

        reduce_func : func
            Function to apply to reduce values.
            Default is mean.

        Returns
        ---------
        pandas.DataFrame
            Chromatin accessibility values reduced per gene.
        """
        msg = "Analysis object lacks 'gene_annotation' dataframe."
        hint = " Call 'analysis.get_peak_gene_annotation' to have region-gene associations."
        if not hasattr(self, "gene_annotation"):
            _LOGGER.error(msg + hint)
            raise AssertionError(msg)

        matrix = getattr(self, matrix).copy()

        # if isinstance(matrix.columns, pd.MultiIndex):
        #     matrix.columns = matrix.columns.get_level_values("sample_name")

        g = self.gene_annotation['gene_name'].str.split(",").apply(pd.Series).stack()
        g.index = g.index.droplevel(1)
        g.name = "gene_name"

        matrix2 = matrix.join(g).drop("gene_name", axis=1)
        matrix2.index = matrix.join(g).reset_index().set_index(['index', 'gene_name']).index
        matrix2.columns = matrix.columns
        matrix3 = matrix2.groupby(level=['gene_name']).apply(reduce_func)

        return matrix3.loc[:, ~matrix3.isnull().all()]

    def get_gene_level_changes(self, differential_results=None, reduce_func=np.mean):
        """
        Redcuce changes in regulatory elements to gene-level by aggregating across regulatory elements.
        Requires a 'gene_annotation' attribute to be set containing a mapping between
        the index of `matrix` and genes (produced from `get_peak_gene_annotation`).

        Parameters
        ----------
        differential_results : pandas.DataFrame
            Matrix with differential results to use.
            Default is a 'differential_results' attribute of self.

        reduce_func : func
            Function to apply to reduce values. Default is mean

        Returns
        ---------
        pandas.DataFrame
            Changes in chromatin accessibility (log2FoldChanges) reduced per gene.
        """
        msg = "Analysis object lacks 'gene_annotation' dataframe."
        hint = " Call 'analysis.get_peak_gene_annotation' to have region-gene associations."
        if not hasattr(self, "gene_annotation"):
            _LOGGER.error(msg + hint)
            raise AssertionError(msg)

        if differential_results is None:
            differential_results = self.differential_results

        g = self.gene_annotation['gene_name'].str.split(",").apply(pd.Series).stack()
        g.index = g.index.droplevel(1)
        g.name = "gene_name"

        dr2 = differential_results.join(g).drop("gene_name", axis=1)
        dr2.index = differential_results.join(g).reset_index().set_index(['index', 'gene_name']).index
        dr2.columns = differential_results.columns
        dr3 = dr2.reset_index().groupby(['gene_name', 'comparison_name']).apply(reduce_func)

        return dr3.loc[:, ~dr3.isnull().all()]

    def plot_peak_characteristics(self, samples=None, by_attribute=None, genome_space=3e9):
        """
        Several diagnostic plots on the peak set and the Sample's signal on them.

        Provides plots with Samples grouped `by_attribute` if given (a string or a list of strings).
        """
        def get_sample_reads(bam_file):
            return pysam.AlignmentFile(bam_file).count()

        def get_peak_number(bed_file):
            return len(open(bed_file, "r").read().split("\n"))

        def get_total_open_chromatin(bed_file):
            peaks = pd.read_csv(bed_file, sep="\t", header=None)
            return (peaks.iloc[:, 2] - peaks.iloc[:, 1]).sum()

        def get_peak_lengths(bed_file):
            peaks = pd.read_csv(bed_file, sep="\t", header=None)
            return (peaks.iloc[:, 2] - peaks.iloc[:, 1])

        def get_peak_chroms(bed_file):
            peaks = pd.read_csv(bed_file, sep="\t", header=None)
            return peaks.iloc[:, 0].value_counts()

        if samples is None:
            samples = self.samples

        # TODO: add parallelization with parmap
        stats = pd.DataFrame([
            map(get_sample_reads, [s.aligned_filtered_bam for s in samples]),
            map(get_peak_number, [s.peaks for s in samples]),
            map(get_total_open_chromatin, [s.peaks for s in samples])],
            index=["reads_used", "peak_number", "open_chromatin"],
            columns=[s.name for s in samples]).T

        stats["peaks_norm"] = (stats["peak_number"] / stats["reads_used"]) * 1e3
        stats["open_chromatin_norm"] = (stats["open_chromatin"] / stats["reads_used"])
        stats.to_csv(os.path.join(self.results_dir, "{}.open_chromatin_space.csv".format(self.name)), index=True)
        stats = pd.read_csv(os.path.join(self.results_dir, "{}.open_chromatin_space.csv".format(self.name)), index_col=0)

        # median lengths per sample (split-apply-combine)
        if by_attribute is not None:
            stats = pd.merge(stats, stats.groupby(by_attribute)['open_chromatin'].median().to_frame(name='group_open_chromatin').reset_index())
            stats = pd.merge(stats, stats.groupby(by_attribute)['open_chromatin_norm'].median().to_frame(name='group_open_chromatin_norm').reset_index())

        # plot
        stats = stats.sort_values("open_chromatin_norm")
        fig, axis = plt.subplots(2, 1, figsize=(4 * 2, 6 * 1))
        sns.barplot(x="index", y="open_chromatin", data=stats.reset_index(), palette="summer", ax=axis[0])
        sns.barplot(x="index", y="open_chromatin_norm", data=stats.reset_index(), palette="summer", ax=axis[1])
        axis[0].set_ylabel("Total open chromatin space (bp)")
        axis[1].set_ylabel("Total open chromatin space (normalized)")
        axis[0].set_xticklabels(axis[0].get_xticklabels(), rotation=45, ha="right")
        axis[1].set_xticklabels(axis[1].get_xticklabels(), rotation=45, ha="right")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.total_open_chromatin_space.per_sample.svg".format(self.name)), bbox_inches="tight")

        if by_attribute is not None:
            # median lengths per group (split-apply-combine)
            stats = pd.merge(stats, stats.groupby(by_attribute)['open_chromatin'].median().to_frame(name='group_open_chromatin').reset_index())
            stats = stats.sort_values("group_open_chromatin")

            fig, axis = plt.subplots(2, 1, figsize=(4 * 2, 6 * 1))
            stats = stats.sort_values("group_open_chromatin")
            sns.barplot(x=by_attribute, y="open_chromatin", data=stats.reset_index(), palette="summer", ax=axis[0])
            sns.stripplot(x=by_attribute, y="open_chromatin", data=stats.reset_index(), palette="summer", ax=axis[0])
            stats = stats.sort_values("group_open_chromatin_norm")
            sns.barplot(x=by_attribute, y="open_chromatin_norm", data=stats.reset_index(), palette="summer", ax=axis[1])
            sns.stripplot(x=by_attribute, y="open_chromatin_norm", data=stats.reset_index(), palette="summer", ax=axis[1])
            axis[0].axhline(stats.groupby(by_attribute)['open_chromatin'].median()["WT"], color="black", linestyle="--")
            axis[1].axhline(stats.groupby(by_attribute)['open_chromatin_norm'].median()["WT"], color="black", linestyle="--")
            axis[0].set_ylabel("Total open chromatin space (bp)")
            axis[1].set_ylabel("Total open chromatin space (normalized)")
            axis[0].set_xticklabels(axis[0].get_xticklabels(), rotation=45, ha="right")
            axis[1].set_xticklabels(axis[1].get_xticklabels(), rotation=45, ha="right")
            sns.despine(fig)
            fig.savefig(os.path.join(self.results_dir, "{}.total_open_chromatin_space.per_{}.svg".format(self.name, by_attribute)), bbox_inches="tight")

        # plot distribution of peak lengths
        sample_peak_lengths = map(get_peak_lengths, [s.peaks for s in samples])
        lengths = pd.melt(pd.DataFrame(sample_peak_lengths, index=[s.name for s in samples]).T, value_name='peak_length', var_name="sample_name").dropna()

        # median lengths per sample (split-apply-combine)
        lengths = pd.merge(lengths, lengths.groupby('sample_name')['peak_length'].median().to_frame(name='mean_peak_length').reset_index())

        lengths = lengths.sort_values("mean_peak_length")
        fig, axis = plt.subplots(2, 1, figsize=(8 * 1, 4 * 2))
        sns.boxplot(x="sample_name", y="peak_length", data=lengths, palette="summer", ax=axis[0], showfliers=False)
        axis[0].set_ylabel("Peak length (bp)")
        axis[0].set_xticklabels(axis[0].get_xticklabels(), visible=False)
        sns.boxplot(x="sample_name", y="peak_length", data=lengths, palette="summer", ax=axis[1], showfliers=False)
        axis[1].set_yscale("log")
        axis[1].set_ylabel("Peak length (bp)")
        axis[1].set_xticklabels(axis[1].get_xticklabels(), rotation=45, ha="right")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.peak_lengths.per_sample.svg".format(self.name)), bbox_inches="tight")

        if by_attribute is not None:
            # median lengths per group (split-apply-combine)
            lengths = pd.merge(lengths, lengths.groupby(by_attribute)['peak_length'].median().to_frame(name='group_mean_peak_length').reset_index())
            lengths = lengths.sort_values("group_mean_peak_length")
            fig, axis = plt.subplots(2, 1, figsize=(8 * 1, 4 * 2))
            sns.boxplot(x=by_attribute, y="peak_length", data=lengths, palette="summer", ax=axis[0], showfliers=False)
            axis[0].set_ylabel("Peak length (bp)")
            axis[0].set_xticklabels(axis[0].get_xticklabels(), visible=False)
            sns.boxplot(x=by_attribute, y="peak_length", data=lengths, palette="summer", ax=axis[1], showfliers=False)
            axis[1].set_yscale("log")
            axis[1].set_ylabel("Peak length (bp)")
            axis[1].set_xticklabels(axis[1].get_xticklabels(), rotation=45, ha="right")
            sns.despine(fig)
            fig.savefig(os.path.join(self.results_dir, "{}.peak_lengths.per_{}.svg".format(self.name, by_attribute)), bbox_inches="tight")

        # peaks per chromosome per sample
        chroms = pd.DataFrame(map(get_peak_chroms, [s.peaks for s in samples]), index=[s.name for s in samples]).fillna(0).T
        chroms_norm = (chroms / chroms.sum(axis=0)) * 100
        chroms_norm = chroms_norm.ix[["chr{}".format(i) for i in range(1, 23) + ['X', 'Y', 'M']]]

        fig, axis = plt.subplots(1, 1, figsize=(8 * 1, 8 * 1))
        sns.heatmap(
            chroms_norm, square=True, cmap="summer",
            xticklabels=True, yticklabels=True, ax=axis)
        axis.set_xticklabels(axis.get_xticklabels(), rotation=90, ha="right")
        axis.set_yticklabels(axis.get_yticklabels(), rotation=0, ha="right")
        fig.savefig(os.path.join(self.results_dir, "{}.peak_location.per_sample.svg".format(self.name)), bbox_inches="tight")

        # Peak set across samples:
        # interval lengths
        fig, axis = plt.subplots()
        sns.distplot([interval.length for interval in self.sites if interval.length < 2000], bins=300, kde=False, ax=axis)
        axis.set_xlabel("peak width (bp)")
        axis.set_ylabel("frequency")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.lengths.svg".format(self.name)), bbox_inches="tight")

        # plot support
        fig, axis = plt.subplots()
        sns.distplot(self.support["support"], bins=40, ax=axis)
        axis.set_ylabel("frequency")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.support.svg".format(self.name)), bbox_inches="tight")

        # Plot distance to nearest TSS
        fig, axis = plt.subplots(2, 2, figsize=(4 * 2, 4 * 2), sharex=False, sharey=False)
        sns.distplot([x for x in self.closest_tss_distances if x < 1e5], bins=1000, kde=False, hist=True, ax=axis[0][0])
        sns.distplot([x for x in self.closest_tss_distances if x < 1e5], bins=1000, kde=False, hist=True, ax=axis[0][1])
        sns.distplot([x for x in self.closest_tss_distances if x < 1e6], bins=1000, kde=True, hist=False, ax=axis[1][0])
        sns.distplot([x for x in self.closest_tss_distances if x < 1e6], bins=1000, kde=True, hist=False, ax=axis[1][1])
        for ax in axis.flat:
            ax.set_xlabel("distance to nearest TSS (bp)")
            ax.set_ylabel("frequency")
        axis[0][1].set_yscale("log")
        axis[1][1].set_yscale("log")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.tss_distance.svg".format(self.name)), bbox_inches="tight")

        # Plot genomic regions
        # these are just long lists with genomic regions
        all_region_annotation = self.region_annotation['genomic_region'].apply(lambda x: pd.Series(x.split(","))).stack().reset_index(level=[1], drop=True)
        all_region_annotation.name = "foreground"
        all_region_annotation_b = self.region_annotation_b['genomic_region'].apply(lambda x: pd.Series(x.split(","))).stack().reset_index(level=[1], drop=True)
        all_region_annotation_b.name = "background"

        # count region frequency
        data = all_region_annotation.value_counts().sort_values(ascending=False)
        background = all_region_annotation_b.value_counts().sort_values(ascending=False)
        data = data.to_frame().join(background)
        data["fold_change"] = np.log2(data['foreground'] / data['background'])
        data.index.name = "region"

        # plot also % of genome space "used"
        self.region_annotation["length"] = self.region_annotation["end"] - self.region_annotation["start"]

        s = self.region_annotation.join(all_region_annotation).groupby("foreground")["length"].sum()
        s.name = "size"
        s = (s / pd.Series(genome_space)).dropna() * 100
        s.name = "percent_space"
        data = data.join(s)

        # plot together
        g = sns.FacetGrid(data=pd.melt(data.reset_index(), id_vars='region'), col="variable", col_wrap=2, sharex=True, sharey=False, size=4)
        g.map(sns.barplot, "region", "value")
        sns.despine(fig)
        g.savefig(os.path.join(self.results_dir, "{}.genomic_regions.svg".format(self.name)), bbox_inches="tight")

        # # Plot chromatin states
        all_chrom_state_annotation = self.chrom_state_annotation['chromatin_state'].apply(lambda x: pd.Series(x.split(","))).stack().reset_index(level=[1], drop=True)
        all_chrom_state_annotation.name = "foreground"
        all_chrom_state_annotation_b = self.chrom_state_annotation_b['chromatin_state'].apply(lambda x: pd.Series(x.split(","))).stack().reset_index(level=[1], drop=True)
        all_chrom_state_annotation_b.name = "background"

        # count region frequency
        data = all_chrom_state_annotation.value_counts().sort_values(ascending=False)
        background = all_chrom_state_annotation_b.value_counts().sort_values(ascending=False)
        data = data.to_frame().join(background)
        data["fold_change"] = np.log2(data['foreground'] / data['background'])
        data.index.name = "region"

        # plot also % of genome space "used"
        self.chrom_state_annotation["length"] = self.chrom_state_annotation["end"] - self.chrom_state_annotation["start"]

        s = self.chrom_state_annotation.join(all_chrom_state_annotation).groupby("foreground")["length"].sum()
        s.name = "size"
        s = (s / pd.Series(genome_space)).dropna() * 100
        s.name = "percent_space"
        data = data.join(s)

        # plot together
        g = sns.FacetGrid(data=pd.melt(data.reset_index(), id_vars='region'), col="variable", col_wrap=2, sharex=True, sharey=False, size=4)
        g.map(sns.barplot, "region", "value")
        sns.despine(fig)
        g.savefig(os.path.join(self.results_dir, "{}.chromatin_states.svg".format(self.name)), bbox_inches="tight")

        # distribution of count attributes
        data = self.coverage_annotated.copy()

        fig, axis = plt.subplots(1)
        sns.distplot(data["mean"], rug=False, ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.mean.distplot.svg".format(self.name)), bbox_inches="tight")

        fig, axis = plt.subplots(1)
        sns.distplot(data["qv2"], rug=False, ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.qv2.distplot.svg".format(self.name)), bbox_inches="tight")

        fig, axis = plt.subplots(1)
        sns.distplot(data["dispersion"], rug=False, ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.dispersion.distplot.svg".format(self.name)), bbox_inches="tight")

        fig, axis = plt.subplots(1)
        sns.distplot(data["support"], rug=False, ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "{}.support.distplot.svg".format(self.name)), bbox_inches="tight")

        # joint
        for metric in ["support", "variance", "std_deviation", "dispersion", "qv2", "amplitude"]:
            p = data[(data["mean"] > 0) & (data[metric] < np.percentile(data[metric], 99) * 3)]
            g = sns.jointplot(p["mean"], p[metric], s=2, alpha=0.2, rasterized=True)
            sns.despine(g.fig)
            g.fig.savefig(os.path.join(self.results_dir, "{}.mean_{}.svg".format(self.name, metric)), bbox_inches="tight", dpi=300)

    def plot_raw_coverage(self, samples=None, by_attribute=None):
        """
        Diagnostic plots on the Sample's signal.
        Provides plots with Samples grouped `by_attribute` if given (a string or a list of strings).

        Parameters
        ----------
        samples : list
            List of peppy.Samples objects to use for plotting.

        by_attribute : str, optional
            Attribute of samples to group by.
            Values will be aggregated across samples by that attribute.
        """
        # TODO: get matrix as input, move to graphics
        if samples is None:
            samples = self.samples

        if by_attribute is None:
            cov = pd.melt(
                np.log2(1 + self.coverage[[s.name for s in samples]]),
                var_name="Sample name", value_name="Raw counts (log2)"
            )
            fig, axis = plt.subplots(1, 1, figsize=(6, 1 * 4))
            sns.violinplot(
                "Raw counts (log2)", "Sample name",
                orient="horizontal", palette="tab20", data=cov, ax=axis)
            sns.despine(fig)
            fig.savefig(
                os.path.join(self.results_dir, self.name + ".raw_counts.violinplot.svg"),
                bbox_inches="tight")
        else:
            attrs = set([getattr(s, by_attribute) for s in samples])
            fig, axis = plt.subplots(len(attrs), 1, figsize=(8, len(attrs) * 6))
            for i, attr in enumerate(attrs):
                _LOGGER.info(attr)
                cov = pd.melt(
                    np.log2(1 + self.coverage[[s.name for s in samples if getattr(s, by_attribute) == attr]]),
                    var_name="Sample name", value_name="Raw counts (log2)"
                )
                sns.violinplot(
                    "Raw counts (log2)", "Sample name",
                    orient="horizontal", palette="tab20", data=cov, ax=axis[i])
                axis[i].set_title(attr)
                axis[i].set_xticklabels(axis[i].get_xticklabels(), rotation=90)
            sns.despine(fig)
            fig.savefig(
                os.path.join(self.results_dir, self.name + ".raw_counts.violinplot.by_{}.svg"
                             .format(by_attribute)),
                bbox_inches="tight")

    def plot_coverage(self):

        # TODO: add plots for overal genome
        # TODO: add raw counts too

        data = self.accessibility.copy()
        # (rewrite to avoid putting them there in the first place)
        variables = ['gene_name', 'genomic_region', 'chromatin_state']

        for variable in variables:
            d = data[variable].str.split(',').apply(pd.Series).stack()  # separate comma-delimited fields
            d.index = d.index.droplevel(1)  # returned a multiindex Series, so get rid of second index level (first is from original row)
            data = data.drop([variable], axis=1)  # drop original column so there are no conflicts
            d.name = variable
            data = data.join(d)  # joins on index

        variables = [
            'chrom', 'start', 'end',
            'ensembl_transcript_id', 'distance', 'ensembl_gene_id', 'support',
            'mean', 'variance', 'std_deviation', 'dispersion', 'qv2',
            'amplitude', 'gene_name', 'genomic_region', 'chromatin_state']
        # Plot
        data_melted = pd.melt(
            data,
            id_vars=variables, var_name="sample", value_name="norm_counts")

        # transform dispersion
        data_melted['dispersion'] = np.log2(1 + data_melted['dispersion'])

        # Together in same violin plot
        fig, axis = plt.subplots(1)
        sns.violinplot("genomic_region", "norm_counts", data=data_melted, ax=axis)
        fig.savefig(os.path.join(self.results_dir, self.name + ".norm_counts.per_genomic_region.violinplot.svg"), bbox_inches="tight")

        # dispersion
        fig, axis = plt.subplots(1)
        sns.violinplot("genomic_region", "dispersion", data=data_melted, ax=axis)
        fig.savefig(os.path.join(self.results_dir, self.name + ".norm_counts.dispersion.per_genomic_region.violinplot.svg"), bbox_inches="tight")

        # dispersion
        fig, axis = plt.subplots(1)
        sns.violinplot("genomic_region", "qv2", data=data_melted, ax=axis)
        fig.savefig(os.path.join(self.results_dir, self.name + ".norm_counts.qv2.per_genomic_region.violinplot.svg"), bbox_inches="tight")

        fig, axis = plt.subplots(1)
        sns.violinplot("chromatin_state", "norm_counts", data=data_melted, ax=axis)
        fig.savefig(os.path.join(self.results_dir, self.name + ".norm_counts.chromatin_state.violinplot.svg"), bbox_inches="tight")

        fig, axis = plt.subplots(1)
        sns.violinplot("chromatin_state", "dispersion", data=data_melted, ax=axis)
        fig.savefig(os.path.join(self.results_dir, self.name + ".norm_counts.dispersion.chromatin_state.violinplot.svg"), bbox_inches="tight")

        fig, axis = plt.subplots(1)
        sns.violinplot("chromatin_state", "qv2", data=data_melted, ax=axis)
        fig.savefig(os.path.join(self.results_dir, self.name + ".norm_counts.qv2.chromatin_state.violinplot.svg"), bbox_inches="tight")

        # separated by variable in one grid
        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "mean", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, self.name + ".norm_counts.mean.per_genomic_region.distplot.svg"), bbox_inches="tight")

        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "dispersion", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, self.name + ".norm_counts.dispersion.per_genomic_region.distplot.svg"), bbox_inches="tight")

        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "qv2", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, self.name + ".norm_counts.qv2.per_genomic_region.distplot.svg"), bbox_inches="tight")

        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "support", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, self.name + ".norm_counts.support.per_genomic_region.distplot.svg"), bbox_inches="tight")

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "mean", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, self.name + ".norm_counts.mean.chromatin_state.distplot.svg"), bbox_inches="tight")

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "dispersion", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, self.name + ".norm_counts.dispersion.chromatin_state.distplot.svg"), bbox_inches="tight")

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "qv2", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, self.name + ".norm_counts.qv2.chromatin_state.distplot.svg"), bbox_inches="tight")

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "support", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, self.name + ".norm_counts.support.chromatin_state.distplot.svg"), bbox_inches="tight")
        plt.close("all")

    def plot_variance(self, samples):

        g = sns.jointplot('mean', "dispersion", data=self.accessibility, kind="kde")
        g.fig.savefig(os.path.join(self.results_dir, self.name + ".norm_counts_per_sample.dispersion.svg"), bbox_inches="tight")

        g = sns.jointplot('mean', "qv2", data=self.accessibility)
        g.fig.savefig(os.path.join(self.results_dir, self.name + ".norm_counts_per_sample.qv2_vs_mean.svg"), bbox_inches="tight")

        g = sns.jointplot('support', "qv2", data=self.accessibility)
        g.fig.savefig(os.path.join(self.results_dir, self.name + ".norm_counts_per_sample.support_vs_qv2.svg"), bbox_inches="tight")

        # Filter out regions which the maximum across all samples is below a treshold
        filtered = self.accessibility[self.accessibility[[sample.name for sample in samples]].max(axis=1) > 3]

        sns.jointplot('mean', "dispersion", data=filtered)
        plt.savefig(os.path.join(self.results_dir, self.name + ".norm_counts_per_sample.dispersion.filtered.svg"), bbox_inches="tight")
        plt.close('all')
        sns.jointplot('mean', "qv2", data=filtered)
        plt.savefig(os.path.join(self.results_dir, self.name + ".norm_counts_per_sample.support_vs_qv2.filtered.svg"), bbox_inches="tight")

    def region_context_enrichment(
            self, regions,
            steps=['genomic_region', 'chromatin_state'],
            background="region_set",
            prefix="region_type_enrichment", output_dir="{results_dir}"):
        """
        Characterize a subset of the regions (e.g. differential regions) in terms of their genomic context.

        Parameters
        ----------
        regions : {list, pandas.DataFrame, pandas.Index}
            Subset of regions of interest to analysis.
            Must be a subset of the universe (i.e. `sites` attribute).

        steps : list, optional
            Steps of enrichment to perform.
            Defaults to all available: ['genomic_region', 'chromatin_state']

        background : str, optional
            Which set to consider as backgroud.
            Options are:
                region_set: the consensus region_set of the analysis
                genome: a randomized set of size as region_set across the genome

        prefix : str, optional
            Prefix for saved files.
            Default is `region_type_enrichment`.

        output_dir : str, optional
            Directory to write results to.

        Returns
        ---------
        pandas.DataFrame
            Enrichment results
        """
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

        options = ['region_set', "genome"]
        if background not in options:
            msg = "Option `background` must be one of '{}'.".format("', '".join(options))
            raise ValueError(msg)

        # compare genomic regions and chromatin_states
        enr = list()
        msg = "'{}' step selected, but analysis does not have '{}'."
        msg2 = "'genome' selected, but analysis does not have '{}'."
        for step, matrix, matrix_b in [
            ("genomic_region", "region_annotation_mapping", "region_annotation_b_mapping"),
            ("chromatin_state", "chrom_state_annotation_mapping", "chrom_state_annotation_b_mapping"),
        ]:
            if step not in steps:
                continue
            if not hasattr(self, matrix):
                _LOGGER.warning(msg.format(step, matrix))
                continue

            _LOGGER.debug("Getting enrichment of regions in '{}'.".format(step))

            # Count foreground occurences
            annot = getattr(self, matrix)
            res = annot.loc[set(regions), step].value_counts().to_frame('foreground')

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
                msg3 = ("Foreground regions contains type of {} not in background: {}"
                        .format(step, "', '".join(m)))
                msg3 += " Continuing without those."
                _LOGGER.warning(msg3)
            res = res.reindex(res_b.index)

            # # join
            res = res.join(res_b, how="outer").fillna(0).astype(int)

            # Calculate log fold enrichment:
            # # normalize to total:
            res.loc[:, 'foreground_fraction'] = res['foreground'] / res['foreground'].sum()
            res.loc[:, 'universe_fraction'] = res['universe'] / res['universe'].sum()
            res.loc[:, 'log2_fold_change'] = np.log2(res['foreground_fraction'] / res['universe_fraction'])
            # Calculate overlap p-value:
            for feature in res['foreground'].index:
                a = res.loc[feature, 'foreground']
                b = res.loc[:, 'foreground'].drop(feature).sum()
                c = annot.loc[(~annot.index.isin(regions)), step].value_counts()[feature]
                d = annot.loc[(~annot.index.isin(regions)), step].value_counts().drop(feature).sum()
                res.loc[feature, 'odds_ratio'], res.loc[feature, 'p_value'] = fisher_exact(
                    [[a, c], [b, d]], alternative="two-sided")
            res.loc[:, 'log2_odds_ratio'] = np.log2(res['odds_ratio'])
            res.loc[:, '-log10(p-value)'] = log_pvalues(res['p_value'])
            res.loc[:, 'region_type'] = step
            # Append
            enr.append(res)

        # save
        enr = pd.concat(enr)
        enr.index.name = "region"
        enr.to_csv(os.path.join(output_dir, prefix + ".csv"), index=True)
        return enr

    def characterize_regions_function(
            self, differential, output_dir, prefix, universe_file=None,
            run=True, genome=None,
            steps=['region', 'lola', 'meme', 'homer', 'enrichr']):
        """
        Performs a range of functional enrichments of a set of regions given in `differential`
        (a dataframe which is typically a subset of an annotated coverage dataframe).
        Will extract regions, their underlying sequence, associated genes, perform enrichment of
        genomic regions, chromatin states against a background, motif enrichment,
        location overlap analysis (LOLA), and gene set enrichment (using the Enrichr API).

        This requires several programs and R libraries:
            - MEME suite (AME)
            - HOMER suite (findMotifsGenome.pl)
            - LOLA (R library)
        Additionally, some genome-specific databases are needed to run these programs.

        Parameters
        ----------
        differential : pandas.DataFrame
            Results of differential analysis for a given comparison of interest.

        output_dir : str
            Directory to output results to.

        prefix : str
            Prefix to use for output files.

        universe_file : str, optional
            Path to BED file with set of universe regions where differential were selected from.
            Default is analysis.sites.

        run : bool, optional
            Whether to run enrichment commands now or to simply prepare the input files for it.
            Default is True.

        genome : str, optional
            Genome assembly of analysis.
            Default is analysis' genome assembly.

        steps : list, optional
            Which steps of the analysis to perform.
            Default is all: ['region', 'lola', 'meme', 'homer', 'enrichr']
        """
        from ngs_toolkit.general import get_genome_reference
        import warnings
        # use all sites as universe
        if universe_file is None:
            try:
                universe_file = getattr(self, "sites").fn
                _LOGGER.info("Using default background region set from 'analysis.sites': {}'."
                             .format(universe_file))
            except AttributeError as e:
                _LOGGER.error("Background region set 'analysis.sites' is not set! Cannot run LOLA!")
                raise e

        if genome is None:
            genome = self.genome

        # make output dirs
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # save to bed
        bed_file = os.path.join(output_dir, "{}_regions.bed".format(prefix))
        differential[['chrom', 'start', 'end']].to_csv(
            bed_file, sep="\t", header=False, index=False)
        # save as tsv
        tsv_file = os.path.join(output_dir, "{}_regions.tsv".format(prefix))
        differential[['chrom', 'start', 'end']].reset_index().to_csv(
            tsv_file, sep="\t", header=False, index=False)

        # export gene names
        clean_gene = differential['gene_name'].str.split(",").apply(pd.Series, 1).stack().drop_duplicates()
        clean_gene = clean_gene[~clean_gene.isin(['.', 'nan', ''])]
        clean_gene.to_csv(
                os.path.join(output_dir, "{}_genes.symbols.txt".format(prefix)),
                header=False, index=False)
        if "ensembl_gene_id" in differential.columns:
            # export ensembl gene names
            clean = differential['ensembl_gene_id'].str.split(",").apply(pd.Series, 1).stack().drop_duplicates()
            clean.to_csv(
                os.path.join(output_dir, "{}_genes.ensembl.txt".format(prefix)),
                header=False, index=False)

        # export gene symbols with scaled absolute fold change
        if "log2FoldChange" in differential.columns:
            differential["score"] = standard_score(abs(differential["log2FoldChange"]))
            differential["abs_fc"] = abs(differential["log2FoldChange"])

            d = differential[['gene_name', 'score']].sort_values('score', ascending=False)

            # split gene names from score if a reg.element was assigned to more than one gene
            a = d['gene_name'].str.split(",").apply(pd.Series, 1).stack()
            a.index = a.index.droplevel(1)
            a.name = 'gene_name'
            d = d[['score']].join(a)
            # reduce various ranks to mean per gene
            d = d.groupby('gene_name').mean().reset_index()
            d.to_csv(
                os.path.join(output_dir, "{}_genes.symbols.score.csv".format(prefix)),
                index=False)

        # get fasta file with sequence underlying region
        if ('meme' in steps) or ('homer' in steps):
            hint = " Will not do motif enrichment analysis."
            fasta_file = os.path.join(output_dir, "{}_regions.fa".format(prefix))
            okay = False
            for f in ['fasta', '2bit']:
                with warnings.catch_warnings():
                    try:
                        genome_file = get_genome_reference(self.organism, file_format=f, overwrite=False)
                        okay = True
                    except ValueError:
                        pass
            if not okay:
                reason = "Could not get genome sequence file in either FASTA or 2bit format."
                _LOGGER.warning(reason + hint)
            else:
                try:
                    bed_to_fasta(input_bed=bed_file, output_fasta=fasta_file, genome_file=genome_file)
                except EnvironmentError:
                    reason = "Could not get FASTA sequence."
                    _LOGGER.warning(reason + hint)

        if not run:
            return

        if 'region' in steps:
            self.region_context_enrichment(
                differential, output_dir=output_dir)

        # MEME
        if "meme" in steps:
            _LOGGER.info("Running MEME-AME for '{}'".format(prefix))
            omap = {"hg38": "human", "hg19": "human", "mm10": "mouse"}
            meme_ame(fasta_file, output_dir, organism=omap[genome])

        # HOMER
        if "homer" in steps:
            _LOGGER.info("Running HOMER for '{}'".format(prefix))
            homer_motifs(bed_file, output_dir, genome=genome)

        # LOLA
        if 'lola' in steps:
            _LOGGER.info("Running LOLA for '{}'".format(prefix))
            try:
                lola(bed_file, universe_file, output_dir, genome=genome)
            except:
                _LOGGER.error("LOLA analysis for '{}' failed!".format(prefix))

        # Enrichr
        if 'enrichr' in steps:
            _LOGGER.info("Running Enrichr for '{}'".format(prefix))
            results = enrichr(clean_gene.to_frame(name="gene_name"))
            results.to_csv(
                os.path.join(output_dir, "{}_regions.enrichr.csv".format(prefix)),
                index=False, encoding='utf-8')
