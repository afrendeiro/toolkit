#!/usr/bin/env python


import os

import numpy as np
import pandas as pd
import pybedtools

from ngs_toolkit.atacseq import ATACSeqAnalysis
from . import _LOGGER


class ChIPSeqAnalysis(ATACSeqAnalysis):
    """
    Class to model analysis of ChIP-seq data.
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
            **kwargs):
        super(ChIPSeqAnalysis, self).__init__(
            name=name,
            data_dir=data_dir,
            results_dir=results_dir,
            pickle_file=pickle_file,
            samples=samples,
            prj=prj,
            from_pickle=from_pickle,
            **kwargs)

        self.data_type = self.__data_type__ = "ChIP-seq"
        self._var_names = "region"
        self._quantity = "binding"
        self._norm_units = "RPM"
        self._raw_matrix_name = "coverage"
        self._norm_matrix_name = "coverage_qnorm"  # or coverage_gc_corrected
        # TODO: have normalization method reset the above value, with info to user
        # or just overwrite
        self._annot_matrix_name = "binding"

    def call_peaks_from_comparisons(
            self, comparison_table, output_dir="{results_dir}/chipseq_peaks",
            permissive=True, overwrite=True, as_jobs=True):
        """
        Call peaks for ChIP-seq samples using an annotation of which samples belong in each comparison and which
        samples represent signal or background.

        :param pandas.DataFrame comparison_table: Comparison table with the following required columns:
            "comparison_name", "sample_name", "comparison_side", "sample_group".
        :param str output_dir: Parent directory where peaks will be created. Will be created if does not exist.
        :param bool permissive: If incomplete/incoherent comparisons should be skipped or an error should be thrown.
        :raises ValueError: Will be raised if not `permissive` and incomplete/incoherent comparisons are detected.
        """
        from tqdm import tqdm

        req_columns = ["comparison_name", "sample_name", "comparison_side", "sample_group"]
        assert all([col in comparison_table.columns for col in req_columns]), "Comparison table is missing some of the following columns: '{}'.".format(",".join(req_columns))

        # Complement default `output_dir`
        if "{results_dir}" in output_dir:
            output_dir = os.path.abspath(output_dir.format(results_dir=self.results_dir))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # For each comparison
        comps = comparison_table['comparison_name'].drop_duplicates().sort_values()
        for comparison in tqdm(comps, total=len(comps), desc="Comparison"):
            # If there aren't two sides to each comparison, skip it or throw error
            if len(set(comparison_table[
                (comparison_table["comparison_name"] == comparison)
            ]["comparison_side"].tolist())) != 2:
                error = "Comparison '{}' does not contain two sides.".format(comparison)
                if permissive:
                    _LOGGER.warn(error)
                    continue
                else:
                    raise ValueError(error)

            # Get the sample names of samples in each side
            pos_names = comparison_table[
                (comparison_table["comparison_name"] == comparison) &
                (comparison_table["comparison_side"] == 1)
            ]["sample_name"].tolist()
            neg_names = comparison_table[
                (comparison_table["comparison_name"] == comparison) &
                (comparison_table["comparison_side"] < 1)
            ]["sample_name"].tolist()

            # Now get the actual samples
            signal_samples = [s for s in self.samples if s.name in pos_names]
            control_samples = [s for s in self.samples if s.name in neg_names]

            if len(signal_samples) == 0 or len(control_samples) == 0:
                error = "Comparison side for '{}' comparison does not contain samples.".format(comparison)
                if permissive:
                    _LOGGER.warn(error)
                    continue
                else:
                    raise ValueError(error)

            _LOGGER.info("Doing comparison '{}' with positive samples '{}' and background samples '{}'".format(
                comparison, [s.name for s in signal_samples], [s.name for s in control_samples]
            ))
            # Call peaks
            cmds = list()
            if overwrite:
                cmds += [macs2_call_chipseq_peak(
                    signal_samples, control_samples, output_dir=output_dir, name=comparison, as_job=as_jobs)]
                cmds += [homer_call_chipseq_peak_job(
                    signal_samples, control_samples, output_dir=output_dir, name=comparison, as_job=as_jobs)]
            else:
                if not os.path.exists(os.path.join(output_dir, comparison, comparison + "_peaks.narrowPeak")):
                    cmds += [macs2_call_chipseq_peak(
                        signal_samples, control_samples, output_dir=output_dir, name=comparison, as_job=as_jobs)]
                if not os.path.exists(os.path.join(output_dir, comparison, comparison + "_homer_peaks.narrowPeak")):
                    cmds += [homer_call_chipseq_peak_job(
                        signal_samples, control_samples, output_dir=output_dir, name=comparison, as_job=as_jobs)]
                else:
                    _LOGGER.warn("Peak files for comparison '{}' already exist. Skipping.".format(comparison))

            if not as_jobs:
                for cmd in cmds:
                    _LOGGER.info("Calling peaks for comparison '{}' with command: '{}'.\n".format(comparison, cmd))
                    os.system(cmd)

    def filter_peaks(self, comparison_table, filter_bed="blacklist.mm10_liftOver.bed", peaks_dir="{results_dir}/chipseq_peaks"):
        """
        Filter BED file for entries that overlap another BED file.

        :param str input_bed: BED file to filter.
        :param str filter_bed: BED file with entries to filter from input_bed.
        :param str output_bed: Output BED file.
        """
        import pybedtools
        from ngs_toolkit.chipseq import homer_peaks_to_bed

        if filter_bed == "blacklist.mm10_liftOver.bed":
            _LOGGER.warn("Using blacklist features of mm10 genome!")

        # Complement default `peaks_dir`
        if "{results_dir}" in peaks_dir:
            peaks_dir = os.path.abspath(peaks_dir.format(results_dir=self.results_dir))
        if not os.path.exists(peaks_dir):
            os.makedirs(peaks_dir)

        def _filter(input_bed, filter_bed, output_bed):
            """
            Filter BED file for entries that overlap another BED file.

            :param str input_bed: BED file to filter.
            :param str filter_bed: BED file with entries to filter from input_bed.
            :param str output_bed: Output BED file.
            """
            (
                pybedtools.BedTool(input_bed)
                .intersect(pybedtools.BedTool(filter_bed), v=True)
                .saveas(output_bed))

        for comparison in comparison_table['comparison_name'].unique():
            # MACS2
            prefix = os.path.join(peaks_dir, comparison, comparison + "_peaks")
            _filter(prefix + ".narrowPeak", filter_bed, prefix + ".filtered.bed")
            # HOMER
            prefix = os.path.join(peaks_dir, comparison, comparison + "_homer_peaks")
            homer_peaks_to_bed(prefix + ".factor.narrowPeak", prefix + ".factor.bed")
            _filter(prefix + ".factor.bed", filter_bed, prefix + ".factor.filtered.bed")
            homer_peaks_to_bed(prefix + ".histone.narrowPeak", prefix + ".histone.bed")
            _filter(prefix + ".histone.bed", filter_bed, prefix + ".histone.filtered.bed")

    def summarize_peaks_from_comparisons(self, comparison_table, output_dir="{results_dir}/chipseq_peaks", filtered=True, permissive=True):
        """
        Call peaks for ChIP-seq samples using an annotation of which samples belong in each comparison and which
        samples represent signal or background.

        :param pandas.DataFrame comparison_table: Comparison table with the following required columns:
            "comparison_name", "sample_name", "comparison_side", "sample_group".
        :param str output_dir: Parent directory where peaks will be created. Will be created if does not exist.
        :param bool permissive: If incomplete/incoherent comparisons should be skipped or an error should be thrown.
        :raises ValueError: Will be raised if not `permissive` and incomplete/incoherent comparisons are detected.
        """
        req_columns = ["comparison_name", "sample_name", "comparison_side", "sample_group"]
        assert all([col in comparison_table.columns for col in req_columns]), "Comparison table is missing some of the following columns: '{}'.".format(",".join(req_columns))

        # Complement default `output_dir`
        if "{results_dir}" in output_dir:
            output_dir = os.path.abspath(output_dir.format(results_dir=self.results_dir))

        # For each comparison, count called peaks
        peak_counts = pd.DataFrame()
        for comparison in comparison_table['comparison_name'].drop_duplicates().sort_values():
            _LOGGER.info(comparison)
            ending = ".filtered.bed" if filtered else ".narrowPeak"
            for peak_type, file in [
                    ("macs", os.path.join(output_dir, comparison, comparison + "_peaks" + ending)),
                    ("homer_factor", os.path.join(output_dir, comparison, comparison + "_homer_peaks.factor" + ending)),
                    ("homer_histone", os.path.join(output_dir, comparison, comparison + "_homer_peaks.histone" + ending))]:
                error = "Peak files for comparison '{}' with '{}' parameters don't exist.".format(comparison, peak_type)

                if "homer" in peak_type and not filtered:
                    try:
                        homer_peaks_to_bed(file, file.replace("narrowPeak", "bed"))
                    except IOError:
                        if permissive:
                            _LOGGER.warn(error)
                            peak_counts = peak_counts.append(
                                pd.Series([comparison, peak_type, np.nan]), ignore_index=True)
                            continue
                        else:
                            raise
                    except pd.errors.EmptyDataError:
                        peak_counts = peak_counts.append(
                                pd.Series([comparison, peak_type, 0.]), ignore_index=True)

                    file = file.replace("narrowPeak", "bed")
                try:
                    df = pd.read_csv(file, sep="\t")
                except IOError:
                    if permissive:
                        _LOGGER.warn(error)
                        peak_counts = peak_counts.append(
                            pd.Series([comparison, peak_type, np.nan]), ignore_index=True)
                        continue
                    else:
                        raise
                except pd.errors.EmptyDataError:
                    peak_counts = peak_counts.append(
                                pd.Series([comparison, peak_type, 0.]), ignore_index=True)

                peak_counts = peak_counts.append(pd.Series([comparison, peak_type, df.shape[0]]), ignore_index=True)
        peak_counts.columns = ['comparison_name', 'peak_type', 'peak_counts']

        return peak_counts  # .fillna(0)

    def get_consensus_sites(
            self, comparison_table, peak_dir="{results_dir}/chipseq_peaks", region_type="peaks", extension=250,
            blacklist_bed="wgEncodeDacMapabilityConsensusExcludable.bed"):
        """
        Get consensus (union) of enriched sites (peaks) across all comparisons.
        If `region_type` --> "summits, regions used will be peak summits which will be extended by `extension`
        before union. Otherwise sample peaks will be used with no modification.

        `blacklist_bed` is a 3 column BED file with genomic positions to exclude from consensus peak set.
        """
        import re
        from tqdm import tqdm
        import pybedtools

        # Complement default `peak_dir`
        if "{results_dir}" in peak_dir:
            peak_dir = os.path.abspath(peak_dir.format(results_dir=self.results_dir))

        first = True
        comps = comparison_table["comparison_name"].drop_duplicates()
        for comparison in tqdm(comps, total=len(comps), desc="Comparison"):
            peak_files = [
                os.path.join(peak_dir, comparison, comparison + "_peaks.narrowPeak"),
                os.path.join(peak_dir, comparison, comparison + "_homer_peaks.factor.bed"),
                os.path.join(peak_dir, comparison, comparison + "_homer_peaks.histone.bed")]
            for peak_file in peak_files:
                genome = comparison_table.loc[comparison_table["comparison_name"] == comparison, "comparison_genome"].drop_duplicates().squeeze()

                assert type(genome) is str, "Could not determine genome of comparison '{}'.".format(comparison)

                # Get peaks
                if region_type == "summits":
                    try:
                        f = re.sub("_peaks.narrowPeak", "_summits.bed", peak_file)
                        peaks = pybedtools.BedTool(f).slop(b=extension, genome=genome)
                    except ValueError:
                        _LOGGER.warn("Summits for comparison {} ({}) not found!".format(comparison, f))
                        continue
                else:
                    try:
                        peaks = pybedtools.BedTool(peak_file)
                    except ValueError:
                        _LOGGER.warn("Peaks for comparison {} ({}) not found!".format(comparison, peak_file))
                        continue
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
        blacklist = pybedtools.BedTool(os.path.join(self.data_dir, "external", blacklist_bed))
        # remove chrM peaks and save
        sites.intersect(v=True, b=blacklist).filter(lambda x: x.chrom != 'chrM').saveas(os.path.join(self.results_dir, self.name + "_peak_set.bed"))

        # Read up again
        self.sites = pybedtools.BedTool(os.path.join(self.results_dir, self.name + "_peak_set.bed"))

    def set_consensus_sites(self, bed_file, overwrite=True):
        """
        Set consensus (union) sites across samples.
        Will be stored in a `sites` attribute.
        """
        self.sites = pybedtools.BedTool(bed_file)
        if overwrite:
            self.sites.saveas(os.path.join(self.results_dir, self.name + "_peak_set.bed"))

    def calculate_peak_support(self, comparison_table, peak_dir="{results_dir}/chipseq_peaks"):
        """
        Calculate a measure of support for each region in peak set
        (i.e. ratio of samples containing a peak overlapping region in union set of peaks).
        """
        from tqdm import tqdm

        if "{results_dir}" in peak_dir:
            peak_dir = os.path.abspath(peak_dir.format(results_dir=self.results_dir))

        # get index
        index = self.sites.to_dataframe()
        index = index['chrom'] + ":" + index['start'].astype(str) + "-" + index['end'].astype(str)

        # calculate support (number of samples overlaping each merged peak)
        comps = comparison_table["comparison_name"].drop_duplicates()
        support = pd.DataFrame(index=index)
        for i, comparison in tqdm(enumerate(comps), total=comps.shape[0], desc="Comparison"):
            peak_files = [
                ("MACS", os.path.join(peak_dir, comparison, comparison + "_peaks.narrowPeak")),
                ("HOMER_factor", os.path.join(peak_dir, comparison, comparison + "_homer_peaks.factor.bed")),
                ("HOMER_histone", os.path.join(peak_dir, comparison, comparison + "_homer_peaks.histone.bed"))]
            for peak_type, peak_file in peak_files:
                try:
                    sample_support = self.sites.intersect(peak_file, wa=True, c=True).to_dataframe()
                except:
                    continue
                sample_support.index = index
                support[(comparison, peak_type)] = sample_support.iloc[:, 3]

        # Make multiindex labeling comparisons and peak type
        support.columns = pd.MultiIndex.from_tuples(support.columns, names=["comparison", "peak_type"])
        support.to_csv(os.path.join(self.results_dir, self.name + "_peaks.binary_overlap_support.csv"), index=True)

        # get % of total consensus regions found per sample
        # m = (
        #     pd.melt(support, ["chrom", "start", "end"], var_name="sample_name")
        #     .groupby("sample_name")
        #     .apply(lambda x: len(x[x["value"] == 1])))

        # divide sum (of unique overlaps) by total to get support value between 0 and 1
        support["support"] = support.astype(bool).astype(int).sum(axis=1) / float(support.shape[1])

        # save
        support.to_csv(os.path.join(self.results_dir, self.name + "_peaks.support.csv"), index=True)

        self.support = support

    def get_supported_peaks(self, comparisons):
        """
        Get mask of sites with 0 support in the given samples.
        Requires support matrix produced by `ngs_toolkit.atacseq.ATACSeqAnalysis.calculate_peak_support`.

        :param list samples: Iterable of peppy.Sample objects to restrict to.
        :returns pd.Series: Boolean Pandas Series with sites with at least one of the \
                            given samples having a peak called.
        """
        return self.support[[c for c in comparisons]].sum(1) != 0


def macs2_call_chipseq_peak(signal_samples, control_samples, output_dir, name, as_job=True):
    """
    Call ChIP-seq peaks with MACS2 in a slurm job.

    :param list signal_samples: Signal Sample objects.
    :param list control_samples: Background Sample objects.
    :param list output_dir: Parent directory where MACS2 outputs will be stored.
    :param str name: Name of the MACS2 comparison being performed.
    :param bool as_job: Whether to submit a SLURM job or to return a string with the runnable.
    """
    output_path = os.path.join(output_dir, name)
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    runnable = (
        """macs2 callpeak -t {0} -c {1} -n {2} --outdir {3}"""
        .format(
            " ".join([s.filtered for s in signal_samples]),
            " ".join([s.filtered for s in control_samples]), name, output_path))

    if as_job:
        from pypiper.ngstk import NGSTk
        import textwrap
        tk = NGSTk()
        job_name = "macs2_{}".format(name)
        cmd = tk.slurm_header(
            job_name,
            os.path.join(output_path, job_name + ".log"),
            cpus_per_task=4)
        cmd += "\t\t" + runnable
        cmd += "\t\t" + tk.slurm_footer() + "\n"
        job_file = os.path.join(output_path, name + ".macs2.sh")
        with open(job_file, "w") as handle:
            handle.write(textwrap.dedent(cmd))
        tk.slurm_submit_job(job_file)
    else:
        return runnable


def homer_call_chipseq_peak_job(signal_samples, control_samples, output_dir, name, as_job=True):
    """
    Call ChIP-seq peaks with MACS2 in a slurm job.

    :param list signal_samples: Signal Sample objects.
    :param list control_samples: Background Sample objects.
    :param list output_dir: Parent directory where MACS2 outputs will be stored.
    :param str name: Name of the MACS2 comparison being performed.
    """
    output_path = os.path.join(output_dir, name)
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    # make tag directory for the signal and background samples separately
    signal_tag_directory = os.path.join(output_dir, "homer_tag_dir_" + name + "_signal")
    fs = " ".join([s.filtered for s in signal_samples])
    runnable = """makeTagDirectory {0} {1}\n""".format(signal_tag_directory, fs)
    background_tag_directory = os.path.join(output_dir, "homer_tag_dir_" + name + "_background")
    fs = " ".join([s.filtered for s in control_samples])
    runnable += """makeTagDirectory {0} {1}\n""".format(background_tag_directory, fs)

    # call peaks
    output_file = os.path.join(output_dir, name, name + "_homer_peaks.factor.narrowPeak")
    runnable += """findPeaks {signal} -style factor -o {output_file} -i {background}\n""".format(
        output_file=output_file, background=background_tag_directory, signal=signal_tag_directory)
    output_file = os.path.join(output_dir, name, name + "_homer_peaks.histone.narrowPeak")
    runnable += """findPeaks {signal} -style histone -o {output_file} -i {background}\n""".format(
        output_file=output_file, background=background_tag_directory, signal=signal_tag_directory)

    if as_job:
        from pypiper.ngstk import NGSTk
        import textwrap
        tk = NGSTk()
        job_name = "homer_findPeaks_{}".format(name)
        cmd = tk.slurm_header(
            job_name,
            os.path.join(output_path, job_name + ".log"),
            cpus_per_task=4)
        cmd += runnable.replace("\n", "\t\t\n")
        cmd += "\t\t" + tk.slurm_footer() + "\n"
        job_file = os.path.join(output_path, name + ".homer.sh")
        with open(job_file, "w") as handle:
            handle.write(textwrap.dedent(cmd))
        tk.slurm_submit_job(job_file)
    else:
        return runnable


def homer_peaks_to_bed(homer_peaks, output_bed):
    """
    Convert HOMER peak calls to BED format.
    The fifth column (score) is the -log10(p-value) of the peak.

    :param str homer_peaks: HOMER output with peak calls.
    :param str output_bed: Output BED file.
    """
    df = pd.read_csv(homer_peaks, sep="\t", comment="#", header=None)
    df['-log_pvalue'] = (-np.log10(df.iloc[:, -2])).replace(pd.np.inf, 1000)
    df['name'] = df[1] + ":" + df[2].astype(str) + "-" + df[3].astype(str)

    (
        df[[1, 2, 3, 'name', '-log_pvalue']]
        .sort_values([1, 2, 3], ascending=True)
        .to_csv(output_bed, header=False, index=False, sep="\t"))
