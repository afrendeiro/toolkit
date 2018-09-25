#!/usr/bin/env python


import os
import pickle
from collections import Counter

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybedtools
import seaborn as sns

from .general import Analysis

# Set settings
pd.set_option("date_dayfirst", True)


class ATACSeqAnalysis(Analysis):
    """
    Class to model analysis of ATAC-seq data.
    Inherits from the `ngs_toolkit.general.Analysis` class.

    :param name: Name to give analysis object. Default is `atacseq_analysis`.
    :param list samples: Iterable of peppy.Sample objects use in analysis.
                         If not provided (`None` is passed) and `prj` is, it will \
                         default to all samples in the `prj` object (`samples` attribute).
    :param peppy.Project prj: Project to tie analysis to.
    :param str data_dir: Directory containing relevant data for analysis. Default is `data`.
    :param str results_dir: Directory to output relevant analysis results. Default is `results`.
    :param str pickle_file: File path to use to save serialized object in `pickle` format.
    :param bool from_pickle: If the analysis should be loaded from an \
                             existing pickle object. Default is `False.
    :returns: Self

    Remaining keyword arguments will be passed to parent class `ngs_toolkit.general.Analysis`.

    :Example:

    .. code-block:: python

        import os
        from peppy import Project
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
            **kwargs):
        super(ATACSeqAnalysis, self).__init__(
            name=name,
            data_dir=data_dir,
            results_dir=results_dir,
            pickle_file=pickle_file,
            samples=samples,
            prj=prj,
            from_pickle=from_pickle,
            **kwargs)

        self.data_type = "ATAC-seq"

    def load_data(self, output_mapping=None, only_these_keys=None, permissive=True, n_header_vars=5):
        """
        Load the output files of the major functions of the Analysis.

        :param dict output_mapping: Dictionary with "attribute name": "path prefix" to load the files.
        :param list | None only_these_keys: Iterable of analysis attributes to load up.
                                            Possible attributes:
                                                "sites", "support", "coverage", "coverage_qnorm",
                                                "nuc", "coverage_gc_corrected", "gene_annotation",
                                                "region_annotation", "region_annotation_b",
                                                "chrom_state_annotation", "chrom_state_annotation_b",
                                                "coverage_annotated", "accessibility".
        :param permissive bool: Whether an error should be thrown if reading a file causes IOError.
        :raises IOError: If not `permissive` and any of the files are not readable.
        :var pybedtools.BedTool sites: Sets a `sites` variable.
        :var pandas.DataFrame support: Sets a `support` variable.
        :var pandas.DataFrame coverage: Sets a `coverage` variable.
        :var pandas.DataFrame coverage_qnorm: Sets a `coverage_qnorm` variable.
        :var pandas.DataFrame nuc: Sets a `nuc` variable.
        :var pandas.DataFrame coverage_gc_corrected: Sets a `coverage_gc_corrected` variable.
        :var pandas.DataFrame gene_annotation: Sets a `gene_annotation` variable.
        :var pandas.DataFrame region_annotation: Sets a `region_annotation` variable.
        :var pandas.DataFrame region_annotation_b: Sets a `region_annotation_b` variable.
        :var pandas.DataFrame chrom_state_annotation: Sets a `chrom_state_annotation` variable.
        :var pandas.DataFrame chrom_state_annotation_b: Sets a `chrom_state_annotation_b` variable.
        :var pandas.DataFrame coverage_annotated: Sets a `coverage_annotated` variable.
        :var pandas.DataFrame accessibility: Sets a `accessibility` variable.
        :returns: `None`

        """
        if output_mapping is None:
            # TODO: get default mapping by having functions declare what they output
            # perhaps also with a dict of kwargs to pass to pandas.read_csv
            output_mapping = {
                "support": "_peaks.support.csv",
                "coverage": "_peaks.raw_coverage.csv",
                "coverage_qnorm": "_peaks.coverage_qnorm.csv",
                "nuc": "_peaks.gccontent_length.csv",
                "coverage_gc_corrected": "_peaks.coverage_gc_corrected.csv",
                "gene_annotation": "_peaks.gene_annotation.csv",
                "region_annotation": "_peaks.region_annotation.csv",
                "region_annotation_b": "_peaks.region_annotation_background.csv",
                "chrom_state_annotation": "_peaks.chromatin_state.csv",
                "chrom_state_annotation_b": "_peaks.chromatin_state_background.csv",
                "coverage_annotated": "_peaks.coverage_qnorm.annotated.csv",
            }
        if only_these_keys is None:
            only_these_keys = [
                "sites", "support", "coverage", "coverage_qnorm",
                "nuc", "coverage_gc_corrected", "gene_annotation",
                "region_annotation", "region_annotation_b",
                "chrom_state_annotation", "chrom_state_annotation_b",
                "coverage_annotated", "accessibility"]

        output_mapping = {k: v for k, v in output_mapping.items() if k in only_these_keys}

        for name, suffix in output_mapping.items():
            file = os.path.join(self.results_dir, self.name + suffix)
            if name in only_these_keys:
                print("Loading '{}' analysis attribute.".format(name))
                try:
                    setattr(self, name, pd.read_csv(file, index_col=0))
                except IOError as e:
                    if not permissive:
                        raise e
                    else:
                        print(e)

        # Special cases
        if "sites" in only_these_keys:
            file = os.path.join(self.results_dir, self.name + "_peak_set.bed")
            try:
                setattr(self, "sites", pybedtools.BedTool(file))
            except IOError as e:
                if not permissive:
                    raise e
                else:
                    print(e)

        if "accessibility" in only_these_keys:
            file = os.path.join(self.results_dir, self.name + ".accessibility.annotated_metadata.csv")
            try:
                # TODO: add configuration for custom number of header lines
                # TODO: or think of smart way to infer
                setattr(self, "accessibility", pd.read_csv(file, index_col=0, header=range(n_header_vars)))
            except IOError as e:
                if not permissive:
                    raise e
                else:
                    print(e)

    def get_consensus_sites(
            self, samples=None, region_type="summits", extension=250,
            blacklist_bed="wgEncodeDacMapabilityConsensusExcludable.bed",
            permissive=False):
        """
        Get consensus (union) of enriched sites (peaks) across samples.

        :param list samples: Iterable of peppy.Sample objects to restrict \
                             to. Must have a `peaks` attribute set. \
                             If not provided (`None` is passed) if will default to all samples \
                             in the analysis (`samples` attribute).
        :param str region_type: One of `summits` or `peaks`. The type of region \
                                to use to create the consensus region set. \
                                If `summits`, peak summits will be extended by `extension` \
                                before union. Otherwise sample peaks will be used with \
                                no modification.
        :param int extension: Amount to extend peaks summits by in both directions.
        :param str blacklist_bed: A (3 column) BED file with genomic positions to \
                                  exclude from consensus peak set.
        :raises ValueError: If not `permissive` and either the peak or summit file \
                            of a sample is not readable.
        :var pybedtools.BedTool sites: Sets a `sites` variable with consensus peak set.
        :returns: `None`
        """
        import re
        from tqdm import tqdm

        if samples is None:
            samples = self.samples

        for i, sample in tqdm(enumerate(samples), total=(len(samples))):
            # print(sample.name)
            if region_type == "summits":
                try:
                    f = re.sub("_peaks.narrowPeak", "_summits.bed", sample.peaks)
                    peaks = pybedtools.BedTool(f).slop(b=extension, genome=sample.genome)
                except ValueError as e:
                    if permissive:
                        print("Summits for sample {} ({}) not found!".format(sample, f))
                        continue
                    else:
                        raise e
            else:
                try:
                    peaks = pybedtools.BedTool(sample.peaks)
                except ValueError as e:
                    if permissive:
                        print("Peaks for sample {} ({}) not found!".format(sample, sample.peaks))
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

        # Merge overlaping peaks across samples
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
        Set consensus (union) sites across samples given a BED file.

        :param str bed_file: BED file to use as consensus sites.
        :param book overwrite: Whether a possibly existing file with a consensus peak set \
                               for this analysis should be overwritten in disk.
        :var pybedtools.BedTool sites: Sets a `sites` variable with consensus peak set.

        """
        self.sites = pybedtools.BedTool(bed_file)
        if overwrite:
            self.sites.saveas(os.path.join(self.results_dir, self.name + "_peak_set.bed"))

    def calculate_peak_support(self, samples=None, region_type="summits"):
        """
        Count number of called peaks per sample in the consensus region set.
        In addition calculate a measure of peak support (or ubiquitouness) by \
        observing the ratio of samples containing a peak overlapping each region.

        :param list samples: Iterable of peppy.Sample objects to restrict \
                             to. Must have a `peaks` attribute set. \
                             If not provided (`None` is passed) if will default to all samples \
                             in the analysis (`samples` attribute).
        :param str region_type: One of `summits` or `peaks`. The type of region \
                                to use to create the consensus region set. \
                                If `summits`, peak summits will be extended by `extension` \
                                before union. Otherwise sample peaks will be used with \
                                no modification.
        """
        import re
        from tqdm import tqdm

        if samples is None:
            samples = self.samples

        # calculate support (number of samples overlaping each merged peak)
        for i, sample in tqdm(enumerate(samples), total=len(samples)):
            # print(sample.name)
            if region_type == "summits":
                peaks = re.sub("_peaks.narrowPeak", "_summits.bed", sample.peaks)
            else:
                peaks = sample.peaks

            if i == 0:
                support = self.sites.intersect(peaks, wa=True, c=True)
            else:
                support = support.intersect(peaks, wa=True, c=True)

        try:
            support = support.to_dataframe()
        except:
            support.saveas("_tmp.peaks.bed")
            support = pd.read_csv("_tmp.peaks.bed", sep="\t", header=None)

        support.columns = ["chrom", "start", "end"] + [sample.name for sample in samples]
        support.index = support['chrom'] + ":" + support['start'].astype(str) + "-" + support['end'].astype(str)
        support.to_csv(os.path.join(self.results_dir, self.name + "_peaks.binary_overlap_support.csv"), index=True)

        # get % of total consensus regions found per sample
        # m = (
        #     pd.melt(support, ["chrom", "start", "end"], var_name="sample_name")
        #     .groupby("sample_name")
        #     .apply(lambda x: len(x[x["value"] == 1])))

        # divide sum (of unique overlaps) by total to get support value between 0 and 1
        support["support"] = support[[sample.name for sample in samples]].apply(lambda x: sum([i if i <= 1 else 1 for i in x]) / float(len(self.samples)), axis=1)
        # save
        support.to_csv(os.path.join(self.results_dir, self.name + "_peaks.support.csv"), index=True)

        self.support = support

    def get_supported_peaks(self, samples):
        """
        Get mask of sites with 0 support in the given samples.
        Requires support matrix produced by `ngs_toolkit.atacseq.ATACSeqAnalysis.calculate_peak_support`.

        :param list samples: Iterable of peppy.Sample objects to restrict to.
        :returns pd.Series: Boolean Pandas Series with sites with at least one of the \
                            given samples having a peak called.
        """
        if samples is None:
            samples = self.samples
        return self.support[[s.name for s in samples]].sum(1) != 0

    def measure_coverage(self, samples=None, sites=None, assign=True, save=True, output_file=None):
        """
        Measure read coverage (counts) of each sample in each region in consensus sites.
        Will try to use parallel computing using the `parmap` library.

        :param list samples: Iterable of peppy.Sample objects to restrict
                             to. Must have a `filtered` attribute set.
                             If not provided (`None` is passed) if will default to all samples
                             in the analysis (`samples` attribute).
        :param pybedtools.BedTool sites: Sites in the genome to quantify.
                                         If `None` the object's `sites` attribute will be used.
        :param pybedtools.BedTool assign: Whether to assign the matrix to an attribute of self
                                          named `coverage`.
        :param pybedtools.BedTool save: Whether to save to disk the coverage matrix with filename
                                        `output_file`.
        :param str output_file: A path to a CSV file with coverage output.
                                Default is `self.results_dir/self.name + "_peaks.raw_coverage.csv"`.
        :raises ValueError: If not `permissive` and either the peak or summit file
                            of a sample is not readable.
        :var pd.DataFrame coverage: Sets a `coverage` variable with DataFrame with read counts
                                    of shape (n_sites, m_samples).
        :returns pd.DataFrame: Pandas DataFrame with read counts of shape (n_sites, m_samples).
        """
        import multiprocessing
        import parmap
        from ngs_toolkit.general import count_reads_in_intervals

        if samples is None:
            samples = self.samples

        missing = [s for s in samples if not os.path.exists(s.filtered)]
        if len(missing) > 0:
            print("Samples have missing BAM file: {}".format("\n".join([s.name for s in missing])))
            samples = [s for s in samples if s not in missing]

        if sites is None:
            sites = self.sites
            if output_file is None:
                default = True
            else:
                default = False

        # Count reads with pysam
        # make strings with intervals
        if type(sites) is pybedtools.bedtool.BedTool:
            sites_str = [str(i.chrom) + ":" + str(i.start) + "-" + str(i.stop) for i in self.sites]
        elif type(sites) is pd.core.frame.DataFrame:
            sites_str = (sites.iloc[:, 0] + ":" + sites.iloc[:, 1].astype(str) + "-" + sites.iloc[:, 2].astype(str)).astype(str).tolist()
        elif type(sites) is str:
            sites_str = [str(i.chrom) + ":" + str(i.start) + "-" + str(i.stop) for i in pybedtools.bedtool.BedTool(sites)]

        # count, create dataframe
        coverage = pd.DataFrame(
            map(
                lambda x:
                    pd.Series(x),
                    parmap.map(
                        count_reads_in_intervals,
                        [sample.filtered for sample in samples],
                        sites_str,
                        parallel=True
                    )
            ),
            index=[sample.name for sample in samples]
        ).T

        # Add interval description to df
        ints = map(
            lambda x: (
                x.split(":")[0],
                x.split(":")[1].split("-")[0],
                x.split(":")[1].split("-")[1]
            ),
            coverage.index
        )
        coverage["chrom"] = [x[0] for x in ints]
        coverage["start"] = [int(x[1]) for x in ints]
        coverage["end"] = [int(x[2]) for x in ints]

        if assign:
            self.coverage = coverage
        if save:
            if output_file is not None:
                coverage.to_csv(output_file, index=True)
            else:
                self.coverage.to_csv(os.path.join(self.results_dir, self.name + "_peaks.raw_coverage.csv"), index=True)

        return coverage

    def get_matrix(self, matrix=None, samples=None, matrix_name="coverage"):
        """
        Return a matrix that is an attribute of self for the requested samples.

        :param pandas.DataFrame matrix: Pandas DataFrame.
        :param list samples: Iterable of peppy.Sample objects to restrict \
                             matrix to. \
                             If not provided (`None` is passed) the matrix will not be subsetted.
        :param str matrix_name: Name of the matrix that is an attribute of the object.
        :returns pandas.DataFrame: Requested DataFrame.
        """
        # default to matrix to be normalized
        if matrix is None:
            to_norm = getattr(self, matrix_name)
        else:
            to_norm = matrix
        # default to all samples in self with matching names in matrix
        if samples is None:
            to_norm = to_norm[[s.name for s in self.samples]]
        else:
            to_norm = to_norm[[s.name for s in samples]]

        return to_norm

    def normalize_coverage_rpm(
            self, matrix=None, samples=None, mult_factor=1e6,
            save=True, assign=True):
        """
        Normalization of matrix of (n_features, n_samples) by total in each sample.

        :param str matrix: Attribute name of matrix to normalize.
        :param list samples: Iterable of peppy.Sample objects to restrict \
                             matrix to. \
                             If not provided (`None` is passed) the matrix will not be subsetted.
        :param float mult_factor: A constant to multiply values for.
        :param bool save: Whether to write normalized DataFrame to disk.
        :param bool assign: Whether to assign the normalized DataFrame to an attribute ``.

        """
        to_norm = self.get_matrix(matrix=matrix, samples=samples, matrix_name="coverage")
        # apply normalization over total
        coverage_rpm = np.log2(((1 + to_norm) / (1 + to_norm).sum()) * mult_factor)
        coverage_rpm = coverage_rpm.join(self.coverage[['chrom', 'start', 'end']])
        if save:
            coverage_rpm.to_csv(os.path.join(self.results_dir, self.name + "_peaks.coverage_rpm.csv"), index=True)
        if assign:
            self.coverage_rpm = coverage_rpm

        return coverage_rpm

    def normalize_coverage_quantiles(
            self, matrix=None, samples=None, implementation="Python",
            log_transform=True, log_constant=0.001, save=True, assign=True):
        """
        Quantile normalization of matrix of (n_features, n_samples).

        :param str matrix: Attribute name of matrix to normalize.
        :param list samples: Iterable of peppy.Sample objects to restrict \
                             matrix to. \
                             If not provided (`None` is passed) the matrix will not be subsetted.
        :param str implementation: One of `"R"` or `"Python"`. Dictates which implementation \
                                   is to be used. The R implementation comes from the \
                                   `preprocessCore` package, and the Python one is from \
                                   here https://github.com/ShawnLYU/Quantile_Normalize.
        :param bool log_transform: Whether to log transform values or not.
        :param float log_constant: A constant to add before log transformation.
        :param bool save: Whether to write normalized DataFrame to disk.
        :param bool assign: Whether to assign the normalized DataFrame to an attribute ``.

        """
        if matrix is None:
            to_norm = self.get_matrix(matrix=matrix, samples=samples, matrix_name="coverage")
        else:
            to_norm = matrix

        if implementation == "R":
            from ngs_toolkit.general import normalize_quantiles_r
            coverage_qnorm = pd.DataFrame(
                normalize_quantiles_r(to_norm.values),
                index=to_norm.index,
                columns=to_norm.columns)
        elif implementation == "Python":
            from ngs_toolkit.general import normalize_quantiles_p
            coverage_qnorm = normalize_quantiles_p(to_norm)
        else:
            raise ValueError("Implementation of quantile normalization must be one of 'R' of 'Python'")

        # Log2 transform
        if log_transform:
            coverage_qnorm = np.log2(log_constant + coverage_qnorm)

        # Add back postition columns
        coverage_qnorm = coverage_qnorm.join(self.coverage[['chrom', 'start', 'end']])
        if save:
            coverage_qnorm.to_csv(os.path.join(self.results_dir, self.name + "_peaks.coverage_qnorm.csv"), index=True)
        if assign:
            self.coverage_qnorm = coverage_qnorm

        return coverage_qnorm

    def get_peak_gccontent_length(
            self, bed_file=None, genome="hg19",
            fasta_file="/home/arendeiro/resources/genomes/{g}/{g}.fa"):
        """
        Get length and GC content of features in region set.

        :param str bed_file: A BED file with regions to calculate GC content on. \
                             Must be a 3-column BED! \
                             If not provided the calculation will be for the analysis \
                             `sites` attribute.
        :param str genome: Genome assembly.
        :param str fasta_file: Fasta file of `genome`.
        :var nuc: DataFrame with nucleotide content and length of each region.
        :returns pandas.DataFrame: DataFrame with nucleotide content and length of each region.

        """
        if bed_file is None:
            sites = self.sites
        else:
            sites = pybedtools.BedTool(bed_file)

        nuc = sites.nucleotide_content(fi=fasta_file.format(g=genome)).to_dataframe(comment="#")[["score", "blockStarts"]]
        nuc.columns = ["gc_content", "length"]
        nuc.index = [str(i.chrom) + ":" + str(i.start) + "-" + str(i.stop) for i in sites]

        # get only the sites matching the coverage (not overlapping blacklist)
        self.nuc = nuc.ix[self.coverage.index]

        self.nuc.to_csv(os.path.join(self.results_dir, self.name + "_peaks.gccontent_length.csv"), index=True)

        return self.nuc

    def normalize_gc_content(self, matrix=None, samples=None, save=True, assign=True):
        """
        Quantile normalization of matrix of (n_features, n_samples) followed by GC content correction by regression.

        Requires the R package "cqn" to be installed:
            >>> source('http://bioconductor.org/biocLite.R')
            >>> biocLite('cqn')

        :param str matrix: Attribute name of matrix to normalize.
        :param list samples: Iterable of peppy.Sample objects to restrict \
                             matrix to. \
                             If not provided (`None` is passed) the matrix will not be subsetted.
        :param bool save: Whether to write normalized DataFrame to disk.
        :param bool assign: Whether to assign the normalized DataFrame to an attribute ``.

        """

        """
        # Run manually:
        library("cqn")
        gc = read.csv("gccontent_length.csv", sep=",", row.names=1)
        cov = read.csv("coverage_qnorm.csv", sep=",", row.names=1)
        cov2 = cov[, 1:(length(cov) - 3)]

        cqn_out = cqn(cov2, x=gc$gc_content, lengths=gc$length)

        y = cqn_out$y +cqn_out$offset
        y2 = cbind(y, cov[, c("chrom", "start", "end")])
        write.csv(y2, "coverage_gc_corrected.csv", sep=",")

        # Fix R's stupid colnames replacement
        sed -i 's/ATAC.seq_/ATAC-seq_/g' coverage_gc_corrected.csv
        """
        def cqn(cov, gc_content, lengths):
            # install R package
            # source('http://bioconductor.org/biocLite.R')
            # biocLite('cqn')
            from rpy2 import robjects
            import rpy2.robjects.pandas2ri
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
            .join(self.coverage[['chrom', 'start', 'end']])
        )

        if save:
            coverage_gc_corrected.to_csv(os.path.join(self.results_dir, self.name + "_peaks.coverage_gc_corrected.csv"), index=True)
        if assign:
            self.coverage_gc_corrected = coverage_gc_corrected

        return coverage_gc_corrected

    def normalize(self, method="quantile", matrix=None, samples=None, save=True, assign=True):
        """
        Normalization of matrix of (n_features, n_samples).

        :param str method: Normalization method to apply. One of:
                           `total`: Reads per total normalization (RPM).
                           `quantile`: Quantile normalization and log2 transformation.
                           `gc_content`: Quantile normalization followed by GC content \
                                         correction by regression (cqn R package) \
                                         and log2 transformation.
        :param str matrix: Attribute name of matrix to normalize.
        :param list samples: Iterable of peppy.Sample objects to restrict \
                             matrix to. \
                             If not provided (`None` is passed) the matrix will not be subsetted.
        :param bool save: Whether to write normalized DataFrame to disk.
        :param bool assign: Whether to assign the normalized DataFrame to an attribute \
                            (see variables below for each respective normalization type).
        :var pd.DataFrame coverage_rpm: If `assign` is True, `coverage_rpm` is a pandas \
                                        DataFrame normalized with RPM.
        :var pd.DataFrame coverage_qnorm: If `assign` is True, `coverage_qnorm` \
                                          is a pandas DataFrame normalized with \
                                          quantile normalization.
        :var pd.DataFrame coverage_gc_corrected: If `assign` is True, `coverage_gc_corrected` \
                                                 is a pandas DataFrame normalized with GC correction \
                                                 and quantile normalization.
        :returns pd.DataFrame: Normalized pandas DataFrame.

        """
        if method == "quantile":
            return self.normalize_coverage_quantiles(matrix=matrix, samples=samples, save=save, assign=assign)
        elif method == "total":
            return self.normalize_coverage_rpm(matrix=matrix, samples=samples, save=save, assign=assign)
        elif method == "gc_content":
            return self.normalize_gc_content(matrix=matrix, samples=samples, save=save, assign=assign)

    def get_peak_gene_annotation(self, tss_file="refseq.refflat.tss.bed"):
        """
        Annotates peaks with closest gene.
        Needs `tss_file` file in `data_dir`/external.
        A pickle file with each feature's distance to the nearest gene is also saved.

        :param str tss_file: A 5 column BED file with TSS positions (cols 1-3), \
                             gene names (col 4) and TSS strand (col 5) exclusively! \
                             For backward compatibility, it assumes this file is \
                             in the `data_dir`/external directory.
                             This will change in the future.
        :var pd.DataFrame gene_annotation: A pandas DataFrame containing the genome
                                             annotations of the region features.
        :returns: None

        # TODO: change file to absolute path.
        # TODO: add options to limit association based on maximum distance.

        """
        # create bedtool with hg19 TSS positions
        tss = pybedtools.BedTool(os.path.join(self.data_dir, "external", tss_file))
        # get closest TSS of each cll peak
        gene_annotation = self.sites.closest(tss, d=True).to_dataframe()
        gene_annotation = gene_annotation[['chrom', 'start', 'end'] + gene_annotation.columns[-3:].tolist()]  # TODO: check this
        gene_annotation.columns = ['chrom', 'start', 'end', 'gene_name', "strand", 'distance']

        # aggregate annotation per peak, concatenate various genes (comma-separated)
        self.gene_annotation = gene_annotation.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(set([str(i) for i in x]))).reset_index()
        self.gene_annotation.index = self.gene_annotation['chrom'] + ":" + self.gene_annotation['start'].astype(str) + "-" + self.gene_annotation['end'].astype(str)

        # save to disk
        self.gene_annotation.to_csv(os.path.join(self.results_dir, self.name + "_peaks.gene_annotation.csv"), index=True)

        # save distances to all TSSs (for plotting)
        self.closest_tss_distances = gene_annotation['distance'].tolist()
        pickle.dump(self.closest_tss_distances, open(os.path.join(self.results_dir, self.name + "_peaks.closest_tss_distances.pickle"), 'wb'))

    def get_peak_genomic_location(
            self, region_files=[
                "ensembl_genes.bed", "ensembl_tss2kb.bed",
                "ensembl_utr5.bed", "ensembl_exons.bed", "ensembl_introns.bed", "ensembl_utr3.bed"],
            genome="hg19"):
        """
        Annotates peaks with its type of genomic location.

        :param list region_files: List of 3 column BED files with positions (cols 1-3)
                                  of the features to be annotated. \
                                  For backward compatibility, it assumes these files are \
                                  in the `data_dir`/external directory.
                                  This will change in the future
        :param str genome: Genome assembly.
        :var pd.DataFrame region_annotation: A pandas DataFrame containing the genome
                                                  annotations of the region features.
        :var pd.DataFrame region_annotation_b: A pandas DataFrame containing the genome
                                               annotations of the genome background.
        :returns: None

        # TODO: change region files to absolute paths.

        """

        # create background
        # shuffle regions in genome to create background (keep them in the same chromossome)
        background = self.sites.shuffle(genome=genome, chrom=True)

        for i, region in enumerate(region_files):
            region_name = region.replace(".bed", "").replace("ensembl_", "")
            r = pybedtools.BedTool(os.path.join(self.data_dir, "external", region))
            if region_name == "genes":
                region_name = "intergenic"
                df = self.sites.intersect(r, wa=True, f=0.2, v=True).to_dataframe()
                dfb = background.intersect(r, wa=True, f=0.2, v=True).to_dataframe()
            else:
                df = self.sites.intersect(r, wa=True, u=True, f=0.2).to_dataframe()
                dfb = background.intersect(r, wa=True, u=True, f=0.2).to_dataframe()
            df['genomic_region'] = region_name
            dfb['genomic_region'] = region_name
            if i == 0:
                region_annotation = df
                region_annotation_b = dfb
            else:
                region_annotation = pd.concat([region_annotation, df])
                region_annotation_b = pd.concat([region_annotation_b, dfb])

        # sort
        region_annotation.sort_values(['chrom', 'start', 'end'], inplace=True)
        region_annotation_b.sort_values(['chrom', 'start', 'end'], inplace=True)
        # remove duplicates (there shouldn't be anyway)
        region_annotation = region_annotation.reset_index(drop=True).drop_duplicates()
        region_annotation_b = region_annotation_b.reset_index(drop=True).drop_duplicates()
        # join various regions per peak
        self.region_annotation = region_annotation.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(set([str(i) for i in x]))).reset_index()
        self.region_annotation_b = region_annotation_b.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(set([str(i) for i in x]))).reset_index()

        # save to disk
        self.region_annotation.to_csv(os.path.join(self.results_dir, self.name + "_peaks.region_annotation.csv"), index=False)
        self.region_annotation_b.to_csv(os.path.join(self.results_dir, self.name + "_peaks.region_annotation_background.csv"), index=False)

    def get_peak_chromatin_state(self, chrom_state_file="E032_15_coreMarks_mnemonics", genome="hg19"):
        """
        Annotates peaks with known chromatin states.

        :param str chrom_state_file: A 4 column BED file with positions (cols 1-3) \
                                     and names (col 4) of the features to be annotated.
        :param str genome: Genome assembly.
        :var pd.DataFrame chrom_state_annotation: A pandas DataFrame containing the chromatin
                                                  annotations of the region features.
        :var pd.DataFrame chrom_state_annotation_b: A pandas DataFrame containing the chromatin
                                                    annotations of the genome background.
        :returns: None

        """
        # create bedtool with chromatin states
        cell_states = pybedtools.BedTool(chrom_state_file)

        # create background
        # shuffle regions in genome to create background (keep them in the same chromossome)
        background = self.sites.shuffle(genome=genome, chrom=True)

        # intersect with cll peaks, to create annotation, get original peaks
        chrom_state_annotation = self.sites.intersect(cell_states, wa=True, wb=True, f=0.2).to_dataframe()[['chrom', 'start', 'end', 'thickStart']]
        chrom_state_annotation_b = background.intersect(cell_states, wa=True, wb=True, f=0.2).to_dataframe()[['chrom', 'start', 'end', 'thickStart']]

        # aggregate annotation per peak, concatenate various annotations (comma-separated)
        self.chrom_state_annotation = chrom_state_annotation.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(x)).reset_index()
        self.chrom_state_annotation.columns = ['chrom', 'start', 'end', 'chromatin_state']

        self.chrom_state_annotation_b = chrom_state_annotation_b.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(x)).reset_index()
        self.chrom_state_annotation_b.columns = ['chrom', 'start', 'end', 'chromatin_state']

        # save to disk
        self.chrom_state_annotation.to_csv(os.path.join(self.results_dir, self.name + "_peaks.chromatin_state.csv"), index=False)
        self.chrom_state_annotation_b.to_csv(os.path.join(self.results_dir, self.name + "_peaks.chromatin_state_background.csv"), index=False)

    def annotate(self, samples=None, quant_matrix="coverage_gc_corrected", permissive=True):
        """
        Annotates analysis regions by aggregating region-wise annotations \
        (region, chromatin state and gene annotations - if present) and \
        calculates numeric metrics across samples (support, mean, variance, \
        deviation, dispersion, amplitude).

        The numeric matrix to be used is specified in `quant_matrix`.

        :param list samples: Iterable of peppy.Sample objects to restrict \
                             matrix to. \
                             If not provided (`None` is passed) the matrix will not be subsetted. \
                             Calculated metrics will be restricted to these samples.
        :param str quant_matrix: Attribute name of matrix to annotate.
        :param str permissive: Whether DataFrames that do not exist should be simply skipped \
                               or an error will be thrown.
        :raises AttributeError: If not `permissive` a required DataFrame does not exist as \
                                an object attribute. 
        :var pd.DataFrame coverage_annotated: A pandas DataFrame containing annotations of
                                              the region features.
        :returns: None

        """
        if samples is None:
            samples = self.samples
        if quant_matrix is None:
            quant_matrix = self.coverage_gc_corrected
        else:
            quant_matrix = getattr(self, quant_matrix)

        next_matrix = quant_matrix
        # add closest gene
        if hasattr(self, "gene_annotation"):
            self.coverage_annotated = pd.merge(
                next_matrix,
                self.gene_annotation, on=['chrom', 'start', 'end'], how="left")
            next_matrix = self.coverage_annotated
        else:
            if not permissive:
                raise AttributeError("`gene_annotation` attribute does not exist.")
        # add genomic location
        if hasattr(self, "region_annotation"):
            self.coverage_annotated = pd.merge(
                next_matrix,
                self.region_annotation[['chrom', 'start', 'end', 'genomic_region']], on=['chrom', 'start', 'end'], how="left")
            next_matrix = self.coverage_annotated
        else:
            if not permissive:
                raise AttributeError("`region_annotation` attribute does not exist.")
        # add chromatin state
        if hasattr(self, "chrom_state_annotation"):
            self.coverage_annotated = pd.merge(
                next_matrix,
                self.chrom_state_annotation[['chrom', 'start', 'end', 'chromatin_state']], on=['chrom', 'start', 'end'], how="left")
            next_matrix = self.coverage_annotated
        else:
            if not permissive:
                raise AttributeError("`chrom_state_annotation` attribute does not exist.")
        # add support
        if hasattr(self, "support"):
            self.coverage_annotated = pd.merge(
                next_matrix,
                self.support[['chrom', 'start', 'end', 'support']], on=['chrom', 'start', 'end'], how="left")
            next_matrix = self.coverage_annotated
        else:
            if not permissive:
                raise AttributeError("`support` attribute does not exist.")

        if not hasattr(self, "coverage_annotated"):
            self.coverage_annotated = next_matrix
        # calculate mean coverage
        self.coverage_annotated['mean'] = self.coverage_annotated[[s.name for s in samples]].mean(axis=1)
        # calculate coverage variance
        self.coverage_annotated['variance'] = self.coverage_annotated[[s.name for s in samples]].var(axis=1)
        # calculate std deviation (sqrt(variance))
        self.coverage_annotated['std_deviation'] = np.sqrt(self.coverage_annotated['variance'])
        # calculate dispersion (variance / mean)
        self.coverage_annotated['dispersion'] = self.coverage_annotated['variance'] / self.coverage_annotated['mean']
        # calculate qv2 (std / mean) ** 2
        self.coverage_annotated['qv2'] = (self.coverage_annotated['std_deviation'] / self.coverage_annotated['mean']) ** 2

        # calculate "amplitude" (max - min)
        self.coverage_annotated['amplitude'] = (
            self.coverage_annotated[[s.name for s in samples]].max(axis=1) -
            self.coverage_annotated[[s.name for s in samples]].min(axis=1)
        )

        # Pair indexes
        assert self.coverage.shape[0] == self.coverage_annotated.shape[0]
        self.coverage_annotated.index = self.coverage.index

        # Save
        self.coverage_annotated.to_csv(os.path.join(self.results_dir, self.name + "_peaks.coverage_qnorm.annotated.csv"), index=True)

    def get_gene_level_accessibility(self, matrix="accessibility", reduce_func=np.mean):
        """
        Get gene-level measurements of coverage.
        Requires a 'gene_annotation' attribute to be set containing a mapping between the index of `matrix` and genes
        (produced from `get_peak_gene_annotation`).

        :param str matrix: Quantification matrix to be used (e.g. 'coverage' or 'accessibility')
        :param def reduce_func: Function to apply to reduce values. Default is mean
        """
        assert hasattr(self, "gene_annotation"), """Analysis object lacks "gene_annotation" dataframe."""

        matrix = getattr(self, matrix).copy()

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
        Requires a 'gene_annotation' attribute to be set containing a mapping between the index of `matrix` and genes
        (produced from `get_peak_gene_annotation`).
        :param pandas.DataFrame differential_results: Matrix with differential results to use.
            Default is a 'differential_results' attribute of self.
        :param def reduce_func: Function to apply to reduce values. Default is mean
        """
        assert hasattr(self, "gene_annotation"), """Analysis object lacks "gene_annotation" dataframe."""

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
            import pysam
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
            map(get_sample_reads, [s.filtered for s in samples]),
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
        """

        if samples is None:
            samples = self.samples

        if by_attribute is None:
            cov = pd.melt(
                np.log2(1 + self.coverage[[s.name for s in samples]]),
                var_name="sample_name", value_name="counts"
            )
            fig, axis = plt.subplots(1, 1, figsize=(6, 1 * 4))
            sns.violinplot("sample_name", "counts", data=cov, ax=axis)
            sns.despine(fig)
            fig.savefig(os.path.join(self.results_dir, self.name + ".raw_counts.violinplot.svg"), bbox_inches="tight")

        else:
            attrs = set([getattr(s, by_attribute) for s in samples])
            fig, axis = plt.subplots(len(attrs), 1, figsize=(8, len(attrs) * 6))
            for i, attr in enumerate(attrs):
                print(attr)
                cov = pd.melt(
                    np.log2(1 + self.coverage[[s.name for s in samples if getattr(s, by_attribute) == attr]]),
                    var_name="sample_name", value_name="counts"
                )
                sns.violinplot("sample_name", "counts", data=cov, ax=axis[i])
                axis[i].set_title(attr)
                axis[i].set_xticklabels(axis[i].get_xticklabels(), rotation=90)
            sns.despine(fig)
            fig.savefig(os.path.join(self.results_dir, self.name + ".raw_counts.violinplot.by_{}.svg".format(by_attribute)), bbox_inches="tight")

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


    def unsupervised(
            self, args, **kwargs):
        """
        ATACSeqAnalysis.unsupervised is provided for backward compatibility only and will be removed in the future.
        Please use ngs_toolkit.general.unsupervised_analysis(ATACSeqAnalysis) in the future.
        """
        from ngs_toolkit.general import unsupervised_analysis

        print(PendingDeprecationWarning(
            "ATACSeqAnalysis.unsupervised is provided for backward compatibility "
            "only and will be removed. Please use "
            "ngs_toolkit.general.unsupervised_analysis(ATACSeqAnalysis) "
            "in the future."))

        unsupervised_analysis(args, **kwargs)


    def unsupervised_enrichment(self, samples, variables=["IL10_status", "subset", "replicate", "batch"]):
        """
        """
        from sklearn.decomposition import PCA
        import itertools
        from statsmodels.sandbox.stats.multicomp import multipletests

        def jackstraw(data, pcs, n_vars, n_iter=100):
            """
            """
            import rpy2.robjects as robj
            from rpy2.robjects.vectors import IntVector
            from rpy2.robjects import pandas2ri
            pandas2ri.activate()

            run = robj.r("""
                run = function(data, pcs, n_vars, B) {
                    library(jackstraw)
                    out <- jackstraw.PCA(data, r1=pcs, r=n_vars, B=B)
                    return(out$p.value)
                }
            """)
            # save to disk just in case
            data.to_csv("_tmp_matrix.jackstraw.csv", index=True)

            if type(pcs) is not int:
                pcs = IntVector(pcs)
            return run(data.values, pcs, n_vars, n_iter)

        def lola(wd="~/projects/breg/results/plots/"):
            """
            Performs location overlap analysis (LOLA) on bedfiles with regions sets.
            """
            import rpy2.robjects as robj

            run = robj.r("""
                function(wd) {
                    setwd(wd)

                    library("LOLA")

                    dbPath1 = "/data/groups/lab_bock/shared/resources/regions/LOLACore/mm10/"
                    dbPath2 = "/data/groups/lab_bock/shared/resources/regions/customRegionDB/mm10/"
                    regionDB = loadRegionDB(c(dbPath1, dbPath2))

                    for (topdir in sample(list.dirs(".", recursive=FALSE))){
                        if (!"allEnrichments.txt" %in% list.files(topdir, recursive=FALSE)){
                            print(topdir)
                            userSet <- LOLA::readBed(paste0(basename(topdir), "/", basename(topdir), ".bed"))
                            userUniverse  <- LOLA::readBed(paste0(basename(topdir), "/", "universe_sites.bed"))

                            lolaResults = runLOLA(list(userSet), userUniverse, regionDB, cores=8)
                            writeCombinedEnrichment(lolaResults, outFolder=topdir, includeSplits=FALSE)
                        }
                    }
                }
            """)

            # convert the pandas dataframe to an R dataframe
            run()

        def vertical_line(x, **kwargs):
            plt.axvline(x.mean(), **kwargs)

        # Get accessibility matrix excluding sex chroms
        X = self.coverage_qnorm_annotated[[s.name for s in samples]]
        X = X.ix[X.index[~X.index.str.contains("chrX|chrY")]]

        # Now perform association analysis (jackstraw)
        n_vars = 3
        max_sig = 1000
        alpha = 0.01
        pcs = range(1, n_vars + 1)
        pcs += list(itertools.combinations(pcs, 2))

        p_values = pd.DataFrame(index=X.index)
        for pc in pcs:
            print(pc)
            out = jackstraw(X, pc, n_vars, 100).flatten()
            if type(pc) is int:
                p_values[pc] = out
            else:
                p_values["+".join([str(x) for x in pc])] = out
        q_values = p_values.apply(lambda x: multipletests(x, method="fdr_bh")[1])
        p_values.to_csv(os.path.join(self.results_dir, "PCA.PC_pvalues.csv"), index=True)

        p_values = pd.read_csv(os.path.join(self.results_dir, "PCA.PC_pvalues.csv"), index_col=0)

        # Get enrichments of each PC-regions
        for pc in p_values.columns:
            p = p_values[pc].sort_values()
            sig = p[p < alpha].index

            # Cap to a maximum number of regions
            if len(sig) > max_sig:
                sig = p.head(max_sig).index

            # Run LOLA
            # save to bed
            output_folder = os.path.join(self.results_dir, "PCA.PC{}_regions".format(pc))
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
            universe_file = os.path.join(output_folder, "universe_sites.bed")
            self.sites.saveas(universe_file)
            bed_file = os.path.join(output_folder, "PCA.PC{}_regions.bed".format(pc))
            self.coverage[['chrom', 'start', 'end']].ix[sig].to_csv(bed_file, sep="\t", header=False, index=False)

        lola("/home/arendeiro/projects/breg/results/plots/")

        lola_enrichments = pd.DataFrame()
        enrichr_enrichments = pd.DataFrame()
        for pc in p_values.columns:
            # # read lola
            # output_folder = os.path.join(self.results_dir, "PCA.PC{}_regions".format(pc))
            # lol = pd.read_csv(os.path.join(output_folder, "allEnrichments.txt"), sep="\t")
            # lol["PC"] = pc
            # lola_enrichments = lola_enrichments.append(lol)

            # Get genes, run enrichr
            p = p_values[pc].sort_values()
            sig = p[p < alpha].index

            # Cap to a maximum number of regions
            if len(sig) > max_sig:
                sig = p.head(max_sig).index

            sig_genes = self.coverage_qnorm_annotated['gene_name'].ix[sig]
            sig_genes = [x for g in sig_genes.dropna().astype(str).tolist() for x in g.split(',')]

            enr = enrichr(pd.DataFrame(sig_genes, columns=["gene_name"]))
            enr["PC"] = pc
            enrichr_enrichments = enrichr_enrichments.append(enr)
        enrichr_enrichments.to_csv(os.path.join(self.results_dir, "PCA.enrichr.csv"), index=False, encoding="utf-8")
        enrichr_enrichments = pd.read_csv(os.path.join(self.results_dir, "PCA.enrichr.csv"))
        lola_enrichments.to_csv(os.path.join(self.results_dir, "PCA.lola.csv"), index=False, encoding="utf-8")

        # Plots

        # p-value distributions
        g = sns.FacetGrid(data=pd.melt(-np.log10(p_values), var_name="PC", value_name="-log10(p-value)"), col="PC", col_wrap=5)
        g.map(sns.distplot, "-log10(p-value)", kde=False)
        g.map(plt.axvline, x=-np.log10(alpha), linestyle="--")
        g.add_legend()
        g.fig.savefig(os.path.join(self.results_dir, "PCA.PC_pvalues.distplot.svg"), bbox_inches="tight")

        # Volcano plots (loading vs p-value)
        # get PCA loadings
        pca = PCA()
        pca.fit(X.T)
        loadings = pd.DataFrame(pca.components_.T, index=X.index, columns=range(1, X.shape[1] + 1))
        loadings.to_csv(os.path.join(self.results_dir, "PCA.loadings.csv"), index=True, encoding="utf-8")

        melted_loadings = pd.melt(loadings.reset_index(), var_name="PC", value_name="loading", id_vars=["index"])
        melted_loadings["PC"] = melted_loadings["PC"].astype(str)
        melted_loadings = melted_loadings.set_index(["index", "PC"])
        melted_p_values = pd.melt((-np.log10(p_values)).reset_index(), var_name="PC", value_name="-log10(p-value)", id_vars=["index"]).set_index(["index", "PC"])
        melted = melted_loadings.join(melted_p_values)

        g = sns.FacetGrid(data=melted.dropna().reset_index(), col="PC", col_wrap=5, sharey=False, sharex=False)
        g.map(plt.scatter, "loading", "-log10(p-value)", s=2, alpha=0.5)
        g.map(plt.axhline, y=-np.log10(alpha), linestyle="--")
        g.add_legend()
        g.fig.savefig(os.path.join(self.results_dir, "PCA.PC_pvalues_vs_loading.scatter.png"), bbox_inches="tight", dpi=300)

        # Plot enrichments
        # LOLA
        # take top n per PC
        import string
        lola_enrichments["set_id"] = lola_enrichments[
            ["collection", "description", "cellType", "tissue", "antibody", "treatment"]].astype(str).apply(string.join, axis=1)

        top = lola_enrichments.set_index('set_id').groupby("PC")['pValueLog'].nlargest(50)
        top_ids = top.index.get_level_values('set_id').unique()

        pivot = pd.pivot_table(
            lola_enrichments,
            index="set_id", columns="PC", values="pValueLog").fillna(0)
        pivot.index = pivot.index.str.replace(" nan", "").str.replace("blueprint blueprint", "blueprint").str.replace("None", "")
        top_ids = top_ids.str.replace(" nan", "").str.replace("blueprint blueprint", "blueprint").str.replace("None", "")

        g = sns.clustermap(
            pivot.ix[top_ids],
            cbar_kws={"label": "Enrichment: -log10(p-value)"}, cmap="Spectral_r",
            col_cluster=True)
        for tick in g.ax_heatmap.get_xticklabels():
            tick.set_rotation(90)
        for tick in g.ax_heatmap.get_yticklabels():
            tick.set_rotation(0)
        g.fig.savefig(os.path.join(self.results_dir, "PCA.PC_pvalues.lola_enrichments.svg"), bbox_inches="tight", dpi=300)

        g = sns.clustermap(
            pivot.ix[top_ids],
            cbar_kws={"label": "Enrichment: p-value z-score"},
            col_cluster=True, z_score=0)
        for tick in g.ax_heatmap.get_xticklabels():
            tick.set_rotation(90)
        for tick in g.ax_heatmap.get_yticklabels():
            tick.set_rotation(0)
        g.fig.savefig(os.path.join(self.results_dir, "PCA.PC_pvalues.lola_enrichments.z_score.svg"), bbox_inches="tight", dpi=300)

        # Enrichr
        for gene_set_library in enrichr_enrichments["gene_set_library"].drop_duplicates():
            enr = enrichr_enrichments[enrichr_enrichments["gene_set_library"] == gene_set_library]

            top = enr.set_index('description').groupby("PC")['p_value'].nsmallest(20)
            top_ids = top.index.get_level_values('description').unique()

            pivot = pd.pivot_table(enr, index="description", columns="PC", values="p_value").fillna(1)
            pivot.index = pivot.index.str.extract("(.*)[,\_\(].*", expand=False).str.replace("_Homo sapiens", "")
            top_ids = top_ids.str.extract("(.*)[,\_\(].*", expand=False).str.replace("_Homo sapiens", "")

            g = sns.clustermap(
                -np.log10(pivot.ix[top_ids]), cmap='BuGn',
                cbar_kws={"label": "Enrichment: -log10(p-value)"}, col_cluster=True, figsize=(6, 15))
            for tick in g.ax_heatmap.get_xticklabels():
                tick.set_rotation(90)
            for tick in g.ax_heatmap.get_yticklabels():
                tick.set_rotation(0)
            g.fig.savefig(os.path.join(self.results_dir, "PCA.PC_pvalues.enrichr_enrichments.{}.svg".format(gene_set_library)), bbox_inches="tight", dpi=300)

            g = sns.clustermap(
                -np.log10(pivot.ix[top_ids]),
                cbar_kws={"label": "Enrichment: p-value z-score"}, col_cluster=True, z_score=0, figsize=(6, 15))
            for tick in g.ax_heatmap.get_xticklabels():
                tick.set_rotation(90)
            for tick in g.ax_heatmap.get_yticklabels():
                tick.set_rotation(0)
            g.fig.savefig(os.path.join(self.results_dir, "PCA.PC_pvalues.enrichr_enrichments.{}.z_score.svg".format(gene_set_library)), bbox_inches="tight", dpi=300)


def characterize_regions_structure(analysis, df, prefix, output_dir, universe_df=None):
    # use all sites as universe
    if universe_df is None:
        universe_df = analysis.coverage_annotated

    # make output dirs
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # compare genomic regions and chromatin_states
    enrichments = pd.DataFrame()
    for i, var in enumerate(['genomic_region', 'chromatin_state']):
        # prepare:
        # separate comma-delimited fields:
        df_count = Counter(df[var].str.split(',').apply(pd.Series).stack().tolist())
        df_universe_count = Counter(universe_df[var].str.split(',').apply(pd.Series).stack().tolist())

        # divide by total:
        df_count = {k: v / float(len(df)) for k, v in df_count.items()}
        df_universe_count = {k: v / float(len(universe_df)) for k, v in df_universe_count.items()}

        # join data, sort by subset data
        both = pd.DataFrame([df_count, df_universe_count], index=['subset', 'all']).T
        both = both.sort_values("subset")
        both['region'] = both.index
        data = pd.melt(both, var_name="set", id_vars=['region']).replace(np.nan, 0)

        # sort for same order
        data.sort_values('region', inplace=True)

        # g = sns.FacetGrid(col="region", data=data, col_wrap=3, sharey=True)
        # g.map(sns.barplot, "set", "value")
        # plt.savefig(os.path.join(output_dir, "%s_regions.%s.svg" % (prefix, var)), bbox_inches="tight")

        fc = pd.DataFrame(np.log2(both['subset'] / both['all']), columns=['value'])
        fc['variable'] = var

        # append
        enrichments = enrichments.append(fc)

    # save
    enrichments.to_csv(os.path.join(output_dir, "%s_regions.region_enrichment.csv" % prefix), index=True)


def characterize_regions_function(analysis, df, output_dir, prefix, universe_file=None, run=True, genome="hg19"):
    """
    Performs a range of functional enrichments of a set of regions given in `df`
    (a dataframe with 'chrom', 'start', 'end', 'gene_name', 'ensebl_gene_id' columns - typically the coverage_annotated dataframe).
    Will extract regions, their underlying sequence, associated genes, perform motif enrichment,
    location overlap analysis (LOLA), and gene set enrichment (using the Enrichr API).

    This requires several programs and R libraries:
        - MEME suite (AME)
        - HOMER suite (findMotifsGenome.pl)
        - LOLA (R library)

    Additionally, some genome-specific databases are needed to run these programs.
    """
    from ngs_toolkit.general import bed_to_fasta, meme_ame, homer_motifs, lola, enrichr, standard_score

    # use all sites as universe
    if universe_file is None:
        universe_file = analysis.sites

    # make output dirs
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # save to bed
    bed_file = os.path.join(output_dir, "{}_regions.bed".format(prefix))
    df[['chrom', 'start', 'end']].to_csv(bed_file, sep="\t", header=False, index=False)
    # save as tsv
    tsv_file = os.path.join(output_dir, "{}_regions.tsv".format(prefix))
    df[['chrom', 'start', 'end']].reset_index().to_csv(tsv_file, sep="\t", header=False, index=False)

    # export gene names
    (df['gene_name']
    .str.split(",")
    .apply(pd.Series, 1)
    .stack()
    .drop_duplicates()
    .to_csv(
        os.path.join(output_dir, "{}_genes.symbols.txt".format(prefix)),
        index=False))
    if "ensembl_gene_id" in df.columns:
        # export ensembl gene names
        (df['ensembl_gene_id']
        .str.split(",")
        .apply(pd.Series, 1)
        .stack()
        .drop_duplicates()
        .to_csv(
            os.path.join(output_dir, "{}_genes.ensembl.txt".format(prefix)),
            index=False))

    # export gene symbols with scaled absolute fold change
    if "log2FoldChange" in df.columns:
        df["score"] = standard_score(abs(df["log2FoldChange"]))
        df["abs_fc"] = abs(df["log2FoldChange"])

        d = df[['gene_name', 'score']].sort_values('score', ascending=False)

        # split gene names from score if a reg.element was assigned to more than one gene
        a = (
            d['gene_name']
            .str.split(",")
            .apply(pd.Series, 1)
            .stack()
        )
        a.index = a.index.droplevel(1)
        a.name = 'gene_name'
        d = d[['score']].join(a)
        # reduce various ranks to mean per gene
        d = d.groupby('gene_name').mean().reset_index()
        d.to_csv(os.path.join(output_dir, "{}_genes.symbols.score.csv".format(prefix)), index=False)

    # Motifs
    # de novo motif finding - enrichment
    fasta_file = os.path.join(output_dir, "{}_regions.fa".format(prefix))
    bed_to_fasta(bed_file, fasta_file, genome=genome)

    if not run:
        return

    omap = {"hg38": "human", "hg19": "human", "mm10": "mouse"}
    meme_ame(fasta_file, output_dir, organism=omap[genome])
    homer_motifs(bed_file, output_dir, genome=genome)

    # Lola
    try:
        lola(bed_file, universe_file, output_dir, genome=genome)
    except:
        print("LOLA analysis for {} failed!".format(prefix))

    # Enrichr
    results = enrichr(df[['chrom', 'start', 'end', "gene_name"]])
    results.to_csv(os.path.join(output_dir, "{}_regions.enrichr.csv".format(prefix)), index=False, encoding='utf-8')


def metagene_plot(bams, labels, output_prefix, region="genebody", genome="hg19"):
    from pypiper import NGSTk
    import textwrap
    import os
    tk = NGSTk()

    job_name = output_prefix
    job_file = output_prefix + ".sh"
    job_log = output_prefix + ".log"

    # write ngsplot config file to disk
    config_file = os.path.join(os.environ['TMPDIR'], "ngsplot_config.txt")
    with open(config_file, "w") as handle:
        for i in range(len(bams)):
            handle.write("\t".join([bams[i], "-1", labels[i]]) + "\n")

    cmd = tk.slurm_header(job_name, job_log, queue="mediumq", time="1-10:00:00", mem_per_cpu=8000, cpus_per_task=8)

    # build plot command
    if region == "genebody":
        cmd += """xvfb-run ngs.plot.r -G {0} -R {1} -C {2} -O {3} -L 3000 -GO km\n""".format(genome, region, config_file, output_prefix)
    elif region == "tss":
        cmd += """xvfb-run ngs.plot.r -G {0} -R {1} -C {2} -O {3} -L 3000 -FL 300\n""".format(genome, region, config_file, output_prefix)

    cmd += tk.slurm_footer()

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurm_submit_job(job_file)


def global_changes(samples, trait="knockout"):
    import glob
    import re

    output_dir = os.path.join("data", "merged")
    # sel_samples = [s for s in samples if not pd.isnull(getattr(s, trait))]
    # groups = sorted(list(set([getattr(s, trait) for s in sel_samples])))
    groups = [os.path.basename(re.sub(".merged.sorted.bam", "", x)) for x in glob.glob(output_dir + "/*.merged.sorted.bam")]

    for region in ["genebody", "tss"]:
        print(metagene_plot(
            [os.path.abspath(os.path.join(output_dir, group + ".merged.sorted.bam")) for group in groups],
            groups,
            os.path.abspath(os.path.join(output_dir, "%s.metaplot" % region)),
            region=region
        ))


def nucleosome_changes(analysis, samples):
    # select only ATAC-seq samples
    df = analysis.prj.sheet.df[analysis.prj.sheet.df["library"] == "ATAC-seq"]

    groups = list()
    for attrs, index in df.groupby(["library", "cell_line", "knockout", "clone"]).groups.items():
        name = "_".join([a for a in attrs if not pd.isnull(a)])
        groups.append(name)
    groups = sorted(groups)

    # nucleosomes per sample
    nucpos_metrics = {
        "z-score": 4,
        "nucleosome_occupancy_estimate": 5,
        "lower_bound_for_nucleosome_occupancy_estimate": 6,
        "upper_bound_for_nucleosome_occupancy_estimat": 7,
        "log_likelihood_ratio": 8,
        "normalized_nucleoatac_signal": 9,
        "cross-correlation_signal_value_before_normalization": 10,
        "number_of_potentially_nucleosome-sized_fragments": 11,
        "number_of_fragments_smaller_than_nucleosome_sized": 12,
        "fuzziness": 13
    }
    # nucleosome-free regions per sample
    nfrpos_metrics = {
        "mean_occupancy": 4,
        "minimum_upper_bound_occupancy": 5,
        "insertion_density": 6,
        "bias_density": 7,
    }
    for data_type, metrics in [("nucpos", nucpos_metrics), ("nfrpos", nfrpos_metrics)]:
        figs = list()
        axizes = list()
        for metric in metrics:
            fig, axis = plt.subplots(5, 6, figsize=(6 * 4, 5 * 4), sharex=True, sharey=True)
            figs.append(fig)
            axizes.append(axis.flatten())

        counts = pd.Series()
        for i, group in enumerate(groups):
            print(data_type, group)
            s = pd.read_csv(
                os.path.join("results", "nucleoatac", group, group + ".{}.bed.gz".format(data_type)),
                sep="\t", header=None)
            counts[group] = s.shape[0]

            for j, (metric, col) in enumerate(metrics.items()):
                print(data_type, group, metric, col)
                sns.distplot(s[col - 1].dropna(), hist=False, ax=axizes[j][i])
                axizes[j][i].set_title(group)

        for i, metric in enumerate(metrics):
            sns.despine(figs[i])
            figs[i].savefig(os.path.join("results", "nucleoatac", "plots", "{}.{}.svg".format(data_type, metric)), bbox_inches="tight")

        fig, axis = plt.subplots(1, 1, figsize=(1 * 4, 1 * 4))
        sns.barplot(counts, counts.index, orient="horizontal", ax=axis, color=sns.color_palette("colorblind")[0])
        axis.set_yticklabels(axis.get_yticklabels(), rotation=0)
        axis.set_xlabel("Calls")
        sns.despine(fig)
        fig.savefig(os.path.join("results", "nucleoatac", "plots", "{}.count_per_sample.svg".format(data_type)), bbox_inches="tight")

    # fragment distribution
    for data_type in ["InsertionProfile", "InsertSizes", "fragmentsizes"]:
        fig, axis = plt.subplots(5, 6, figsize=(6 * 4, 5 * 4))
        axis = axis.flatten()

        data = pd.DataFrame()
        for i, group in enumerate(groups):
            print(data_type, group)
            s = pd.read_csv(
                os.path.join("results", "nucleoatac", group, group + ".{}.txt".format(data_type)),
                sep="\t", header=None, squeeze=True, skiprows=5 if data_type == "fragmentsizes" else 0)
            if data_type == "InsertionProfile":
                a = (len(s.index) - 1) / 2.
                s.index = np.arange(-a, a + 1)
            if data_type == "fragmentsizes":
                s = s.squeeze()
            data[group] = s

            axis[i].plot(s)
            axis[i].set_title(group)
        sns.despine(fig)
        fig.savefig(os.path.join("results", "nucleoatac", "plots", "{}.svg".format(data_type)), bbox_inches="tight")

        norm_data = data.apply(lambda x: x / data['ATAC-seq_HAP1_WT_C631'])
        if data_type == "fragmentsizes":
            norm_data = norm_data.loc[50:, :]
        fig, axis = plt.subplots(5, 6, figsize=(6 * 4, 5 * 4), sharey=True)
        axis = axis.flatten()
        for i, group in enumerate(groups):
            print(data_type, group)
            axis[i].plot(norm_data[group])
            axis[i].set_title(group)
        sns.despine(fig)
        fig.savefig(os.path.join("results", "nucleoatac", "plots", "{}.WT_norm.svg".format(data_type)), bbox_inches="tight")

    # Vplots and
    v_min, v_max = (105, 251)
    # Vplots over WT
    for data_type in ["VMat"]:
        fig, axis = plt.subplots(5, 6, figsize=(6 * 4, 5 * 4))
        fig2, axis2 = plt.subplots(5, 6, figsize=(6 * 4, 5 * 4), sharey=True)
        axis = axis.flatten()
        axis2 = axis2.flatten()

        group = 'ATAC-seq_HAP1_WT_C631'
        wt = pd.read_csv(
            os.path.join("results", "nucleoatac", group, group + ".{}".format(data_type)),
            sep="\t", header=None, skiprows=7)
        wt.index = np.arange(v_min, v_max)
        a = (len(wt.columns) - 1) / 2.
        wt.columns = np.arange(-a, a + 1)
        wt = wt.loc[0:300, :]

        for i, group in enumerate(groups):
            print(data_type, group)
            m = pd.read_csv(
                os.path.join("results", "nucleoatac", group, group + ".{}".format(data_type)),
                sep="\t", header=None, skiprows=7)
            m.index = np.arange(v_min, v_max)
            a = (len(m.columns) - 1) / 2.
            m.columns = np.arange(-a, a + 1)
            m = m.loc[0:300, :]

            n = m / wt

            axis[i].imshow(m.sort_index(ascending=False))
            axis[i].set_title(group)
            axis2[i].imshow(n.sort_index(ascending=False))
            axis2[i].set_title(group)
        sns.despine(fig)
        sns.despine(fig2)
        fig.savefig(os.path.join("results", "nucleoatac", "plots", "{}.svg".format(data_type)), bbox_inches="tight")
        fig2.savefig(os.path.join("results", "nucleoatac", "plots", "{}.WT_norm.svg".format(data_type)), bbox_inches="tight")


def investigate_nucleosome_positions(self, samples, cluster=True):
    df = self.prj.sheet.df[self.prj.sheet.df["library"] == "ATAC-seq"]
    groups = list()
    for attrs, index in df.groupby(["library", "cell_line", "knockout", "clone"]).groups.items():
        name = "_".join([a for a in attrs if not pd.isnull(a)])
        groups.append(name)
    groups = sorted(groups)

    def get_differential(diff, conditions_for, directions_for):
        # filter conditions
        d = diff[diff['comparison'].isin(conditions_for)]

        # filter directions
        d = pd.concat([d[f(d['log2FoldChange'], x)] for f, x in directions_for])

        # make bed format
        df = pd.Series(d.index.str.split(":")).apply(lambda x: pd.Series([x[0]] + x[1].split("-"))).drop_duplicates().reset_index(drop=True)
        df[0] = df[0].astype(str)
        df[1] = df[1].astype(np.int64)
        df[2] = df[2].astype(np.int64)
        return df

    def center_window(bedtool, width=1000):
        chroms = ["chr" + str(x) for x in range(1, 23)] + ["chrX", "chrY"]

        sites = list()
        for site in bedtool:
            if site.chrom in chroms:
                mid = site.start + ((site.end - site.start) / 2)
                sites.append(pybedtools.Interval(site.chrom, mid - (width / 2), (mid + 1) + (width / 2)))
        return pybedtools.BedTool(sites)

    def center_series(series, width=1000):
        mid = series[1] + ((series[2] - series[1]) / 2)
        return pd.Series([series[0], mid - (width / 2), (mid + 1) + (width / 2)])

    # Regions to look at
    regions_pickle = os.path.join(self.results_dir, "nucleoatac", "all_types_of_regions.pickle")
    if os.path.exists():
        regions = pickle.load(open(regions_pickle, "rb"))
    else:
        regions = dict()
        # All accessible sites
        out = os.path.join(self.results_dir, "nucleoatac", "all_sites.bed")
        center_window(self.sites).to_dataframe()[['chrom', 'start', 'end']].to_csv(out, index=False, header=None, sep="\t")
        regions["all_sites"] = out

        # get bed file of promoter proximal sites
        promoter_sites = os.path.join(self.results_dir, "nucleoatac", "promoter_proximal.bed")
        self.coverage_annotated[self.coverage_annotated['distance'].astype(int) < 5000][['chrom', 'start', 'end']].to_csv(
            promoter_sites, index=False, header=None, sep="\t")
        regions["promoter"] = promoter_sites

        # All accessible sites - that are distal
        out = os.path.join(self.results_dir, "nucleoatac", "distal_sites.bed")
        center_window(self.sites.intersect(pybedtools.BedTool(promoter_sites), v=True, wa=True)).to_dataframe()[['chrom', 'start', 'end']].to_csv(out, index=False, header=None, sep="\t")
        regions["distal_sites"] = out

        # Differential sites
        diff = pd.read_csv(os.path.join("results", "deseq_knockout", "deseq_knockout.knockout.csv"), index_col=0)
        diff = diff[(diff["padj"] < 0.01) & (abs(diff["log2FoldChange"]) > 1.)]

        # Loosing accessibility with ARID1A/SMARCA4 KO
        less_ARID1ASMARCA4 = get_differential(
            diff,
            ['ARID1A-WT', 'SMARCA4-WT'],
            [(np.less_equal, -1), (np.less_equal, -1)]).apply(center_series, axis=1)
        out = os.path.join(self.results_dir, "nucleoatac", "diff_sites.less_ARID1ASMARCA4.bed")
        less_ARID1ASMARCA4.to_csv(out, index=False, header=None, sep="\t")
        regions["diff_sites.less_ARID1ASMARCA4"] = out

        # Gaining accessibility with ARID1A/SMARCA4 KO
        more_ARID1ASMARCA4 = get_differential(
            diff,
            ['ARID1A-WT', 'SMARCA4-WT'],
            [(np.greater_equal, 1), (np.greater_equal, 1)]).apply(center_series, axis=1)
        out = os.path.join(self.results_dir, "nucleoatac", "diff_sites.more_ARID1ASMARCA4.bed")
        more_ARID1ASMARCA4.to_csv(out, index=False, header=None, sep="\t")
        regions["diff_sites.more_ARID1ASMARCA4"] = out

        # TFBSs
        tfs = ["CTCF", "BCL", "SMARC", "POU5F1", "SOX2", "NANOG", "TEAD4"]
        for tf in tfs:
            tf_bed = pybedtools.BedTool("/home/arendeiro/resources/genomes/hg19/motifs/TFs/{}.true.bed".format(tf))
            out = os.path.join(self.results_dir, "nucleoatac", "tfbs.%s.bed" % tf)
            center_window(tf_bed.intersect(self.sites, wa=True)).to_dataframe()[['chrom', 'start', 'end']].to_csv(out, index=False, header=None, sep="\t")
            regions["tfbs.%s" % tf] = out

        pickle.dump(regions, open(regions_pickle, "wb"))

    # Launch jobs
    for group in groups:
        output_dir = os.path.join(self.results_dir, "nucleoatac", group)

        # Signals to measure in regions
        signal_files = [
            ("signal", os.path.join(self.data_dir, "merged", group + ".merged.sorted.bam")),
            ("nucleosome", os.path.join(self.data_dir, "merged", group + ".nucleosome_reads.bam")),
            ("nucleosome_free", os.path.join(self.data_dir, "merged", group + ".nucleosome_free_reads.bam")),
            ("nucleoatac", os.path.join("results", "nucleoatac", group, group + ".nucleoatac_signal.smooth.bedgraph.gz")),
            ("dyads", os.path.join("results", "nucleoatac", group, group + ".nucpos.bed.gz"))
        ]
        for region_name, bed_file in regions.items():
            for label, signal_file in signal_files:
                print(group, region_name, label)
                # run job
                run_coverage_job(bed_file, signal_file, label, ".".join([group, region_name, label]), output_dir, window_size=2001)
                # run vplot
                if label == "signal":
                    run_vplot_job(bed_file, signal_file, ".".join([group, region_name, label]), output_dir)

    # Collect signals
    signals = pd.DataFrame(columns=['group', 'region', 'label'])
    # signals = pd.read_csv(os.path.join(self.results_dir, "nucleoatac", "collected_coverage.csv"))

    import re
    for group in [g for g in groups if any([re.match(".*%s.*" % x, g) for x in ["C631", "HAP1_ARID1", "HAP1_SMARCA"]])]:
        output_dir = os.path.join(self.results_dir, "nucleoatac", group)
        signal_files = [
            # ("signal", os.path.join(self.data_dir, "merged", group + ".merged.sorted.bam")),
            # ("nucleosome", os.path.join(self.data_dir, "merged", group + ".nucleosome_reads.bam")),
            # ("nucleosome_free", os.path.join(self.data_dir, "merged", group + ".nucleosome_free_reads.bam")),
            ("nucleoatac", os.path.join("results", "nucleoatac", group, group + ".nucleoatac_signal.smooth.bedgraph.gz")),
            # ("dyads", os.path.join("results", "nucleoatac", group, group + ".nucpos.bed.gz"))
        ]
        for region_name, bed_file in [regions.items()[-1]]:
            for label, signal_file in signal_files:
                # Skip already done
                if len(signals[
                        (signals["group"] == group) &
                        (signals["region"] == region_name) &
                        (signals["label"] == label)
                ]) > 0:
                    print("Continuing", group, region_name, label)
                    continue

                print(group, region_name, label)
                df = pd.read_csv(os.path.join(output_dir, "{}.coverage_matrix.csv".format(".".join([group, region_name, label, label]))), index_col=0)
                df = df.mean(0).reset_index(name="value").rename(columns={"index": "distance"})
                df["group"] = group
                df["region"] = region_name
                df["label"] = label
                signals = signals.append(df, ignore_index=True)
    signals.to_csv(os.path.join(self.results_dir, "nucleoatac", "collected_coverage.csv"), index=False)

    signals = pd.read_csv(os.path.join(self.results_dir, "nucleoatac", "collected_coverage.csv"))

    signals = signals[(signals["distance"] > -400) & (signals["distance"] < 400)]

    # plot together
    region_order = sorted(signals["region"].unique(), reverse=True)
    group_order = sorted(signals["group"].unique(), reverse=True)
    label_order = sorted(signals["label"].unique(), reverse=True)

    # raw
    g = sns.FacetGrid(signals, hue="group", col="region", row="label", hue_order=group_order, row_order=label_order, col_order=region_order, sharex=False, sharey=False)
    g.map(plt.plot, "distance", "value")
    g.add_legend()
    g.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.raw_mean_coverage.svg"), bbox_inches="tight")

    # normalized
    signals["norm_values"] = signals.groupby(["region", "label", "group"])["value"].apply(lambda x: (x - x.mean()) / x.std())

    g = sns.FacetGrid(signals, hue="group", col="region", row="label", hue_order=group_order, row_order=label_order, col_order=region_order, sharex=False, sharey=False)
    g.map(plt.plot, "distance", "norm_values")
    g.add_legend()
    g.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.norm_mean_coverage.svg"), bbox_inches="tight")

    # normalized smoothed
    signals["norm_smooth_values"] = signals.groupby(["region", "label", "group"])["value"].apply(lambda x: pd.rolling_window(((x - x.mean()) / x.std()), 10))

    g = sns.FacetGrid(signals, hue="group", col="region", row="label", hue_order=group_order, row_order=label_order, col_order=region_order, sharex=False, sharey=False)
    g.map(plt.plot, "distance", "norm_smooth_values")
    g.add_legend()
    g.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.norm_mean_coverage.smooth.svg"), bbox_inches="tight")

    #
    # specific regions/samples
    specific_signals = signals[
        (signals["group"].str.contains("ARID1|SMARCA|C631")) &
        (~signals["group"].str.contains("OV90|GFP")) &

        # (signals["region"].str.contains("diff_sites")) &

        (signals["label"] == "nucleoatac")
    ]

    region_order = sorted(specific_signals["region"].unique(), reverse=True)
    group_order = sorted(specific_signals["group"].unique(), reverse=True)
    label_order = sorted(specific_signals["label"].unique(), reverse=True)

    g = sns.FacetGrid(specific_signals, hue="group", col="region", col_wrap=4, hue_order=group_order, row_order=label_order, col_order=region_order, sharex=False, sharey=False)
    g.map(plt.plot, "distance", "value")
    g.add_legend()
    g.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.specific.extended.svg"), bbox_inches="tight")

    # normalized (centered), zoom in center
    specific_signals["norm_values"] = specific_signals.groupby(["region", "label", "group"])["value"].apply(lambda x: (x - x.mean()) / x.std())
    specific_signals = specific_signals[
        (specific_signals["distance"] < 200) &
        (specific_signals["distance"] > -200)]
    g = sns.FacetGrid(specific_signals, hue="group", col="region", row="label", hue_order=group_order, row_order=label_order, col_order=region_order, sharex=False, sharey=False)
    g.map(plt.plot, "distance", "norm_values")
    g.add_legend()
    g.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.specific.norm.svg"), bbox_inches="tight")

    #

    # Violinplots of nucleosome occupancies

    # Heatmap of nucleosome occupancies
    # Collect signals
    sel_groups = [x for x in groups if "OV90" not in x and "GFP" not in x and ("ARID1" in x or "SMARCA" in x or "WT" in x)]
    regions = pickle.load(open(regions_pickle, "rb"))
    sel_regions = {k: v for k, v in regions.items() if "diff" in k}

    # get parameters based on WT accessibility
    region_order = dict()
    region_vmax = dict()
    region_norm_vmax = dict()
    output_dir = os.path.join(self.results_dir, "nucleoatac", "ATAC-seq_HAP1_WT_C631")
    for region_name in sel_regions.keys():
        df = pd.read_csv(os.path.join(output_dir, "{}.coverage_matrix.csv".format(".".join(["ATAC-seq_HAP1_WT_C631", region_name, "nucleoatac", "nucleoatac"]))), index_col=0)
        # vmax
        region_vmax[region_name] = np.percentile(df, 95)
        region_norm_vmax[region_name] = np.percentile((df - df.mean(0)) / df.std(0), 95)
        # region order
        region_order[region_name] = df.sum(axis=1).sort_values().index  # sorted by mean
        # region_order[region_name] = g.dendrogram_row.dendrogram

    # plot all
    fig, axis = plt.subplots(len(sel_regions), len(sel_groups), figsize=(len(sel_groups) * 4, len(sel_regions) * 4))
    fig2, axis2 = plt.subplots(len(sel_regions), len(sel_groups), figsize=(len(sel_groups) * 4, len(sel_regions) * 4))
    for j, group in enumerate(sorted(sel_groups, reverse=True)):
        output_dir = os.path.join(self.results_dir, "nucleoatac", group)
        signal_files = [
            ("nucleoatac", os.path.join("results", "nucleoatac", group, group + ".nucleoatac_signal.smooth.bedgraph.gz"))
        ]
        for i, (region_name, bed_file) in enumerate(sel_regions.items()):
            for label, signal_file in signal_files:
                print(group, region_name, label)
                df = pd.read_csv(os.path.join(output_dir, "{}.coverage_matrix.csv".format(".".join([group, region_name, label, label]))), index_col=0)

                d = df.ix[region_order[region_name]]
                axis[i, j].imshow(
                    d,
                    norm=None, cmap="inferno", vmax=region_vmax[region_name], extent=[-500, 500, 0, 10], aspect="auto")  # aspect=100
                d_norm = (d - d.mean(0)) / d.std(0)
                axis2[i, j].imshow(
                    d_norm,
                    norm=None, cmap="inferno", vmax=region_norm_vmax[region_name], extent=[-500, 500, 0, 10], aspect="auto")  # aspect=100
                for ax in [axis, axis2]:
                    ax[i, j].set_title(group)
                    ax[i, j].set_xlabel("distance")
                    ax[i, j].set_ylabel(region_name)
    sns.despine(fig, top=True, right=True, left=True, bottom=True)
    sns.despine(fig2, top=True, right=True, left=True, bottom=True)
    fig.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.specific.heatmap.png"), bbox_inches="tight", dpi=300)
    fig.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.specific.heatmap.svg"), bbox_inches="tight")
    fig2.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.specific.heatmap.centered.png"), bbox_inches="tight", dpi=300)
    fig2.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "collected_coverage.specific.heatmap.centered.svg"), bbox_inches="tight")


def phasograms(self, samples, max_dist=10000, rolling_window=50, plotting_window=(0, 500)):
    from ngs_toolkit.general import detect_peaks

    df = self.prj.sheet.df[self.prj.sheet.df["library"] == "ATAC-seq"]
    groups = list()
    for attrs, index in df.groupby(["library", "cell_line", "knockout", "clone"]).groups.items():
        name = "_".join([a for a in attrs if not pd.isnull(a)])
        groups.append(name)
    groups = sorted(groups)

    def difference_matrix(a):
        x = np.reshape(a, (len(a), 1))
        return x - x.transpose()

    distances = dict()

    for group in groups:
        print(group)
        # Get dyad calls from nucleoatac
        df = pd.read_csv(os.path.join(self.results_dir, "nucleoatac", group, group + ".nucpos.bed.gz"), sep="\t", header=None)

        # separate by chromosome (groupby?)
        # count pairwise distance
        dists = list()
        for chrom in df[0].unique():
            d = abs(difference_matrix(df[df[0] == chrom][1]))
            dd = d[(d < max_dist) & (d != 0)]
            dists += dd.tolist()
        distances[group] = dists

    pickle.dump(distances, open(os.path.join(self.results_dir, "nucleoatac", "phasogram.distances.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
    distances = pickle.load(open(os.path.join(self.results_dir, "nucleoatac", "phasogram.distances.pickle"), "rb"))

    # Plot distances between dyads
    from scipy.ndimage.filters import gaussian_filter1d
    n_rows = n_cols = int(np.ceil(np.sqrt(len(groups))))
    n_rows -= 1
    fig, axis = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(n_cols * 3, n_rows * 2))
    fig2, axis2 = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(n_cols * 3, n_rows * 2))
    axis = axis.flatten()
    axis2 = axis2.flatten()
    for i, group in enumerate(groups):
        # Count frequency of dyad distances
        x = pd.Series(distances[group])
        y = x.value_counts().sort_index()
        y = y.ix[range(plotting_window[0], plotting_window[1])]
        y /= y.sum()

        # Find peaks
        y2 = pd.Series(gaussian_filter1d(y, 5), index=y.index)
        peak_indices = detect_peaks(y2.values, mpd=73.5)[:3]
        print(group, y2.iloc[peak_indices].index)

        # Plot distribution and peaks
        axis[i].plot(y.index, y, color="black", alpha=0.6, linewidth=0.5)
        axis[i].plot(y2.index, y2, color=sns.color_palette("colorblind")[0], linewidth=1)
        axis[i].scatter(y2.iloc[peak_indices].index, y2.iloc[peak_indices], s=25, color="orange")
        for peak in y2.iloc[peak_indices].index:
            axis[i].axvline(peak, color="black", linestyle="--")
        axis[i].set_title(group)

        # Transform into distances between nucleosomes
        # Plot distribution and peaks
        axis2[i].plot(y.index - 147, y, color="black", alpha=0.6, linewidth=0.5)
        axis2[i].plot(y2.index - 147, y2, color=sns.color_palette("colorblind")[0], linewidth=1)
        axis2[i].scatter(y2.iloc[peak_indices].index - 147, y2.iloc[peak_indices], s=25, color="orange")
        for peak in y2.iloc[peak_indices].index:
            axis2[i].axvline(peak - 147, color="black", linestyle="--")
        axis2[i].set_title(group)
    sns.despine(fig)
    sns.despine(fig2)
    fig.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "phasograms.dyad_distances.peaks.svg"), bbox_inches="tight")
    fig2.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "phasograms.nucleosome_distances.peaks.svg"), bbox_inches="tight")

    # Get NFR per knockout
    lengths = dict()

    for group in groups:
        print(group)
        # Get NFR calls from nucleoatac
        df = pd.read_csv(os.path.join(self.results_dir, "nucleoatac", group, group + ".nfrpos.bed.gz"), sep="\t", header=None)
        # Get lengths
        lengths[group] = (df[2] - df[1]).tolist()

    pickle.dump(lengths, open(os.path.join(self.results_dir, "nucleoatac", "nfr.lengths.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)
    lengths = pickle.load(open(os.path.join(self.results_dir, "nucleoatac", "nfr.lengths.pickle"), "rb"))

    # plot NFR lengths
    from scipy.ndimage.filters import gaussian_filter1d
    n_rows = n_cols = int(np.ceil(np.sqrt(len(groups))))
    n_rows -= 1
    fig, axis = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(n_cols * 3, n_rows * 2))
    axis = axis.flatten()
    for i, group in enumerate(groups):
        # Count NFR lengths
        x = pd.Series(lengths[group])
        y = x.value_counts().sort_index()
        y = y.ix[range(plotting_window[0], 300)]
        y /= y.sum()

        # Find peaks
        y2 = pd.Series(gaussian_filter1d(y, 5), index=y.index)
        peak_indices = [detect_peaks(y2.values, mpd=73.5)[0]]
        print(group, y2.iloc[peak_indices].index)

        # Plot distribution and peaks
        axis[i].plot(y.index, y, color="black", alpha=0.6, linewidth=0.5)
        axis[i].plot(y2.index, y2, color=sns.color_palette("colorblind")[0], linewidth=1)
        axis[i].scatter(y2.iloc[peak_indices].index, y2.iloc[peak_indices], s=25, color="orange")
        for peak in y2.iloc[peak_indices].index:
            axis[i].axvline(peak, color="black", linestyle="--")
        axis[i].set_title(group)

    sns.despine(fig)
    fig.savefig(os.path.join(self.results_dir, "nucleoatac", "plots", "phasograms.nfr_lengths.peaks.svg"), bbox_inches="tight")


def run_coverage_job(bed_file, bam_file, coverage_type, name, output_dir, window_size=1001):
    from pypiper import NGSTk
    tk = NGSTk()
    job_file = os.path.join(output_dir, "%s.run_enrichment.sh" % name)
    log_file = os.path.join(output_dir, "%s.run_enrichment.log" % name)

    cmd = """#!/bin/bash
#SBATCH --partition=shortq
#SBATCH --ntasks=1
#SBATCH --time=12:00:00

#SBATCH --cpus-per-task=2
#SBATCH --mem=8000
#SBATCH --nodes=1

#SBATCH --job-name=baf-kubicek-run_enrichment_{}
#SBATCH --output={}

#SBATCH --mail-type=end
#SBATCH --mail-user=

# Start running the job
hostname
date

cd /home/arendeiro/baf-kubicek/

python /home/arendeiro/jobs/run_profiles.py \
--bed-file {} \
--bam-file {} \
--coverage-type {} \
--window-size {} \
--name {} \
--output-dir {}

date
""".format(
        name,
        log_file,
        bed_file,
        bam_file,
        coverage_type,
        window_size,
        name,
        output_dir
    )

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(cmd)

    tk.slurm_submit_job(job_file)


def run_vplot_job(bed_file, bam_file, name, output_dir):
    from pypiper import NGSTk
    tk = NGSTk()
    job_file = os.path.join(output_dir, "%s.run_enrichment.sh" % name)
    log_file = os.path.join(output_dir, "%s.run_enrichment.log" % name)

    cmd = """#!/bin/bash
#SBATCH --partition=mediumq
#SBATCH --ntasks=1
#SBATCH --time=1-12:00:00

#SBATCH --cpus-per-task=8
#SBATCH --mem=24000
#SBATCH --nodes=1

#SBATCH --job-name=baf-kubicek-run_enrichment_{name}
#SBATCH --output={log}

#SBATCH --mail-type=end
#SBATCH --mail-user=

# Start running the job
hostname
date

\t\t# Remove everything to do with your python and env.  Even reset your home dir
\t\tunset PYTHONPATH
\t\tunset PYTHON_HOME
\t\tmodule purge
\t\tmodule load python/2.7.6
\t\tmodule load slurm
\t\tmodule load gcc/4.8.2

\t\tENV_DIR=/scratch/users/arendeiro/nucleoenv
\t\texport HOME=$ENV_DIR/home

\t\t# Activate your virtual env
\t\texport VIRTUALENVWRAPPER_PYTHON=/cm/shared/apps/python/2.7.6/bin/python
\t\tsource $ENV_DIR/bin/activate

\t\t# Prepare to install new python packages
\t\texport PATH=$ENV_DIR/install/bin:$PATH
\t\texport PYTHONPATH=$ENV_DIR/install/lib/python2.7/site-packages


cd /home/arendeiro/baf-kubicek/

pyatac vplot \
--bed {peaks} \
--bam {bam} \
--out {output_prefix}.vplot \
--cores 8 \
--lower 30 \
--upper 1000 \
--flank 500 \
--scale \
--plot_extra

pyatac sizes \
--bam {bam} \
--bed {peaks} \
--out {output_prefix}.sizes  \
--lower 30 \
--upper 1000

pyatac bias \
--fasta ~/resources/genomes/hg19/hg19.fa \
--bed {peaks} \
--out {output_prefix}.bias \
--cores 8

pyatac bias_vplot \
--bed {peaks} \
--bg {output_prefix}.bias.Scores.bedgraph.gz \
--fasta ~/resources/genomes/hg19/hg19.fa \
--sizes {output_prefix}.sizes.fragmentsizes.txt \
--out {output_prefix}.bias_vplot \
--cores 8 \
--lower 30 \
--upper 1000 \
--flank 500 \
--scale \
--plot_extra

date
""".format(
        name=name,
        log=log_file,
        peaks=bed_file,
        bam=bam_file,
        output_prefix=os.path.join(output_dir, name))

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(cmd)

    tk.slurm_submit_job(job_file)


def piq_prepare_motifs(
        motifs_file="data/external/jaspar_human_motifs.txt",
        output_dir="footprinting",
        piq_source_dir="/home/arendeiro/workspace/piq-single/",
        motif_numbers=None):
    """
    Prepare motifs for footprinting with PIQ.
    This is typically done only once per genome.

    `motifs_file` is a JASPAR-format file with PWMs of motifs.
    `output_dir` is a root directory to use for the output of a footprint analysis.
    `piq_source_dir` is the root directory of the PIQ code.
    `motif_numbers` is a range of integers enumerating the TF motifs from
    the `motifs_file` to prepare. This is by default, all in the file.
    """
    import textwrap
    from pypiper import NGSTk
    tk = NGSTk()

    output_dir = os.path.abspath(output_dir)
    motifs_dir = os.path.join(output_dir, "motifs")
    jobs_dir = os.path.join(motifs_dir, "jobs")

    for path in [output_dir, motifs_dir, jobs_dir]:
        if not os.path.exists(path):
            os.makedirs(path)

    if motif_numbers is None:
        n_motifs = open(motifs_file, "r").read().count(">")
        motif_numbers = range(1, n_motifs + 1)

    for motif in motif_numbers:
        # skip if exists
        a = os.path.exists(os.path.join(motifs_dir, "{}.pwmout.RData".format(motif)))
        b = os.path.exists(os.path.join(motifs_dir, "{}.pwmout.rc.RData".format(motif)))
        if a and b:
            continue

        # prepare job
        log_file = os.path.join(jobs_dir, "piq_preparemotifs.motif{}.slurm.log".format(motif))
        job_file = os.path.join(jobs_dir, "piq_preparemotifs.motif{}.slurm.sh".format(motif))
        cmd = tk.slurm_header(
            "piq_preparemotifs.{}".format(motif), log_file,
            cpus_per_task=1, queue="shortq", mem_per_cpu=8000)

        # change to PIQ dir (required in PIQ because of hard-coded links in code)
        cmd += """cd {}\n\t\t""".format(piq_source_dir)

        # Actual PIQ command
        cmd += """Rscript {}/pwmmatch.exact.r""".format(piq_source_dir)
        cmd += " {}/common.r".format(piq_source_dir)
        cmd += " {}".format(os.path.abspath(motifs_file))
        cmd += " " + str(motif)
        cmd += """ {}\n""".format(motifs_dir)

        # write job to file
        with open(job_file, 'w') as handle:
            handle.writelines(textwrap.dedent(cmd))

        tk.slurm_submit_job(job_file)


def piq_prepare_bams(
        bam_files,
        group_name,
        output_dir="footprinting",
        piq_source_dir="/home/arendeiro/workspace/piq-single/"):
    """
    Prepare single or group of BAM files for footprinting.

    `bam_files` is a string path to a BAM file or a list of the same.
    `group_name` is a label matching the one given to `piq_prepare_bams`.
    `output_dir` is a root directory to use for the output of a footprint analysis.
    `piq_source_dir` is the root directory of the PIQ code.

    This will launch jobs for each motif.
    """
    import textwrap
    from pypiper import NGSTk
    tk = NGSTk()

    if type(bam_files) is str:
        bam_files = [bam_files]

    output_dir = os.path.abspath(output_dir)
    cache_dir = os.path.join(output_dir, "piq_cache")
    jobs_dir = os.path.join(cache_dir, "jobs")

    for path in [output_dir, cache_dir, jobs_dir]:
        if not os.path.exists(path):
            os.makedirs(path)

    merged_bam = os.path.join(cache_dir, "{}.merged.bam".format(group_name))
    merged_cache = os.path.join(cache_dir, "{}.merged.RData".format(group_name))
    log_file = os.path.join(jobs_dir, "piq_preparebams.{}.slurm.log".format(group_name))
    job_file = os.path.join(jobs_dir, "piq_preparebams.{}.slurm.sh".format(group_name))

    # Build job
    cmd = tk.slurm_header(
        "piq_preparebams.{}".format(group_name), log_file,
        cpus_per_task=4, queue="mediumq", time="11:00:00", mem_per_cpu=8000
    )

    # merge all bam files
    if not os.path.exists(merged_bam):
        cmd += """sambamba merge -t 12 {0} {1}\n\t\t""".format(merged_bam, " ".join(bam_files))

    # change to PIQ dir (required in PIQ because of hard-coded links in code)
    cmd += """cd {}\n""".format(piq_source_dir)

    # Actual PIQ command
    cmd += """\t\tRscript {0}/bam2rdata.r {0}/common.r {1} {2}\n""".format(
        piq_source_dir, merged_cache, merged_bam)

    # slurm footer
    cmd += tk.slurm_footer()

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurm_submit_job(job_file)


def footprint(
        group_name, motif_numbers,
        output_dir="footprinting",
        piq_source_dir="/home/arendeiro/workspace/piq-single/",
        total_job_limit=250,
        min_time_between_jobs=5,
        refresh_time=10):
    """
    Perform TF footprinting with PIQ.

    `group_name` is a label matching the one given to `piq_prepare_bams`.
    `motif_numbers` is a range of integers enumerating the TF motifs from
    the motif file given to `piq_prepare_motifs` to footprint.
    `output_dir` is a root directory to use for the input/output of a footprint analysis.
    `piq_source_dir` is the root directory of the PIQ code.

    This will launch jobs in `min_time_between_jobs` (seconds) intervals
    until the whole job queue has `total_job_limit` jobs
    and retry after `refresh_time` (seconds).
    """
    import os
    import pybedtools
    import pandas as pd
    import subprocess
    import time
    import textwrap
    from pypiper import NGSTk

    tk = NGSTk()

    output_dir = os.path.abspath(output_dir)
    motifs_dir = os.path.join(output_dir, "motifs/")  # slash is important - PIQ reasons
    cache_dir = os.path.join(output_dir, "piq_cache")
    foots_dir = os.path.join(output_dir, "footprint_calls")
    tmp_dir = os.path.join(foots_dir, "tmp_output")
    jobs_dir = os.path.join(foots_dir, "jobs")
    merged_cache = os.path.join(cache_dir, "{}.merged.RData".format(group_name))

    for path in [output_dir, motifs_dir, cache_dir, foots_dir, tmp_dir, jobs_dir]:
        if not os.path.exists(path):
            os.makedirs(path)

    for motif_number in motif_numbers:
        if not os.path.exists(os.path.join(motifs_dir, "{}.pwmout.RData".format(motif_number))):
            print("PIQ file for motif {} does not exist".format(motif_number))
            continue

        t_dir = os.path.join(tmp_dir, group_name)
        o_dir = os.path.join(foots_dir, group_name)
        for path in [t_dir, o_dir]:
            if not os.path.exists(path):
                os.mkdir(path)
        log_file = os.path.join(jobs_dir, "piq_footprinting.{}.motif{}.slurm.log".format(group_name, motif_number))
        job_file = os.path.join(jobs_dir, "piq_footprinting.{}.motif{}.slurm.sh".format(group_name, motif_number))

        # prepare slurm job header
        cmd = tk.slurm_header(
            "piq_footprinting.{}.motif{}".format(group_name, motif_number),
            log_file, cpus_per_task=1, queue="shortq", mem_per_cpu=8000)

        # change to PIQ dir (required in PIQ because of hard-coded links in code)
        cmd += """cd {}\n""".format(piq_source_dir)

        # Actual PIQ command
        # footprint
        cmd += """\t\tRscript {0}/pertf.r {0}/common.r {1} {2} {3} {4} {5}\n""".format(
            piq_source_dir, motifs_dir, t_dir, o_dir, merged_cache, motif_number)

        # slurm footer
        cmd += tk.slurm_footer()

        # write job to file
        with open(job_file, 'w') as handle:
            handle.writelines(textwrap.dedent(cmd))

        # submit jobs slowly and keeping total numbe of jobs below a certain ammount
        submit = (len(subprocess.check_output("squeue").split("\n")) - 1) < total_job_limit
        while not submit:
            time.sleep(refresh_time)
            # print("Waited {}. Checking again...".format(refresh_time))
            submit = (len(subprocess.check_output("squeue").split("\n")) - 1) < total_job_limit

        tk.slurm_submit_job(job_file)
        print("Submitted job of group {}, motif {}.".format(group_name, motif_number))
        time.sleep(min_time_between_jobs)


def tfbs_to_gene(
        bed_file,
        tss_file="data/external/ensembl.tss.bed",
        promoter_and_genesbody_file="data/external/ensembl.promoter_and_genesbody.bed"
    ):
    """
    Assign a TF binding site (output of PIQ footprinting) to a gene.
    TFBS are assigned to a gene if they overlap with their promoter or genebody,
    else to the nearest TSS (regardless of the distance).
    This distance can be used later to weight the interaction e.g. in combination
    with the confidence measures of the footprint (binding).

    `bed_file` is a PIQ output 8-column BED file:
    chrom, start, end, pwm, shape, strand, score, purity.

    `tss_file` and `promoter_and_genesbody_file` are 6-column BED files:
    chrom, start, end, gene_id, transcript_id, strand.

    Returns TFBS-gene assignment matrix.
    """
    import os
    import pybedtools
    import pandas as pd

    # read in gene body + promoter info
    promoter_and_genesbody = pybedtools.BedTool(promoter_and_genesbody_file)
    # read in TSS info
    tss = pybedtools.BedTool(tss_file)
    # columns
    columns = ["chrom", "start", "end", "pwm", "shape", "strand", "score", "purity",
               "chrom_gene", "start_gene", "end_gene", "gene", "transcript", "strand_gene"]

    # Assign TFBS to gene if they overlap with gene body or promoter (5kb around TSS -> 2.5kb upstream)
    gene_assignments = (
        pybedtools.BedTool(os.path.join(bed_file))
        .intersect(promoter_and_genesbody, wa=True, wb=True)
        .to_dataframe(names=columns))

    # For the remaining TFBSs, assign TFBS to closest TSS regardless of distance
    # (distance is not so important for assignment because distance is a penalyzing effect during TF-gene interaction score calculation)
    # 1. get genes not assigned previously
    all_ = (
        pybedtools.BedTool(os.path.join(bed_file))
        .to_dataframe(names=columns[:8]))

    merged = pd.merge(all_, gene_assignments, how="left")
    remaining = merged[merged['gene'].isnull()]

    # 2. assign to nearest
    closest_tss = (
        pybedtools.BedTool(
            remaining.iloc[:, range(all_.shape[1])]
            .to_string(index=False, header=False)
            .replace(" ", "\t"),
            from_string=True,
        )
        .closest(tss, d=True)
        .to_dataframe(names=columns + ['distance']))

    # put the two together
    gene_assignments = pd.concat([gene_assignments, closest_tss])

    # set overlapping distance to 0
    gene_assignments.loc[gene_assignments['distance'].isnull(), 'distance'] = 0

    return gene_assignments


def piq_to_network(
        group_name,
        motif_numbers,
        peak_universe_file,
        output_dir="footprinting"):
    """
    Parse PIQ output, filter footprints and score TF-gene interactions.
    Returns matrix with score of each TF regulating each gene.

    `group_name` is a label matching the one given to `piq_prepare_bams`.
    `motif_numbers` is a range of integers enumerating the TF motifs from
    the motif file given to `piq_prepare_motifs` to footprint.
    `peak_universe_file` is a 3-column BED file with high confidence genomic locations
    to filter interactions for (generally a peak set).
    `output_dir` is a root directory to use for the input/output of a footprint analysis.

    Records in `peak_universe_file` are recommended to be enlarged (e.g. up to 500bp in each direction).
    """
    import re
    import os
    import pybedtools
    import pandas as pd

    output_dir = os.path.abspath(output_dir)
    foots_dir = os.path.join(output_dir, "footprint_calls")
    group_foot_dir = os.path.join(foots_dir, group_name)

    for path in [output_dir, foots_dir, group_foot_dir]:
        if not os.path.exists(path):
            os.makedirs(path)

    # list results_dir
    files = os.listdir(group_foot_dir)

    if len(files) == 0:
        print("There are not footprint calls for group '{}' in '{}'".format(group_name, group_foot_dir))

    # use universe set of ATAC-seq peaks to filter data
    all_peaks = pybedtools.BedTool(peak_universe_file)

    # dataframe to store TFBS assignment to genes
    assignments = pd.DataFrame()

    # dataframe to store TF->gene interactions
    interactions = pd.DataFrame()

    # dataframe to store stats about the TFBS and the interactions
    stats = pd.DataFrame()

    # loop through motifs/TFs, filter and establish relationship between TF and gene
    for motif in motif_numbers:
        print("Gathering footprint calls of motif '{}' for group '{}'".format(motif, group_name))
        # get both forward and reverse complement PIQ output files
        result_files = list()
        for f in files:
            m = re.match(r'%i-.*-calls\.csv$' % motif, f)
            if hasattr(m, "string"):
                result_files.append(m.string)

        # make bed file from it
        # concatenate files (forward and reverse complement are treated differently by PIQ)
        for i, result_file in enumerate(result_files):
            df = pd.read_csv(os.path.join(group_foot_dir, result_file), index_col=0)
            df.rename(columns={"coord": "start"}, inplace=True)
            # fix coordinates
            if "RC-calls.csv" not in result_file:
                df["end"] = df["start"] + 1
                df['strand'] = "+"
            else:
                df["end"] = df["start"]
                df["start"] = df["start"] - 1
                df['strand'] = "-"
            # concatenate
            if i == 0:
                df2 = df
            else:
                df2 = pd.concat([df, df2])

        # add total TFBS to stats
        stats.loc[motif, "TFBS"] = len(df2)
        stats.loc[motif, "TFBS_+"] = len(df2[df2['strand'] == "+"])
        stats.loc[motif, "TFBS_-"] = len(df2[df2['strand'] == "-"])

        # Filter for purity
        footprints = df2[df2["purity"] > 0.7]
        stats.loc[motif, "pur0.7"] = len(footprints)

        # If less than 500 significant interactions, ignore TF
        if footprints.shape[0] < 500:
            continue

        # filter for motifs overlapping CLL peaks
        footprints = (
            pybedtools.BedTool(
                (
                    footprints[['chr', 'start', 'end', 'pwm', 'shape', 'strand', 'score', 'purity']]
                    .to_string(header=None, index=False)
                    .replace(' ', '\t')
                ),
                from_string=True)
            .intersect(all_peaks, wa=True)
            .to_dataframe()
        )
        footprints.columns = ['chr', 'start', 'end', 'pwm', 'shape', 'strand', 'score', 'purity']
        footprints.to_csv(os.path.join("tmp.bed"), sep="\t", index=False, header=False)
        stats.loc[motif, "overlap_universe"] = footprints.shape[0]

        # assign TFBS to gene
        gene_assignments = tfbs_to_gene(os.path.join("tmp.bed"))
        gene_assignments.loc[:, "TF"] = motif
        stats.loc[motif, "gene_overlap_count"] = len(gene_assignments[gene_assignments['distance'] == 0])
        stats.loc[motif, "dist_gene_mean"] = gene_assignments['distance'].mean()
        stats.loc[motif, "dist_gene_median"] = gene_assignments['distance'].median()
        stats.loc[motif, "dist_gene_std"] = gene_assignments['distance'].std()

        # Get weighted values
        # weigh with footprint purity and distance to tss
        gene_assignments['interaction_score'] = gene_assignments.apply(lambda x: 2 * (x['purity'] - 0.5) * 10 ** -(x['distance'] / 1000000.), axis=1)
        # sum scores for each gene
        scores = gene_assignments.groupby(['gene'])['interaction_score'].apply(sum).reset_index()
        scores.loc[:, "TF"] = motif

        # filter out potentially not assigned bindings
        scores = scores[scores["gene"] != "."]

        # add mean score for each gene
        stats.loc[motif, "score_gene_mean"] = scores['interaction_score'].mean()
        stats.loc[motif, "score_gene_std"] = scores['interaction_score'].std()
        stats.loc[motif, "TF"] = motif

        # save
        scores.to_csv(os.path.join(group_foot_dir, "scores.motif{}.csv".format(motif)), index=False)
        gene_assignments.to_csv(os.path.join(group_foot_dir, "gene_assignments.motif{}.csv".format(motif)), index=False)

        # add to dataframe with all TF-gene interactions
        interactions = interactions.append(scores, ignore_index=True)
        assignments = assignments.append(gene_assignments, ignore_index=True)

    # save
    assignments.to_csv(os.path.join(group_foot_dir, "assignments.all_motifs.csv".format(motif)), index=False)
    interactions.to_csv(os.path.join(group_foot_dir, "interactions.all_motifs.csv".format(motif)), index=False)
    stats.to_csv(os.path.join(group_foot_dir, "stats.all_motifs.csv".format(motif)), index=False)

    return (assignments, interactions, stats)


def differential_interactions(
        group_name1, group_name2,
        output_dir="footprinting",
        motifs_mapping="data/external/jaspar_human_motifs.id_mapping.txt"):
    """
    Compare TF-gene interactions between two sets of groups (e.g. KO and WT)
    and visualize differences.
    """
    import re
    from scipy.stats import mannwhitneyu
    from statsmodels.sandbox.stats.multicomp import multipletests

    def add_xy_line(axis):
        lims = [np.min([axis.get_xlim(), axis.get_ylim()]), np.max([axis.get_xlim(), axis.get_ylim()])]
        axis.plot(lims, lims, '--', alpha=0.5, zorder=0)
        axis.set_aspect('equal')
        axis.set_xlim(lims)
        axis.set_ylim(lims)

    if group_name1 == group_name2:
        print("The two groups are the same! Skipping...")
        return

    comparison_name = "{}-{}".format(group_name1, group_name2)

    output_dir = os.path.abspath(output_dir)
    foots_dir = os.path.join(output_dir, "footprint_calls")
    diff_dir = os.path.join(output_dir, "differential_calls")

    for path in [output_dir, foots_dir, diff_dir]:
        if not os.path.exists(path):
            os.makedirs(path)

    f1 = pd.read_csv(os.path.join(foots_dir, group_name1, "interactions.all_motifs.csv"))
    f2 = pd.read_csv(os.path.join(foots_dir, group_name2, "interactions.all_motifs.csv"))

    # Global changes in TFs
    # are the TFs binding differently globally on average?
    tf_stats = f1.groupby("TF")['interaction_score'].mean().to_frame(name=group_name1)
    tf_stats = tf_stats.join(
        f2.groupby("TF")['interaction_score'].mean().to_frame(name=group_name2).squeeze(),
        how="outer")

    # calculate some tf_stats
    for tf in tf_stats.index:
        a = f1.loc[f1['TF'] == tf, "interaction_score"]
        b = f2.loc[f2['TF'] == tf, "interaction_score"]
        tf_stats.loc[tf, "group1_log_interactions"] = np.log((a.shape[0] / float(f1.shape[0])) * 1e3)
        tf_stats.loc[tf, "group2_log_interactions"] = np.log((b.shape[0] / float(f2.shape[0])) * 1e3)
        tf_stats.loc[tf, "p_value"] = mannwhitneyu(a, b)[1]
    tf_stats.loc[:, "q_value"] = multipletests(tf_stats.loc[:, "p_value"], method="fdr_bh")[1]
    tf_stats.loc[:, "log_fold_change"] = np.log2(tf_stats.loc[:, group_name1] / tf_stats.loc[:, group_name2])
    tf_stats.loc[:, "a"] = (1 / 2.) * np.log2(tf_stats.loc[:, group_name1] * tf_stats.loc[:, group_name2])

    # annotate TF ids with names
    annot = pd.read_table(motifs_mapping, header=None, names=['tf', "jaspar_id", "tf_name"]).set_index('tf')
    tf_stats = tf_stats.join(annot).set_index("tf_name")

    # save summary
    tf_stats.to_csv(os.path.join(diff_dir, "tf_differential_binding.{}.csv".format(comparison_name)))

    # Scatter
    fig, axis = plt.subplots(1, figsize=(4, 4))
    axis.scatter(
        tf_stats.loc[:, "group1_log_interactions"],
        tf_stats.loc[:, "group2_log_interactions"],
        alpha=0.5,
        rasterized=True)
    add_xy_line(axis)
    axis.set_xlabel("{} interactions (log)".format(group_name1))
    axis.set_ylabel("{} interactions (log)".format(group_name2))
    sns.despine(fig)
    fig.savefig(
        os.path.join(diff_dir, "tf_differential_binding.{}.scatter.log_interactions.svg".format(comparison_name)),
        bbox_inches="tight", dpi=300)

    fig, axis = plt.subplots(1, figsize=(4, 4))
    axis.scatter(
        tf_stats.loc[:, group_name1],
        tf_stats.loc[:, group_name2],
        alpha=0.5,
        rasterized=True)
    add_xy_line(axis)
    axis.set_xlabel("{} interactions mean".format(group_name1))
    axis.set_ylabel("{} interactions mean".format(group_name2))
    sns.despine(fig)
    fig.savefig(
        os.path.join(diff_dir, "tf_differential_binding.{}.scatter.mean_interactions.svg".format(comparison_name)),
        bbox_inches="tight", dpi=300)

    # MA
    fig, axis = plt.subplots(1, figsize=(4, 4))
    axis.scatter(
        tf_stats.loc[:, "a"],
        tf_stats.loc[:, "log_fold_change"],
        alpha=0.5,
        rasterized=True)
    axis.axhline(0, linestyle="--", alpha=0.5, zorder=0)
    axis.set_xlabel("Intensity (A)".format(group_name1))
    axis.set_ylabel("Log2 fold change ({} / {})".format(group_name1, group_name2))
    sns.despine(fig)
    fig.savefig(
        os.path.join(diff_dir, "tf_differential_binding.{}.ma.svg".format(comparison_name)),
        bbox_inches="tight", dpi=300)

    # Volcano
    fig, axis = plt.subplots(1, figsize=(4, 4))
    axis.scatter(
        tf_stats.loc[:, "log_fold_change"],
        -np.log(tf_stats.loc[:, "p_value"]),
        alpha=0.5,
        rasterized=True)
    axis.axvline(0, linestyle="--", alpha=0.5, zorder=0)
    axis.set_xlabel("Log2 fold change ({} / {})".format(group_name1, group_name2))
    axis.set_ylabel("-log p-value")
    sns.despine(fig)
    fig.savefig(
        os.path.join(diff_dir, "tf_differential_binding.{}.volcano.svg".format(comparison_name)),
        bbox_inches="tight", dpi=300)
