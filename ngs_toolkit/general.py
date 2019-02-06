#!/usr/bin/env python


import glob
import json
import os
import subprocess
import textwrap
import warnings

import distutils.spawn
import matplotlib.pyplot as plt
from ngs_toolkit import _CONFIG, _LOGGER
from ngs_toolkit.graphics import savefig
from ngs_toolkit.utils import download_gzip_file, download_file, r2pandas_df
import numpy as np
import pandas as pd
import patsy
from pybedtools import BedTool
import pybedtools
from pypiper.ngstk import NGSTk
import pysam
import requests
import rpy2
from rpy2.rinterface import RRuntimeError, RRuntimeWarning
from rpy2.robjects import numpy2ri, pandas2ri
import rpy2.robjects as robjects
from scipy import stats
from scipy.linalg import lstsq
from scipy.stats import gaussian_kde
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.sandbox.stats.multicomp import multipletests
from tqdm import tqdm


def get_genome_reference(
        organism, genome_assembly=None, output_dir=None,
        genome_provider="UCSC", file_format="2bit", dry_run=False, overwrite=True):
    """
    Get genome FASTA/2bit file.
    Saves results to disk and returns path to file.

    Attributes:

    :param str organism:
        Organism to get annotation for. Currently supported: "human" and "mouse".

    :param str,optional output_dir:
        Directory to write output to.
        Defaults to current directory

    :param str,optional genome_provider:
        Which genome provider to use. One of 'UCSC' or 'Ensembl'.

    :param str,optional file_format:
        File format to get. One of 'fasta' or '2bit'.

    :param bool,optional dry_run:
        Whether to not download and just return path to file.

    :param bool,optional overwrite:
        Whether existing files should be overwritten by new ones.
        Otherwise they will be kept and no action is made.
        Defaults to True.

    :returns str|tuple:
        If not ``dry_run``, path to genome FASTA/2bit file,
        otherwise tuple of URL of reference genome and path to file.
    """
    def index_fasta(fasta):
        """
        # The idea is to use a hidden method of bedtools
        # to create an index (to skip having e.g. samtools as dependency)
        # and use the bedtool nuc command to to do it.
        # This actually fails to get nucleotide content every time due to this 'bug':
        # https://github.com/daler/pybedtools/issues/147
        # but nonetheless creates an index :whatever:
        """
        bed = pd.DataFrame([['chr1', 1, 2]])
        try:
            bed = pybedtools.BedTool().from_dataframe(bed)
            bed.nucleotide_content(fi=fasta)
        except pybedtools.helpers.BEDToolsError:
            pass

    def twobit_to_fasta(genome_file):
        if distutils.spawn.find_executable("twoBitToFa") is not None:
            args = ["twoBitToFa"] + [genome_file, genome_file.replace(".2bit", ".fa")]
            subprocess.check_output(args, shell=False)
            index_fasta(genome_file.replace(".2bit", ".fa"))

    if output_dir is None:
        output_dir = os.path.join(os.path.abspath(os.path.curdir), "reference")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    opts = ['UCSC', 'Ensembl']
    if genome_provider not in opts:
        msg = "`genome_provider` attribute must be one of '{}'.".format(", ".join(opts))
        _LOGGER.error(msg)
        raise ValueError(msg)

    opts = ['fasta', '2bit']
    if file_format not in opts:
        msg = "`file_format` attribute must be one of '{}'.".format(", ".join(opts))
        _LOGGER.error(msg)
        raise ValueError(msg)

    if (genome_provider == "Ensembl") and (file_format == '2bit'):
        msg = "Ensembl does not provide 2bit files."
        hint = " Use for example 'faToTwoBit' to convert the FASTA file."
        _LOGGER.error(msg + hint)
        raise ValueError(msg)

    if genome_provider == "UCSC":
        organisms = {
            "human": "hg19", "hsapiens": "hg19", "homo_sapiens": "hg19",
            "mouse": "mm10", "mmusculus": "mm10", "mus_musculus": "mm10",
            "yeast": "sacCer3", "scerevisiae": "sacCer3", "saccharomyces_cerevisiae": "sacCer3",
        }
        base_link = "http://hgdownload.cse.ucsc.edu/goldenPath/{assembly}/bigZips/{assembly}"
        base_link += '.fa.gz' if file_format == 'fasta' else '.2bit'
        if genome_assembly is None:
            genome_assembly = organisms[organism]
        url = base_link.format(assembly=genome_assembly)

    elif genome_provider == "Ensembl":
        organisms = {
            "human": {"long_species": "homo_sapiens", "version": "grch37", "release": "75"},
            "hsapiens": {"long_species": "homo_sapiens", "version": "grch37", "release": "75"},
            "homo_sapiens": {"long_species": "homo_sapiens", "version": "grch37", "release": "75"},
            "mouse": {"long_species": "mus_musculus", "version": "grcm38", "release": "94"},
            "mmusculus": {"long_species": "mus_musculus", "version": "grcm38", "release": "94"},
            "mus_musculus": {"long_species": "mus_musculus", "version": "grcm38", "release": "94"},
            "yeast": {"long_species": "saccharomyces_cerevisiae", "version": "R64", "release": "94"},
            "scerevisiae": {"long_species": "saccharomyces_cerevisiae", "version": "R64", "release": "94"},
            "saccharomyces_cerevisiae": {"long_species": "saccharomyces_cerevisiae", "version": "R64", "release": "94"}
        }
        if genome_assembly is None:
            genome_assembly = organisms[organism]['version'].replace("grc", "GRC")
        base_link = "ftp://ftp.ensembl.org/pub/release-{release}/fasta/{long_organism}/dna/"
        base_link += "{Clong_organism}.{assembly}."
        base_link += "{}.".format(organisms[organism]['release']) if genome_assembly.endswith("37") else ""
        base_link += "{sequence_type}.{id_type}.{id}fa.gz".format(
            sequence_type="dna", id_type="primary_assembly", id="")
        url = base_link.format(release=organisms[organism]['release'],
                               long_organism=organisms[organism]['long_species'],
                               Clong_organism=organisms[organism]['long_species'].capitalize(),
                               assembly=genome_assembly)

    if (genome_provider == "UCSC") and (file_format == 'fasta') and (genome_assembly != 'hg38'):
        msg = "UCSC does not provide FASTA files for the {} assembly.".format(genome_assembly)
        hint = " Download a 2bit file and use for example 'TwoBitToFa' to convert."
        _LOGGER.error(msg + hint)
        raise ValueError(msg)

    genome_file = os.path.join(output_dir, "{}.{}.{}".format(
        organism, genome_assembly, 'fa.gz' if file_format == 'fasta' else '2bit'))

    if os.path.exists(genome_file) and (not overwrite):
        msg = "Genome file already exists and 'overwrite' is set to False."
        hint = " Returning existing file: {}".format(genome_file)
        _LOGGER.warn(msg + hint)
        # even so, if 2bit and FASTA not there try to get fasta
        if file_format == "2bit":
            if not os.path.exists(genome_file.replace(".2bit", ".fa")):
                twobit_to_fasta(genome_file)
        return genome_file

    # create .fai index for fasta file
    if file_format == 'fasta':
        if not dry_run:
            download_gzip_file(url, genome_file)
            index_fasta(genome_file)
            return genome_file
    else:
        if not dry_run:
            download_file(url, genome_file)
            twobit_to_fasta(genome_file)
            return genome_file

    return (url, genome_file)


def get_blacklist_annotations(
        organism, genome_assembly=None, output_dir=None, overwrite=True):
    """
    Get annotations of blacklisted genomic regions for a given organism/genome assembly.
    Saves results to disk and returns a path to a BED file.

    Attributes:

    :param str organism:
        Organism to get annotation for. Currently supported: "human" and "mouse".

    :param str,optional genome_assembly:
        Ensembl assembly/version to use.
       Default for "human" is "hg19/grch37" and for "mouse" is "mm10/grcm38".

    :param str,optional output_dir:
        Directory to write output to.
        Defaults to "reference" in current directory.

    :param bool,optional overwrite:
        Whether existing files should be overwritten by new ones.
        Otherwise they will be kept and no action is made.
        Defaults to True.

    :returns str:
        Path to blacklist BED file
    """
    if output_dir is None:
        output_dir = os.path.join(os.path.abspath(os.path.curdir), "reference")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    organisms = {
        "human": "hg19",
        "mouse": "mm10"}
    if genome_assembly is None:
        genome_assembly = organisms[organism]
        _LOGGER.warn("Genome assembly not selected. Using assembly '{}' for '{}'."
                     .format(genome_assembly, organism))

    output_file = os.path.join(output_dir, "{}.{}.blacklist.bed"
                               .format(organism, genome_assembly))
    if os.path.exists(output_file) and (not overwrite):
        msg = "Annotation file already exists and 'overwrite' is set to False."
        hint = " Returning existing annotation file: {}".format(output_file)
        _LOGGER.warn(msg + hint)
        return output_file

    url = "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/"
    if genome_assembly not in ['hg19']:
        url += "{0}-{1}/{0}.blacklist.bed.gz".format(genome_assembly, organism)
    else:
        url += "{0}-human/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz".format(genome_assembly)

    download_gzip_file(url, output_file)
    return output_file


def get_tss_annotations(
        organism, genome_assembly=None, save=True, output_dir=None, chr_prefix=True,
        gene_types=['protein_coding', 'processed_transcript', 'lincRNA', 'antisense'],
        overwrite=True):
    """
    Get annotations of TSS for a given organism/genome assembly.
    This is a simple approach using Biomart's API querying the Ensembl database.
    Saves results to disk and returns a dataframe.

    Attributes:

    :param str organism:
        Organism to get annotation for.
        Currently supported: "human" and "mouse".

    :param str,optional genome_assembly:
        Ensembl assembly/version to use.
        Default for "human" is "grch37" and for "mouse" is "grcm38".

    :param bool,optional save:
        Whether to save to disk under ``output_dir``.
        Defaults to True.

    :param str,optional output_dir:
        Directory to write output to.
        Defaults to "reference" in current directory.

    :param bool,optional chr_prefix:
        Whether chromosome names should have the "chr" prefix.
        Defaults to True

    :param list,optional gene_types:
        Subset of transcript biotypes to keep.
      See here the available biotypes https://www.ensembl.org/Help/Faq?id=468
      Defaults to 'protein_coding', 'processed_transcript', 'lincRNA', 'antisense'.

    :param bool,optional overwrite:
        Whether existing files should be overwritten by new ones.
        Otherwise they will be kept and no action is made.
        Defaults to True.

    :returns pandas.DataFrame:
        DataFrame with genome annotations
    """
    organisms = {
        "human": {"species": "hsapiens", "ensembl_version": "grch37"},
        "mouse": {"species": "mmusculus", "ensembl_version": "grcm38"},
        "yeast": {"species": "scerevisiae", "ensembl_version": "R64"}
    }

    if genome_assembly is None:
        genome_assembly = organisms[organism]['ensembl_version']
    if genome_assembly == "hg19":
        genome_assembly = "grch37"
    if genome_assembly == "hg38":
        genome_assembly = "grch38"
    if genome_assembly == "mm10":
        genome_assembly = "grcm38"
    if genome_assembly == "sacCer3":
        genome_assembly = "R64"

    output_file = os.path.join(output_dir, "{}.{}.gene_annotation.tss.bed"
                               .format(organism, genome_assembly))
    if os.path.exists(output_file) and (not overwrite):
        msg = "Annotation file already exists and 'overwrite' is set to False."
        hint = " Returning existing annotation from file: {}".format(output_file)
        _LOGGER.warn(msg + hint)
        annot = pd.read_csv(output_file, header=None, sep="\t")
        return annot

    if save:
        if output_dir is None:
            output_dir = os.path.join(os.path.abspath(os.path.curdir), "reference")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    attributes = [
        "chromosome_name", "start_position",
        "ensembl_gene_id", "external_gene_name",
        "strand", "transcript_biotype"]
    res = query_biomart(
        attributes=attributes,
        species=organisms[organism]['species'],
        ensembl_version=genome_assembly)

    if gene_types is None:
        res = res.drop(['transcript_biotype'], axis=1).drop_duplicates()
    else:
        res = res.loc[res['transcript_biotype'].isin(gene_types), :]
    if chr_prefix:
        res.loc[:, 'chromosome_name'] = "chr" + res['chromosome_name']
    res.loc[:, 'start_position'] = res['start_position'].astype(int)
    res.loc[:, 'strand'] = res['strand'].replace("1", "+").replace("-1", "-")
    res.loc[:, 'end'] = res.apply(
        lambda x: x['start_position'] + 1 if x['strand'] == "+" else x['start_position'],
        axis=1)
    res.loc[:, 'start_position'] = res.apply(
        lambda x: x['start_position'] - 1 if x['strand'] == "-" else x['start_position'],
        axis=1)

    # drop same gene duplicates if starting in same position but just annotated with
    # different biotypes
    res = (
        res[attributes[:2] + ["end"] + attributes[2:]]
        .sort_values(by=res.columns.tolist(), axis=0)
        .drop_duplicates(subset=res.columns[:3], keep="last"))

    # make real BED format
    res.loc[:, 'fill'] = '.'
    cols = [
        'chromosome_name', 'start_position', 'end', 'external_gene_name',
        'fill', 'strand', 'ensembl_gene_id', 'transcript_biotype']
    res = res.loc[:, cols]
    res.columns = [
        'chr', 'start', 'end', 'gene_name', 'score', 'strand',
        'ensembl_gene_id', 'transcript_biotype']

    # save
    if save:
        res.to_csv(output_file, index=False, header=False, sep="\t")

        output_file = os.path.join(output_dir, "{}.{}.gene_annotation.protein_coding.tss.bed"
                                   .format(organism, genome_assembly))
        res[res['transcript_biotype'] == "protein_coding"].drop(attributes[-1], axis=1).to_csv(
            output_file, index=False, header=False, sep="\t")

    return res


def get_genomic_context(
        organism, genome_assembly=None, save=True, output_dir=None, chr_prefix=True,
        region_subset=['promoter', 'exon', '5utr', '3utr', 'intron', 'genebody', 'intergenic'],
        gene_types=['protein_coding', 'processed_transcript', 'lincRNA', 'antisense'],
        promoter_width=3000, overwrite=True):
    """
    Get annotations of TSS for a given organism/genome assembly.
    This is a simple approach using Biomart's API querying the Ensembl database.
    Saves results to disk and returns a dataframe.

    The API call to BioMart can take a bit, so the function should take ~4 min for a human genome.

    Attributes:

    :param str organism:
        Organism to get annotation for. Currently supported: "human" and "mouse".

    :param str,optional genome_assembly:
        Ensembl assembly/version to use.
        Default for "human" is "grch37" and for "mouse" is "grcm38".

    :param bool,optional save:
        Whether to save to disk under ``output_dir``.
        Defaults to True.

    :param str,optional output_dir:
        Directory to write output to.
        Defaults to "reference" in current directory.

    :param bool,optional chr_prefix:
        Whether chromosome names should have the "chr" prefix. Defaults to True

    :param list,optional gene_types:
        Subset of transcript biotypes to keep.
        See here the available biotypes https://www.ensembl.org/Help/Faq?id=468
        Defaults to 'protein_coding', 'processed_transcript', 'lincRNA', 'antisense'.

    :param bool,optional overwrite:
        Whether existing files should be overwritten by new ones.
        Otherwise they will be kept and no action is made.
        Defaults to True.

    :returns pandas.DataFrame:
        DataFrame with genome annotations
    """
    organisms = {
        "human": {"species": "hsapiens", "ensembl_version": "grch37", "ucsc_version": "hg19"},
        "mouse": {"species": "mmusculus", "ensembl_version": "grcm38", "ucsc_version": "mm10"},
        "yeast": {"species": "scerevisiae", "ensembl_version": "R64"}
    }

    if genome_assembly is None:
        genome_assembly = organisms[organism]['ucsc_version']
        ensembl_genome_assembly = organisms[organism]['ensembl_version']
    elif genome_assembly == "hg19":
        ensembl_genome_assembly = "grch37"
    elif genome_assembly == "hg38":
        ensembl_genome_assembly = "grch38"
    elif genome_assembly == "mm10":
        ensembl_genome_assembly = "grcm38"
    elif genome_assembly == "sacCer3":
        ensembl_genome_assembly = "R64"
    else:
        _LOGGER.warning()
        ensembl_genome_assembly = genome_assembly

    output_file = os.path.join(output_dir, "{}.{}.genomic_context.bed"
                               .format(organism, ensembl_genome_assembly))
    if os.path.exists(output_file) and (not overwrite):
        msg = "Annotation file already exists and 'overwrite' is set to False."
        hint = " Returning existing annotation from file: {}".format(output_file)
        _LOGGER.warn(msg + hint)
        annot = pd.read_csv(output_file, header=None, sep="\t")
        return annot

    if save:
        if output_dir is None:
            output_dir = os.path.join(os.path.abspath(os.path.curdir), "reference")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    attrs = [
        "ensembl_gene_id", "chromosome_name", "start_position", "end_position",
        "exon_chrom_start", "exon_chrom_end",
        "5_utr_start", "5_utr_end", "3_utr_start", "3_utr_end",
        "strand", "external_gene_name", "gene_biotype"]
    res = query_biomart(
        attributes=attrs,
        species=organisms[organism]['species'],
        ensembl_version=ensembl_genome_assembly)

    if chr_prefix:
        res['chromosome_name'] = "chr" + res['chromosome_name']

    # Keep only desired gene types
    res = res.loc[res['gene_biotype'].isin(gene_types), :]
    # Get strand
    res.loc[:, 'strand'] = res['strand'].replace("1", "+").replace("-1", "-")
    # remove exons which are also UTRs

    bed_cols = ['chromosome_name', 'start_position', 'end_position', 'strand']
    # # Promoters = genebody:start +/- N
    promoter = res[bed_cols].drop_duplicates()
    for col in promoter.columns:
        if ("start" in col) or ("end" in col):
            promoter.loc[:, col] = promoter.loc[:, col].astype(int)

    promoter.loc[promoter['strand'] == "+", 'start_position'] = (
        promoter.loc[promoter['strand'] == "+", 'start_position'] - (promoter_width / 2))
    promoter.loc[promoter['strand'] == "-", 'start_position'] = (
        promoter.loc[promoter['strand'] == "-", 'end_position'] - (promoter_width / 2))  # + 1
    promoter.loc[promoter['strand'] == "+", 'end_position'] = (
        promoter.loc[:, 'start_position'] + (promoter_width / 2))
    for col in bed_cols[1:3]:
        promoter.loc[:, col] = promoter.loc[:, col].astype(int)

    # # # fix where promoter start < 0
    promoter.loc[promoter['start_position'] < 0, 'start_position'] = 0
    # # # remove end < start
    promoter = promoter.loc[~(promoter['start_position'] > promoter['end_position'])]

    # # genebody = start->end + promoter
    gb = res[bed_cols].drop_duplicates()
    for col in gb.columns:
        if ("start" in col) or ("end" in col):
            gb.loc[:, col] = gb.loc[:, col].astype(int)
    gb = gb.append(promoter)
    # for col in bed_cols[1:3]:
    #     gb.loc[:, col] = gb.loc[:, col].astype(int)
    gb = gb.sort_values(gb.columns.tolist())
    genebody_bed = BedTool.from_dataframe(gb)
    genebody_bed = genebody_bed.sort().merge()
    genebody = genebody_bed.to_dataframe()

    # # Exon
    exon = (
        res[["chromosome_name", "exon_chrom_start", "exon_chrom_end"]]
        .drop_duplicates().dropna())
    # # 5utr
    utr5 = (
        res[["chromosome_name", "5_utr_start", "5_utr_end"]]
        .drop_duplicates().dropna())
    # # 3utr
    utr3 = (
        res[["chromosome_name", "3_utr_start", "3_utr_end"]]
        .drop_duplicates().dropna())
    for d in [exon, utr5, utr3]:
        for col in d.columns:
            if ("start" in col) or ("end" in col):
                d.loc[:, col] = d.loc[:, col].astype(int)

    # # Introns = genebody - (promoter + exon + 5utr + 3utr)
    promoter = promoter.drop(['strand'], axis=1)
    exon.columns = utr3.columns = utr5.columns = promoter.columns = bed_cols[:3]
    others = promoter.append(exon).append(utr5).append(utr3)
    for col in others.columns[1:]:
        others.loc[:, col] = others.loc[:, col].astype(int)
    intron = genebody_bed.subtract(BedTool.from_dataframe(others))
    intron = intron.to_dataframe()

    # # Intergenic = ChromSizes - genebody
    cs = pybedtools.get_chromsizes_from_ucsc(genome_assembly)
    cs = pd.DataFrame(cs).T.reset_index()
    if not chr_prefix:
        cs['index'] = cs['index'].str.subtract("chr", "")
    cs = BedTool.from_dataframe(cs)
    intergenic = cs.subtract(genebody_bed)
    intergenic = intergenic.to_dataframe()

    # join all
    features = [
        "promoter", "genebody", "exon", "intron",
        "utr5", "utr3", "intergenic"]
    annot = list()
    for name in features:
        cur = locals()[name]
        cur = cur.iloc[:, list(range(3))]
        cur.columns = bed_cols[:3]
        cur.loc[:, "region_type"] = name
        if save:
            cur.sort_values(cur.columns.tolist()).drop('region_type', axis=1).to_csv(
                os.path.join(output_dir, "{}.{}.genomic_context.{}.bed"
                             .format(organism, ensembl_genome_assembly, name)),
                index=False, header=False, sep="\t")
        annot.append(cur)

    annot = pd.concat(annot)
    annot = annot.sort_values(annot.columns.tolist())

    # save
    if save:
        annot.to_csv(output_file, index=False, header=False, sep="\t")
    return annot


def count_reads_in_intervals(bam, intervals):
    """
    Count total number of reads in a iterable holding strings
    representing genomic intervals of the form ``"chrom:start-end"``.

    Attributes:

    :param str bam:
        Path to BAM file.

    :param list intervals:
        List of strings with genomic coordinates in format
        ``"chrom:start-end"``.

    :returns dict:
        Dict of read counts for each interval.
    """
    counts = dict()

    bam = pysam.AlignmentFile(bam, mode='rb')

    chroms = ["chr" + str(x) for x in range(1, 23)] + ["chrX", "chrY"]

    for interval in intervals:
        if interval.split(":")[0] not in chroms:
            continue
        counts[interval] = bam.count(region=interval.split("|")[0])
    bam.close()

    return counts


def normalize_quantiles_r(array):
    """
    Quantile normalization with a R implementation.
    Requires the "rpy2" library and the R library "preprocessCore".

    Requires the R package "cqn" to be installed:
        >>> source('http://bioconductor.org/biocLite.R')
        >>> biocLite('preprocessCore')

    Attributes:

    :param numpy.array array:
        Numeric array to normalize.

    :returns numpy.array:
        Normalized numeric array.
    """
    warnings.filterwarnings("ignore", category=RRuntimeWarning)
    rpy2.robjects.numpy2ri.activate()

    robjects.r('require("preprocessCore")')
    normq = robjects.r('normalize.quantiles')
    return np.array(normq(array))


def normalize_quantiles_p(df_input):
    """
    Quantile normalization with a ure Python implementation.
    Code from https://github.com/ShawnLYU/Quantile_Normalize.

    Attributes:

    :param pandas.DataFrame df_input:
        Dataframe to normalize.

    :returns numpy.array:
        Normalized numeric array.
    """
    df = df_input.copy()
    # compute rank
    dic = {}
    for col in df:
        dic.update({col: sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis=1).tolist()
    # sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df


def deseq_analysis(
        count_matrix, experiment_matrix, comparison_table, formula,
        output_dir, output_prefix,
        overwrite=True, alpha=0.05, independent_filtering=False,
        create_subdirectories=True):
    """
    Perform differential comparison analysis with DESeq2.

    .. note::
        Do not include hyphens ("-") in any of the samples or groups names!
        R freaks out with this.

    # TODO: fix hyphens in names issue

    Attributes:

    :param pandas.DataFrame count_matrix:
        Data frame of shape (samples, variables) with raw read counts.

    :param pandas.DataFrame experiment_matrix:
        Data frame with columns "sample_name" and any other variables used in the `formula`.

    :param pandas.DataFrame comparison_table:
        Data frame with columns "comparison_name", "sample_group" and sample_name".

    :param str formula:
        Formula to test in R/patsy notation. Usually something like "~ batch + group".

    :param str output_dir:
        Output directory for produced files.

    :param str output_prefix:
        Prefix to add to produced files.

    :param bool,optional overwrite:
        Whether files existing should be overwritten. Defaults to True.

    :param number,optional alpha:
        Significance level to reject null hypothesis.
        This in practice has no effect as results for all features will be returned.
        Defaults to 0.05.

    :param bool create_subdirectories:
        Whether to create subdirectories for the result of each comparison.

    :returns pandas.DataFrame: Data frame with results, statistics for each feature.
    """
    numpy2ri.activate()
    pandas2ri.activate()
    warnings.filterwarnings("ignore", category=RRuntimeWarning)

    robjects.r('require("DESeq2")')
    _as_formula = robjects.r('as.formula')
    _DESeqDataSetFromMatrix = robjects.r('DESeqDataSetFromMatrix')
    _DESeq = robjects.r('DESeq')
    _results = robjects.r('results')
    _as_data_frame = robjects.r('as.data.frame')

    # order experiment and count matrices in same way
    experiment_matrix = experiment_matrix.set_index("sample_name").loc[count_matrix.columns, :]

    # save the matrices just in case
    count_matrix.to_csv(os.path.join(output_dir, output_prefix + ".count_matrix.tsv"), sep="\t")
    experiment_matrix.to_csv(os.path.join(output_dir, output_prefix + ".experiment_matrix.tsv"), sep="\t")
    comparison_table.to_csv(os.path.join(output_dir, output_prefix + ".comparison_table.tsv"), sep="\t")

    # Rename samples to avoid R errors with sample names containing symbols
    count_matrix.columns = ["S{}".format(i) for i in range(len(count_matrix.columns))]
    experiment_matrix.index = ["S{}".format(i) for i in range(len(experiment_matrix.index))]
    # Replace hyphens with underscores
    experiment_matrix = experiment_matrix.replace("-", "_")
    comparison_table = comparison_table.replace("-", "_")

    # Run DESeq analysis
    dds = _DESeqDataSetFromMatrix(
        countData=count_matrix.astype(int),
        colData=experiment_matrix,
        design=_as_formula(formula))
    dds = _DESeq(dds, parallel=True)
    # _save(dds, file=os.path.join(output_dir, output_prefix + ".deseq_dds_object.Rdata"))

    results = pd.DataFrame()
    comps = comparison_table["comparison_name"].drop_duplicates().sort_values()
    for comp in tqdm(comps, total=len(comps), desc="Comparison"):
        if create_subdirectories:
            if not os.path.exists(os.path.join(output_dir, comp)):
                os.makedirs(os.path.join(output_dir, comp))
            out_file = os.path.join(output_dir, comp, output_prefix + ".deseq_result.{}.csv".format(comp))
        else:
            out_file = os.path.join(output_dir, output_prefix + ".deseq_result.{}.csv".format(comp))

        if not overwrite and os.path.exists(out_file):
            continue
        _LOGGER.info("Doing comparison '{}'".format(comp))
        a = comparison_table.loc[
            (comparison_table["comparison_name"] == comp) &
            (comparison_table["comparison_side"] >= 1),
            "sample_group"].drop_duplicates().squeeze()
        b = comparison_table.loc[
            (comparison_table["comparison_name"] == comp) &
            (comparison_table["comparison_side"] <= 0),
            "sample_group"].drop_duplicates().squeeze()

        contrast = np.array(["sample_group", a, b])
        try:
            res = _as_data_frame(
                _results(dds, contrast=contrast, alpha=alpha,
                         independentFiltering=independent_filtering, parallel=True))
        except RRuntimeError as e:
            _LOGGER.warning("DESeq2 group contrast '{}' didn't work!".format(contrast))
            _LOGGER.error(e)
            contrast = ["sample_group" + a, "sample_group" + b]
            try:
                res = _as_data_frame(
                    _results(dds, contrast=contrast, alpha=alpha,
                             independentFiltering=independent_filtering, parallel=True))
                _LOGGER.warning("DESeq2 group contrast '{}' did work now!".format(", ".join(contrast)))
            except RRuntimeError as e2:
                _LOGGER.warning("DESeq2 group contrast '{}' didn't work either!".format(", ".join(contrast)))
                _LOGGER.error(e2)
                raise e2

        if isinstance(res, rpy2.robjects.vectors.DataFrame):
            # convert to pandas dataframe
            res = r2pandas_df(res)

        res.loc[:, "comparison_name"] = comp
        res.index.name = "index"

        # save
        res.sort_values('pvalue').to_csv(out_file)
        # append
        results = results.append(res.reset_index(), ignore_index=True)

    # save all
    results.index.name = "index"
    results.to_csv(os.path.join(output_dir, output_prefix + ".deseq_result.all_comparisons.csv"), index=False)

    # return
    return results


def least_squares_fit(
        quant_matrix, design_matrix, test_model,
        null_model="~ 1", standardize_data=True,
        multiple_correction_method="fdr_bh"):
    """
    Fit a least squares model with only categorical predictors.
    Computes p-values by comparing the log likelihood ratio of the chosen model to a `null_model`.

    Attributes:

    :param pandas.DataFrame quant_matrix:
        A Data frame of shape (samples, variables).

    :param pandas.DataFrame design_matrix:
        A Data frame of shape (samples, variables) with all the variables in `test_model`.

    :param str test_model:
        Model design to test in R/patsy notation.

    :param str,optional null_model:
        Null model design in R/patsy notation. Defaults to "~ 1".

    :param bool,optional standardize_data:
        Whether data should be standardized prior to fitting. Defaults to True.

    :param str,optional multiple_correction_method:
        Method to use for multiple test correction.
        See statsmodels.sandbox.stats.multicomp.multipletests. Defaults to "fdr_bh".

    :returns pandas.DataFrame:
        Statistics of model fitting and comparison between models for each feature.

    :Example:

    quant_matrix = np.random.random(10000000).reshape(100, 100000)
    P = np.concatenate([[0] * 50, [1] * 50])  # dependent variable
    Q = np.concatenate([[0] * 25, [1] * 25] + [[0] * 25, [1] * 25])  # covariate
    design_matrix = pd.DataFrame([P, Q], index=["P", "Q"]).T
    quant_matrix = quant_matrix.T * ((1 + design_matrix.sum(axis=1)) * 4).values
    quant_matrix = pd.DataFrame(quant_matrix.T)
    test_model = "~ Q + P"
    null_model = "~ Q"
    res = least_squares_fit(quant_matrix, design_matrix, test_model, null_model)
    res.head()

    """
    if standardize_data:
        norm = StandardScaler()
        quant_matrix = pd.DataFrame(
            norm.fit_transform(quant_matrix),
            index=quant_matrix.index, columns=quant_matrix.columns)

    A1 = patsy.dmatrix(test_model, design_matrix)
    betas1, residuals1, _, _ = lstsq(A1, quant_matrix)

    A0 = patsy.dmatrix(null_model, design_matrix)
    betas0, residuals0, _, _ = lstsq(A0, quant_matrix)

    results = pd.DataFrame(betas1.T, columns=A1.design_info.column_names, index=quant_matrix.columns)

    # Calculate the log-likelihood ratios
    n = float(quant_matrix.shape[0])
    results['model_residuals'] = residuals1
    results['null_residuals'] = residuals0
    results['model_log_likelihood'] = (-n / 2.) * np.log(2 * np.pi) - n / 2. * np.log(results['model_residuals'] / n) - n / 2.
    results['null_log_likelihood'] = (-n / 2.) * np.log(2 * np.pi) - n / 2. * np.log(results['null_residuals'] / n) - n / 2.

    results['log_likelihood_ratio'] = results['model_log_likelihood'] - results['null_log_likelihood']
    results['D_statistic'] = 2 * results['log_likelihood_ratio']
    results['p_value'] = stats.chi2.sf(results['log_likelihood_ratio'], df=betas1.shape[0] - betas0.shape[0])
    results['q_value'] = multipletests(results['p_value'], method=multiple_correction_method)[1]

    if not standardize_data:
        results["mean"] = quant_matrix.mean(axis=0)

    return results


def differential_from_bivariate_fit(
        comparison_table, matrix,
        output_dir, output_prefix,
        n_bins=250, multiple_correction_method="fdr_bh",
        plot=True, palette="colorblind", make_values_positive=False):
    """
    Perform differential analysis using a bivariate gaussian fit
    on the relationship between mean and fold-change for each comparison.

    Attributes:

    :param pandas.DataFrame comparison_table:
        Dataframe with 'comparison_name', 'comparison_side' and 'sample_name', 'sample_group' columns.

    :param pandas.DataFrame matrix:
        Matrix of `n_features, n_samples` with normalized, log-transformed values to perform analysis on.

    :param str output_dir:
        Output directory

    :param str output_prefix:
        Prefix for outputs.

    :param int n_bins:
        Number of bins of mean values along which to standardize fold-changes.

    :param str multiple_correction_method:
        Multiple correction method from `statsmodels.sandbox.stats.multicomp.multipletests`.

    :param bool plot:
        Whether to generate plots.

    :param str palette:
        Color palette to use. This can be any matplotlib palette and is passed to `sns.color_palette`.

    :param bool make_values_positive:
        Whether to transform `matrix` to have minimum value 0. Default False.

    :returns pandas.DataFrame:
        Results of fitting and comparison between groups for each feature.
    """
    comparisons = comparison_table['comparison_name'].drop_duplicates().sort_values()
    if plot:
        fig, axis = plt.subplots(
            2, len(comparisons),
            figsize=(4 * len(comparisons), 4 * 2), sharex=True, sharey='row')

    if make_values_positive:
        matrix = matrix + abs(matrix.min().min())

    results = pd.DataFrame()
    for i, comparison in enumerate(comparisons):
        _LOGGER.info("Doing comparison '{}'".format(comparison))
        out_file = os.path.join(output_dir, output_prefix + ".fit_result.{}.csv".format(comparison))

        sa = comparison_table.loc[
            (comparison_table['comparison_name'] == comparison) &
            (comparison_table['comparison_side'] == 1),
            "sample_name"]
        ga = comparison_table.loc[
            (comparison_table['comparison_name'] == comparison) &
            (comparison_table['comparison_side'] == 1),
            "sample_group"].squeeze()
        sb = comparison_table.loc[
            (comparison_table['comparison_name'] == comparison) &
            (comparison_table['comparison_side'] == 0),
            "sample_name"]
        gb = comparison_table.loc[
            (comparison_table['comparison_name'] == comparison) &
            (comparison_table['comparison_side'] == 0),
            "sample_group"].squeeze()
        a = matrix.loc[:, sa].mean(axis=1)
        a.name = ga
        b = matrix.loc[:, sb].mean(axis=1)
        b.name = gb

        # assemble stats
        res = a.to_frame()
        res = res.join(b)
        res['global_mean'] = matrix.mean(axis=1)
        res['global_std'] = matrix.std(axis=1)
        res['comparison_mean'] = res.mean(axis=1)
        res['comparison_std'] = res.mean(axis=1)
        res['log2FoldChange'] = np.log2(res[ga] / res[gb])
        res['comparison_name'] = comparison

        # standardize fold change
        bounds = np.linspace(0, res['comparison_mean'].max(), n_bins)
        for (start, end) in zip(bounds[:-2], bounds[1: -1]):
            r = res.loc[
                (res['comparison_mean'] > start) &
                (res['comparison_mean'] < end)].index
            v = res.loc[r, 'log2FoldChange']
            res.loc[r, 'norm_log2FoldChange'] = (v - np.nanmean(v)) / np.nanstd(v)

        # let's try a bivariate gaussian kernel
        # separately for positive and negative to avoid biases in center of mass
        kernel = gaussian_kde(res.loc[res['norm_log2FoldChange'] > 0, [
                              "comparison_mean", "norm_log2FoldChange"]].T.values)
        res.loc[res['norm_log2FoldChange'] > 0, "density"] = kernel(
            res.loc[res['norm_log2FoldChange'] > 0, ["comparison_mean", "norm_log2FoldChange"]].T.values)
        kernel = gaussian_kde(res.loc[res['norm_log2FoldChange'] <= 0, [
                              "comparison_mean", "norm_log2FoldChange"]].T.values)
        res.loc[res['norm_log2FoldChange'] <= 0, "density"] = kernel(
            res.loc[res['norm_log2FoldChange'] <= 0, ["comparison_mean", "norm_log2FoldChange"]].T.values)

        # Let's calculate something like an empirical p-value on the density
        res['pvalue'] = (res['density'] - res['density'].min()) / \
            (res['density'].max() - res['density'].min())
        res['padj'] = multipletests(res['pvalue'].fillna(1), method=multiple_correction_method)[1]

        res['direction'] = (res['norm_log2FoldChange'] >= 0).astype(int).replace(0, -1)
        res.to_csv(out_file)

        if plot:
            axis[0, i].scatter(res["comparison_mean"],
                               res["log2FoldChange"],
                               alpha=0.2, s=5, color=sns.color_palette(palette)[0], rasterized=True)
            axis[0, i].axhline(0, color="black", linestyle="--")
            axis[1, i].scatter(res["comparison_mean"],
                               res["norm_log2FoldChange"],
                               alpha=0.2, s=5, color=sns.color_palette(palette)[0], rasterized=True)
            diff = res.loc[(res['pvalue'] < 0.05) & (res['norm_log2FoldChange'].abs() >= 2), :].index
            axis[0, i].scatter(res.loc[diff, "comparison_mean"],
                               res.loc[diff, "log2FoldChange"],
                               alpha=0.2, s=5, color=sns.color_palette(palette)[1], rasterized=True)
            axis[1, i].scatter(res.loc[diff, "comparison_mean"],
                               res.loc[diff, "norm_log2FoldChange"],
                               alpha=0.2, s=5, color=sns.color_palette(palette)[1], rasterized=True)
            axis[1, i].axhline(0, color="black", linestyle="--")
            axis[0, i].set_title(comparison + "\n" + ga + " vs " + gb)
            axis[1, i].set_xlabel("Comparison mean")
            if i == 0:
                axis[0, i].set_ylabel("log2 fold-change")
                axis[1, i].set_ylabel("Norm(log2 fold-change)")

        results = results.append(res.reset_index(), ignore_index=True)

    # save figure
    savefig(fig, os.path.join(output_dir, output_prefix + ".deseq_result.all_comparisons.scatter.svg"))

    # save all
    results = results.set_index("index")
    results.to_csv(os.path.join(output_dir, output_prefix + ".deseq_result.all_comparisons.csv"), index=True)

    return results


# def independent_filtering(df, alpha=0.05, n_quantiles=100):
#     """
#     """
#     req_columns = ["pvalue", "baseMean"]
#     assert all([x in df.columns for x in req_columns])
#     # compute quantiles accross mean and pvalue distributions
#     stats = pd.DataFrame()
#     p = (np.arange(n_quantiles) / float(n_quantiles)) * 100
#     p = np.append(p, 100.)
#     for start, end in zip(p, p[1:]):
#         m = np.log2(1 + df['baseMean'])
#         i = df.loc[
#             (m >= np.percentile(m, start)) &
#             (m <= np.percentile(m, end)), :].index
#         stats.loc[start, "n"] = i.shape[0]
#         stats.loc[start, "mean"] = df.loc[i, "baseMean"].mean()
#         stats.loc[start, "mean_p"] = df.loc[i, "pvalue"].mean()
#         stats.loc[start, "n_sig_p"] = (df.loc[i, "pvalue"] < alpha).sum()
#     # plot
#     fig, axis = plt.subplots(1, 2, figsize=(4 * 2, 4 * 1))
#     axis[0].scatter(stats.index, stats.loc[:, "n_sig_p"])
#     axis[1].scatter(stats.index, -np.log10(stats.loc[:, "mean_p"]))
#     # choose inflection point
#     return


# def fit_curve():
#     def func(x, a, b, c):
#         return a * np.exp(-b * x) + c
#     plt.plot(stats.index, -np.log10(stats['mean_p']), 'b-', label='data')
#     popt, pcov = curve_fit(func, stats.index, -np.log10(stats['mean_p']))
#     plt.plot(stats.index, func(stats.index, *popt), 'r-', label='fit')
#     popt, pcov = curve_fit(func, stats.index, stats['mean_p'], bounds=(0, [3., 2., 1.]))
#     plt.plot(stats.index, func(stats.index, *popt), 'g--', label='fit-with-bounds')
#     plt.xlabel('x')
#     plt.ylabel('y')
#     plt.legend()
#     plt.show()


def lola(bed_files, universe_file, output_folder, output_prefixes=None, genome="hg19", cpus=8):
    """
    Perform location overlap analysis (LOLA).

    If bed_files is a list with more than one element, use ``output_prefixes`` to pass a list of
    prefixes to label the output files for each input BED file.

    Files will be created in ``output_folder`` mimicking the output that the R function
    LOLA::writeCombinedEnrichment writes.

    Requires the R package "LOLA" to be installed:
        >>> source('http://bioconductor.org/biocLite.R')
        >>> biocLite('LOLA')

    Attributes:

    :param str,list bed_files:
        A string path to a BED file or a list of paths.

    :param str universe_file:
        A path to a BED file representing the universe from where the BED file(s) come from.

    :param str output_folder:
        Output folder for resulting files.

    :param list,optional output_prefixes:
        A list of strings with prefixes to be used in case ``bed_files`` is a list.
        Defaults to None

    :param str,optional genome:
        Genome assembly from which the BED files come from.
        This is used to get the LOLA databases from the ngs_toolkit._CONFIG parameters.
        Defaults to "hg19".

    :param int,optional cpus:
        Number of CPUs/threads to use.
        Defaults to 8
    """
    numpy2ri.activate()
    pandas2ri.activate()
    warnings.filterwarnings("ignore", category=RRuntimeWarning)

    _LOGGER.info("Loading LOLA R library thorugh rpy2.")
    robjects.r('require("LOLA")')
    _loadRegionDB = robjects.r('LOLA::loadRegionDB')
    _readBed = robjects.r('LOLA::readBed')
    _runLOLA = robjects.r('LOLA::runLOLA')
    # _writeCombinedEnrichment = robjects.r('LOLA::writeCombinedEnrichment')

    # Get region databases from config
    _LOGGER.info("Getting LOLA databases for genome '{}' from configuration.".format(genome))

    msg = "LOLA database values in configuration could not be found or understood. "
    msg += "Please add a list of value(s) to this section 'resources:lola:region_databases:{}'. ".format(genome)
    msg += "For an example, see https://github.com/afrendeiro/toolkit/tree/master/ngs_toolkit/config/example.yaml"
    try:
        databases = _CONFIG['resources']['lola']['region_databases'][genome]
    except KeyError:
        _LOGGER.error(msg)
        return

    if not isinstance(databases, list):
        if isinstance(databases, str):
            databases = list(databases)
        else:
            _LOGGER.error(msg)
            return

    if len(databases) < 1:
        _LOGGER.error(msg)
        return

    if isinstance(bed_files, str):
        bed_files = [bed_files]
    if output_prefixes is None:
        if len(bed_files) > 1:
            msg = "Running more than one BED file at once while only specifying `output_folder` argument"
            msg += " will cause output files to be named in the form '{output_folder}/{region_database}.{input_file}.tsv'."
            msg += " To prevent this behaviour, pass a list of arguments to `output_prefixes`."
            _LOGGER.warning(msg)
            output_prefixes = [r.replace(os.path.sep, "__").replace(".bed", ".") for r in bed_files]
        else:
            output_prefixes = ["."]

    _LOGGER.info("Reading up universe file '{}'.".format(universe_file))
    universe = _readBed(universe_file)
    _LOGGER.info("Loading region set databases.")
    regionDB = _loadRegionDB(np.array(databases))
    for suffix, bed_file in zip(output_prefixes, bed_files):
        _LOGGER.info("Reading up BED file '{}'.".format(bed_file))
        user_set = _readBed(bed_file)
        _LOGGER.info("Running LOLA testing for file '{}'.".format(bed_file))
        _lola_results = _runLOLA(user_set, universe, regionDB, cores=cpus)
        _LOGGER.info("Converting results from R to Python")
        lola_results = r2pandas_df(_lola_results)
        _LOGGER.info("Saving all results for file '{}'.".format(bed_file))
        lola_results.to_csv(
            os.path.join(output_folder, "allEnrichments" + suffix + "tsv"), index=False, sep="\t")
        for region_set in lola_results['collection'].drop_duplicates():
            _LOGGER.info("Saving results for collection '{}' only.".format(region_set))
            lola_results[lola_results['collection'] == region_set].to_csv(
                os.path.join(output_folder, "col_" + region_set + suffix + "tsv"), index=False, sep="\t")


def bed_to_fasta(bed_file, fasta_file, genome="hg19", genome_2bit=None):
    if genome_2bit is None:
        # Get region databases from config
        _LOGGER.info("Getting 2bit reference genome for genome '{}' from configuration.".format(genome))

        msg = "Reference genome in 2bit format value in configuration could not be found or understood. "
        msg += "Please add a list of value(s) to this section 'resources:genomes:2bit:{}'. ".format(genome)
        msg += "For an example, see https://github.com/afrendeiro/toolkit/tree/master/ngs_toolkit/config/example.yaml"
        try:
            genome_2bit = _CONFIG['resources']['genomes']['2bit'][genome]
        except KeyError:
            _LOGGER.error(msg)
            return

        if not isinstance(genome_2bit, str):
            _LOGGER.error(msg)
            return

    if not os.path.exists(genome_2bit):
        _LOGGER.error("Reference genome in 2bit does not exist or can't be open: '{}'".format(genome_2bit))

    # write name column
    bed = pd.read_csv(bed_file, sep='\t', header=None)
    bed['name'] = bed[0] + ":" + bed[1].astype(str) + "-" + bed[2].astype(str)
    bed[1] = bed[1].astype(int)
    bed[2] = bed[2].astype(int)
    bed.to_csv(bed_file + ".tmp.bed", sep='\t', header=None, index=False)

    # do enrichment
    cmd = "twoBitToFa {0} -bed={1} {2}".format(genome_2bit, bed_file + ".tmp.bed", fasta_file)

    subprocess.call(cmd.split(" "))
    # subprocess.call("rm %s" % bed_file + ".tmp.bed")


def meme_ame(
        input_fasta,
        output_dir,
        background_fasta=None,
        organism="human",
        motif_database_file=None):
    if motif_database_file is None:
        # Get region databases from config
        _LOGGER.info("Getting 2bit reference genome for genome '{}' from configuration.".format(organism))

        msg = "Reference genome in 2bit format value in configuration could not be found or understood. "
        msg += "Please add a list of value(s) to this section 'resources:meme:motif_databases:{}'. ".format(organism)
        msg += "For an example, see https://github.com/afrendeiro/toolkit/tree/master/ngs_toolkit/config/example.yaml"
        try:
            motif_database_file = _CONFIG['resources']['meme']['motif_databases'][organism]
        except KeyError:
            _LOGGER.error(msg)
            return

        if not isinstance(motif_database_file, str):
            _LOGGER.error(msg)
            return

    # shuffle input in no background is provided
    if background_fasta is None:
        shuffled = input_fasta + ".shuffled"
        cmd = """
        fasta-dinucleotide-shuffle -c 1 -f {0} > {1}
        """.format(input_fasta, shuffled)
        subprocess.call(cmd.split(" "))

    cmd = """
    ame --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 0.05 \\
    --control {0} -o {1} {2} {3}
    """.format(
        background_fasta if background_fasta is not None else shuffled,
        output_dir, input_fasta, motif_database_file)
    subprocess.call(cmd.split(" "))
    # subprocess.call("rm {}".format(shuffled).split(" "))


def homer_motifs(bed_file, output_dir, genome="hg19"):
    cmd = "findMotifsGenome.pl {bed} {genome}r {out_dir} \
    -size 1000 -h -p 2 -len 8,10,12,14 -noknown".format(
        bed=bed_file, genome=genome, out_dir=output_dir)
    subprocess.call(cmd.split(" "))


def homer_combine_motifs(
        comparison_dirs, output_dir,
        region_prefix="differential_analysis",
        reduce_threshold=0.6, match_threshold=10, info_value=0.6,
        p_value_threshold=1e-25, fold_enrichment=None,
        cpus=8, run=True, as_jobs=True, genome="hg19",
        motif_database=None, known_vertebrates_TFs_only=False):
    """
    Create consensus of de novo discovered motifs from HOMER

    Attributes:

    :param list comparison_dirs:
        Iterable of comparison directories where homer was run.
        Should contain a "homerMotifs.all.motifs" file.

    :param str output_dir:
        Output directory.

    :param number,optional p_value_threshold:
        Threshold for inclusion of a motif in the consensus set.
        Defaults to 1e-5

    :param number,optional cpus:
        Number of available CPUS/threads for multithread processing.
        Defaults to 8

    :param bool,optional run:
        Whether to run enrichment of each comparison in the consensus motifs.
        Default is True

    :param bool,optional as_jobs:
        Whether to run enrichment as a cluster job.
        Default is True

    :param str genome:
        Genome assembly of the data.
        Default is 'hg19'.

    :param bool known_vertebrates_TFs_only:
        Deprecated. Pass a given motif_database to `motif_database` directly.

    :param str motif_database:
        Motif database to restrict motif matching too.

    :returns str:
        If `run` is `False`, returns path to consensus motif file. Otherwise `None`.
    """
    if known_vertebrates_TFs_only:
        _LOGGER.warning("WARNING! `known_vertebrates_TFs_only` option is deprecated!" +
                        "Pass a given motif_database to `motif_database` directly.")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # concatenate files
    out_file = os.path.join(output_dir, "homerMotifs.combined.motifs")
    with open(out_file, 'a') as outfile:
        for dir_ in comparison_dirs:
            with open(os.path.join(dir_, "homerMotifs.all.motifs"), "r") as infile:
                for line in infile:
                    outfile.write(line)

    # Filter and get motif consensus
    extra = ""
    if motif_database:
        extra = " -known {}".format(motif_database)
    if fold_enrichment is None:
        fold_enrichment = ""
    else:
        fold_enrichment = " -F " + str(fold_enrichment)
    subprocess.call("compareMotifs.pl {} {} -reduceThresh {} -matchThresh {}{} -pvalue {} -info {}{} -nofacts -cpu {}"
                    .format(out_file, output_dir, reduce_threshold, match_threshold,
                            extra, p_value_threshold, info_value, fold_enrichment, cpus).split(" "))

    # concatenate consensus motif files
    files = glob.glob(os.path.join(output_dir, "homerResults/*motif"))
    combined_motifs = os.path.join(output_dir, "homerMotifs.filtered.motifs")
    with open(combined_motifs, 'w') as outfile:
        for f in files:
            if ("similar" in f) or ("RV" in f):
                continue
            _LOGGER.info(f)
            with open(f, "r") as infile:
                for line in infile:
                    outfile.write(line)

    if run:
        for dir_ in comparison_dirs:
            # delete previous results if existing
            results = os.path.join(dir_, "knownResults.txt")
            if os.path.exists(results):
                _LOGGER.warning("Deleting previously existing results file: '{}'".format(results))
                os.remove(results)
            # prepare enrichment command with consensus set
            cmd = (
                "findMotifsGenome.pl {bed} {genome}r {dir} -p {cpus} -nomotif -mknown {motif_file}"
                .format(bed=os.path.join(dir_, region_prefix + "_regions.bed"),
                        genome=genome, cpus=cpus, dir=dir_,
                        motif_file=combined_motifs))
            # run
            if as_jobs:
                subprocess.call("sbatch -J homer.{d} -o {dir}.homer.log -p shortq -c 8 --mem 20000"
                                .format(d=os.path.basename(dir_), dir=dir_)
                                .split(" ")
                                + ['--wrap', cmd])
            else:
                subprocess.call(cmd.split(" "))


def enrichr(dataframe, gene_set_libraries=None, kind="genes"):
    """
    Use Enrichr on a list of genes (currently only genes supported through the API).

    :param str dataframe:
        DataFrame with column "gene_name".

    :param list,optional gene_set_libraries:
        Gene set libraries to use.
        Defaults to values in initial configuration file.
        To see them, do: ``ngs_toolkit._CONFIG['resources']['enrichr']['gene_set_libraries']``

    :param str,optional kind:
        Kind of input. Right now, only "genes" is supported.
        Defaults to "genes"

    :returns pandas.DataFrame:
        Results of enrichment analysis

    :raises: Exception
    """
    ENRICHR_ADD = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    ENRICHR_RETRIEVE = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'

    if gene_set_libraries is None:
        # Get region databases from config
        _LOGGER.info("Getting Enrichr gene set libraries from configuration.")

        msg = "Enrichr gene set libraries value in configuration could not be found or understood. "
        msg += "Please add a list of value(s) to this section 'resources:mem:motif_databases'. "
        msg += "For an example, see https://github.com/afrendeiro/toolkit/tree/master/ngs_toolkit/config/example.yaml"
        try:
            gene_set_libraries = _CONFIG['resources']['enrichr']['gene_set_libraries']
        except KeyError:
            _LOGGER.error(msg)
            return

        if not isinstance(gene_set_libraries, list):
            if isinstance(gene_set_libraries, str):
                gene_set_libraries = list(gene_set_libraries)
            else:
                _LOGGER.error(msg)
                return

        if len(gene_set_libraries) < 1:
            _LOGGER.error(msg)
            return

    results = pd.DataFrame()
    for gene_set_library in tqdm(gene_set_libraries, total=len(gene_set_libraries), desc="Gene set library"):
        _LOGGER.info("Using enricher on %s gene set library." % gene_set_library)

        if kind == "genes":
            # Build payload with bed file
            attr = "\n".join(dataframe["gene_name"].dropna().tolist())
        elif kind == "regions":
            # Build payload with bed file
            attr = "\n".join(dataframe[['chrom', 'start', 'end']].apply(
                lambda x: "\t".join([str(i) for i in x]), axis=1).tolist())

        payload = {
            'list': (None, attr),
            'description': (None, gene_set_library)
        }
        # Request adding gene set
        response = requests.post(ENRICHR_ADD, files=payload)
        if not response.ok:
            raise Exception('Error analyzing gene list')

        # Track gene set ID
        user_list_id = json.loads(response.text)['userListId']

        # Request enriched sets in gene set
        response = requests.get(
            ENRICHR_RETRIEVE + query_string % (user_list_id, gene_set_library)
        )
        if not response.ok:
            raise Exception('Error fetching enrichment results')

        # Get enriched sets in gene set
        res = json.loads(response.text)
        # If there's no enrichemnt, continue
        if len(res) < 0:
            continue

        # Put in dataframe
        res = pd.DataFrame([pd.Series(s) for s in res[gene_set_library]])
        if res.shape[0] == 0:
            continue
        cols = ["rank", "description", "p_value", "z_score",
                "combined_score", "genes", "adjusted_p_value"]
        if len(res.columns) == 7:
            res.columns = cols
        elif len(res.columns) == 9:
            res.columns = cols + ["old_p_value", "old_adjusted_p_value"]

        # Remember gene set library used
        res["gene_set_library"] = gene_set_library

        # Append to master dataframe
        results = results.append(res, ignore_index=True)

    return results


def run_enrichment_jobs(
        analysis_name, results_dir, genome,
        background_bed="results/{PROJECT_NAME}_peak_set.bed",
        steps=['lola', 'meme', 'homer', 'enrichr']):
    """
    Submit enrichment jobs for a specifc analysis.
    """
    cmds = list()

    # LOLA
    if 'lola' in steps:
        cmds += ["""for F in `find {results_dir} -name "*_regions.bed"`; do
DIR=`dirname $F`
if [ ! -f ${{DIR}}/allEnrichments.tsv ]; then
echo $DIR $F
sbatch -J lola.$F -o $F.lola.log -p shortq -c 8 --mem 24000 \
--wrap "Rscript ~/jobs/run_LOLA.R $F {background_bed} {GENOME}"
fi
done""".format(results_dir=results_dir, background_bed=background_bed, GENOME=genome)]

    # AME
    if 'meme' in steps:
        dbs = {
            "human": "~/resources/motifs/motif_databases/HUMAN/HOCOMOCOv10.meme",
            "mouse": "~/resources/motifs/motif_databases/MOUSE/uniprobe_mouse.meme"}
        omap = {"hg38": "human", "hg19": "human", "mm10": "mouse"}

        cmds += ["""for F in `find {results_dir} -name "*_regions.fa"`; do
DIR=`dirname $F`
if [ ! -f ${{DIR}}/ame.html ]; then
echo $DIR $F
sbatch -J "meme_ame.${{F}}" -o "${{F}}.meme_ame.log" -p shortq -c 1 --mem 4000 \
--wrap "fasta-dinucleotide-shuffle -c 1 -f "${{F}}" > "${{F}}".shuffled.fa; \
ame --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 0.05 \
--control "${{F}}".shuffled.fa -o "${{DIR}}" "${{F}}" {motifs}"
fi
done""".format(results_dir=results_dir, motifs=dbs[omap[genome]])]

    # HOMER
    if 'homer' in steps:
        cmds += ["""for F in `find {results_dir} -name "*_regions.bed"`; do
DIR=`dirname $F`
if [ ! -f ${{DIR}}/homerResults.html ]; then
echo $DIR $F
sbatch -J "homer.${{F}}" -o "${{F}}.homer.log" -p shortq -c 8 --mem 20000 \
--wrap "findMotifsGenome.pl ${{F}} {GENOME}r ${{DIR}} -size 1000 -h -p 2 -len 8,10,12,14 -noknown"
fi
done""".format(results_dir=results_dir, GENOME=genome)]

    # Enrichr
    if 'enrichr' in steps:
        cmds += ["""for F in `find {results_dir} -name "*.gene_symbols.txt"`; do
if [ ! -f ${{F/gene_symbols.txt/enrichr.csv}} ]; then
echo $F
sbatch -J enrichr.$F -o $F.enrichr.log -p shortq -c 1 --mem 4000 \
--wrap "python -u ~/jobs/run_Enrichr.py --input-file "$F" --output-file "${{F/gene_symbols.txt/enrichr.csv}}" "
fi
done""".format(results_dir=results_dir)]
        cmds += ["""for F in `find {results_dir} -name "*_genes.symbols.txt"`; do
if [ ! -f ${{F/symbols.txt/enrichr.csv}} ]; then
echo $F
sbatch -J enrichr.$F -o $F.enrichr.log -p shortq -c 1 --mem 4000 \
--wrap "python -u ~/jobs/run_Enrichr.py --input-file "$F" --output-file "${{F/symbols.txt/enrichr.csv}}" "
fi
done""".format(results_dir=results_dir)]

    for cmd in cmds:
        subprocess.call(cmd.split(" "))


def project_to_geo(
        project, output_dir="geo_submission",
        samples=None, distributed=False, dry_run=False):
    """
    Prepare raw sequencing files for submission to GEO.
    Files will be copied or generated in a new directory `output_dir`.
    It will get the raw BAM file(s) of each sample, and in case of
    ATAC-seq/ChIP-seq samples, the bigWig and peak files. If multiple BAM
    files exist for each sample, all will be copied and sequencially named
    with the "fileN" suffix, where "N" is the file number.

    For each copied file a md5sum will be calculated.

    A pandas DataFrame with info on the sample's files and md5sums will be returned.

    :param peppy.Project project:
        A peppy Project object to process.

    :param str,optional output_dir:
        Directory to create output. Will be created/overwriten if existing.
        Defaults to "geo_submission".

    :param list,optional samples:
        List of peppy.Sample objects in project to restrict to.
        Defaults to all samples in project.

    :param bool,optional distributed:
        Whether processing should be distributed as jobs in a computing cluster for each sample.
        Currently available implementation supports a SLURM cluster only.
        Defaults to False.

    :param bool,optional dry_run:
        Whether copy/execution/submisison commands should be not be run to test.
        Default is False.

    :returns pandas.DataFrame:
        Annotation of samples and their BAM, BigWig, narrowPeak files and respective md5sums.
    """
    output_dir = os.path.abspath(output_dir)
    if samples is None:
        samples = project.samples
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    cmds = list()
    annot = pd.DataFrame(index=pd.Index([], name="sample_name"))
    for sample in samples:
        various = len(sample.data_path.split(" ")) > 1
        cmd = ""
        for i, file in enumerate(sample.data_path.split(" ")):
            suffix = ".file{}".format(i) if various else ""
            # Copy raw file
            bam_file = os.path.join(output_dir, sample.name + "{}.bam".format(suffix))
            cmd += "cp {} {}; ".format(file, bam_file)
            cmd += "chmod 644 {}; ".format(bam_file)
            annot.loc[sample.name, "bam_file{}".format(i)] = bam_file

            # Copy or generate md5sum
            md5_file = bam_file + ".md5"
            if os.path.exists(file + ".md5"):
                cmd += "cp {} {}; ".format(file + ".md5", md5_file)
            else:
                b = os.path.basename(file)
                cmd += "md5sum {} > {}; ".format(os.path.join(output_dir, b), md5_file)
            cmd += "chmod 644 {}; ".format(md5_file)
            annot.loc[sample.name, "bam_file{}_md5sum".format(i)] = md5_file

        # Copy bigWig files
        if sample.library in ["ATAC-seq", "ChIP-seq"]:
            if hasattr(sample, "bigwig"):
                bigwig_file = os.path.join(output_dir, sample.name + ".bigWig")
                cmd += "cp {} {}; ".format(sample.bigwig, bigwig_file)
                cmd += "chmod 644 {}; ".format(bigwig_file)
                annot.loc[sample.name, "bigwig_file"] = bigwig_file

                # Copy or generate md5sum
                md5_file = bigwig_file + ".md5"
                if os.path.exists(sample.bigwig + ".md5"):
                    cmd += "cp {} {}; ".format(sample.bigwig + ".md5", md5_file)
                else:
                    b = os.path.basename(sample.bigwig)
                    cmd += "md5sum {} > {}; ".format(os.path.join(output_dir, b), md5_file)
                cmd += "chmod 644 {}; ".format(md5_file)
                annot.loc[sample.name, "bigwig_file_md5sum"] = md5_file
            else:
                _LOGGER.warning("'{}' sample '{}' does not have a 'bigwig'"
                                .format(sample.library, sample.name) +
                                " attribute set. Skipping bigWig file.")
        # Copy peaks
        if sample.library == "ATAC-seq":
            if hasattr(sample, "peaks"):
                peaks_file = os.path.join(output_dir, sample.name + ".peaks.narrowPeak")
                cmd += "cp {} {}; ".format(sample.peaks, peaks_file)
                cmd += "chmod 644 {}; ".format(peaks_file)
                annot.loc[sample.name, "peaks_file"] = peaks_file

                # Copy or generate md5sum
                md5_file = peaks_file + ".md5"
                if os.path.exists(sample.peaks + ".md5"):
                    cmd += "cp {} {}; ".format(sample.peaks + ".md5", md5_file)
                else:
                    b = os.path.basename(sample.peaks)
                    cmd += "md5sum {} > {}; ".format(peaks_file, md5_file)
                cmd += "chmod 644 {}; ".format(md5_file)
                annot.loc[sample.name, "peaks_file_md5sum"] = md5_file
            else:
                _LOGGER.warning("'{}' sample '{}' does not have a 'peaks' attribute set."
                                .format(sample.library, sample.name) +
                                " Skipping peaks file.")
        if distributed:
            tk = NGSTk()

            job_name = "project_to_geo.{}".format(sample.name)
            log_file = os.path.join(output_dir, job_name + ".log")
            job_file = os.path.join(output_dir, job_name + ".sh")

            job = textwrap.dedent(tk.slurm_header(
                job_name=job_name, output=log_file,
                cpus_per_task=1, mem_per_cpu=8000))
            job += "\n" + "\n".join(cmd.split("; ")) + "\n"
            job += textwrap.dedent(tk.slurm_footer())

            with open(job_file, "w") as handle:
                handle.write(textwrap.dedent(job))
            if not dry_run:
                tk.slurm_submit_job(job_file)
        else:
            cmds.append(cmd)

    if not distributed and not dry_run:
        for i, cmd in enumerate(cmds):
            _LOGGER.info(i, cmd)
            subprocess.call(cmd.split(" "))

    return annot


def rename_sample_files(
        annotation_mapping,
        old_sample_name_column="old_sample_name",
        new_sample_name_column="new_sample_name",
        tmp_prefix="rename_sample_files",
        results_dir="results_pipeline",
        dry_run=False):
    """
    Rename existing directories with pipeline outputs for samples based on mapping of
    old/new sample names.
    All files within the directory with the old sample name will be renamed recursively.
    Old and new sample names can overlap - this procedure will handle these cases correctly
    by a 2-step process with temporary sample names with prefix `prefix`.

    NEEDS TESTING!

    :param pandas.DataFrame annotation_mapping:
        DataFrame with mapping of old (column "previous_sample_name") vs new ("new_sample_name") sample names.

    :param str,optional old_sample_name_column:
        Name of column with old sample names.
        Defaults to "old_sample_name"

    :param str,optional new_sample_name_column:
        Name of column with new sample names.
        Defaults to "new_sample_name"

    :param str,optional tmp_prefix:
        Prefix for temporary files to avoid overlap between old and new names.
        Defaults to "rename_sample_files"

    :param str,optional results_dir:
        Pipeline output directory containing sample output directories.
        Defaults to "results_pipeline"

    :param bool,optional dry_run:
        Whether to print commands instead of running them. Defaults to False

    :returns: None
    """
    cmds = list()
    # 1) move to tmp name
    for i, series in annotation_mapping.iterrows():
        o = series[old_sample_name_column]
        t = "{}-{}".format(tmp_prefix, i)
        cmds.append("# Moving old sample '{}' to '{}'.".format(
            o, t))

        # directory
        cmds.append("mv {} {}".format(
            os.path.join(results_dir, o), os.path.join(results_dir, t)))
        # further directories
        cmds.append("find {} -type d -exec rename {} {} {{}} \\;".format(
            os.path.join(results_dir, t), o, t))
        # files
        cmds.append("find {} -type f -exec rename {} {} {{}} \\;".format(
            os.path.join(results_dir, t), o, t))

    # 2) move to new name
    for i, series in annotation_mapping.iterrows():
        t = "{}-{}".format(tmp_prefix, i)
        n = series[new_sample_name_column]
        cmds.append("# Moving old sample '{}' to '{}'.".format(
            t, n))

        # directory
        cmds.append("mv {} {}".format(
            os.path.join(results_dir, t), os.path.join(results_dir, n)))
        # further directories
        cmds.append("find {} -type d -exec rename {} {} {{}} \\;".format(
            os.path.join(results_dir, n), t, n))
        # files
        cmds.append("find {} -type f -exec rename {} {} {{}} \\;".format(
            os.path.join(results_dir, n), t, n))

    if dry_run:
        _LOGGER.info("\n".join(cmds))
    else:
        for i, cmd in enumerate(cmds):
            _LOGGER.info(i, cmd)
            if cmd.startswith("#"):
                continue
            try:
                r = subprocess.call(cmd.split(" "))
            except OSError as e:
                raise e
            if r != 0:
                raise OSError("Command '{}' failed.".format(cmd))


def query_biomart(
        attributes=None,
        species="hsapiens", ensembl_version="grch37"):
    """
    Query Biomart (https://www.ensembl.org/biomart/martview/).

    Query Biomart for gene attributes. Returns pandas dataframe with query results.
    If a certain field contains commas, it will attemp to return dataframe but it might fail.

    :param list,optional attributes:
        List of ensembl atrributes to query.
        Defaults to ["ensembl_gene_id", "external_gene_name", "hgnc_id", "hgnc_symbol"].

    :param str,optional species:
        Ensembl string of species to query. Must be vertebrate.
        Defaults to "hsapiens".

    :param str,optional ensembl_version:
        Ensembl version to query. Currently "grch37", "grch38" and "grcm38" are tested.
        Defaults to "grch37".

    :returns pandas.DataFrame:
        Dataframe with required attributes for each entry.
    """
    supported = ['grch37', 'grch38', 'grcm38']
    if ensembl_version not in supported:
        msg = "Ensembl version might not be supported."
        msg += " Tested versions are '{}'.".format("','".join(supported))
        msg += " Will try anyway."
        _LOGGER.warning(msg)

    if attributes is None:
        attributes = ["ensembl_gene_id", "external_gene_name", "hgnc_id", "hgnc_symbol"]

    # Build request XML
    ens_ver = "" if ensembl_version.endswith("38") else ensembl_version + "."
    url_query = "".join([
        """http://{}ensembl.org/biomart/martservice?query=""".format(ens_ver),
        """<?xml version="1.0" encoding="UTF-8"?>""",
        """<!DOCTYPE Query>""",
        """<Query  virtualSchemaName="default" formatter="CSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6" >""",
        """<Dataset name="{}_gene_ensembl" interface="default" >""".format(species)] +
        ["""<Attribute name="{}" />""".format(attr) for attr in attributes] +
        ["""</Dataset>""",
         """</Query>"""])
    req = requests.get(url_query, stream=True)
    if not req.ok:
        msg = "Request to Biomart API was not successful."
        _LOGGER.error(msg)
        raise ValueError(msg)
    content = list(req.iter_lines())

    if (len(content) == 1) and (content[0].startswith("Query ERROR")):
        msg = "Request to Biomart API was not successful. Check your input.\n{}".format(content[0])
        _LOGGER.error(msg)
        raise ValueError(msg)

    if isinstance(content[0], bytes):
        content = [x.decode("utf-8") for x in content]

    try:
        mapping = pd.DataFrame([x.strip().split(",") for x in content], columns=attributes)
    except AssertionError as e:
        msg = "Could not simply return dataframe with results."
        msg += " It is likely this is because one of the requested attributes has commas as is quoted."
        msg += " Or because querying an organism not present in vertebrate database."
        msg += " Will try to replace these fields and parse."
        _LOGGER.warning(msg)
        # input probably contains commas inside fields
        c = pd.Series([x.strip() for x in content])
        # well's assume that fields with commas have been quoted
        # get text inside double quotes
        cc = c.str.extract("\"(.*)\"")
        if cc.shape[1] > 1:
            _LOGGER.error("Attempt to fix failed.")
            raise e
        cc = cc.squeeze().str.replace(",", "")
        # now get text until first quote and concatenate with clean text
        f = c.str.extract("(.*),\"").fillna("").squeeze() + "," + cc
        # now get back together with instances that didn't have quotes
        c.update(f)
        try:
            mapping = pd.DataFrame([x.strip().split(",") for x in c], columns=attributes)
        except AssertionError:
            _LOGGER.error("Attempt to fix failed.")
            raise e
        _LOGGER.info("Attempt to fix successful.")

    return mapping.replace("", np.nan)


def subtract_principal_component(
        X, pc=1, norm=False, plot=True, plot_name="PCA_based_batch_correction.svg", pcs_to_plot=6):
    """
    Given a matrix (n_samples, n_variables), remove `pc` (1-based) from matrix.
    """
    pc -= 1

    # All regions
    if norm:
        X = StandardScaler().fit_transform(X)

    # PCA
    pca = PCA()
    X_hat = pca.fit_transform(X)

    # Remove PC
    X2 = X - np.outer(X_hat[:, pc], pca.components_[pc, :])

    # plot
    if plot:
        X2_hat = pca.fit_transform(X2)
        fig, axis = plt.subplots(pcs_to_plot, 2, figsize=(4 * 2, 4 * pcs_to_plot))
        for pc in range(pcs_to_plot):
            # before
            for j, sample in enumerate(X.index):
                axis[pc, 0].scatter(X_hat[j, pc], X_hat[j, pc + 1], s=50, rasterized=True)
            axis[pc, 0].set_xlabel("PC{}".format(pc + 1))
            axis[pc, 0].set_ylabel("PC{}".format(pc + 2))
            # after
            for j, sample in enumerate(X2.index):
                axis[pc, 1].scatter(X2_hat[j, pc], X2_hat[j, pc + 1], s=35, alpha=0.8, rasterized=True)
            axis[pc, 1].set_xlabel("PC{}".format(pc + 1))
            axis[pc, 1].set_ylabel("PC{}".format(pc + 2))
        fig.savefig(plot_name)

    return X2


def subtract_principal_component_by_attribute(df, attributes, pc=1):
    """
    Given a matrix (n_samples, n_variables), remove `pc` (1-based) from matrix.
    """
    pc -= 1

    X2 = pd.DataFrame(index=df.index, columns=df.columns)
    for attr in attributes:
        _LOGGER.info(attr)
        sel = df.index[df.index.str.contains(attr)]
        X = df.loc[sel, :]

        # PCA
        pca = PCA()
        X_hat = pca.fit_transform(X)

        # Remove PC
        X2.loc[sel, :] = X - np.outer(X_hat[:, pc], pca.components_[pc, :])
    for sample in df.index:
        if X2.loc[sample, :].isnull().all():
            X2.loc[sample, :] = df.loc[sample, :]
    return X2


def fix_batch_effect_limma(matrix, batch_variable="batch", covariates=None):
    """
    Fix batch effect in matrix using limma.

    Requires the R package "limma" to be installed:
        >>> source('http://bioconductor.org/biocLite.R')
        >>> biocLite('limma')

    :param matrix: DataFrame with multiindex for potential covariate annotations
    :type matrix: [type]
    :param formula: [description], defaults to "~knockout"
    :type formula: str,optional
    :returns: [description]
    :rtype: {[type]}
    """
    warnings.filterwarnings("ignore", category=RRuntimeWarning)
    numpy2ri.activate()
    pandas2ri.activate()

    robjects.r('require("limma")')
    _removeBatchEffect = robjects.r('removeBatchEffect')

    if covariates is None:
        covariates = []

    if len(covariates) > 0:
        cov = patsy.dmatrix("~{} - 1".format(" + ".join(covariates)), matrix.columns.to_frame())
        fixed = _removeBatchEffect(
            x=matrix.values,
            batch=matrix.columns.get_level_values(batch_variable),
            design=cov)
    else:
        fixed = _removeBatchEffect(
            x=matrix.values,
            batch=matrix.columns.get_level_values(batch_variable))
    fixed = pd.DataFrame(
        np.asarray(fixed),
        index=matrix.index,
        columns=matrix.columns)
    return fixed
    # fixed.to_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.limma_fixed.csv"))
