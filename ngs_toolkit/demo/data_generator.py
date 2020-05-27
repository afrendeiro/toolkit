#!/usr/bin/env python

"""
A module dedicated to the generation of Analysis, Projects and their data.
"""

import os
import string
import tempfile

import numpy as np
import pandas as pd
import pybedtools
import yaml

from ngs_toolkit.general import query_biomart
from ngs_toolkit.utils import location_index_to_bed
from ngs_toolkit.project_manager import create_project as init_proj


REGION_BASED_DATA_TYPES = [
    "ATAC-seq",
    "ChIP-seq",
]  # CNV technically is too but it does not need a `sites` attr
DEFAULT_CNV_RESOLUTIONS = ["1000kb", "100kb", "10kb"]


def generate_count_matrix(
    n_factors=1,
    n_replicates=4,
    n_features=1000,
    intercept_mean=4,
    intercept_std=2,
    coefficient_stds=0.4,
    size_factors=None,
    size_factors_std=0.1,
    dispersion_function=None,
):
    """
    Generate count matrix for groups of samples by sampling from a
    negative binomial distribution.
    """
    import patsy

    if isinstance(coefficient_stds, (int, float)):
        coefficient_stds = [coefficient_stds] * n_factors

    if dispersion_function is None:
        dispersion_function = _disp

    # Build sample vs factors table
    dcat = pd.DataFrame(patsy.demo_data(*(list(string.ascii_lowercase[:n_factors]))))
    dcat.columns = dcat.columns.str.upper()
    for col in dcat.columns:
        dcat[col] = dcat[col].str.upper()
    if n_replicates > 1:
        dcat = (
            pd.concat([dcat for _ in range(int(np.ceil(n_replicates / 2)))])
            .sort_values(dcat.columns.tolist())
            .reset_index(drop=True)
        )
    dcat.index = ["S{}_{}".format(str(i + 1).zfill(2), dcat.loc[i, :].sum()) for i in dcat.index]
    m_samples = dcat.shape[0]

    # make model design table
    design = np.asarray(
        patsy.dmatrix("~ 1 + " + " + ".join(string.ascii_uppercase[:n_factors]), dcat)
    )

    # get means
    beta = np.asarray(
        [np.random.normal(intercept_mean, intercept_std, n_features)]
        + [np.random.normal(0, std, n_features) for std in coefficient_stds]
    ).T

    if size_factors is None:
        size_factors = np.random.normal(1, size_factors_std, (m_samples, 1))

    mean = (2 ** (design @ beta.T) * size_factors).T

    # now sample counts
    dispersion = (1 / dispersion_function(2 ** (beta[:, 1:]))).mean(1).reshape(-1, 1)
    dnum = pd.DataFrame(
        np.random.negative_binomial(n=mean, p=dispersion, size=mean.shape), columns=dcat.index
    )
    dcat.index.name = dnum.columns.name = "sample_name"
    return dnum, dcat


def generate_data(
    n_factors=1,
    n_replicates=4,
    n_features=1000,
    coefficient_stds=0.4,
    data_type="ATAC-seq",
    genome_assembly="hg38",
    **kwargs
):
    """
    Creates real-looking data dependent on the data type.

    Parameters
    ----------
    n_factors : :obj:`int`, optional
        Number of factors influencing variance between groups.
        For each factor there will be two groups of samples.

        Defaults to 1.
    n_replicates : :obj:`int`, optional
        Number of replicates per group.

        Defaults to 4.
    n_features : :obj:`int`, optional
        Number of features (i.e. genes, regions) in matrix.

        Defaults to 1000.
    coefficient_stds : {:obj:`int`, :obj:`list`}, optional
        Standard deviation of the coefficients between groups.
        If a list, must match the number of ``n_factors``.

        Defaults to 1.
    data_type : :obj:`bool`, optional
        Data type of the project. Must be one of the ``ngs_toolkit`` classes.

        Default is "ATAC-seq"
    genome_assembly : :obj:`bool`, optional
        Genome assembly of the project.

        Default is "hg38"
    **kwargs : :obj:`dict`
        Additional keyword arguments will be passed to
        :func:`ngs_toolkit.demo.data_generator.generate_count_matrix`.

    Returns
    -------
    :obj:`tuple`
        A tuple of :class:`pandas.DataFrame` objects with numeric
        and categorical data respectively.
    """
    dnum, dcat = generate_count_matrix(
        n_factors=n_factors, n_features=n_features, coefficient_stds=coefficient_stds, **kwargs
    )

    # add random location indexes
    if data_type in ["ATAC-seq", "ChIP-seq"]:
        dnum.index = get_random_genomic_locations(n_features, genome_assembly=genome_assembly)
    elif data_type in ["RNA-seq"]:
        dnum.index = get_random_genes(n_features, genome_assembly=genome_assembly)
    elif data_type in ["CNV"]:
        # choose last value of last factor as "background"
        factor = dcat.columns[-1]
        value = dcat[factor].unique()[-1]
        background_samples = dcat[factor] == value
        # get log2 ratios
        dnum += 1
        background = dnum.loc[:, background_samples].median(1)
        signal = dnum.loc[:, ~background_samples].median(1)
        dnum = np.log2(dnum.T / background).T

        # now sort feature by group difference
        dnum = dnum.reindex((signal - background).sort_values().index)

        # and assign contiguous locations in the genome
        dnum_res = dict()
        for resolution in DEFAULT_CNV_RESOLUTIONS:
            x = dnum.copy()
            x.index = get_genomic_bins(
                n_features, genome_assembly=genome_assembly, resolution=resolution
            )
            dnum_res[resolution] = x
        dnum = dnum_res
    elif data_type in ["CRISPR"]:
        dnum.index = get_random_grnas(n_features, genome_assembly=genome_assembly)
    return dnum, dcat


def generate_project(
    output_dir=None,
    project_name="test_project",
    organism="human",
    genome_assembly="hg38",
    data_type="ATAC-seq",
    n_factors=1,
    only_metadata=False,
    sample_input_files=False,
    initialize=True,
    **kwargs
):
    """
    Creates a real-looking PEP-based project with respective input files
    and quantification matrix.

    Parameters
    ----------
    output_dir : :obj:`str`, optional
        Directory to write files to.

        Defaults to a temporary location in the user's ``${TMPDIR}``.
    project_name : :obj:`bool`, optional
        Name for the project.

        Default is "test_project".
    organism : :obj:`bool`, optional
        Organism of the project.

        Default is "human"
    genome_assembly : :obj:`bool`, optional
        Genome assembly of the project.

        Default is "hg38"
    data_type : :obj:`bool`, optional
        Data type of the project. Must be one of the ``ngs_toolkit`` classes.

        Default is "ATAC-seq"
    only_metadata : obj:`bool`, optional
        Whether to only generate metadata for the project or
        input files in addition.

        Default is :obj:`False`.
    sample_input_files : obj:`bool`, optional
        Whether the input files for the respective data type should be produced.

        This would be BAM and peak files for ATAC-seq or BAM files for RNA-seq.

        Default is :obj:`True`.
    initialize : obj:`bool`, optional
        Whether the project should be initialized into an Analysis object
        for the respective ``data_type`` or simply return the path to a PEP
        configuration file.

        Default is :obj:`True`.
    **kwargs : :obj:`dict`
        Additional keyword arguments will be passed to
        :func:`ngs_toolkit.demo.data_generator.generate_data`.

    Returns
    -------
    {:class:`ngs_toolkit.analysis.Analysis`, :obj:`str`}
        The Analysis object for the project or a path to its PEP configuration
        file.
    """
    from ngs_toolkit import _LOGGER

    if output_dir is None:
        output_dir = tempfile.mkdtemp()

    # Create project with projectmanager
    init_proj(
        project_name,
        genome_assemblies={organism: genome_assembly},
        overwrite=True,
        root_projects_dir=output_dir,
    )

    # Generate random data
    dnum, dcat = generate_data(
        n_factors=n_factors, genome_assembly=genome_assembly, data_type=data_type, **kwargs
    )

    # add additional sample info
    dcat["protocol"] = data_type
    dcat["organism"] = organism
    dcat["data_source"] = "local"
    # now save it
    dcat.to_csv(os.path.join(output_dir, project_name, "metadata", "annotation.csv"))

    # Make comparison table
    comp_table_file = os.path.join(output_dir, project_name, "metadata", "comparison_table.csv")
    ct = pd.DataFrame()
    factors = list(string.ascii_uppercase[:n_factors])
    for factor in factors:
        for side, f in [(1, "2"), (0, "1")]:
            ct2 = dcat.query("{} == '{}'".format(factor, factor + f)).index.to_frame()
            ct2["comparison_side"] = side
            ct2["comparison_name"] = "Factor_" + factor + "_" + "2vs1"
            ct2["sample_group"] = "Factor_" + factor + f
            ct = ct.append(ct2)
    ct["comparison_type"] = "differential"
    ct["data_type"] = data_type
    ct["comparison_genome"] = genome_assembly
    ct.to_csv(comp_table_file, index=False)

    # add the sample_attributes and group_attributes depending on the number of factors
    config_file = os.path.join(output_dir, project_name, "metadata", "project_config.yaml")
    config = yaml.safe_load(open(config_file, "r"))
    factors = list(string.ascii_uppercase[:n_factors])
    config["sample_attributes"] = ["sample_name"] + factors
    config["group_attributes"] = factors
    config["comparison_table"] = comp_table_file
    yaml.safe_dump(config, open(config_file, "w"))

    # prepare dirs
    dirs = [os.path.join(output_dir, project_name, "results")]
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    if not only_metadata:
        # Add consensus region set
        if data_type in REGION_BASED_DATA_TYPES:
            bed = location_index_to_bed(dnum.index)
            bed.to_csv(
                os.path.join(output_dir, project_name, "results", project_name + ".peak_set.bed"),
                index=False,
                sep="\t",
                header=False,
            )
        if isinstance(dnum, pd.DataFrame):
            dnum.to_csv(
                os.path.join(output_dir, project_name, "results", project_name + ".matrix_raw.csv")
            )
        elif isinstance(dnum, dict):
            for res, d in dnum.items():
                d.to_csv(
                    os.path.join(
                        output_dir,
                        project_name,
                        "results",
                        project_name + "." + res + ".matrix_raw.csv",
                    )
                )

    # Here we are raising the logger level to omit messages during
    # sample file creation which would have been too much clutter,
    # especially since I make use of `load_data` with permissive=True.
    # In dev, I still want to see them.
    # if DEV:
    # _LOGGER.debug("Shutting logger for now")
    prev_level = _LOGGER.getEffectiveLevel()
    _LOGGER.setLevel("ERROR")
    an = initialize_analysis_of_data_type(data_type, config_file)
    an.load_data(permissive=True)
    if sample_input_files:
        generate_sample_input_files(an, dnum)
    # if DEV:
    _LOGGER.setLevel(prev_level)
    # _LOGGER.debug("Reactivated logger")
    if not initialize:
        return config_file
    return an


def generate_projects(
    output_path=None,
    project_prefix_name="demo-project",
    data_types=["ATAC-seq", "ChIP-seq", "CNV", "RNA-seq"],
    organisms=["human", "mouse"],
    genome_assemblies=["hg38", "mm10"],
    n_factors=[1, 2, 5],
    n_features=[100, 1000, 10000],
    n_replicates=[1, 3, 5],
    **kwargs
):
    """
    Create a list of Projects given ranges of parameters, which will be passed
    to :func:`ngs_toolkit.demo.data_generator.generate_project`.
    """
    ans = list()
    for data_type in data_types:
        for organism, genome_assembly in zip(organisms, genome_assemblies):
            for factors in n_factors:
                for features in n_features:
                    for replicates in n_replicates:
                        project_name = "_".join(
                            str(x)
                            for x in [
                                project_prefix_name,
                                data_type,
                                genome_assembly,
                                factors,
                                features,
                                replicates,
                            ]
                        )

                        an = generate_project(
                            output_dir=output_path,
                            project_name=project_name,
                            organism=organism,
                            genome_assembly=genome_assembly,
                            data_type=data_type,
                            n_factors=factors,
                            n_replicates=replicates,
                            n_features=features,
                            **kwargs
                        )
                        ans.append(an)
    return ans


def generate_bam_file(
    count_vector, output_bam, genome_assembly="hg38", chrom_sizes_file=None, index=True
):
    """Generate BAM file containing reads matching the counts in a vector of features"""
    s = location_index_to_bed(count_vector.index)

    # get reads per region
    i = [i for i, c in count_vector.iteritems() for _ in range(c)]
    s = s.reindex(i)

    # shorten/enlarge by a random fraction; name reads
    d = s["end"] - s["start"]
    s = s.assign(
        start=(s["start"] + d * np.random.uniform(-0.2, 0.2, s.shape[0])).astype(int),
        end=(s["end"] + d * np.random.uniform(-0.2, 0.2, s.shape[0])).astype(int),
        name=["{}_read_{}".format(count_vector.name, i) for i in range(s.shape[0])],
    )

    s = pybedtools.BedTool.from_dataframe(s).truncate_to_chrom(genome=genome_assembly).sort()
    # get a file with chromosome sizes (usually not needed but only for bedToBam)
    if chrom_sizes_file is None:
        chrom_sizes_file = tempfile.NamedTemporaryFile().name
        pybedtools.get_chromsizes_from_ucsc(genome=genome_assembly, saveas=chrom_sizes_file)
    s.to_bam(g=chrom_sizes_file).saveas(output_bam)

    if index:
        import pysam

        pysam.index(output_bam)


def generate_peak_file(peak_set, output_peak, genome_assembly="hg38", summits=False):
    """Generate peak files containing regions from a fraction of a given set of features"""
    if not isinstance(peak_set, pybedtools.BedTool):
        peak_set = pybedtools.BedTool(peak_set)

    s = peak_set.to_dataframe()

    # choose a random but non-empty fraction of sites to keep
    while True:
        s2 = s.sample(frac=np.random.uniform())
        if not s2.empty:
            break
    s = pybedtools.BedTool.from_dataframe(s2)
    # shorten/enlarge sites by a random fraction
    s = s.slop(
        l=np.random.uniform(-0.2, 0.2),
        r=np.random.uniform(-0.2, 0.2),
        pct=True,
        genome=genome_assembly,
    )

    if summits:
        # get middle basepair
        s = s.to_dataframe()
        mid = ((s["end"] - s["start"]) / 2).astype(int)
        s.loc[:, "start"] += mid
        s.loc[:, "end"] -= mid
        s = pybedtools.BedTool.from_dataframe(s)

    s = s.sort()
    s.saveas(output_peak)


def generate_log2_profiles(sample_vector, background_vector, output_file):
    """Generate a file with read count profile of CNVs for the sample and background"""
    # raise NotImplementedError

    s = location_index_to_bed(sample_vector.index)
    s["log2." + sample_vector.name + ".trimmed.bowtie2.filtered.bam"] = sample_vector
    s["log2.background.trimmed.bowtie2.filtered.bam"] = background_vector
    s.index.name = "Feature"
    with open(output_file, "w") as handle:
        handle.write("# This is a comment made by the CNV caller program\n")
        handle.write("# This is a another\n")
        s.to_csv(handle, sep="\t")


def generate_sample_input_files(analysis, matrix):
    """Generate input files (BAM, peaks) for a sample depending on its data type."""
    if analysis.data_type in REGION_BASED_DATA_TYPES:
        chrom_sizes_file = tempfile.NamedTemporaryFile().name
        pybedtools.get_chromsizes_from_ucsc(genome=analysis.genome, saveas=chrom_sizes_file)

        if not hasattr(analysis, "sites"):
            analysis.load_data(only_these_keys=["sites"], permissive=True)
        if not hasattr(analysis, "sites"):
            raise AttributeError("Need a consensus peak set to generate sample input files.")

    for sample in analysis.samples:
        if hasattr(sample, "aligned_filtered_bam"):
            if sample.aligned_filtered_bam is not None:
                d = os.path.dirname(sample.aligned_filtered_bam)
                os.makedirs(d, exist_ok=True)
                generate_bam_file(
                    matrix.loc[:, sample.name],
                    sample.aligned_filtered_bam,
                    genome_assembly=analysis.genome,
                    chrom_sizes_file=chrom_sizes_file,
                )
        if hasattr(sample, "peaks"):
            if sample.peaks is not None:
                d = os.path.dirname(sample.peaks)
                os.makedirs(d, exist_ok=True)
                generate_peak_file(
                    analysis.sites, sample.peaks, summits=False, genome_assembly=analysis.genome
                )
        if hasattr(sample, "summits"):
            if sample.summits is not None:
                d = os.path.dirname(sample.summits)
                os.makedirs(d, exist_ok=True)
                generate_peak_file(
                    analysis.sites, sample.summits, summits=True, genome_assembly=analysis.genome
                )

        if hasattr(sample, "log2_read_counts"):
            if sample.log2_read_counts is not None:
                for res, file in sample.log2_read_counts.items():
                    os.makedirs(os.path.dirname(file), exist_ok=True)
                    generate_log2_profiles(
                        (2 ** matrix[res].loc[:, sample.name]).astype(int),
                        (2 ** matrix[res].loc[:, sample.name]).astype(
                            int
                        ),  # this should be the background vector
                        file,
                    )


def initialize_analysis_of_data_type(data_type, pep_config, *args, **kwargs):
    """Initialize an Analysis object from a PEP config with the appropriate ``data_type``."""
    from ngs_toolkit import Analysis

    sub = Analysis.__subclasses__()
    sub += [sc for x in sub for sc in x.__subclasses__()]

    m = {t._data_type: t for t in sub}
    return m[data_type](from_pep=pep_config, *args, **kwargs)


def get_random_genomic_locations(
    n_regions, width_mean=500, width_std=400, min_width=300, genome_assembly="hg38"
):
    """Get `n_regions`` number of random genomic locations respecting the boundaries of the ``genome_assembly``"""
    from ngs_toolkit.utils import bed_to_index

    # weight chroms by their size, excluding others
    csizes = {
        k: v[-1] for k, v in dict(pybedtools.chromsizes(genome_assembly)).items() if "_" not in k
    }
    gsize = sum(csizes.values())
    csizes = {k: v / gsize for k, v in csizes.items()}
    chrom = pd.Series(
        np.random.choice(a=list(csizes.keys()), size=n_regions, p=list(csizes.values()))
    )
    start = np.array([0] * n_regions)
    end = np.absolute(np.random.normal(width_mean, width_std, n_regions)).astype(int)
    df = pd.DataFrame([chrom.tolist(), start.tolist(), end.tolist()]).T
    df.loc[(df[2] - df[1]) < min_width, 2] += min_width
    bed = (
        pybedtools.BedTool.from_dataframe(df)
        .shuffle(genome=genome_assembly, chromFirst=True, noOverlapping=True, chrom=True)
        .sort()
        .to_dataframe()
    )
    return bed_to_index(bed)


def get_random_genes(n_genes, genome_assembly="hg38"):
    """Get ``n_genes`` number of random genes from the set of genes of the ``genome_assembly``"""

    m = {"hg19": "grch37", "hg38": "grch38", "mm10": "grcm38"}
    o = {"hg19": "hsapiens", "hg38": "hsapiens", "mm10": "mmusculus"}

    g = (
        query_biomart(
            attributes=["external_gene_name"],
            ensembl_version=m[genome_assembly],
            species=o[genome_assembly],
        )
        .squeeze()
        .drop_duplicates()
    )
    return g.sample(n=n_genes, replace=False).sort_values()


def get_random_grnas(n_genes, genome_assembly="hg38", n_grnas_per_gene=4):
    """Get ``n_genes`` number of random genes from the set of genes of the ``genome_assembly``"""
    import itertools

    gs = int(n_genes / n_grnas_per_gene)
    genes = get_random_genes(gs, genome_assembly=genome_assembly)
    grna_ns = itertools.chain(*[range(n_grnas_per_gene) for _ in range(gs)])

    return [a + "__" + str(b) for a, b in zip(genes.repeat(n_grnas_per_gene), list(grna_ns))]


def get_genomic_bins(n_bins, genome_assembly="hg38", resolution=None):
    """Get a ``size`` number of random genomic bins respecting the boundaries of the ``genome_assembly``"""
    from ngs_toolkit.utils import bed_to_index

    bed = pybedtools.BedTool.from_dataframe(
        pd.DataFrame(dict(pybedtools.chromsizes(genome_assembly))).T.reset_index()
    )
    w = bed.makewindows(
        genome=genome_assembly, w=sum([i.length for i in bed]) / n_bins
    ).to_dataframe()
    if resolution is not None:
        if isinstance(resolution, str):
            resolution = int(resolution.replace("kb", "000"))
        w["end"] = w["start"] + resolution
    return bed_to_index(w.head(n_bins))


def _disp(x):
    return 4 / x + 0.1
