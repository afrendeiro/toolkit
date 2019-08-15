#!/usr/bin/env python


import os
import string

from ngs_toolkit.general import query_biomart
from ngs_toolkit.utils import location_index_to_bed
from ngs_toolkit.project_manager import create_project
import numpy as np
import pandas as pd
import patsy
import pybedtools
import yaml


class RandomDataGenerator(object):
    def generate_random_data(
        self,
        n_factors=2,
        n_variables=20000,
        n_replicates=5,
        distribution="negative_binomial",
        group_fold_differences=5,
        fraction_of_different=0.2,
        data_type="ATAC-seq",
        genome_assembly="hg19",
    ):
        if not isinstance(group_fold_differences, list):
            # _LOGGER.warning("Assuming same fold difference for all factors between samples")
            group_fold_differences = [group_fold_differences] * n_factors

        n_samples = n_factors * n_replicates
        s = list((string.ascii_lowercase[15:] * n_variables)[:n_variables])
        s = [i + str(j) for i, j in zip(s, range(n_variables))]
        d = patsy.demo_data(
            *(list(string.ascii_lowercase[:n_factors]) + s),
            nlevels=2,
            min_rows=n_samples
        )
        d = pd.DataFrame(d)
        dcat = d.loc[:, d.dtypes != np.float]
        dnum = d.loc[:, d.dtypes == np.float]

        # add sample names
        names = ["S{}_{}".format(i + 1, dcat.loc[i, :].sum()) for i in dcat.index]
        dcat.index = dnum.index = names
        dcat.index.name = dnum.index.name = "sample_name"

        # add variation to groups
        for i, factor in enumerate(dcat.columns):
            af = dnum.columns.to_series().sample(frac=fraction_of_different).tolist()
            levels = dcat.loc[:, factor].unique()

            if len(levels) == 2:
                for level, f in zip(levels, [1, -1]):
                    diff = (
                        dnum.loc[:, af].mean(axis=0)
                        * (group_fold_differences[i] / 2)
                        * f
                    )  # * np.absolute(np.random.normal(0, 0.1)))
                    dnum.loc[
                        dcat.loc[:, factor] == level, dnum.columns.isin(af)
                    ] += diff
            else:
                raise NotImplementedError

        # add intersect from distribution
        dist = getattr(np.random, distribution)
        dnum = dnum + dist(1, 0.1, n_variables)

        # add overdispersion across all features
        mean = dnum.mean(axis=0)
        space = np.linspace(mean.max(), mean.min(), 100)
        step = space[0] - space[1]
        for i, lim in enumerate(space):
            cur = mean[(mean <= lim + step) & (mean > lim)].index
            dnum.loc[:, cur] *= np.random.normal(0, i + 1, len(cur))

        # make non-negative
        if dnum.min().min() < 0:
            dnum -= dnum.min().min()

        # make integer
        dnum = dnum.astype(int)

        # add random location indexes
        if data_type in ["ATAC-seq", "ChIP-seq"]:
            dnum.columns = self.get_random_genomic_locations(
                n_variables, genome_assembly=genome_assembly
            )
        if data_type in ["RNA-seq"]:
            dnum.columns = self.get_random_genes(
                n_variables, genome_assembly=genome_assembly
            )
        if data_type in ["CNV"]:
            from ngs_toolkit.utils import z_score

            dnum.columns = self.get_genomic_bins(
                n_variables, genome_assembly=genome_assembly
            )
            dnum = z_score(dnum)

        return dnum.T, dcat

    @staticmethod
    def get_random_genomic_locations(
        size, width_mean=500, width_std=400, min_width=300, genome_assembly="hg19"
    ):
        from ngs_toolkit.utils import bed_to_index

        chrom = ["chr1"] * size
        start = np.array([0] * size)
        end = np.absolute(np.random.normal(width_mean, width_std, size)).astype(int)
        df = pd.DataFrame([chrom, start.tolist(), end.tolist()]).T
        df.loc[(df[2] - df[1]) < min_width, 2] += min_width
        bed = (
            pybedtools.BedTool.from_dataframe(df)
            .shuffle(genome=genome_assembly)
            .sort()
            .to_dataframe()
        )
        return bed_to_index(bed)

    @staticmethod
    def get_random_genes(size, genome_assembly="hg19"):
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
        return pd.Series(np.random.choice(g, size, replace=False)).sort_values()

    @staticmethod
    def get_genomic_bins(n_bins, distribution="normal", genome_assembly="hg19"):
        from ngs_toolkit.utils import bed_to_index

        bed = pybedtools.BedTool.from_dataframe(
            pd.DataFrame(dict(pybedtools.chromsizes(genome_assembly))).T.reset_index()
        )
        w = bed.makewindows(
            genome=genome_assembly, w=sum([i.length for i in bed]) / n_bins
        ).to_dataframe()
        return bed_to_index(w.head(n_bins))


def generate_project(
    output_dir="tests",
    project_name="test_project",
    organism="human",
    genome_assembly="hg19",
    data_type="ATAC-seq",
    only_metadata=False,
    sample_input_files=False,
    **kwargs
):
    output_dir = os.path.abspath(output_dir)

    # Create project with projectmanager
    create_project(
        project_name,
        genome_assemblies={organism: genome_assembly},
        overwrite=True,
        root_projects_dir=output_dir,
    )

    # Generate random data
    g = RandomDataGenerator()
    n, c = g.generate_random_data(
        genome_assembly=genome_assembly, data_type=data_type, **kwargs
    )

    # add additional sample info
    c["protocol"] = data_type
    c["organism"] = organism
    # now save it
    c.to_csv(os.path.join(output_dir, project_name, "metadata", "annotation.csv"))

    # Make comparison table
    if "n_factors" in kwargs:
        table_file = os.path.join(
            output_dir, project_name, "metadata", "comparison_table.csv"
        )
        ct = pd.DataFrame()
        for factor in list(string.ascii_lowercase[: kwargs["n_factors"]]):
            for side, f in [(1, "2"), (0, "1")]:
                ct2 = c.loc[c[factor] == factor + f].index.to_frame()
                ct2["comparison_side"] = side
                ct2["comparison_name"] = "Factor_" + factor + "_" + "2vs1"
                ct2["sample_group"] = "Factor_" + factor + f
                ct = ct.append(ct2)
        ct["comparison_type"] = "differential"
        ct["data_type"] = data_type
        ct["comparison_genome"] = genome_assembly
        ct.to_csv(table_file, index=False)

    # add the sample_attributes and group_attributes depending on the number of factors
    if "n_factors" in kwargs:
        config_file = os.path.join(
            output_dir, project_name, "metadata", "project_config.yaml"
        )
        config = yaml.safe_load(open(config_file, "r"))
        factors = list(string.ascii_lowercase[: kwargs["n_factors"]])
        config["sample_attributes"] = ["sample_name"] + factors
        config["group_attributes"] = factors
        config["metadata"]["comparison_table"] = table_file
        yaml.safe_dump(config, open(config_file, "w"))

    # prepare dirs
    dirs = [os.path.join(output_dir, project_name, "results")]
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    if not only_metadata:
        if data_type == "ATAC-seq":
            bed = location_index_to_bed(n.index)
            bed.to_csv(
                os.path.join(
                    output_dir, project_name, "results", project_name + ".peak_set.bed"
                ),
                index=False,
                sep="\t",
                header=False,
            )
        n.to_csv(
            os.path.join(
                output_dir, project_name, "results", project_name + ".matrix_raw.csv"
            )
        )


def generate_bam_file(peak_set, output_bam, genome_assembly="hg38", chrom_sizes_file=None, index=True):
    import tempfile
    if not isinstance(peak_set, pybedtools.BedTool):
        peak_set = pybedtools.BedTool(peak_set)

    s = peak_set.to_dataframe()

    # choose a random number of reads per region
    r = np.random.negative_binomial(1000, 0.99, s.shape[0])
    i = [i for i in r for _ in range(i)]
    s = s.reindex(i)

    # shorten/enlarge by a random fraction; name reads
    d = s['end'] - s['start']
    s = s.assign(
        start=(s['start'] + d * np.random.uniform(-0.2, 0.2, s.shape[0])).astype(int),
        end=(s['end'] + d * np.random.uniform(-0.2, 0.2, s.shape[0])).astype(int),
        name=["read_{}".format(i) for i in range(s.shape[0])])

    s = pybedtools.BedTool.from_dataframe(s)
    s = s.truncate_to_chrom(genome=genome_assembly)
    s = s.sort()
    # get a file with chromosome sizes (usually not needed but only for bedToBam)
    if chrom_sizes_file is None:
        chrom_sizes_file = tempfile.NamedTemporaryFile().name
        pybedtools.get_chromsizes_from_ucsc(genome=genome_assembly, saveas=chrom_sizes_file)
    b = s.to_bam(g=chrom_sizes_file)

    b.saveas(output_bam)

    if index:
        import pysam
        pysam.index(output_bam)


def generate_peak_file(peak_set, output_peak, genome_assembly="hg38", summits=False):
    if not isinstance(peak_set, pybedtools.BedTool):
        peak_set = pybedtools.BedTool(peak_set)

    s = peak_set.to_dataframe()

    # choose a random fraction of sites to keep
    s = pybedtools.BedTool.from_dataframe(
            s.sample(frac=np.random.uniform()))
    # shorten/enlarge sites by a random fraction
    s = s.slop(
        l=np.random.uniform(-0.2, 0.2),
        r=np.random.uniform(-0.2, 0.2),
        pct=True, genome=genome_assembly)

    if summits:
        # get middle basepair
        s = s.to_dataframe()
        mid = ((s['end'] - s['start']) / 2).astype(int)
        s.loc[:, 'start'] += mid
        s.loc[:, 'end'] -= (mid - 1)
        s = pybedtools.BedTool.from_dataframe(s)

    s = s.sort()
    s.saveas(output_peak)


def add_sample_input_files(analysis):
    import tempfile

    chrom_sizes_file = tempfile.NamedTemporaryFile().name
    pybedtools.get_chromsizes_from_ucsc(genome=analysis.genome, saveas=chrom_sizes_file)

    for sample in analysis.samples:
        if hasattr(sample, "aligned_filtered_bam"):
            d = os.path.dirname(sample.aligned_filtered_bam)
            if not os.path.exists(d):
                os.makedirs(d)
            generate_bam_file(
                analysis.sites, sample.aligned_filtered_bam,
                genome_assembly=analysis.genome,
                chrom_sizes_file=chrom_sizes_file)
        if hasattr(sample, "peaks"):
            d = os.path.dirname(sample.peaks)
            if not os.path.exists(d):
                os.makedirs(d)
            generate_peak_file(
                analysis.sites, sample.peaks, summits=False,
                genome_assembly=analysis.genome)
        if hasattr(sample, "summits"):
            d = os.path.dirname(sample.summits)
            if not os.path.exists(d):
                os.makedirs(d)
            generate_peak_file(
                analysis.sites, sample.summits, summits=True,
                genome_assembly=analysis.genome)
