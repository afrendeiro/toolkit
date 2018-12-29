#!/usr/bin/env python


class RandomDataGenerator(object):
    def __init__(self):
        pass

    def generate_random_data(
            self,
            n_factors=2, n_variables=100000, n_replicates=5,
            distribution="negative_binomial", group_fold_differences=2, variation=1,
            data_type="ATAC-seq", genome_assembly="hg19"):
        import string
        import patsy
        import pandas as pd
        import numpy as np

        if not isinstance(group_fold_differences, list):
            # _LOGGER.warning("Assuming same fold difference for all factors between samples")
            group_fold_differences = [group_fold_differences] * n_factors

        n_samples = n_factors * n_replicates
        s = list((string.ascii_lowercase[15:] * n_variables)[:n_variables])
        s = [i + str(j) for i, j in zip(s, range(n_variables))]
        d = patsy.demo_data(
            * (list(string.ascii_lowercase[:n_factors]) + s), nlevels=2, min_rows=n_samples)
        d = pd.DataFrame(d)
        dcat = d.loc[:, d.dtypes != np.float]
        dnum = d.loc[:, d.dtypes == np.float]

        # add sample names
        names = ["S{}_{}".format(i + 1, dcat.loc[i, :].sum()) for i in dcat.index]
        dcat.index = dnum.index = names
        dcat.index.name = dnum.index.name = "sample_name"

        # add variation to groups
        for i, factor in enumerate(dcat.columns):
            levels = dcat.loc[:, factor].unique()

            if len(levels) == 2:
                for level, f in zip(levels, [1, -1]):
                    diff = dnum.mean(axis=0) * (group_fold_differences[i] / 2) * f
                    dnum.loc[dcat.loc[:, factor] == level, :] += diff
            else:
                raise NotImplementedError

        # add intersect from distribution
        dist = getattr(np.random, distribution)
        dnum = (dnum + dist(1, 0.1, n_variables))

        # make non-negative
        if dnum.min().min() < 0:
            dnum -= dnum.min().min()

        # make integer
        dnum = dnum.astype(int)

        # add random location indexes
        if data_type in ["ATAC-seq", "ChIP-seq"]:
            dnum.columns = self.get_random_genomic_locations(n_variables, genome_assembly=genome_assembly)
        if data_type in ["RNA-seq"]:
            dnum.columns = self.get_random_genes(n_variables, genome_assembly=genome_assembly)

        return dnum.T, dcat

    @staticmethod
    def get_random_genomic_locations(
            size,
            width_mean=500, width_std=400, min_width=300, distribution="normal",
            genome_assembly="hg19"):
        import pybedtools
        import numpy as np
        import pandas as pd

        chrom = ['chr1'] * size
        start = np.array([0] * size)
        end = np.absolute(np.random.normal(width_mean, width_std, size)).astype(int)
        df = pd.DataFrame([chrom, start.tolist(), end.tolist()]).T
        df.loc[(df[2] - df[1]) < min_width, 2] += min_width
        bed = pybedtools.BedTool.from_dataframe(df).shuffle(genome=genome_assembly).sort().to_dataframe()
        return bed['chrom'] + ":" + bed['start'].astype(str) + "-" + bed['end'].astype(str)

    @staticmethod
    def get_random_genes(
            size,
            genome_assembly="hg19"):
        from ngs_toolkit.general import query_biomart
        import numpy as np
        import pandas as pd

        m = {"hg19": "grch37", "hg38": "grch38", "mm10": "grcm38"}

        g = query_biomart(
            attributes=['external_gene_name'],
            ensembl_version=m[genome_assembly]).squeeze()
        return pd.Series(np.random.choice(g, size)).sort_values()


def generate_project(
        output_dir="tests", project_name="test_project",
        organism="human", genome_assembly="hg19", data_type="ATAC-seq",
        only_metadata=False, **kwargs):
    from ngs_toolkit.project_manager import create_project
    import os
    import yaml

    # Create project with projectmanager
    create_project(
        project_name, output_dir,
        genome_assemblies={organism: genome_assembly},
        overwrite=True)

    # Generate random data
    g = RandomDataGenerator()
    n, c = g.generate_random_data(genome_assembly=genome_assembly, **kwargs)

    # add additional sample info
    c['protocol'] = data_type
    c['organism'] = organism
    # now save it
    c.to_csv(os.path.join(output_dir, project_name, "metadata", "annotation.csv"))

    # add the sample_attributes and group_attributes depending on the number of factors
    if "n_factors" in kwargs:
        import string
        config_file = os.path.join(output_dir, project_name, "metadata", "project_config.yaml")
        config = yaml.safe_load(
            open(config_file, "r"))
        factors = list(string.ascii_lowercase[:kwargs['n_factors']])
        config['sample_attributes'] = ['sample_name'] + factors
        config['group_attributes'] = factors
        yaml.safe_dump(config, open(config_file, "w"))

    # prepare dirs
    dirs = [os.path.join(output_dir, project_name, "results")]
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)

    if not only_metadata:
        if data_type == "ATAC-seq":
            from ngs_toolkit.general import location_index_to_bed
            bed = location_index_to_bed(n.index)
            bed.to_csv(os.path.join(
                output_dir, project_name, "results", project_name + "_peak_set.bed"), index=False, sep="\t", header=False)
            n.to_csv(os.path.join(
                output_dir, project_name, "results", project_name + "_peaks.raw_coverage.csv"))
        elif data_type == "RNA-seq":
            n.to_csv(os.path.join(
                output_dir, project_name, "results", project_name + ".expression_counts.gene_level.csv"))
