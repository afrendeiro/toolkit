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
            dnum.columns = self.get_random_genomic_locations(n_variables, genome_assembly=genome_assembly)

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
        bed = pybedtools.BedTool.from_dataframe(df).shuffle(genome=genome_assembly).to_dataframe()
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
        return pd.Series(np.random.choice(g, size))


def generate_project(
        output_dir="tests", project_name="test_project",
        organism="human", genome_assembly="hg19", data_type="ATAC-seq",
        only_metadata=False, **kwargs):
    from ngs_toolkit.project_manager import create_project
    import os

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


def test_project_manager(tmp_path):
    from ngs_toolkit import _CONFIG
    from ngs_toolkit.project_manager import create_project
    import os
    import pandas as pd
    import shutil

    project_name = "test_project"
    annotation_vars = [
        "sample_name", "toggle", "pass_qc", "protocol", "library",
        "cell_line", "cell_type", "condition",
        "experimental_batch", "experiment_name", "replicate",
        "organism", "flowcell", "lane", "BSF_name", "data_source"]

    genome_assemblies = {k: v for x in _CONFIG["default_genome_assemblies"] for k, v in x.items()}
    create_project(
        project_name, tmp_path, genome_assemblies=genome_assemblies, overwrite=True)

    expected_files = [
        os.path.join(tmp_path, project_name, ".git"),
        os.path.join(tmp_path, project_name, "metadata"),
        os.path.join(tmp_path, project_name, "metadata", "project_config.yaml"),
        os.path.join(tmp_path, project_name, "metadata", "annotation.csv"),
        os.path.join(tmp_path, project_name, "metadata", "merge_table.csv"),
        os.path.join(tmp_path, project_name, "metadata", "comparison_table.csv"),
    ]
    for f in expected_files:
        assert os.path.exists(f)

    df = pd.read_csv(os.path.join(tmp_path, project_name, "metadata", "annotation.csv"))
    assert df.shape == (0, len(annotation_vars))
    assert all(c in df.columns for c in annotation_vars)

    shutil.rmtree(tmp_path)


def test_analysis_creation(tmp_path):
    from ngs_toolkit.general import Analysis
    import os
    from peppy import Project
    import yaml
    import shutil

    name = "test_analysis"

    a = Analysis(name=name)
    assert a.__repr__() == "Analysis object named '{}'.".format(name)
    assert "samples" not in a.__repr__()

    # Let's make several "reallish" test projects
    project_prefix_name = "test-project"
    data_types = ["ATAC-seq", "RNA-seq", "ChIP-seq"]  # "CNV"
    genome_assemblies = [("human", "hg19"), ("human", "hg38"), ("mouse", "mm10")]

    params = {
        "ATAC-seq": {
            "n_factors": [1, 2, 3],
            "n_variables": [100, 1000, 10000],
            "n_replicates": [1, 2, 5],
        },
        "ChIP-seq": {
            "n_factors": [1, 2, 3],
            "n_variables": [100, 1000, 10000],
            "n_replicates": [1, 2, 5],
        },
        "RNA-seq": {
            "n_factors": [1, 2, 3],
            "n_variables": [100, 1000, 25000],
            "n_replicates": [1, 2, 5],
        },
    }

    for data_type in data_types:
        n_factors = params[data_type]['n_factors'][0]
        n_variables = params[data_type]['n_variables'][0]
        n_replicates = params[data_type]['n_replicates'][0]
        for organism, genome_assembly in genome_assemblies:

            project_name = "{}_{}_{}_{}_{}_{}".format(
                project_prefix_name, data_type, genome_assembly,
                n_factors, n_variables, n_replicates
            )
            n_samples = (n_factors * n_replicates) + n_factors

            generate_project(
                output_dir=tmp_path,
                project_name=project_name, genome_assembly=genome_assembly, data_type=data_type,
                n_factors=n_factors, n_replicates=n_replicates, n_variables=n_variables)

            # first edit the defaul path to the annotation sheet
            config = os.path.join(
                tmp_path, project_name, "metadata", "project_config.yaml")
            c = yaml.safe_load(open(config, 'r'))
            c['metadata']['sample_annotation'] = os.path.abspath(
                os.path.join(tmp_path, project_name, "metadata", "annotation.csv"))
            c['metadata']['comparison_table'] = os.path.abspath(
                os.path.join(tmp_path, project_name, "metadata", "comparison_table.csv"))
            yaml.safe_dump(c, open(config, "w"))

            # project and associated analysis
            prj = Project(config)
            a = Analysis(name=project_name, prj=prj)
            assert a.__repr__() == (
                "Analysis object named '{}' with {} samples of genome '{}'."
                .format(project_name, n_samples, genome_assembly))
            assert len(prj.samples) == len(a.samples)
            assert all([x == y for x, y in zip(prj.samples, a.samples)])

            shutil.rmtree(tmp_path)


def test_analysis_serialization(tmp_path):
    from ngs_toolkit.general import Analysis
    import shutil
    import os
    import numpy as np

    pickle_file = os.path.join(tmp_path, "pickle")
    a = Analysis(pickle_file=pickle_file)
    assert not os.path.exists(pickle_file)
    a.to_pickle()
    assert os.path.exists(pickle_file)
    assert os.stat(pickle_file).st_size > 0

    previous_size = os.stat(pickle_file).st_size
    a.random = np.random.random((100, 100))
    a.to_pickle()
    new_size = os.stat(pickle_file).st_size
    assert new_size > previous_size

    shutil.rmtree(tmp_path)


def test_analysis_loading(tmp_path):
    from ngs_toolkit.general import Analysis
    import shutil
    import os

    pickle_file = os.path.join(tmp_path, "pickle")
    a = Analysis(pickle_file=pickle_file)
    a.secret = "I've existed before"
    a.to_pickle()

    a2 = Analysis(pickle_file=pickle_file, from_pickle=True)
    assert a2.secret == "I've existed before"

    a3 = Analysis()
    a3.update(pickle_file)
    assert a3.secret == "I've existed before"

    a4 = Analysis(pickle_file=pickle_file).from_pickle()
    assert a4.secret == "I've existed before"

    shutil.rmtree(tmp_path)
