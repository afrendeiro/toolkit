#!/usr/bin/env python


import os
import pandas as pd


class RandomDataGenerator(object):
    def __init__(self):
        pass

    def generate_random_data(
            self,
            n_factors=2, n_variables=100000, n_replicates=5,
            distribution="negative_binomial", group_fold_differences=2, variation=1,
            data_type="ATAC-seq", genome="hg19"):
        import string
        import patsy
        import numpy as np

        if not isinstance(group_fold_differences, list):
            # _LOGGER.warning("Assuming same fold difference for all factors between samples")
            group_fold_differences = [group_fold_differences] * n_factors

        n_samples = n_factors * n_replicates
        s = list((string.ascii_lowercase[15:] * n_variables)[:n_variables])
        s = [i + str(j) for i, j in zip(s, range(n_variables))]
        d = patsy.demo_data(
            *string.ascii_lowercase[:n_factors], *s, nlevels=2, min_rows=n_samples)
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
            dnum.columns = self.get_random_genomic_locations(n_variables)

        return dnum.T, dcat

    @staticmethod
    def get_random_genomic_locations(
            size,
            width_mean=500, width_std=400, min_width=300, distribution="normal",
            genome="hg19"):
        import pybedtools
        import numpy as np

        chrom = ['chr1'] * size
        start = np.array([0] * size)
        end = np.absolute(np.random.normal(width_mean, width_std, size)).astype(int)
        df = pd.DataFrame([chrom, start.tolist(), end.tolist()]).T
        df.loc[(df[2] - df[1]) < min_width, 2] += min_width
        bed = pybedtools.BedTool.from_dataframe(df).shuffle(genome=genome).to_dataframe()
        return bed['chrom'] + ":" + bed['start'].astype(str) + "-" + bed['end'].astype(str)


class GenerateProject(object):
    def __init__(
            self,
            output_dir="tests", project_name="test_project",
            genome="hg19", data_type="ATAC-seq",
            only_metadata=False, **kwargs):
        from ngs_toolkit.project_manager import create_project

        # Create project with projectmanager
        create_project(
            project_name, output_dir, overwrite=True,
            username="arendeiro", email="{username}@cemm.oeaw.ac.at",
            url="http://biomedical-sequencing.at/bocklab/{username}/{project_name}")

        # Generate random data
        g = RandomDataGenerator()
        n, c = g.generate_random_data(**kwargs)

        # add additional sample info
        c['protocol'] = data_type
        if genome.startswith("hg"):
            c['organism'] = "human"
        elif genome.startswith("mm"):
            c['organism'] = "mouse"
        else:
            c['organism'] = "human"

        # now save it
        c.to_csv(os.path.join(output_dir, project_name, "metadata", "annotation.csv"))

        dirs = [os.path.join(output_dir, project_name, "results")]
        for d in dirs:
            if not os.path.exists(d):
                os.makedirs(d)

        if not only_metadata:
            if data_type == "ATAC-seq":
                n.to_csv(os.path.join(output_dir, project_name, "results", project_name + "_peaks.raw_coverage.csv"))
            elif data_type == "RNA-seq":
                n.to_csv(os.path.join(output_dir, project_name, "results", project_name + ".expression_counts.gene_level.csv"))


def test_version_matches():
    from ngs_toolkit import __version__ as installed_version
    import pkgutil

    file_version = pkgutil.get_data('ngs_toolkit', "_version.py").decode().strip().split(" = ")[1].replace("\"", "")
    assert installed_version == file_version


def test_all_requirements_are_importable():
    import requests
    import importlib

    package_name = 'ngs-toolkit'
    url = 'https://pypi.python.org/pypi/' + str(package_name) + '/json'
    data = requests.get(url).json()

    replace = {"scikit-learn": "sklearn", "piper": "pypiper"}

    requirements = [x.split(" ")[0] for x in data['info']['requires_dist']]
    # versions = [x.split(" ")[1] if len(x.split(" ")) > 1 else "" for x in data['info']['requires_dist']]

    for req in requirements:
        if req in replace:
            requirements.pop(requirements.index(req))
            requirements.append(replace[req])

    for req in requirements:
        try:
            importlib.import_module(req)
        except ImportError:
            assert False


def test_config_has_all_required_fields():
    from ngs_toolkit import _CONFIG as local_config
    import pkgutil
    import yaml

    def _dicts_same_keys(d1, d2):
        if type(d1) != type(d2):
            return False

        for k in d1.keys():
            if k not in d2:
                return False
            else:
                if type(d1[k]) is dict:
                    return _dicts_same_keys(d1[k], d2[k])
                else:
                    return True

    file_config = pkgutil.get_data('ngs_toolkit', "config/default.yaml").decode().strip()
    file_config = yaml.load(file_config)

    assert _dicts_same_keys(file_config, local_config)


def test_project_manager():
    from ngs_toolkit.project_manager import create_project
    import shutil

    project_name = "test_project"
    root_dir = "tests"
    annotation_vars = [
        "sample_name", "toggle", "pass_qc", "protocol", "library",
        "cell_line", "cell_type", "condition",
        "experimental_batch", "experiment_name", "replicate",
        "organism", "flowcell", "lane", "BSF_name", "data_source"]

    create_project(
        project_name, root_dir, overwrite=True,
        username="arendeiro", email="{username}@cemm.oeaw.ac.at",
        url="http://biomedical-sequencing.at/bocklab/{username}/{project_name}")

    expected_files = [
        os.path.join(root_dir, project_name, ".git"),
        os.path.join(root_dir, project_name, "metadata"),
        os.path.join(root_dir, project_name, "metadata", "project_config.yaml"),
        os.path.join(root_dir, project_name, "metadata", "annotation.csv"),
        os.path.join(root_dir, project_name, "metadata", "merge_table.csv"),
        os.path.join(root_dir, project_name, "metadata", "comparison_table.csv"),
    ]
    for f in expected_files:
        assert os.path.exists(f)

    df = pd.read_csv(os.path.join(root_dir, project_name, "metadata", "annotation.csv"))
    assert df.shape == (0, len(annotation_vars))
    assert all(c in df.columns for c in annotation_vars)

    shutil.rmtree(os.path.join("tests", project_name))


def test_analysis_creation():
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
    project_name = "test_project"
    data_type = "ATAC-seq"
    genome_assembly = "hg19"
    n_factors = 2
    n_variables = 100000
    n_replicates = 5

    n_samples = (n_factors * n_replicates) + n_factors

    GenerateProject(project_name=project_name, genome=genome_assembly, n_factors=n_factors, n_replicates=n_replicates)

    # first edit the defaul path to the annotation sheet
    config = os.path.join("tests", project_name, "metadata", "project_config.yaml")
    c = yaml.safe_load(open(config, 'r'))
    c['metadata']['sample_annotation'] = os.path.abspath(os.path.join("tests", project_name, "metadata", "annotation.csv"))
    c['metadata']['comparison_table'] = os.path.abspath(os.path.join("tests", project_name, "metadata", "comparison_table.csv"))
    yaml.safe_dump(c, open(config, "w"))

    # project and associated analysis
    prj = Project(config)
    a = Analysis(name=project_name, prj=prj)
    assert a.__repr__() == "Analysis object named '{}' with {} samples of genome '{}'.".format(project_name, n_samples, genome_assembly)
    assert len(prj.samples) == len(a.samples)
    assert all([x == y for x, y in zip(prj.samples, a.samples)])

    shutil.rmtree(os.path.join("tests", project_name))
