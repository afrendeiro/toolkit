#!/usr/bin/env python


import glob
import os
import shutil

import numpy as np
import pytest

from ngs_toolkit.analysis import Analysis
from ngs_toolkit.utils import get_this_file_or_timestamped
from .conftest import file_exists, file_not_empty, COMBAT


class TestAnalysis:
    def test_analysis_representation(self):
        name = "test_analysis"

        an = Analysis(name=name)
        assert an.__repr__() == "Analysis '{}'.".format(name)
        assert "samples" not in an.__repr__()

    def test_with_object_as(self):
        name = "test_analysis"

        an = Analysis(name=name)
        with an as _an:
            assert an is _an
            assert an == _an
            assert _an.__repr__() == "Analysis '{}'.".format(name)
            assert "samples" not in _an.__repr__()

    def test_analysis_creation(self, tmp_path):
        from ngs_toolkit.demo.data_generator import generate_project

        tmp_path = str(tmp_path)

        # Let's make several "reallish" test projects
        project_prefix_name = "test-project"
        data_types = ["ATAC-seq", "RNA-seq", "ChIP-seq"]  # "CNV"
        genome_assemblies = [("human", "hg38"), ("mouse", "mm10")]  # ("human", "hg19"),

        params = {
            "ATAC-seq": {
                "n_factors": [1, 2, 3],
                "n_features": [100, 1000, 10000],
                "n_replicates": [1, 2, 5],
                "analysis": "ATACSeqAnalysis",
            },
            "ChIP-seq": {
                "n_factors": [1, 2, 3],
                "n_features": [100, 1000, 10000],
                "n_replicates": [1, 2, 5],
                "analysis": "ChIPSeqAnalysis",
            },
            "RNA-seq": {
                "n_factors": [1, 2, 3],
                "n_features": [100, 1000, 25000],
                "n_replicates": [1, 2, 5],
                "analysis": "RNASeqAnalysis",
            },
        }

        for data_type in data_types:
            n_factors = params[data_type]["n_factors"][0]
            n_features = params[data_type]["n_features"][0]
            n_replicates = params[data_type]["n_replicates"][0]
            for organism, genome_assembly in genome_assemblies:

                project_name = "{}_{}_{}_{}_{}_{}".format(
                    project_prefix_name,
                    data_type,
                    genome_assembly,
                    n_factors,
                    n_features,
                    n_replicates,
                )

                an = generate_project(
                    output_dir=tmp_path,
                    project_name=project_name,
                    organism=organism,
                    genome_assembly=genome_assembly,
                    data_type=data_type,
                    n_factors=n_factors,
                    n_replicates=n_replicates,
                    n_features=n_features,
                    only_metadata=True,
                )
                # n_samples = (n_factors * n_replicates) + n_factors
                # assert an.__repr__() == (
                #     "'{}' analysis '{}' with {} samples of organism '{}' ({}).".format(
                #         data_type, project_name, n_samples, organism, genome_assembly
                #     )
                # )
                # assert len(n_factors * 2 * n_replicates * 2) == len(an.prj.samples) == len(an.samples)
                assert all([x == y for x, y in zip(an.prj.samples, an.samples)])

                shutil.rmtree(tmp_path)

    def test_analysis_serialization(self, tmp_path):

        tmp_path = str(tmp_path)

        pickle_file = os.path.join(tmp_path, "analysis.pickle")
        a = Analysis(pickle_file=pickle_file)
        assert not file_exists(pickle_file)
        a.to_pickle()
        assert file_exists(pickle_file)
        assert file_not_empty(pickle_file)

        previous_size = os.stat(get_this_file_or_timestamped(pickle_file)).st_size
        a.random = np.random.random((100, 100))
        a.to_pickle()
        new_size = os.stat(get_this_file_or_timestamped(pickle_file)).st_size
        assert new_size > previous_size

        previous_size = os.stat(get_this_file_or_timestamped(pickle_file)).st_size
        a.random = np.random.random((100, 100))
        a.to_pickle(timestamp=True)
        assert len(glob.glob(os.path.join(tmp_path, "*.pickle"))) == 2

    def test_analysis_loading(self, tmp_path):
        tmp_path = str(tmp_path)
        pickle_file = os.path.join(tmp_path, "pickle")
        secret = "I've existed before"

        a = Analysis()
        a.pickle_file = pickle_file
        a.secret = secret
        a.to_pickle()

        a2 = Analysis(from_pickle=pickle_file)
        assert a2.secret == secret

        a3 = Analysis()
        a3.update(pickle_file)
        assert a3.secret == secret

        a4 = Analysis()
        a4.pickle_file = pickle_file
        a4 = a4.from_pickle()
        assert a4.secret == secret

        shutil.rmtree(tmp_path)

    def test__overwride_sample_representation(self, atac_analysis):

        prev = atac_analysis.samples[0].__repr__
        Analysis._overwride_sample_representation()
        new = atac_analysis.samples[0].__repr__

        assert prev != new

    def test__check_data_type_is_supported(self):
        assert Analysis._check_data_type_is_supported("ATAC-seq")
        assert Analysis._check_data_type_is_supported("ChIP-seq")
        assert Analysis._check_data_type_is_supported("RNA-seq")
        assert Analysis._check_data_type_is_supported("CNV")
        assert not Analysis._check_data_type_is_supported("Microarray")

    def test__get_data_type(self, atac_analysis):
        assert atac_analysis._get_data_type() == "ATAC-seq"
        assert atac_analysis._get_data_type(data_type="ATAC-seq") == "ATAC-seq"
        assert atac_analysis._get_data_type(data_type="RNA-seq") == "RNA-seq"

        with pytest.raises(ValueError):
            atac_analysis._get_data_type(data_type="Microarray")

        atac_analysis.data_type = None
        with pytest.raises(ValueError):
            atac_analysis._get_data_type()

        del atac_analysis.data_type
        with pytest.raises(AttributeError):
            atac_analysis._get_data_type()

    def test__check_samples_have_file(self, atac_analysis):
        with pytest.raises(AttributeError):
            atac_analysis._check_samples_have_file("NOTEXISTING")

        # assert not atac_analysis._check_samples_have_file("summits")

        assert not atac_analysis._check_samples_have_file("sample_name")

    def test__get_samples_have_file(self, atac_analysis):
        assert not atac_analysis._get_samples_have_file("sample_name")

    def test__get_samples_missing_file(self, atac_analysis):
        with pytest.raises(AttributeError):
            atac_analysis._get_samples_have_file("NOTEXISTING")

        assert not atac_analysis._get_samples_have_file("sample_name")

        assert not atac_analysis._get_samples_have_file("aligned_filtered_bam")

    def test__get_samples_with_input_file(self, atac_analysis):
        with pytest.raises(AttributeError):
            atac_analysis._get_samples_with_input_file("NOTEXISTING")

        with pytest.raises(IOError):
            atac_analysis._get_samples_with_input_file("sample_name")

        with pytest.raises(IOError):
            atac_analysis._get_samples_with_input_file("aligned_filtered_bam")

        assert not atac_analysis._get_samples_with_input_file(
            "aligned_filtered_bam", permissive=True
        )
        assert not atac_analysis._get_samples_with_input_file("peaks", permissive=True)
        assert not atac_analysis._get_samples_with_input_file("summits", permissive=True)

    @pytest.mark.parametrize(
        "env_var,string",
        [
            ("_${USER}_", "_{}_".format(os.environ.get("USER"))),
            # ("_$PATH_", "_{}_".format(os.environ.get("PATH"))),
        ],
    )
    def test__format_string_with_environment_variables(self, env_var, string):
        assert string == Analysis._format_string_with_environment_variables(env_var)

    def test__format_string_with_attributes_simple(self):
        t = Analysis()
        t.a = 1
        t.b = ""
        assert "1" == Analysis._format_string_with_attributes(t, "{a}{b}")

    @pytest.mark.parametrize(
        "env_var,string",
        [("{data_type}", "ATAC-seq"), ("{name}", "test-project_ATAC-seq_human_hg38_1_250_2"),],
    )
    def test__format_string_with_attributes(self, atac_analysis, env_var, string):
        assert string == atac_analysis._format_string_with_attributes(env_var)

    def test_record_output_file(self, atac_analysis):
        assert hasattr(atac_analysis, "output_files")
        assert len(atac_analysis.output_files) == 0
        atac_analysis.record_output_file("a", name="analysis")
        assert hasattr(atac_analysis, "output_files")
        assert len(atac_analysis.output_files) == 1
        assert atac_analysis.output_files[0][0] == "analysis"
        assert atac_analysis.output_files[0][1] == "a"


def test_project_with_subprojects(subproject_config):
    from ngs_toolkit import Analysis

    a = Analysis(from_pep=subproject_config)
    assert len(a.samples) == 0

    a = Analysis(from_pep=subproject_config, amendments=["test_subproject"])
    assert len(a.samples) > 0


@pytest.mark.skipif(not COMBAT, reason="Combat not installed")
def test_remove_factor(atac_analysis_many_factors):
    import pandas as pd

    a = atac_analysis_many_factors
    a.matrix_norm = a.matrix_norm.dropna()

    prefix = os.path.join(a.results_dir, "unsupervised_analysis_{}".format(a.data_type), a.name)
    # inspect
    a.unsupervised_analysis(output_prefix="before", steps=["pca_association"])

    f = prefix + ".before.pca.variable_principle_components_association.csv"
    p = pd.read_csv(get_this_file_or_timestamped(f))

    # extract the name of the factor with highest contribution
    factor = p.iloc[p.query("pc == 1")["p_value"].idxmin()]["attribute"]
    # check if it's significant
    assert p.query("attribute == '{}' and pc < 15".format(factor))["p_value"].min() < 0.05

    # remove factor without regard for the other factors
    m = a.remove_factor_from_matrix(factor=factor, assign=False, save=False)
    a.unsupervised_analysis(matrix=m, output_prefix="after_simple", steps=["pca_association"])

    f = prefix + ".after_simple.pca.variable_principle_components_association.csv"
    p2 = pd.read_csv(get_this_file_or_timestamped(f))
    assert p2.query("attribute == '{}' and pc < 15".format(factor))["p_value"].min() > 0.05

    # remove factor accounting for the other factors
    m = a.remove_factor_from_matrix(
        factor=factor,
        covariates=[x for x in a.group_attributes if x != factor],
        assign=False,
        save=False,
    )
    a.unsupervised_analysis(matrix=m, output_prefix="after_covariates", steps=["pca_association"])

    f = prefix + ".after_covariates.pca.variable_principle_components_association.csv"
    p3 = pd.read_csv(get_this_file_or_timestamped(f))
    assert p3.query("attribute == '{}' and pc < 15".format(factor))["p_value"].min() > 0.05
