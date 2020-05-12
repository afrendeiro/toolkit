#!/usr/bin/env python

import os

import pybedtools
from ngs_toolkit.general import lola

import pytest
from .conftest import file_exists_and_not_empty, CI  # , RPY2


class Test_annotate_samples:
    def test_no_arguments(self, analysis_normalized):
        analysis_normalized.annotate_samples()

    def test_matrix_raw(self, atac_analysis):
        atac_analysis.annotate_samples(matrix="matrix_raw")


def test_get_matrix(atac_analysis):
    import numpy as np
    import pandas as pd

    matrix = atac_analysis.get_matrix(matrix="matrix_raw")
    assert np.array_equal(matrix.values, atac_analysis.matrix_raw.values)
    assert (matrix == atac_analysis.matrix_raw).all().all()
    atac_analysis.dummy = atac_analysis.matrix_raw + 1
    matrix = atac_analysis.get_matrix(matrix="dummy")
    assert (matrix == (atac_analysis.matrix_raw + 1)).all().all()

    matrix = atac_analysis.get_matrix(matrix=atac_analysis.matrix_raw)
    assert np.array_equal(matrix.values, atac_analysis.matrix_raw.values)
    assert (matrix == atac_analysis.matrix_raw).all().all()
    atac_analysis.dummy = atac_analysis.matrix_raw + 1
    matrix = atac_analysis.get_matrix(matrix="dummy")
    assert (matrix == (atac_analysis.matrix_raw + 1)).all().all()

    # sample subssetting
    matrix = atac_analysis.get_matrix(
        matrix="matrix_raw", samples=atac_analysis.samples[:2])
    assert (pd.Series([
        s.name
        for s in atac_analysis.samples[:2]]) == matrix.columns).all()


# +++ get_genome_reference
# index_fasta
# twobit_to_fasta
# +++ get_blacklist_annotations
# +++ get_tss_annotations
# +++ get_genomic_context
# +++ get_chromosome_sizes
# +++ deseq_analysis
# least_squares_fit


def test_differential_from_bivariate_fit(analysis_normalized):
    from ngs_toolkit.general import differential_from_bivariate_fit

    with analysis_normalized as an:
        out_dir = os.path.join(an.results_dir, "diff")
        out_prefix = os.path.join(out_dir, "bivariate_fit")
        differential_from_bivariate_fit(
            an.comparison_table, an.matrix_norm,
            out_dir, out_prefix)

        outs = [
            ".deseq_result.all_comparisons.csv",
            ".deseq_result.all_comparisons.scatter.svg",
            ".fit_result.Factor_A_2vs1.csv"]
        for f in outs:
            assert file_exists_and_not_empty(out_prefix + f)


@pytest.mark.skipif(
    CI,
    reason="LOLA testing is not performed on CI.")
class Test_LOLA():
    def test_lola_function(self, tmp_path):
        bed = pybedtools.example_bedtool('hg38-base.bed')
        univ = bed.slop(l=0, r=10, genome='hg38')
        bed_file = bed.fn
        universe_file = univ.fn
        output_folder = os.path.dirname(tmp_path)
        genome = "hg38"

        lola(bed_file, universe_file, output_folder, genome)

        output_files = [
            "allEnrichments.tsv",
            "col_codex.tsv"]
        for file in output_files:
            assert file_exists_and_not_empty(os.path.join(output_folder, file))

    def test_lola_function_multiple_inputs(self, tmp_path):
        import shutil
        bed = pybedtools.example_bedtool('hg38-base.bed')
        univ = bed.slop(l=0, r=10, genome='hg38')
        bed_file = bed.fn
        shutil.copy(bed_file, "A.bed")
        shutil.copy(bed_file, "B.bed")
        universe_file = univ.fn
        output_folder = os.path.dirname(tmp_path)
        genome = "hg38"

        lola(["A.bed", "B.bed"], universe_file, output_folder, genome)

        output_files = [
            "allEnrichments",
            "col_codex"]
        for file in output_files:
            for i in ['A', 'B']:
                assert file_exists_and_not_empty(
                    os.path.join(output_folder, file + i + ".tsv"))


    def test_lola_through_differential_enrichment(
            self, analysis_with_differential):
        with analysis_with_differential as an:
            an.differential_enrichment(steps=['lola'])

        output_files = [
            "allEnrichments.tsv",
            "col_codex.tsv"]

        for file in output_files:
            for direction in ['up', 'down']:
                assert file_exists_and_not_empty(os.path.join(
                    an.results_dir,
                    "differential_analysis_ATAC-seq/enrichments/Factor_A_2vs1."
                    + direction, file))

    def test_lola_through_differential_enrichment_distributed(
            self, analysis_with_differential):
        with analysis_with_differential as an:
            an.differential_enrichment(steps=['lola'], distributed=True)

        output_files = [
            "allEnrichments.tsv",
            "col_codex.tsv"]

        for file in output_files:
            for direction in ['up', 'down']:
                assert file_exists_and_not_empty(os.path.join(
                    an.results_dir,
                    "differential_analysis_ATAC-seq/enrichments/Factor_A_2vs1."
                    + direction, file))

    # def test_lola__plot_differential_enrichment(self):
    #     pass

# meme_ame

@pytest.mark.skipif(
    CI,
    reason="HOMER testing is not performed on CI.")
class TestHomer():
    def test_homer_function(self, tmp_path):
        from ngs_toolkit.general import homer_motifs

        bed = pybedtools.example_bedtool('hg38-base.bed')
        univ = bed.slop(l=0, r=10, genome='hg38')
        bed_file = bed.fn
        universe_file = univ.fn
        output_dir = os.path.dirname(tmp_path)
        genome_assembly = "hg38"

        homer_motifs(bed_file, output_dir, genome_assembly)
        assert os.path.exists(os.path.join(output_dir, "homerMotifs.all.motifs"))


# homer_combine_motifs
# +++ enrichr
# run_enrichment_jobs


def test_project_to_geo(atac_analysis_with_unmapped_input_files):
    from ngs_toolkit.general import project_to_geo

    with atac_analysis_with_unmapped_input_files as an:
        out_dir = os.path.join(an.root_dir, "geo_submission")
        annot = project_to_geo(
            an.prj,
            output_dir=out_dir, steps=['bam', 'peaks'],
            computing_configuration="default")

        cols = [
            'bam_file0', 'bam_file0_md5sum',
            # 'bigwig_file', 'bigwig_file_md5sum',
            'peaks_file', 'peaks_file_md5sum']
        assert all(annot.columns == cols)

        outs = [
            "project_to_geo.{}.sh",
            "{}.bam",
            "{}.bam.md5",
            # "{}.bigWig",
            # "{}.bigWig.md5",
            "{}.peaks.narrowPeak",
            "{}.peaks.narrowPeak.md5"]
        for sample in an.samples:
            for f in outs:
                assert file_exists_and_not_empty(
                    os.path.join(out_dir, f.format(sample.name)))


def test_rename_sample_files(atac_analysis_with_input_files):
    import pandas as pd
    from ngs_toolkit.general import rename_sample_files

    with atac_analysis_with_input_files as an:

        df = pd.DataFrame(
            [['S01_A1', 'S02_A1'], ['SXX_ZZ', 'SYY_ZZ']],
            index=['old_sample_name', 'new_sample_name']).T

        rename_sample_files(df, results_dir=an.data_dir)

        for sample in ['SXX_ZZ', 'SYY_ZZ']:
            outs = [
                os.path.join("mapped", sample + '.trimmed.bowtie2.filtered.bam'),
                os.path.join("mapped", sample + '.trimmed.bowtie2.filtered.bam.bai'),
                os.path.join("peaks", sample + '_peaks.narrowPeak'),
                os.path.join("peaks", sample + '_summits.bed')]
            for f in outs:
                assert file_exists_and_not_empty(os.path.join(an.data_dir, sample, f))


# +++ query_biomart


def test_subtract_principal_component(analysis_normalized):
    from ngs_toolkit.general import subtract_principal_component

    with analysis_normalized as an:

        plot = os.path.join(an.root_dir, "subtract_plot.svg")
        df = subtract_principal_component(
            an.matrix_norm.T, plot_name=plot).T

        assert df.shape == an.matrix_norm.shape
        assert file_exists_and_not_empty(plot)


# fix_batch_effect_limma
