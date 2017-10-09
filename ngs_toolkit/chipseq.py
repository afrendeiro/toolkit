#!/usr/bin/env python


import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pybedtools
import seaborn as sns

from ngs_toolkit.atacseq import ATACSeqAnalysis


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

    def call_peaks_from_comparisons(self, comparison_table, output_dir="{results_dir}/chipseq_peaks", permissive=True):
        """
        Call peaks for ChIP-seq samples using an annotation of which samples belong in each comparison and which
        samples represent signal or background.

        :param pd.DataFrame comparison_table: Comparison table with the following required columns:
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
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # For each comparison
        for comparison in comparison_table['comparison_name'].drop_duplicates().sort_values():
            # If there aren't two sides to each comparison, skip it or throw error
            if len(set(comparison_table[
                (comparison_table["comparison_name"] == comparison)
            ]["comparison_side"].tolist())) != 2:
                error = "Comparison '{}' does not contain two sides.".format(comparison)
                if permissive:
                    print(error)
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
                (comparison_table["comparison_side"] == 0)
            ]["sample_name"].tolist()

            # Now get the actual samles
            signal_samples = [s for s in self.samples if s.name in pos_names]
            control_samples = [s for s in self.samples if s.name in neg_names]

            print("Doing comparison '{}' with positive samples '{}' and background samples '{}'".format(
                comparison, [s.name for s in signal_samples], [s.name for s in control_samples]
            ))
            # Call peaks
            call_chipseq_peak_job(
                signal_samples, control_samples, output_dir=output_dir, name=comparison)


def call_chipseq_peak_job(signal_samples, control_samples, output_dir, name):
    """
    Call ChIP-seq peaks with MACS2 in a slurm job.

    :param list signal_samples: Signal Sample objects.
    :param list control_samples: Background Sample objects.
    :param list output_dir: Parent directory where MACS2 outputs will be stored.
    :param str name: Name of the MACS2 comparison being performed.
    """
    from pypiper.ngstk import NGSTk
    import textwrap

    tk = NGSTk()

    output_path = os.path.join(output_dir, name)

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    job_name = "macs2_%s" % name

    # Build job script
    # slurm header
    cmd = tk.slurm_header(
        job_name,
        os.path.join(output_path, name + ".log"),
        cpus_per_task=4)

    # load macs2
    cmd += """
\t\t/home/arendeiro/.local/bin/macs2 callpeak -t {0} -c {1} -n {2} --outdir {3}
""".format(" ".join([s.mapped for s in signal_samples]), " ".join([s.mapped for s in control_samples]), name, output_path)

    # Slurm footer
    cmd += "\t\t" + tk.slurm_footer() + "\n"

    # Write job to file
    job_file = os.path.join(output_path, name + ".sh")
    with open(job_file, "w") as handle:
        handle.write(textwrap.dedent(cmd))

    # Submit
    tk.slurm_submit_job(job_file)
