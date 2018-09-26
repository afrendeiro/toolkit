Recipes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`ngs_toolkit` provides scripts to perform routine tasks on NGS data - they are called recipes.

Recipes are distributed with ngs_toolkit and can be seen in the `github repository <https://github.com/afrendeiro/toolkit/tree/master/ngs_toolkit/recipes>`_.

To make it convenient to run the scripts on data from a project, recipes can be run with the command ``projectmanager recipe <recipe_name> <project_config.yaml>``.


ngs_analysis
=============================

This recipe will perform general NGS analysis on 3 data types: ATAC-seq, ChIP-seq and RNA-seq.
For ATAC and ChIP-seq, quantification and annotation of genomic regions will be performed.
Standard analysis appropriate for each data type will proceed with cross-sample normalization, unsupervised analysis and supervised analysis if a ``comparison_table`` is provided.

This recipe uses variables provided in the project configuration file ``project_name``, ``sample_attributes`` and ``group_attributes``.

Here are the command-line arguments to use it in a stand-alone script:

.. code-block:: none

	usage: python ngs_analysis_recipe [-h] [-n NAME] [-o RESULTS_DIR]
	                           [-t {ATAC-seq,RNA-seq,ChIP-seq}] [-q] [-a ALPHA]
	                           [-f ABS_FOLD_CHANGE]
	                           config_file

	positional arguments:
	  config_file           YAML project configuration file.

	optional arguments:
	  -h, --help            show this help message and exit
	  -n NAME, --analysis-name NAME
	                        Name of analysis. Will be the prefix of output_files.
	                        By default it will be the name of the Project given in
	                        the YAML configuration.
	  -o RESULTS_DIR, --results-output RESULTS_DIR
	                        Directory for analysis output files. Default is
	                        'results' under the project roort directory.
	  -t {ATAC-seq,RNA-seq,ChIP-seq}, --data-type {ATAC-seq,RNA-seq,ChIP-seq}
	                        Data type to restrict analysis to. Default is to run
	                        separate analysis for each data type.
	  -q, --pass-qc         Whether only samples with a 'pass_qc' value of '1' in
	                        the annotation sheet should be used.
	  -a ALPHA, --alpha ALPHA
	                        Alpha value of confidence for supervised analysis.
	  -f ABS_FOLD_CHANGE, --fold-change ABS_FOLD_CHANGE
	                        Absolute log2 fold change value for supervised
	                        analysis.


call_peaks
=============================

This recipe will call peaks for samples in a fashion described in a `comparison table <https://ngs-toolkit.readthedocs.io/en/latest/comparison_table.html>`_.

It is capable of parallelizing work in jobs if a SLURM cluster is available. 

Here are the command-line arguments to use it in a stand-alone script:

.. code-block:: none

	usage: call_peaks [-h] [-c COMPARISON_TABLE] [-t] [-j]
	                           [-o RESULTS_DIR]
	                           config_file

	Call peaks recipe.

	positional arguments:
	  config_file           YAML project configuration file.

	optional arguments:
	  -h, --help            show this help message and exit
	  -c COMPARISON_TABLE, --comparison-table COMPARISON_TABLE
	                        Comparison table to use for peak calling. If not
	                        provided will use a filenamed `comparison_table.csv`
	                        in the same directory of the given YAML Project
	                        configuration file.
	  -t, --only-toggle     Whether only comparisons with 'toggle' value of '1' in
	                        the should be performed.
	  -j, --as-jobs         Whether jobs should be created for each sample, or it
	                        should run in serial mode.
	  -o RESULTS_DIR, --results-output RESULTS_DIR
	                        Directory for analysis output files. Default is
	                        'results' under the project roort directory.



region_set_frip
=============================

This recipe will perform fraction of reads in peaks (FRiP) for ATAC-seq or ChIP-seq samples based on a set of regions discovered across all samples in a given project or in an external gold region set.

If the external region set is not given, a region set derived from all samples already exists (e.g. from running the ngs_analysis recipe) the same one will be used, otherwise it will be produced.

Here are the command-line arguments to use it in a stand-alone script:

.. code-block:: none

	usage: region_set_frip [-h] [-n NAME] [-r REGION_SET] [-q] [-j]
	                           [-o RESULTS_DIR]
	                           config_file

	Region set FRiP recipe.

	positional arguments:
	  config_file           YAML project configuration file.

	optional arguments:
	  -h, --help            show this help message and exit
	  -n NAME, --analysis-name NAME
	                        Name of analysis. Will be the prefix of output_files.
	                        By default it will be the name of the Project given in
	                        the YAML configuration.
	  -r REGION_SET, --region-set REGION_SET
	                        BED file with region set derived from several samples
	                        or Oracle region set. If unset, will try to get the
	                        `sites` attribute of an existing analysis object if
	                        existing, otherwise will create a region set from the
	                        peaks of all samples.
	  -q, --pass-qc         Whether only samples with a 'pass_qc' value of '1' in
	                        the annotation sheet should be used.
	  -j, --as-jobs         Whether jobs should be created for each sample, or it
	                        should run in serial mode.
	  -o RESULTS_DIR, --results-output RESULTS_DIR
	                        Directory for analysis output files. Default is
	                        'results' under the project roort directory.


merge_signal
=============================

This recipe will merge signal from various ATAC-seq or ChIP-seq samples given a set of attributes to group samples by.

It produces merged BAM and bigWig files for all signal in the samples but is also capable of producing this for nucleosomal/nucleosomal free signal based on fragment length distribution if data is paired-end sequenced. This signal may optionally be normalized for each group. It is also capable of parallelizing work in jobs if a SLURM cluster is available.

Here are the command-line arguments to use it in a stand-alone script:

.. code-block:: none

	usage: merge_signal [-h] [-a ATTRIBUTES] [-q] [-j] [-n] [--nucleosome]
	                    [--overwrite] [-o OUTPUT_DIR]
	                    config_file

	Merge signal recipe.

	positional arguments:
	  config_file           YAML project configuration file.

	optional arguments:
	  -h, --help            show this help message and exit
	  -a ATTRIBUTES, --attributes ATTRIBUTES
	                        Attributes to merge samples by. By default will use
	                        values in the project config `sample_attributes`.
	  -q, --pass-qc         Whether only samples with a 'pass_qc' value of '1' in
	                        the annotation sheet should be used.
	  -j, --as-jobs         Whether jobs should be created for each sample, or it
	                        should run in serial mode.
	  -n, --normalize       Whether tracks should be normalized to total sequenced
	                        depth.
	  --nucleosome          Whether to produce nucleosome/nucleosome-free signal
	                        files.
	  --overwrite           Whether to overwrite existing files.
	  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
	                        Directory for output files. Default is 'merged' under
	                        the project roort directory.
