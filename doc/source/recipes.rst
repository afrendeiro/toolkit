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

.. code-block:: bash

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

