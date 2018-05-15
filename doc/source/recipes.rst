Recipes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`ngs_toolkit` provides scripts to perform routine tasks on NGS data - recipes.
To make it convinient to run the scripts on data from a project, recipes can be run with the command ``projectmanager recipe <recipe_name> <project_config.yaml>``.


ngs_analysis
=============================

This recipe will perform general NGS analysis on 3 data types: ATAC-seq, ChIP-seq and RNA-seq.
For ATAC and ChIP-seq, quantification and annotation of genomic regions will be performed.
Standard analysis appropriate for each data type will proceed with cross-sample normalization, unsupervised analysis and supervised analysis if a ``comparison_table`` is provided.

This recipe uses variables provided in the project configuration file ``project_name``, ``sample_attributes`` and ``group_attributes``.
