=============================
Changelog
=============================

All notable changes to this project will be documented in this file.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.


[Unreleased]
*****************************

Added
-----------------------------
  - New CRISPR module for the analysis of pooled CRISPR screens
  - Generation of input files for RNA-seq and CNV data types
  - Testing of RNA-seq and CNV modules specific functionality
  - Testing of recipes

Changed
-----------------------------
  - More simplicity and abstraction for functions in main :class:`ngs_toolkit.analysis.Analysis` class.


[0.22.0] - 2020-02-07
*****************************

Changed
-----------------------------

  - Minor changes to allow pandas==1.0.0


[0.21.1] - 2020-02-06
*****************************

Changed
-----------------------------

  - Fix bug preventing plotting of all parts in :func:`ngs_toolkit.analysis.Analysis.plot_differential_enrichment`.


[0.21.0] - 2020-02-06
*****************************

Changed
-----------------------------

  - Remove rpy2 from mandatory requirements
  - Documentation updates
  - Compatibility with new naming of Enrichr gene set libraries


[0.20.0] - 2019-11-29
*****************************

Added
-----------------------------

  - Add :func:`ngs_toolkit.general.get_chromosome_sizes` to get the size of the chromosomes of a genome assembly


Changed
-----------------------------

  - Change default genome assembly of human to hg38/GRCh38
  - Improvements to several recipes
  - Remove cffi from requirements.


[0.19.3] - 2019-10-18
*****************************

Changed
-----------------------------

  - Pin cffi version to fix `bug in rpy2 with specific cffi version <https://bitbucket.org/rpy2/rpy2/issues/591/>`_.


[0.19.2] - 2019-10-13
*****************************

Added
-----------------------------
  
  - New module :class:`ngs_toolkit.demo` which generates random data as PEP-formatted projects
  - New :class:`ngs_toolkit.recipes.generate_project` recipe to generate a new project using CLI
  - New normalization method: variance stabilizing transformation (VST) available
  - Add function to run :class:`ngs_toolkit.recipes.ngs_analysis` recipe from initalized analysis object
  - Add ``distributed`` mode to :func:`ngs_toolkit.atacseq.ATACSeqAnalysis.measure_coverage` using the new coverage recipe
  - New :class:`ngs_toolkit.recipes.coverage` recipe for calculating the coverage of a BAM file in BED regions
  - Docs: usage of ``sphinx-issues`` for connecting to issue tracking and ``sphinx-argparse`` for automatic documentation of CLI of recipes


Changed
-----------------------------

  - Generator of random data now based on proper negative-binomial model
  - Test suite now uses :class:`ngs_toolkit.demo` module.
  - change default of :func:`ngs_toolkit.analysis.Analysis.differential_analysis` to ``filter_support=False``.
  - fix boundary passing bug in subtract_principal_component
  - Adopt `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_ changelog style.


[0.19.1] - 2019-10-09
*****************************

Added
-----------------------------

- Add ``save`` and ``assign`` arguments to :func:`ngs_toolkit.atacseq.ATACSeqAnalysis.get_consensus_sites`.
- New :class:`ngs_toolkit.recipes.coverage`

Changed
-----------------------------
 - Stopped special handling reading of ``matrix_norm`` in :func:`ngs_toolkit.analysis.Analysis.load_data`. This now assumes a non-MultiIndex dataframe; fix :issue:`59`.
 - Change default of :func:`ngs_toolkit.analysis.Analysis.annotate_samples` ``save`` and ``assign`` arguments to :obj:`False`.
 - :func:`ngs_toolkit.analysis.Analysis.remove_factor_from_matrix` now drops features with no variance.
 - Change ``filter_mito_chr`` to ``filter_chroms`` argument of :func:`ngs_toolkit.atacseq.ATACSeqAnalysis.get_consensus_sites` in order to allow filtering arbitrary chromosomes out from consensus sites using a regex pattern. Supports multiple patterns by using the "|" operator.
 - Much more efficient algorithm underlying :func:`ngs_toolkit.atacseq.ATACSeqAnalysis.get_consensus_sites` with speedup times up to 20x.
 - Change method to compute coverage for :func:`ngs_toolkit.atacseq.ATACSeqAnalysis.measure_coverage` with ``distributed=True`` from ``bedtools coverage`` to :class:`ngs_toolkit.recipes.coverage`. This has the following advantages: less reliance on bedtools; outputing a result for each region; same function as in serial mode.
 - :func:`ngs_toolkit.utils.count_reads_in_intervals` now outputs coverage for every chromosome, outputs number of errors to log.
 - Fix bug :issue:`61` - missing time parameter for enrichr job.
 - Pin ``yacman`` version to 0.6.0.

[0.18.1] - 2019-10-03
*****************************

Added
-----------------------------

  - Add ``overwrite`` argument to :func:`ngs_toolkit.analysis.Analysis.measure_coverage`.

Changed
-----------------------------

  - Fix :issue:`60`: building of confusion matrix for Fisher Test in :func:`ngs_toolkit.analysis.Analysis.differential_overlap`.
  - Use ``-sorted`` argument to ``bedtools coverage`` in :func:`ngs_toolkit.analysis.Analysis.measure_coverage` for fast and low-memory computing of coverage of BAM file in BED file when ``distributed=True``.
  - Set ``setuptools_scm``, ``requests``, ``rpy2`` and ``natsort`` versions.
  - Extensive updates to documentation


[0.17.6] - 2019-09-25
*****************************

Added
-----------------------------
  - Using ``setuptools_scm`` for semantic versioning.
  - More testing of DESeq2 functionality.

Changed
-----------------------------
  - Removed ``_version.py`` file due to ``setuptools_scm`` adoption. API does not change though.
  - Fixed continuous deployment in Travis.
  - Dockerfile


[0.17.3] - 2019-09-24
*****************************

Changed
-----------------------------
  - Fixed continuous deployment in Travis.

[0.17.2] - 2019-09-23
*****************************

Changed
-----------------------------
  - Always display ``ngs_toolkit`` version in html report even if ``pip_versions=False``.
  - Pretty README on PyPI by specifying ``long_description_content_type="text/markdown"`` on setup.py.


[0.17.1] - 2019-09-23
*****************************

Added
-----------------------------

 - Continuous deployment through Travis.
 - Gitpod configuration
 - Functionality to test whether ``bedtools`` version is acceptable.
 - :func:`ngs_toolkit.analysis.Analysis.get_sample_annotation` for convinient getting of a current sample annotation based on `sample_attributes` and `group_attributes` given to ``ngs_toolkit``.
 - Add ``deseq_kwargs`` argument to :func:`ngs_toolkit.analysis.Analysis.differential_analysis` to allow passing of arguments to DESeq2 main function.
 - Add functionality to repeat API call to ``BioMart`` in :func:`ngs_toolkit.general.query_biomart` with ``max_api_retries`` argument.

Changed
-----------------------------

  - Switched from custom install of non-Python requirements to Debian packages
  - Required bedtools version is now 2.17.1
  - Abstraction of :func:`ngs_toolkit.decorators.check_organism_genome` and :func:`ngs_toolkit.decorators.check_has_sites` to :func:`ngs_toolkit.decorators.check_has_attributes` which now accepts arguments.
  - Add ``save``, ``assign`` and ``output_prefix`` to :func:`ngs_toolkit.analysis.Analysis.get_matrix_stats`, :func:`ngs_toolkit.analysis.Analysis.annotate_features`, :func:`ngs_toolkit.atacseq.ATACSeqAnalysis.get_peak_gene_annotation` :func:`ngs_toolkit.atacseq.ATACSeqAnalysis.get_peak_genomic_location`, :func:`ngs_toolkit.atacseq.ATACSeqAnalysis.get_peak_chromatin_state`
  - Set default arguments of :func:`ngs_toolkit.analysis.Analysis.annotate_samples` to :obj:`False`.
  - Revamp of :func:`ngs_toolkit.atacseq.ATACSeqAnalysis.plot_peak_characteristics` with specific tests, but backward compatible API call.
  - Avoid ``matplotlib`` version 3.1.1 due to bug manifesting on seaborn. Requirement now set to matplotlib>=2.1.1,<3.1.1.


[0.16.1] - 2019-09-04
*****************************

Changed
-----------------------------
  - Fix bug in setup.py

[0.16] - 2019-09-04
*****************************

Added
-----------------------------

  - Dockerfile

Changed
-----------------------------

  - Fixed bug in general.get_genomic_context due to a bug in the timestamping workflow
  - Various fixes of timestamping and html reporting functionality
  - Distributing tests with library for more portable testing
  - Move rpy2 requirement to mandatory
  - Make test data cases smaller for faster CI

[0.14] - 2019-08-28
*****************************

Added
-----------------------------

  - Add recording of analysis outputs under Analysis.output_files
  - Add timestamping of table and figure Analysis outputs
  - Add HTML report with continuous generation
  - :func:`ngs_toolkit.analysis.Analysis.remove_factor_from_matrix` for Combat removal of batch effects
  - More documentation on installation particularly for setting up in a Conda environment

Changed
-----------------------------

  - Now testing on Ubuntu 18.04 for Python 3.6 and 3.7 only.
  - CNV module updated
  - recipe updated

[0.12] - 2019-08-12
*****************************

Changed
-----------------------------

  - change of unsupervised_analysis API call: homogeneization with remaining functions
  - optional saving of embeddings and loadings of PCA and manifolds in unsupervised_analysis

[0.11] - 2019-08-08
*****************************

Added
-----------------------------

  - support for additional keyword arguments passed to Project initialization when using `from_pep`

Changed
-----------------------------

  - adapt to latest versions of pepkit stack
  - better colouring of sample group levels in get_level_colors

[0.10]
*****************************

Added
-----------------------------

  - ``normalize_by_background`` function to ChIPSeqAnalysis to normalize over background samples


Changed
-----------------------------

  - revamp RNASeqAnalysis
  - adapt ChIPSeqAnalysis to new API
  - fix logger due to accidental deactivation

[0.9]
*****************************

Added
-----------------------------

  - add annotate_matrix to call both above.


Changed
-----------------------------

  - rename annotate to annotate_features and annotate_with_sample_metadata to annotate_samples

[0.8]
*****************************

Changed
-----------------------------

  - usage of the same interpreter running ngs_toolkit to run jobs
  - revamp recipes, usage of recipes for common work functions that run in distributed mode
  - allow import of classes from root of library.

[0.7]
*****************************

Added
-----------------------------

  - implement running of local or distributed jobs using ``divvy``

[0.6]
*****************************

Changed
-----------------------------

  - rename merge_table to sample_subannotation

[0.5]
*****************************

Changed
-----------------------------

  - major API changes
  - implementing only two types of matrices: matrix_raw, matrix_norm
  - unify normalization methods, each overwrites matrix_norm and sets flag with name of method used

[0.2.1] - 2018-12-13
*****************************

  - minor:

    - change default directory for enrichment results
    - add class method to overwrite Sample object representation
    - add configuration to merge_signal recipe
    - add graphics functions
    - add optional requirements for single cell analysis
    - add possibility of different prefixes when collecting enrichments
    - remove requirement of some comparison_table and attributes_to_plot arguments
    - remove obsolete functions
    - more powerful Analysis objects by leveraging on known Project attributes
    - simplify plot of number of differential regions per comparison in plot_differential

  - bug fixes:

    - fix pipy install on Python 3: requirements.txt is now distributed with package
    - update merge_signal recipe - fix bug when grouping samples by only one attribute
    - better error catching
    - fix LOLA output collection saving when running in serial mode
    - fix choice of common p-value color overlay to plot in plot_differential
    - fix creating job in merge_signal recipe
    - fix invalid yaml in configs
    - fix mistake in requirements for peppy
    - fix some security issues

[0.1.6.0] - 2018-12-05
*****************************

  - major:

    - New CNV module
    - many fixes and improvements to run various enrichment analysis in serial
    - add specific attributes to classes - this will be the basis of the new API revamp
    - add support for running specific steps of enrichment analysis
    - better utf8 encoding support to all Enrichr inputs/outputs
    - add support for plotting 1 attribute in unsupervised_analysis
    - add support for limma regression without covariates; more help messages
    - fix bug in plot_differential when plotting scatter with colours per p-value
    - improved general.query_biomart to handle fields with multiple values  
    - update requirements

  - minor:

    - now plotting MA, scatter and volcano plots even if there are no significant variables
    - plot log variance in PCA
    - better docstring styling (in progress)
    - plotting signed p-value heatmap
    - support case when only one feature is differential
    - add option to turn on independnent filtering in DESeq2
    - add y log scale to p-value and fold-change distplots
    - homogeneize range of p-value colouring of scatter, volcano and MA plots across comparisons - new colormap
    - better handling of missing comparisons in general.plot_differential
    - better plotting of plot_differential p-values
    - fix example config to correct paths
    - add verbosity to manager programs
    - reporting more info for plot_differential

[0.1.5.1] - 2018-11-25
*****************************

  - add config file support for better system-independent operation (specially for enrichment analysis)
  - add logger through "logging" library
  - add batch effect correction with limma
  - add GREAT parser
  - add colouring by p-value for plot_differential
  - add set n. of PCs to calculate to PCA
  - add better colorbars
  - add serial processing of peak commands as option for ChIP-seq peak calling


[0.1.4.2] - 2018-10-29
*****************************

  - fix important lack of ngs_toolkit.recipes module in setup.py: recipes should now be correctly added to $PATH
  - fix and add full support to comparison_table in recipes.ngs_analysis
  - add region_set_frip recipe
  - add merge_signal recipe
  - add PEP badge

  - ngs_toolkit.general:

    - fix when general.collect_differential_enrichment reads an empty motif enrichment file
    - delete existing files if existing in general.homer_combine_motifs
    - report unmatched differnetial and total features in general.plot_differential
    - general.collect_differential_analysis now sets index of differential_results dataframe correctly
    - add more manifold methods to general.unsupervised_analysis: Isomap, LocallyLinearEmbedding, SpectralEmbedding, TSNE in addition to MDS (and PCA)
    - add ChIP-seq as a valid data type to general.unsupervised_analysis
    - add standardization to parameters of general.unsupervised_analysis
    - add level labels to group labeling of general.unsupervised_analysis and general.plot_differential
    - new default color palletes
    - add option of matching motifs only to known vertebrate in general.homer_consensus
    - add heatmap plotting of enrichment over background for homer consensus in plot_differential_enrichment
    - change default output_dir of general.unsupervised_analysis
    - add more descriptive labels to tqdm loops;
    - add CPU/mem parametrization of general.differential_analysis when running in job mode
    - quick fix for pypiper.ngstk >= 0.6 compatibility (tabs vs spaces) in general.differential_analysis - needs revision
    - resolve pandas warnings of setting without .loc

  - ngs_toolkit.chipseq:

    - add function to filter_peaks
    - add more descriptive labels to tqdm loops;
    - fix overaping peaks calling job files in chipseq.summarize_peaks_from_comparisons

  - ngs_toolkit.atacseq:

    - add more descriptive labels to tqdm loops;

[0.1.4] - 2018-09-25
*****************************

  - Update to peppy version 0.9.1

  - fixes/improvements:

    - add fold enrichment vs p-value plots in homer_consensus plot_differential_enrichments()
    - add index name to DESeq2 CSV output; fix import on homer_combine_motifs
    - add recipes to command-line entry in setup.py; remove cPickle import; better style
    - add scatterplots to enrichr type of data in plot_differential_enrichment
    - add self.data_type to Analysis objects
    - added "homer_consensus" as one type of possible enrichment in plot_differential_enrichment, related to `issue 21 <https://github.com/afrendeiro/toolkit/issues/21>`_
    - crunch landscape bad score for __init__;
    - default color range of numeric values in get_level_colors to min-max
    - fix `issue 19 <https://github.com/afrendeiro/toolkit/issues/19>`_
    - fix `issue 24 <https://github.com/afrendeiro/toolkit/issues/24>`_; general.project_to_geo file referencing
    - implement homer consensus motifs as in `issue 21 <https://github.com/afrendeiro/toolkit/issues/21>`_; add collectiong and plotting of homer enrichmnts
    - moved annotate_with_sample_metadata to base Analysis class
    - tested qnorm implementations; switched to Python implementation, close `issue 12 <https://github.com/afrendeiro/toolkit/issues/12>`_

  - documentation:

    - docs for the region_set_frip, merge_signal and call_peaks recipes

[0.1.3.6] - 2018-08-05
*****************************

  - add two new recipes:

    - region_set_frip: calculate FRiP in a consensus region set of interest for all samples (rsFRiP)
    - merge_signal: create merged signal data for samples sharing specific attributes. Creates BAM files, bigWig files, and BAM files for nucleosomal and nucleosomal-free reads based on fragment size

  - trackmanager:

    - Fix issue #16: trackmanager output indentation
    - add default attributes to specified in project_config.group_attributes or otherwise to ['sample_name']
    - fix empty subGroups in UCSC trackDb file
    - remove required attributes if no value is found

  - Fix issue #20: len(attributes_to_plot) in general.unsupervised_analysis can be 1 now
  - add Makefile to upload to Pypi
  - update looper template folder of projectmanager
  - add default time to longq in analysis_job task in projectmanager Makefile
  - add unbuferred output to ngs_analysis recipe
  - add atacseq.get_gene_level_changes
  - improve atacseq.get_gene_level_accessibility
  - add 2D support to general.signed_mean

[0.1.3.5.3b] - 2018-06-12
*****************************

  - Fixes:

    - general.deseq_analysis: fix hyphen character conversion; better contrasts for DESeq2

[0.1.3.5.3] - 2018-05-31
*****************************

  - Fixes:

    - projectmanager: fix Makefile creation
    - ngs_analysis recipe: change selection of samples on "pass_qc"; do differential_overlap only when >1 comparison


[0.1.3.5.2] - 2018-05-30
*****************************

  - Fixes:

    - trackmanager: trackHeight attribute typo making tracks have 128 pixels height
    - trackmanager: sampleGroup attribute indentation

  - New:

    - general.plot_differential: center divergent heatmaps on 0 in, add sorted heatmaps
    - General `ngs_analysis` recipe for general analysis of NGS projects.


[0.1.3.5] - 2018-05-15
*****************************

  - New:

    - Extended documentation
    - Command-line interface (CLI) based on sub-commands for ``projectmanager``.
    - Recipes: scripts which ``projectmanager`` can run.
    - General `ngs_analysis` recipe for general analysis of NGS projects.
