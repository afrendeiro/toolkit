Changelog
******************************

- **v1.0** (*2019-06-01*):

  - New release, with new API.

------------

- **v0.10** (development, pre-1.0):

  - revamp RNASeqAnalysis
  - adapt ChIPSeqAnalysis to new API
  - add ``normalize_by_background`` function to ChIPSeqAnalysis to normalize over background samples
  - fix logger due to accidental deactivation

- **v0.9** (development, pre-1.0):

  - rename annotate to annotate_features and annotate_with_sample_metadata to annotate_samples
  - add annotate_matrix to call both above.

- **v0.8** (development, pre-1.0):

  - usage of the same interpreter running ngs_toolkit to run jobs
  - revamp recipes, usage of recipes for common work functions that run in distributed mode
  - allow import of classes from root of library.

- **v0.7** (development, pre-1.0):

  - implement running of local or distributed jobs using ``divvy``

- **v0.6** (development, pre-1.0):

  - rename merge_table to sample_subannotation

- **v0.5** (development, pre-1.0):

  - major API changes
  - implementing only two types of matrices: matrix_raw, matrix_norm
  - unify normalization methods, each overwrites matrix_norm and sets flag with name of method used

------------

- **v0.2.1** (*2018-12-13*):

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

------------

- **v0.1.6.0** (*2018-12-05*):

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

------------

- **v0.1.5.1** (*2018-11-25*):

  - add config file support for better system-independent operation (specially for enrichment analysis)
  - add logger through "logging" library
  - add batch effect correction with limma
  - add GREAT parser
  - add colouring by p-value for plot_differential
  - add set n. of PCs to calculate to PCA
  - add better colorbars
  - add serial processing of peak commands as option for ChIP-seq peak calling

------------


- **v0.1.4.2** (*2018-10-29*):

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

------------

- **v0.1.4** (*2018-09-25*):

  - Update to peppy version v0.9.1

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

------------

- **v0.1.3.6** (*2018-08-05*):

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

------------

- **v0.1.3.5.3b** (*2018-06-12*):

  - Fixes:

    - general.deseq_analysis: fix hyphen character conversion; better contrasts for DESeq2

------------

- **v0.1.3.5.3** (*2018-05-31*):

  - Fixes:

    - projectmanager: fix Makefile creation
    - ngs_analysis recipe: change selection of samples on "pass_qc"; do differential_overlap only when >1 comparison


------------

- **v0.1.3.5.2** (*2018-05-30*):

  - Fixes:

    - trackmanager: trackHeight attribute typo making tracks have 128 pixels height
    - trackmanager: sampleGroup attribute indentation

  - New:

    - general.plot_differential: center divergent heatmaps on 0 in, add sorted heatmaps
    - General `ngs_analysis` recipe for general analysis of NGS projects.


------------

- Major release: **v0.1.3.5** (*2018-05-15*):

  - New:

    - Extended documentation
    - Command-line interface (CLI) based on sub-commands for ``projectmanager``.
    - Recipes: scripts which ``projectmanager`` can run.
    - General `ngs_analysis` recipe for general analysis of NGS projects.
