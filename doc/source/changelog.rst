Changelog
******************************

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

.. |br| raw:: html

   <br />

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

.. |br| raw:: html

   <br />

- **v0.1.3.5.3b** (*2018-06-12*):

  - Fixes:

    - general.deseq_analysis: fix hyphen character conversion; better contrasts for DESeq2

.. |br| raw:: html

   <br />

- **v0.1.3.5.3** (*2018-05-31*):

  - Fixes:

    - projectmanager: fix Makefile creation
    - ngs_analysis recipe: change selection of samples on "pass_qc"; do differential_overlap only when >1 comparison


.. |br| raw:: html

   <br />

- **v0.1.3.5.2** (*2018-05-30*):

  - Fixes:

    - trackmanager: trackHeight attribute typo making tracks have 128 pixels height
    - trackmanager: sampleGroup attribute indentation

  - New:

    - general.plot_differential: center divergent heatmaps on 0 in, add sorted heatmaps
    - General `ngs_analysis` recipe for general analysis of NGS projects.


.. |br| raw:: html

   <br />

- Major release: **v0.1.3.5** (*2018-05-15*):

  - New:

    - Extended documentation
    - Command-line interface (CLI) based on sub-commands for ``projectmanager``.
    - Recipes: scripts which ``projectmanager`` can run.
    - General `ngs_analysis` recipe for general analysis of NGS projects.


------------


- Upcoming releases:

  - New:

    - `ngs_toolkit.utils` to hold small helper functions.
    - Reconstructing genome static files for various genomes through API or script
    - Wrapper function `annotate_regions` in atacseq.ATACAnalysis to run all region annotation functions. Should get external files if needed

  - Changes:

    - Remove requirement to have ``pipelines`` repository installed in order to extend base Sample objects
    - Decoupling of static files from ``data/external``
    - Maintenance of sample attributes as provided by user by means of reading them in as strings (to be improved further)
    - Improved serialization of Sample objects
    - Better hg38 support.
