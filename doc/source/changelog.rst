Changelog
******************************

- **v0.1.3.5.4** (*2018-08-02*):

  - add Makefile to upload to Pypi
  - update looper template folder of projectmanager
  - add unbuferred output to ngs_analysis recipe
  - add atacseq.get_gene_level_changes
  - improve atacseq.get_gene_level_accessibility
  - add 2D support to general.signed_mean
  - Fix issue #20: len(attributes_to_plot) in general.unsupervised_analysis can be 1 now
  - trackmanager:

    - add default attributes to specified in project_config.group_attributes or otherwise to ['sample_name']
    - Fix issue #16: trackmanager output indentation
    - fix empty subGroups in UCSC trackDb file
    - remove required attributes if no value is found

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
