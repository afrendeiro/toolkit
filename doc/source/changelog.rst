Changelog
******************************

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
