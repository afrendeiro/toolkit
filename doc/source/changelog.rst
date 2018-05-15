Changelog
******************************

- Release: **v0.1.3.5** (*2018-06-15*):
  - New
    - Extended documentation
    - Command-line interface (CLI) based on sub-commands for ``projectmanager``.
    - Recipes: scripts which ``projectmanager`` can run.
    - General `ngs_analysis` recipe for general analysis of NGS projects.


- Upcoming:
   - New:
     - `ngs_toolkit.utils` to hold small helper functions.
   - Changes:
     - Removed requirement to have ``pipelines`` repository installed in order to extend base Sample objects
     - Decoupling of static files from ``data/external``.
     - Maintenance of sample attributes as provided by user by means of reading them in as strings (to be improved further)
     - Improved serialization of Sample objects
     - Better hg38 support.
