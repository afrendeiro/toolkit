Changelog
******************************

- Release: **v0.2** (*2018-06-12*):
  - New
    - Command-line interface (CLI) based on sub-commands for ``projectmanager``.
    - Recipes: scripts which ``projectmanager`` can run.
    - General `ngs_analysis` recipe for general analysis of NGS projects.
  - Changed
    - Removed requirement to have ``pipelines`` repository installed in order to extend base Sample objects
    - Maintenance of sample attributes as provided by user by means of reading them in as strings (to be improved further)
    - Improved serialization of Sample objects

- Upcoming:
   - New:
     - `ngs_toolkit.utils` to contain small helper functions.
     - 
   - Changed:
     - Decoupling of static files from ``data/external``.
     - Better hg38 support.
     - 