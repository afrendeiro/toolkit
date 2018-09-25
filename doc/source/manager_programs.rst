Manager programs 
******************************

`ngs_toolkit` comes with two programs that provide a command line interface (CLI):
 - ``projectmanager`` handles the creation of a `looper` project, providing sensible configuration templates and git-enabled tracking of changes.
 - ``trackmanager`` handles the creation of a UCSC trackhub or IGV link for ATAC/ChIP-seq data based on bigWig files created by ``looper`` pipelines.


Here you can see the command-line usage instructions for the main looper command and for each subcommand:


projectmanager
=============================

.. code-block:: none

	usage: projectmanager [-h] {create,recipe} ...

	projectmanager - A project manager.

	positional arguments:
	  {create,recipe}
	    create         Create project.
	    recipe         Run ngs_toolkit recipe for a given project.

	optional arguments:
	  -h, --help       show this help message and exit

	https://github.com/afrendeiro/toolkit



projectmanager::create
-----------------------------

.. code-block:: none

	usage: projectmanager create [-h] [-r ROOT_DIR] [-d] [--overwrite]
	                             project_name

	Create project.

	positional arguments:
	  project_name          Project name.

	optional arguments:
	  -h, --help            show this help message and exit
	  -r ROOT_DIR, --root-dir ROOT_DIR
	                        Root directory to create projects.
	  -d, --dry-run         Don't actually do anything.
	  --overwrite           Don't overwrite any existing directory or file.


projectmanager::recipe
-----------------------------

.. code-block:: none

	usage: projectmanager recipe [-h] recipe_name project_config

	Run recipe.

	positional arguments:
	  recipe_name     Recipe name.
	  project_config  Project config.

	optional arguments:
	  -h, --help      show this help message and exit


trackmanager
=============================

.. code-block:: none

	usage: trackmanager [-h] [-a [ATTRIBUTES]] [-c COLOR_ATTRIBUTE] [-r] [-l]
	                    project_config_file
	
	positional arguments:
	  project_config_file
	
	optional arguments:
	  -h, --help            show this help message and exit
	  -a [ATTRIBUTES], --attrs [ATTRIBUTES]
	                        Sample attributes (annotation sheet columns) to use to
	                        order tracks. Add attributes comma-separated with no
	                        whitespace.
	  -c COLOR_ATTRIBUTE, --color-attr COLOR_ATTRIBUTE
	                        Sample attribute to use to color tracks with. Default
	                        is first attribute passed.
	  -r, --overlay-replicates
	                        Whether replicate samples should be overlaied in same
	                        track. Default=False.
	  -l, --link            Whether bigWig files should be soft-linked to the
	                        track database directory. Default=False.


.. note:: `Copying vs linking bigWig files files in trackmanager`
	
	The intention of trackmanager is to create a hierarchy of files in a HTTP server which can be used by genome browsers.
	This requires files (and their parent directories) to be readable and executable.
	When soft-linking files, they will retain the permission attributes of the original files and this may not be appropriate to serve through a server.
	Be aware that copying or linking these files does not always works (manual movement of files might be required).


.. note:: `Changing permissions of files and directories in bigwig directory`
	
	Trackmanager will try to change the permissions of the bigwig files and their parent directories to allow reading and execution by everyone.
	Be aware that this does not always works (manual permission changes might be required).

