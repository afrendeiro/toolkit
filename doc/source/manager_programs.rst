Manager programs 
******************************

`ngs_toolkit` comes with two programs that provide a command line interface (CLI):
 - ``projectmanager`` handles the creation of a `looper` project, providing sensible configuration templates and git-enabled tracking of changes.
 - ``trackmanager`` handles the creation of a UCSC trackhub or IGV link for ATAC/ChIP-seq data based on bigWig files created by ``looper`` pipelines.

Here you can see the command-line usage instructions for the main looper command and for each subcommand:

``projectmanager --help``
----------------------------------

.. code-block:: none

	usage: projectmanager [-h] [-d] [--overwrite] project_name
	
	projectmanager - A project creator.
	
	positional arguments:
	  project_name   Project name.
	
	optional arguments:
	  -h, --help     show this help message and exit
	  -d, --dry-run  Don't actually do anything. (default: False)
	  --overwrite    Don't overwrite any existing directory or file. (default:
	                 False)
	
	https://github.com/afrendeiro/toolkit

``trackmanager --help``
----------------------------------

.. code-block:: none

	usage: trackmanager [-h] [-a [ATTRIBUTES]] [-c COLOR_ATTRIBUTE] [-r]
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
