Manager programs 
******************************

`ngs_toolkit` comes with two programs that provide a command line interface (CLI):
 - ``projectmanager`` handles the creation of a `looper` project, providing sensible configuration templates and git-enabled tracking of changes.
 - ``trackmanager`` handles the creation of a UCSC trackhub or IGV link for ATAC/ChIP-seq data based on bigWig files created by ``looper`` pipelines.


.. note:: `Copying vs linking bigWig files files in trackmanager`
	
	The intention of trackmanager is to create a hierarchy of files in a HTTP server which can be used by genome browsers.
	This requires files (and their parent directories) to be readable and executable.
	When soft-linking files, they will retain the permission attributes of the original files and this may not be appropriate to serve through a server.
	Be aware that copying or linking these files does not always works (manual movement of files might be required).


.. note:: `Changing permissions of files and directories in bigwig directory`
	
	Trackmanager will try to change the permissions of the bigwig files and their parent directories to allow reading and execution by everyone.
	Be aware that this does not always works (manual permission changes might be required).


Here you can see the command-line usage instructions for the main looper command and for each subcommand:
