Configuration and logging
******************************


Logging
=============================

`ngs_toolkit` will log its operations and errors using the Python standard logging library.

This will happen by default to standard output (sys.stdout) but also to a file in ``$HOME/.ngs_toolkit.log.txt``.

The location of the log file and the level of events to be reported can be customized in the ``ngs_toolkit.setup_logger()`` function.


Configuration
=============================

`ngs_toolkit` allows user configuration to allow usage of static resources across computing environments.

To accomplish this, the user can provide its own configuration in two ways:

* In a YAML file located in ``$HOME/.ngs_toolkit.config.yaml``;
* A user provided file given during interactive runtime passed to ``ngs_toolkit.setup_config()``.

If more than one is given values in the configuration files will be updated in the following order:

1. A minimal configuration file from the package data;
2. The user provided file in ``$HOME/.ngs_toolkit.config.yaml``;
3. The user provided file passed to ``ngs_toolkit.setup_config()``.

To see how to structure the YAML file, see section below.



Example configuration file
-----------------------------


.. note:: `Not all elements are required`
    
    In fact none of it is required, but it is recommended to have a look at the template configuration file and set custom options.


.. code-block:: yaml

    username: user
    email: user@mail.com
    website_root: userwebsite.web.com
    preferences:
      # For the next item, environment variables are formatted if they are of the form ${VAR}
      root_reference_dir: ${USER}/reference_data
      root_projects_dir: ${USER}/projects
      default_genome_assemblies:
        - human: hg38
        - mouse: mm10
      computing_configuration: 'slurm'
    sample_input_files:
      ATAC-seq:
        aligned_filtered_bam: "{data_dir}/{sample_name}/mapped/{sample_name}.bowtie2.filtered.bam"
        peaks: "{data_dir}/{sample_name}/peaks/{sample_name}_peaks.narrowPeak"
        summits: "{data_dir}/{sample_name}/peaks/{sample_name}_summits.narrowPeak"
      ChIP-seq:
        aligned_filtered_bam: "{data_dir}/{sample_name}/mapped/{sample_name}.bowtie2.filtered.bam"
      CNV:
        aligned_filtered_bam: "{data_dir}/{sample_name}/mapped/{sample_name}.bowtie2.filtered.bam"
      RNA-seq:
        aligned_filtered_bam: "{data_dir}/{sample_name}/mapped/{sample_name}.bowtie2.filtered.bam"
        bitseq_counts: "{data_dir}/{sample_name}/quantification/{sample_name}_bitseq.tsv"
