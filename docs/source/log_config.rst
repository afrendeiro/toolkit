Configuration, logging and versioning
*************************************

.. _Configuration:


Configuration
=============================

`ngs_toolkit` uses a YAML configuration file.

While entirely optional, this allows the user to specify preferences, patterns and allows usage across different computing environments.

The user can provide its own configuration in two ways:

* In a YAML file located in ``$HOME/.ngs_toolkit.config.yaml``;
* A user provided file given during interactive runtime passed to ``ngs_toolkit.setup_config()``.

If more than one is given values in the configuration files will be updated in the following order:

1. A minimal configuration file from the package data;
2. The user provided file in ``$HOME/.ngs_toolkit.config.yaml``;
3. The user provided file passed to ``ngs_toolkit.setup_config()``.

To see how to structure the YAML file, see section below.



Example configuration files
-----------------------------

To see all available configuration fields have a look at the default configuration file: https://github.com/afrendeiro/toolkit/blob/master/ngs_toolkit/config/default.yaml#L1

For a full example of a fully configured file have a look at the example configuration file: https://github.com/afrendeiro/toolkit/blob/master/ngs_toolkit/config/example.yaml#L1

However, the configuration file does not need to include all fields. Below is a minimal example of a configuration file.

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
      # Below is the name of the divvy package configuration (http://divvy.databio.org/en/latest/)
      computing_configuration: 'slurm'
    sample_input_files:
      ATAC-seq:
        aligned_filtered_bam: "{data_dir}/{sample_name}/mapped/{sample_name}.bowtie2.filtered.bam"
        peaks: "{data_dir}/{sample_name}/peaks/{sample_name}_peaks.narrowPeak"
        summits: "{data_dir}/{sample_name}/peaks/{sample_name}_summits.narrowPeak"
      ChIP-seq:
        aligned_filtered_bam: "{data_dir}/{sample_name}/mapped/{sample_name}.bowtie2.filtered.bam"
      CNV:
        log2_read_counts: "{data_dir}/{sample_name}/{sample_name}_{resolution}/CNAprofiles/log2_read_counts.igv"
      RNA-seq:
        aligned_filtered_bam: "{data_dir}/{sample_name}/mapped/{sample_name}.bowtie2.filtered.bam"
        bitseq_counts: "{data_dir}/{sample_name}/quantification/{sample_name}_bitseq.tsv"


.. note:: `Not all elements are required`
    
    In fact none of it is required, but it is recommended to have a look at the template configuration file and set custom options.

.. _Logging:

Logging
=============================

`ngs_toolkit` will log its operations and errors using the Python standard logging library.

This will happen by default to standard output (sys.stdout) but also to a file in ``$HOME/.ngs_toolkit.log.txt``.

The location of the log file and the level of events to be reported can be customized in the ``ngs_toolkit.setup_logger()`` function.


.. _Versioning:

Versioning
=============================

`ngs_toolkit` will by default timestamp every output it produces (CSV and figure files).

This behaviour can be controlled independently for tables and figures by setting the respective values of the configuration file:

.. code-block:: yaml

    preferences:
      report:
        timestamp_figures: False
        timestamp_tables: False
