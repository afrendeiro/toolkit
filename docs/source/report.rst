Analysis reports
******************************

.. _Report:

Each analysis object in the `ngs_toolkit` will by default record the outputs it produces (e.g. tables, figures).

This allows the collection of all outputs in a standardized way and the generation of an HTML report.

By default every time a new output is produced, a new report is generated, in a way that analysis progress can be easily monitored in a user-friendly way by simply refreshing the HTML report file. This continuous generation behaviour can be controlled in the configuration file.

The recording behaviour can also be controlled independently for tables and figures by setting the respective values of the configuration file:

.. code-block:: yaml

    preferences:
      report:
        record_figures: True
        record_csv: True
        continuous_generation: True

The report will by default be generated in the root of the project directory, but this can be controlled by manually calling the ``Analysis.generate_report`` function at the user's will.
