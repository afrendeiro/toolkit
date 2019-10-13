Distributed computing
=============================

divvy
-----------------------------

Certain functions in the ``ngs_toolkit`` toolkit can make use of distributed
computing. To achieve this for a variety of computing configurations
it uses the `divvy library <http://divvy.databio.org/>`_.

Divvy provides an abstract way of submitting a job to various job managers by
shipping job templates for each configuration.

When ``divvy`` starts, a configuration is chosen (the ``compute_configuration``
attribute) and that template gets filled with the attributes of the job -
the code to be executed, the resouce requirements and others
(e.g. "cores", "mem", "time" attributes).

To see all supported compute configurations run:

.. code-block:: bash

    divvy list

For more information on how to configure ``divvy``, see its documentation:
http://divvy.databio.org/

To let ``ngs_toolkit`` know which ``divvy`` configuration to use by default,
modify the following section in the ``ngs_toolkit`` configuration file:

.. code-block:: yaml

    preferences:
      # The next item is the default computing configuration to use from divvy.
      # Run "divvy list" to see all options.
      # See more here: http://code.databio.org/divvy/
      computing_configuration: 'slurm'

This will make ``ngs_toolkit`` send jobs to a slurm cluster if wanted.

All functions that allow running a task in a distributed manner have a
``distributed`` keyword argument.

In addition, they also accept additional keyword arguments (`kwargs` in the
function signature) where additional options can be passed.
These options must match fields available to format of the currently selected
``compute_configuration``. 

Sending jobs and collecting output
----------------------------------

Performing a taks in a distributed manner can therefore be as simple as calling 
the desired function with ``distributed=True``. Jobs will be sent to the
job manager of the chosen computing configuration.

However, since the jobs are often run individually for a sample/group of samples,
functions called with ``distributed=True`` may not return the same output as
``distributed=False``.

For that reason, for all such functions, there is a reciprocal function of
identical name as the first prefixed with ``collect``.

.. code-block:: python

    from ngs_toolkit.demo import generate_project
    an = generate_project(sample_input_files=True)
    an.measure_coverage(distributed=True)
    coverage = collect_coverage()

Implementing automatic collection of job outputs in part of future plans.

Example
-----------------------------

The :func:`ngs_toolkit.atacseq.ATACSeqAnalysis.measure_coverage` function has
``distributed`` and ``kwargs`` options.

This provides code portability and allows customization of various aspects of
the jobs:

.. code-block:: python

    from ngs_toolkit.demo import generate_project
    an = generate_project(sample_input_files=True)
    # in serial
    cov1 = an.measure_coverage()
    # as slurm jobs (because the config computing_configuration is set to 'slurm')
    an.measure_coverage(distributed=True)
    cov2 = collect_coverage()
    # confirm the output is the same
    assert (cov2 == cov1).all().all()

.. code-block:: python

    # as slurm jobs to a particular queue and more memory
    an.measure_coverage(distributed=True, partition="longq", mem=24000)
    # here 'partition' and 'mem' are attributes of the slurm divvy template
    # and not magic attributes
