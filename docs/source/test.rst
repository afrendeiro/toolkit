Testing
=============================

To make sure everything is correctly configured, the user is encouraged to test the library prior to use.

In order to do this, install testing requirements and simply run ``pytest``:

.. code-block:: bash

   pip install ngs-toolkit[testing]
   pytest --pyargs ngs_toolkit


Pytest will output summary results (`see for example <https://travis-ci.org/afrendeiro/toolkit/jobs/580167563>`_) and further outputs can be seen in ``${TMPDIR}/pytest-of-${USER}/`` or ``/tmp/pytest-of-${USER}/`` if $TMPDIR is not defined.

