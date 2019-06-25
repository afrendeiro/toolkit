Concepts
******************************

A few notes on the way some of the objects were designed to be used.


Analysis objects
==============================

One easy and recommended way to instantiate ``Analysis`` objects is with a ``PEP Project`` file.
This allows the usage of looper processed samples and project variables through the project configuration file.

These objects hold attributes and functions relevant to the data type under analysis.


Most of the time, functions will assign their outputs to the object itself. To see which
attribute holds it, note the ``Attributes`` section of the respective function documentation.

.. Example:
	from ngs_toolkit import Analysis
	a = Analysis("demo")
	a.normalize()  # see the ``Attributes`` section of the documentation for this function
	a.matrix_norm.head()



..Inheritance example


Workflow
------------------------------

Most functions with specialized functions to a data type will take some input (usually a dataframe),
apply some transformation and assign the result to a variable of the same Analysis object.

To see what variable has been assigned within a given function check the relevant function in the ``API``,
specifically the `Variables` value.

Assignment allows the exchange of information between analysis steps without the user always providing
all required inputs, which would make using such a toolkit quite verbose.


Attributes
----------

Default

data_type, thresholds

Assigned 

matrix_raw, norm_method


Low-level functions - ``utils``
===============================

Functions from Analysis objects are generally pretty high level functions, often performing several
tasks by calling other more general-purpose functions. However, one of the concepts I really wanted
to have is that the user retains as much control as they wish.

They may choose to use the high level functions which generally provide sensible defaults, or retain
more control and build their analysis pipeline from the lower level helper functions.

One example: calling `ATACSeqAnalysis.normalize()` will by default run 3-4 other functions to return
a quantile normalized, GC-corrected, log-transformed output - a fairly complex normalization procedure
but made simple by providing sensible defaults.

A user may easily change the procedure by choosing one of the ~4 types of normalization using keyword
arguments or implement an alternative method which can be plugged in to the next step of the analysis.

In the future the low level functions will be moved to `ngs_toolkit.utils` and the data type-specific
modules will have only classes and functions specific to those data which are usually more high-level.
