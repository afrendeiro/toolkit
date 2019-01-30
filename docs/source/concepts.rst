Concepts
******************************

A few notes on the way some of the objects were designed to be used.


Analysis objects
==============================

Analysis objects are generally instantiated with ``peppy`` ``Project`` objects, which instead build on ``PEP``.
This allows the usage of looper processed samples and project variables through the project configuration file.

These objects will hold attributes and functions relevant to the data type under analysis. Depending on the type of function these functions may assign their outputs to the object or return them.


Assigning vs returning results
------------------------------

Most functions with specialized functions to a data type will take some input (usually a dataframe), apply some transformation and assign the result to a variable of the same Analysis object.

To see what variable has been assigned within a given function check the relevant function in the ``API``, specifically the `Variables` value.

Assignment allows the exchange of information between analysis steps without the user always providing all required inputs, which would make using such a toolkit quite verbose.

Lower level functions are generally not attached to an Analysis class and tend to return their output.


High- vs low-level functions
------------------------------

Functions from Analysis objects are genrally pretty high level functions, often performing several tasks by calling other more general-purpose functions. However, one of the concepts I really wanted to have is that the user retains as much control as they wish.

They may choose to use the high level functions which generally provide sensible defaults, or retain more control and build their analysis pipeline from the lower level helper functions.

One example: calling `ATACSeqAnalysis.normalize()` will by default run 3-4 other functions to return a quantile normalized, GC-corrected, log-transformed output - a fairly complex normalization procedure but made simple by providing sensible defaults.

A user may easily change the procedure by choosing one of the ~4 types of normalization using keyword arguments or implement an alternative method which can be plugged in to the next step of the analysis.

In the future the low level functions will be moved to `ngs_toolkit.utils` and the data type-specific modules will have only classes and functions specific to those data which are usually more high-level.


Cross-data type functions
------------------------------

Some high-level functions that perform tasks common to various analysis types generally work by passing an Analysis object as argument. This is unlike most other high-level functions in that usually the function comes from the object itself.
