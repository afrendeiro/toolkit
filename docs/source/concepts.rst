Concepts
******************************

A few notes on the way some of the library and its objects were designed to be used.

.. _AnalysisObjects:

Analysis objects
==============================

The ``Analysis`` object and its data-type specific dependents are central to the usage of ``ngs-toolkit``. These objects hold attributes and functions relevant to the analysis, such as ``Sample`` objects (and their attributes), Dataframes with numerical values, and others.

.. _LeveragingOnThePEPFormat:

Leveraging on the PEP format
----------------------------

One easy and recommended way to instantiate ``Analysis`` objects is with a ``PEP Project`` file.
This has several advantages:

- Usage of the language-agnostic PEP format to store a project description and interoperability with other tools (see https://github.com/pepkit for other tools);
- Initialization of project-specific variables into the ``Analysis`` object that are derived from the PEP. Examples: analysis samples, genome(s), sample and sample group attributes, sample comparison table.

The example below shows how this works:

.. code-block:: python

	>>> from ngs_toolkit import Analysis
	>>> an = Analysis(from_pep="my_project/metadata/project_config.yaml")
	[INFO] > Setting project's 'sample_attributes' as the analysis 'sample_attributes'.
	[INFO] > Setting project's 'group_attributes' as the analysis 'group_attributes'.
	[INFO] > Setting project's 'comparison_table' as the analysis 'comparison_table'.
	[INFO] > Setting analysis organism as 'mouse'.
	[INFO] > Setting analysis genome as 'mm10'.
	>>> print(an)
	Analysis 'my_project' with 12 samples of organism 'mouse' (mm10).

.. note:: The verbosity of ``ngs-toolkit`` can be controlled

		See the section on `logging <log_config.html#Configuration>`__ to control the verbosity of ``ngs-toolkit``.

.. _ReasonableDefaults:

Reasonable defaults with full customization
-------------------------------------------

Functions in the ``Analysis`` object are aware of these attributes and will 
use them by default, making calling the functions very simple (other overiding 
arguments can be passed though).

In the example below, we will generate a consensus peak set for ATAC-seq 
analyss using the ``get_consensus_sites`` function. This will demonstrate 
several things that "come for free":

.. code-block:: python

	>>> from ngs_toolkit import ATACSeqAnalysis
	>>> an = ATACSeqAnalysis(from_pep="my_project/metadata/project_config.yaml")
	[INFO] > Setting project's 'sample_attributes' as the analysis 'sample_attributes'.
	[INFO] > Setting project's 'group_attributes' as the analysis 'group_attributes'.
	[INFO] > Setting project's 'comparison_table' as the analysis 'comparison_table'.
	[INFO] > Subsetting samples for samples of type 'ATAC-seq'.
	[INFO] > Subsetting comparison_table for comparisons of type 'ATAC-seq'.
	[INFO] > Setting analysis organism as 'mouse'.
	[INFO] > Setting analysis genome as 'mm10'.
	>>> an.get_consensus_sites() 

- even though the PEP project includes samples from several data types (ATAC-, ChIP- and RNA-seq), the current analysis will only consider ATAC-seq samples.
- the necessary files with peak calls for each sample are not specified - ``ngs-toolkit`` knows where to find them;
- a BED file with ENCODE blacklisted regions will not be given, but these regions will be filtered out - ``ngs-toolkit`` will download this and use it. No static files are distributed with the package.
- related to the above, the correct blacklist file is downloaded because the genome assembly for the project is infered from the samples - even though it is not directly specified.

.. _Workflow:

Workflow
------------------------------

Most functions of the ``Analysis`` object will take some input (usually a 
dataframe), apply some transformation and assign the result to a 
variable of the same Analysis object.

To see what variable has been assigned within a given function check the 
relevant function in the `API <api.html#Configuration>`__, specifically the 
`Variables` value. Some functions will assign attributes that are used almost
ubiquitily. See the `common attributes section <concepts.html#CommonAttributes>`__ for some examples.

High-level functions will also often assign their outputs to the object 
itself. To see which attribute holds it, note the ``Attributes`` section of 
the respective function documentation.
Assignment allows the exchange of information between analysis steps without 
the user always providing all required inputs, which would make using such a 
toolkit quite verbose.

The example below illustrates this:

.. code-block:: python

	>>> from ngs_toolkit import ATACSeqAnalysis
	>>> an = ATACSeqAnalysis(from_pep="my_project/metadata/project_config.yaml")
	>>> print(an)
	'ATAC-seq' analysis 'test-project_ATAC-seq_mm10_1_100_1' with 2 samples of organism 'mouse' (mm10).
	>>> an.get_consensus_sites()
	>>> an.measure_coverage()
	>>> print(an.matrix_raw.head())
	                        S1_a1  S2_a2
	region                              
	chr1:42447241-42447627    955   2211
	chr1:44445678-44446750   1939   2122
	chr1:44743959-44744926   1264   1443
	chr1:90513210-90513978   1262   1354
	chr1:93565764-93567191    911    892
	>>> an.normalize()
	>>> print(an.matrix_norm.head())
	region                      S1_a1      S2_a2
	chr1:42447241-42447627  12.681954  13.822151
	chr1:44445678-44446750  13.703582  13.762881
	chr1:44743959-44744926  13.086324  13.206576
	chr1:90513210-90513978  13.084040  13.114743
	chr1:93565764-93567191  12.613915  12.512715

All three ``get_consensus_sites``, ``measure_coverage`` and ``normalize`` build
on the output of each other, but the user doesn't have to specify the input to
any. Changing either the name of the attribute that stores either output or the
location of files outputed is nonetheless easy.

Many functions also have a ``save`` argument which will save the result as a
``CSV`` file.

.. _CommonAttributes:

Common attributes
-----------------

To allow a uniform usage across different data types and analysis types,
a few but important attributes of the ``Analysis`` object and its derivatives
have naming conventions:

- ``data_type``: The type of data of the analysis. Matches the object type.
- ``matrix_raw``: A dataframe of raw, unnormalized values of shape (features, samples)
- ``matrix_norm``: A dataframe of normalized values of shape (features, samples)
- ``quantity``: The name of the units of the values measured. E.g. "expression" for RNA-seq or "accessibility" for ATAC-seq
- ``var_unit_name``: The name of the variables measured. E.g. "gene" for RNA-seq or "region" for ATAC-seq or ChIP-seq
- ``norm_method``: The method used to normalize the ``matrix_norm`` dataframe
- ``thresholds``: A dictionary with keys "log_fold_change" and "p_value" storing thresholds used in the analysis

.. _ComparisonTable:

Comparison table
===============================

``ngs-toolkit`` has functions to perform supervised differntial comparisons 
between groups of samples. The sample groupings are specified in a CSV file called ``comparison_table``.

An example of a typical "case vs control" comparison table is given below:

.. csv-table:: Typical example of comparison_table
   :header: "comparison_name", "comparison_side", "sample_name", "sample_group"
   :widths: 30, 30, 30, 30

	"KOA_vs_WT",	"1",	"ATAC-seq_KOA_r1",	"KO_A"
	"KOA_vs_WT",	"1",	"ATAC-seq_KOA_r2",	"KO_A"
	"KOA_vs_WT",	"0",	"ATAC-seq_WT_r1",	"WT"
	"KOA_vs_WT",	"0",	"ATAC-seq_WT_r2",	"WT"
	"KOB_vs_WT",	"1",	"ATAC-seq_KOB_r1",	"KO_B"
	"KOB_vs_WT",	"1",	"ATAC-seq_KOB_r2",	"KO_B"
	"KOB_vs_WT",	"0",	"ATAC-seq_WT_r1",	"WT"
	"KOB_vs_WT",	"0",	"ATAC-seq_WT_r2",	"WT"


Each row is reserved for a given sample. Samples of the same group (typically 
replicates) should have the same value of "sample_group" and same 
"comparison_side". The group of interest (comparison foreground) should have a 
value of 1 as "comparison_side" and the background a value of 0. Finally, the 
comparison will be labeled with the value of "comparison_name", which should 
be constant for all samples in both foreground and background groups.


For an all-vs-all group comparison, I recommend labeling all background sample groups as a new group in the following manner:

.. csv-table:: "All-vs-all" example of comparison table
   :header: "comparison_name", "comparison_side", "sample_name", "sample_group"
   :widths: 30, 30, 30, 30

	"celltypeA",	"1",	"ATAC-seq_celltypeA_r1",	"ct_A"
	"celltypeA",	"1",	"ATAC-seq_celltypeA_r2",	"ct_A"
	"celltypeA",	"0",	"ATAC-seq_celltypeB_r1",	"ct_A_background"
	"celltypeA",	"0",	"ATAC-seq_celltypeB_r2",	"ct_A_background"
	"celltypeA",	"0",	"ATAC-seq_celltypeC_r1",	"ct_A_background"
	"celltypeA",	"0",	"ATAC-seq_celltypeC_r2",	"ct_A_background"
	"celltypeB",	"1",	"ATAC-seq_celltypeB_r1",	"ct_B"
	"celltypeB",	"1",	"ATAC-seq_celltypeB_r2",	"ct_B"
	"celltypeB",	"0",	"ATAC-seq_celltypeA_r1",	"ct_B_background"
	"celltypeB",	"0",	"ATAC-seq_celltypeA_r2",	"ct_B_background"
	"celltypeB",	"0",	"ATAC-seq_celltypeC_r1",	"ct_B_background"
	"celltypeB",	"0",	"ATAC-seq_celltypeC_r2",	"ct_B_background"


Additional useful columns are `data_type` (to subset comparisons based on type 
of NGS data), `comparison_type` to specify the type of comparison to perform
(e.g. one of 'differential' or 'peaks') and `toggle` for subsetting 
comparisons to perform.


.. note:: **Hyphens and other symbols in comparison_table**
	
	Since differential comparisons are perfomed using DESeq2, R is used 
	(throught the Python-R interface library rpy2).
	ngs_toolkit will create the required tables by DESeq2 which includes names 
	of samples and comparisons as dataframe columns. Unfortunately due to the 
	way R handles column names, these get changed.

	In the future this will be accounted for but for now avoid using hyphens 
	and any other symbols as values for sample names or groups.


.. _LowLevelFunctions:

Low-level functions - ``utils``
===============================

Functions from Analysis objects are generally pretty high level functions, 
often performing several tasks by calling other more general-purpose 
functions. However, one of the concepts I really wanted to have is that the 
user retains as much control as they wish.

They may choose to use the high level functions which generally provide 
sensible defaults, or retain more control and build their analysis pipeline 
from the lower level helper functions.

One example: calling ``ATACSeqAnalysis.normalize()`` will by default run 3-4
other functions to return a quantile normalized, GC-corrected, log-transformed
output - a fairly complex normalization procedure but made simple by providing 
sensible defaults.

A user may easily change the procedure by choosing one of the ~4 types of
normalization using keyword arguments or implement an alternative method which 
can be plugged in to the next step of the analysis.

In the future the low level functions will be moved to `ngs_toolkit.utils` and 
the data type-specific modules will have only classes and functions specific 
to those data which are usually more high-level.
