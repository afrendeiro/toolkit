Comparison table
******************************


`ngs_toolkit` has functions to perform supervised differntial comparisons between groups of samples.
The sample groupings are specified in a CSV file called `comparison_table`.


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


Each row is reserved for a given sample. Samples of the same group (typically replicates) should have the same value of "sample_group" and same "comparison_side". The group of interest (comparison foreground) should have a value of 1 as "comparison_side" and the background a value of 0. Finally, the comparison will be labeled with the value of "comparison_name", which should be constant for all samples in both foreground and background groups.


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


Additional useful columns are `data_type` (to subset comparisons based on type of NGS data), `comparison_type` to specify the type of comparison to perform (e.g. one of 'differential' or 'peaks') and `toggle` for subsetting comparisons to perform.


.. note:: **Hyphens and other symbols in comparison_table**
	
	Since differential comparisons are perfomed using DESeq2, R is used (throught the Python-R interface library rpy2).
	ngs_toolkit will create the required tables by DESeq2 which includes names of samples and comparisons as dataframe columns. Unfortunately due to the way R handles column names, these get changed.

	In the future this will be accounted for but for now avoid using hyphens and any other symbols as values for sample names or groups.
