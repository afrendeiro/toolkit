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
	
	In fact none of it is required, but this allows usage of various tools for enrichment analysis for example.


.. code-block:: yaml


	username: user
	email: user@email.com
	metadata:
	  root_projects_dir: /scratch/lab/shared/projects/
	resources:
	  genomes:
	    2bit:
	      hg19: /scratch/lab/shared/resources/genomes/hg19/hg19.2bit
	      hg38: /scratch/lab/shared/resources/genomes/hg38/hg38.2bit
	      mm10: /scratch/lab/shared/resources/genomes/mm10/mm10.2bit
	    fasta:
	      hg19: /scratch/lab/shared/resources/genomes/hg19/hg19.fa
	      hg38: /scratch/lab/shared/resources/genomes/hg38/hg38.fa
	      mm10: /scratch/lab/shared/resources/genomes/mm10/mm10.fa
	  lola:
	    region_databases:
	      # under each section, there should be a list of items
	      hg19:
	        - /data/groups/lab/shared/resources/regions/LOLACore/hg19/
	        - /data/groups/lab/shared/resources/regions/customRegionDB/hg19/
	      hg38:
	        - /data/groups/lab/shared/resources/regions/LOLACore/hg38/
	        - /data/groups/lab/shared/resources/regions/customRegionDB/hg38/
	      mm10:
	        - /data/groups/lab/shared/resources/regions/LOLACore/mm10/
	        - /data/groups/lab/shared/resources/regions/customRegionDB/mm10/
	  meme:
	    motif_databases:
	      human: /scratch/lab/shared/resources/motifs/motif_databases/HUMAN/HOCOMOCOv10.meme
	      mouse: /scratch/lab/shared/resources/motifs/motif_databases/MOUSE/uniprobe_mouse.meme
	      vertebrate: /home/user/.local/homer_4.8/data/knownTFs/vertebrates/known.motifs
	    motif_id_mapping:
	      mouse: /scratch/lab/shared/resources/motifs/motif_databases/MOUSE/uniprobe_mouse.id_mapping.tsv
	  enrichr:
	    gene_set_libraries:
	      # this should be a list of items
	      - "GO_Biological_Process_2015"
	      - "ChEA_2015"
	      - "KEGG_2016"
	      - "ESCAPE"
	      - "Epigenomics_Roadmap_HM_ChIP-seq"
	      - "ENCODE_TF_ChIP-seq_2015"
	      - "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"
	      - "ENCODE_Histone_Modifications_2015"
	      - "OMIM_Expanded"
	      - "TF-LOF_Expression_from_GEO"
	      - "Single_Gene_Perturbations_from_GEO_down"
	      - "Single_Gene_Perturbations_from_GEO_up"
	      - "Disease_Perturbations_from_GEO_down"
	      - "Disease_Perturbations_from_GEO_up"
	      - "Drug_Perturbations_from_GEO_down"
	      - "Drug_Perturbations_from_GEO_up"
	      - "WikiPathways_2016"
	      - "Reactome_2016"
	      - "BioCarta_2016"
	      - "NCI-Nature_2016"
	executables:
	  twoBitToFa: twoBitToFa
	  fasta-dinucleotide-shuffle: fasta-dinucleotide-shuffle
	  ame: ame
	  findMotifsGenome.pl: findMotifsGenome.pl
	  compareMotifs.pl: compareMotifs.pl