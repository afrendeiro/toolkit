#!/usr/bin/env python

organism_to_species_mapping = {
    "human": "hsapiens",
    "mouse": "mmusculus",
    "yeast": "scerevisiae",
}
organism_to_latest_ensembl_mapping = {
    "human": "grch38",
    "mouse": "grcm38",
    "yeast": "R64",
}
genome_to_organism_mapping = {
    "hg38": "human",
    "hg19": "human",
    "mm10": "mouse"
}
ucsc_to_ensembl_mapping = {
    "hg38": "grch38",
    "hg19": "grch37",
    "mm10": "grcm38",
    "mm9": "grcm37",
}
genome_to_ensembl_mapping = ucsc_to_ensembl_mapping
