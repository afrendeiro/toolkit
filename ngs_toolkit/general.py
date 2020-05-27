#!/usr/bin/env python


import os

import numpy as np
import pandas as pd

from ngs_toolkit import _CONFIG, _LOGGER, MEMORY
from ngs_toolkit.utils import get_this_file_or_timestamped
from ngs_toolkit.exceptions import NetworkError


def get_genome_reference(
    organism,
    genome_assembly=None,
    output_dir=None,
    genome_provider="UCSC",
    file_format="2bit",
    dry_run=False,
    overwrite=True,
):
    """
    Get genome FASTA/2bit file.
    Saves results to disk and returns path to file.

    Parameters
    ----------
    organism : :obj:`str`
        Organism to get annotation for. Currently supported: "human" and "mouse".

    output_dir : :obj:`str`, optional
        Directory to write output to.
        Defaults to current directory

    genome_provider : :obj:`str`, optional
        Which genome provider to use. One of 'UCSC' or 'Ensembl'.

    file_format : :obj:`str`, optional
        File format to get. One of 'fasta' or '2bit'.

    dry_run: :obj:`bool`, optional
        Whether to not download and just return path to file.

    overwrite: :obj:`bool`, optional
        Whether existing files should be overwritten by new ones.
        Otherwise they will be kept and no action is made.
        Defaults to :obj:`True`.

    Returns
    -------
    {str, tuple}
        Path to genome FASTA/2bit file,
        but if `dry_run` tuple of URL of reference genome and path to file.

    Raises
    --------
    ValueError
        If arguments are not in possible options or if desired combination
        is not available.
    """
    import pybedtools
    import subprocess
    import distutils.spawn

    from ngs_toolkit.utils import download_gzip_file, download_file

    def index_fasta(fasta):
        """
        # The idea is to use a hidden method of bedtools
        # to create an index (to skip having e.g. samtools as dependency)
        # and use the bedtool nuc command to to do it.
        # This actually fails to get nucleotide content every time due to this 'bug':
        # https://github.com/daler/pybedtools/issues/147
        # but nonetheless creates an index :whatever:
        """
        bed = pd.DataFrame([["chr1", 1, 2]])
        try:
            bed = pybedtools.BedTool().from_dataframe(bed)
            bed.nucleotide_content(fi=fasta)
        except pybedtools.helpers.BEDToolsError:
            pass

    def twobit_to_fasta(genome_file):
        if distutils.spawn.find_executable("twoBitToFa") is not None:
            args = ["twoBitToFa"] + [genome_file, genome_file.replace(".2bit", ".fa")]
            subprocess.check_output(args, shell=False)
            index_fasta(genome_file.replace(".2bit", ".fa"))

    if output_dir is None:
        output_dir = os.path.join(os.path.abspath(os.path.curdir), "reference")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    opts = ["UCSC", "Ensembl"]
    if genome_provider not in opts:
        msg = "`genome_provider` attribute must be one of '{}'.".format(", ".join(opts))
        _LOGGER.error(msg)
        raise ValueError(msg)

    opts = ["fasta", "2bit"]
    if file_format not in opts:
        msg = "`file_format` attribute must be one of '{}'.".format(", ".join(opts))
        _LOGGER.error(msg)
        raise ValueError(msg)

    if (genome_provider == "Ensembl") and (file_format == "2bit"):
        msg = "Ensembl does not provide 2bit files."
        hint = " Use for example 'faToTwoBit' to convert the FASTA file."
        _LOGGER.error(msg + hint)
        raise ValueError(msg)

    if genome_provider == "UCSC":
        default_organisms = {
            "human": "hg38",
            "hsapiens": "hg38",
            "homo_sapiens": "hg38",
            "mouse": "mm10",
            "mmusculus": "mm10",
            "mus_musculus": "mm10",
            "yeast": "sacCer3",
            "scerevisiae": "sacCer3",
            "saccharomyces_cerevisiae": "sacCer3",
        }
        base_link = "http://hgdownload.cse.ucsc.edu/goldenPath/{assembly}/bigZips/{assembly}"
        base_link += ".fa.gz" if file_format == "fasta" else ".2bit"
        if genome_assembly is None:
            genome_assembly = default_organisms[organism]
        url = base_link.format(assembly=genome_assembly)

    elif genome_provider == "Ensembl":
        default_organisms = {
            "human": {"long_species": "homo_sapiens", "version": "grch38", "release": "75",},
            "hsapiens": {"long_species": "homo_sapiens", "version": "grch38", "release": "75",},
            "homo_sapiens": {"long_species": "homo_sapiens", "version": "grch38", "release": "75",},
            "mouse": {"long_species": "mus_musculus", "version": "grcm38", "release": "94",},
            "mmusculus": {"long_species": "mus_musculus", "version": "grcm38", "release": "94",},
            "mus_musculus": {"long_species": "mus_musculus", "version": "grcm38", "release": "94",},
            "yeast": {
                "long_species": "saccharomyces_cerevisiae",
                "version": "R64",
                "release": "94",
            },
            "scerevisiae": {
                "long_species": "saccharomyces_cerevisiae",
                "version": "R64",
                "release": "94",
            },
            "saccharomyces_cerevisiae": {
                "long_species": "saccharomyces_cerevisiae",
                "version": "R64",
                "release": "94",
            },
        }
        if genome_assembly is None:
            genome_assembly = default_organisms[organism]["version"].replace("grc", "GRC")
        base_link = "ftp://ftp.ensembl.org/pub/release-{release}/fasta/{long_organism}/dna/"
        base_link += "{Clong_organism}.{assembly}."
        base_link += (
            "{}.".format(default_organisms[organism]["release"])
            if genome_assembly.endswith("37")
            else ""
        )
        base_link += "{sequence_type}.{id_type}.{id}fa.gz".format(
            sequence_type="dna", id_type="primary_assembly", id=""
        )
        url = base_link.format(
            release=default_organisms[organism]["release"],
            long_organism=default_organisms[organism]["long_species"],
            Clong_organism=default_organisms[organism]["long_species"].capitalize(),
            assembly=genome_assembly,
        )

    if (genome_provider == "UCSC") and (file_format == "fasta") and (genome_assembly != "hg38"):
        msg = "UCSC does not provide FASTA files for the {} assembly.".format(genome_assembly)
        hint = " Download a 2bit file and use for example 'TwoBitToFa' to convert."
        _LOGGER.error(msg + hint)
        raise ValueError(msg)

    genome_file = os.path.join(
        output_dir,
        "{}.{}.{}".format(organism, genome_assembly, "fa.gz" if file_format == "fasta" else "2bit"),
    )

    if os.path.exists(genome_file) and (not overwrite):
        msg = "Genome file already exists and 'overwrite' is set to False."
        hint = " Returning existing file: {}".format(genome_file)
        _LOGGER.warning(msg + hint)
        # even so, if 2bit and FASTA not there try to get fasta
        if file_format == "2bit":
            if not os.path.exists(genome_file.replace(".2bit", ".fa")):
                twobit_to_fasta(genome_file)
        return genome_file

    # create .fai index for fasta file
    if file_format == "fasta":
        if not dry_run:
            download_gzip_file(url, genome_file)
            index_fasta(genome_file)
            return genome_file
    else:
        if not dry_run:
            download_file(url, genome_file)
            twobit_to_fasta(genome_file)
            return genome_file

    return (url, genome_file)


def get_blacklist_annotations(organism, genome_assembly=None, output_dir=None, overwrite=True):
    """
    Get annotations of blacklisted genomic regions for a given organism/genome assembly.
    Saves results to disk and returns a path to a BED file.

    Parameters
    ----------
    organism : :obj:`str`
        Organism to get annotation for. Currently supported: "human" and "mouse".

    genome_assembly : :obj:`str`, optional
        Ensembl assembly/version to use.
        Default for "human" is "hg38/grch38" and for "mouse" is "mm10/grcm38".

    output_dir : :obj:`str`, optional
        Directory to write output to.
        Defaults to "reference" in current directory.

    overwrite: :obj:`bool`, optional
        Whether existing files should be overwritten by new ones.
        Otherwise they will be kept and no action is made.
        Defaults to :obj:`True`.

    Returns
    -------
    str
        Path to blacklist BED file
    """
    from ngs_toolkit.utils import download_gzip_file

    if output_dir is None:
        output_dir = os.path.join(os.path.abspath(os.path.curdir), "reference")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    default_organisms = {"human": "hg38", "mouse": "mm10"}
    if genome_assembly is None:
        genome_assembly = default_organisms[organism]
        _LOGGER.warning(
            "Genome assembly not selected. Using assembly '{}' for '{}'.".format(
                genome_assembly, organism
            )
        )

    output_file = os.path.join(output_dir, "{}.{}.blacklist.bed".format(organism, genome_assembly))
    if os.path.exists(output_file) and (not overwrite):
        msg = "Annotation file already exists and 'overwrite' is set to False."
        hint = " Returning existing annotation file: {}".format(output_file)
        _LOGGER.warning(msg + hint)
        return output_file

    url = "http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/"
    if genome_assembly not in ["hg19"]:
        url += "{0}-{1}/{0}.blacklist.bed.gz".format(genome_assembly, organism)
    else:
        url += "{0}-human/wgEncodeHg19ConsensusSignalArtifactRegions.bed.gz".format(genome_assembly)

    try:
        download_gzip_file(url, output_file)
    except OSError:
        msg = "Could not download file: {}".format(url)
        _LOGGER.error(msg)
        raise OSError(msg)
    return output_file


def get_tss_annotations(
    organism,
    genome_assembly=None,
    save=True,
    output_dir=None,
    chr_prefix=True,
    gene_types=["protein_coding", "processed_transcript", "lincRNA", "antisense"],
    overwrite=True,
):
    """
    Get annotations of TSS for a given organism/genome assembly.
    This is a simple approach using Biomart's API querying the Ensembl database.
    Saves results to disk and returns a dataframe.

    Parameters
    ----------
    organism : :obj:`str`
        Organism to get annotation for.
        Currently supported: "human" and "mouse".

    genome_assembly : :obj:`str`, optional
        Ensembl assembly/version to use.
        Default for "human" is "grch38" and for "mouse" is "grcm38".

    save: :obj:`bool`, optional
        Whether to save to disk under ``output_dir``.
        Defaults to :obj:`True`.

    output_dir : :obj:`str`, optional
        Directory to write output to.
        Defaults to "reference" in current directory.

    chr_prefix: :obj:`bool`, optional
        Whether chromosome names should have the "chr" prefix.
        Defaults to True

    gene_types : :obj:`list`, optional
        Subset of transcript biotypes to keep.
        See here the available biotypes https://www.ensembl.org/Help/Faq?id=468
        Defaults to 'protein_coding', 'processed_transcript', 'lincRNA', 'antisense'.

    overwrite: :obj:`bool`, optional
        Whether existing files should be overwritten by new ones.
        Otherwise they will be kept and no action is made.
        Defaults to :obj:`True`.

    Returns
    -------
    :obj:`pandas.DataFrame`
        DataFrame with genome annotations
    """
    default_organisms = {
        "human": {"species": "hsapiens", "ensembl_version": "grch38"},
        "mouse": {"species": "mmusculus", "ensembl_version": "grcm38"},
        "yeast": {"species": "scerevisiae", "ensembl_version": "R64"},
    }

    if genome_assembly is None:
        genome_assembly = default_organisms[organism]["ensembl_version"]
    if genome_assembly == "hg19":
        genome_assembly = "grch37"
    if genome_assembly == "hg38":
        genome_assembly = "grch38"
    if genome_assembly == "mm10":
        genome_assembly = "grcm38"
    if genome_assembly == "sacCer3":
        genome_assembly = "R64"

    output_file = os.path.join(
        output_dir, "{}.{}.gene_annotation.tss.bed".format(organism, genome_assembly)
    )
    o = get_this_file_or_timestamped(output_file)
    if os.path.exists(o) and (not overwrite):
        msg = "Annotation file already exists and 'overwrite' is set to False."
        hint = " Returning existing annotation from file: {}".format(o)
        _LOGGER.warning(msg + hint)
        annot = pd.read_csv(o, header=None, sep="\t")
        return annot

    if save:
        if output_dir is None:
            output_dir = os.path.join(os.path.abspath(os.path.curdir), "reference")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    attributes = [
        "chromosome_name",
        "start_position",
        "ensembl_gene_id",
        "external_gene_name",
        "strand",
        "transcript_biotype",
    ]
    res = query_biomart(
        attributes=attributes,
        species=default_organisms[organism]["species"],
        ensembl_version=genome_assembly,
    )

    if gene_types is None:
        res = res.drop(["transcript_biotype"], axis=1).drop_duplicates()
    else:
        res = res.loc[res["transcript_biotype"].isin(gene_types), :]
    if chr_prefix:
        res.loc[:, "chromosome_name"] = "chr" + res["chromosome_name"]
    res.loc[:, "start_position"] = res["start_position"].astype(int)
    res.loc[:, "strand"] = res["strand"].replace("1", "+").replace("-1", "-")
    res.loc[:, "end"] = res.apply(
        lambda x: x["start_position"] + 1 if x["strand"] == "+" else x["start_position"], axis=1,
    )
    res.loc[:, "start_position"] = res.apply(
        lambda x: x["start_position"] - 1 if x["strand"] == "-" else x["start_position"], axis=1,
    )

    # drop same gene duplicates if starting in same position but just annotated with
    # different biotypes
    res = (
        res[attributes[:2] + ["end"] + attributes[2:]]
        .sort_values(by=res.columns.tolist(), axis=0)
        .drop_duplicates(subset=res.columns[:3], keep="last")
    )

    # make real BED format
    res.loc[:, "fill"] = "."
    cols = [
        "chromosome_name",
        "start_position",
        "end",
        "external_gene_name",
        "fill",
        "strand",
        "ensembl_gene_id",
        "transcript_biotype",
    ]
    res = res.loc[:, cols]
    res.columns = [
        "chr",
        "start",
        "end",
        "gene_name",
        "score",
        "strand",
        "ensembl_gene_id",
        "transcript_biotype",
    ]

    # save
    if save:
        res.to_csv(output_file, index=False, header=False, sep="\t")

        output_file = os.path.join(
            output_dir,
            "{}.{}.gene_annotation.protein_coding.tss.bed".format(organism, genome_assembly),
        )
        res[res["transcript_biotype"] == "protein_coding"].drop(attributes[-1], axis=1).to_csv(
            output_file, index=False, header=False, sep="\t"
        )

    return res


def get_genomic_context(
    organism,
    genome_assembly=None,
    save=True,
    output_dir=None,
    chr_prefix=True,
    region_subset=["promoter", "exon", "5utr", "3utr", "intron", "genebody", "intergenic",],
    gene_types=["protein_coding", "processed_transcript", "lincRNA", "antisense"],
    promoter_width=3000,
    overwrite=True,
):
    """
    Get annotations of TSS for a given organism/genome assembly.
    This is a simple approach using Biomart's API querying the Ensembl database.
    Saves results to disk and returns a dataframe.

    The API call to BioMart can take a bit, so the function should take ~4 min for a human genome.

    Parameters
    ----------
    organism : :obj:`str`
        Organism to get annotation for. Currently supported: "human" and "mouse".

    genome_assembly : :obj:`str`, optional
        Ensembl assembly/version to use.
        Default for "human" is "grch38" and for "mouse" is "grcm38".

    save: :obj:`bool`, optional
        Whether to save to disk under ``output_dir``.
        Defaults to :obj:`True`.

    output_dir : :obj:`str`, optional
        Directory to write output to.
        Defaults to "reference" in current directory.

    chr_prefix: :obj:`bool`, optional
        Whether chromosome names should have the "chr" prefix. Defaults to True

    gene_types : :obj:`list`, optional
        Subset of transcript biotypes to keep.
        See here the available biotypes https://www.ensembl.org/Help/Faq?id=468
        Defaults to 'protein_coding', 'processed_transcript', 'lincRNA', 'antisense'.

    overwrite: :obj:`bool`, optional
        Whether existing files should be overwritten by new ones.
        Otherwise they will be kept and no action is made.
        Defaults to :obj:`True`.

    Returns
    -------
    :obj:`pandas.DataFrame`
        DataFrame with genome annotations
    """
    from pybedtools import BedTool
    import pybedtools

    default_organisms = {
        "human": {"species": "hsapiens", "ensembl_version": "grch38", "ucsc_version": "hg38",},
        "mouse": {"species": "mmusculus", "ensembl_version": "grcm38", "ucsc_version": "mm10",},
        "yeast": {"species": "scerevisiae", "ensembl_version": "R64"},
    }

    if genome_assembly is None:
        genome_assembly = default_organisms[organism]["ucsc_version"]
        ensembl_genome_assembly = default_organisms[organism]["ensembl_version"]
    elif genome_assembly == "hg19":
        ensembl_genome_assembly = "grch37"
    elif genome_assembly == "hg38":
        ensembl_genome_assembly = "grch38"
    elif genome_assembly == "mm10":
        ensembl_genome_assembly = "grcm38"
    elif genome_assembly == "sacCer3":
        ensembl_genome_assembly = "R64"
    else:
        _LOGGER.warning()
        ensembl_genome_assembly = genome_assembly

    output_file = os.path.join(
        output_dir, "{}.{}.genomic_context.bed".format(organism, ensembl_genome_assembly)
    )
    o = get_this_file_or_timestamped(output_file)
    if os.path.exists(o) and (not overwrite):
        msg = "Annotation file already exists and 'overwrite' is set to False."
        hint = " Returning existing annotation from file: {}".format(output_file)
        _LOGGER.warning(msg + hint)
        annot = pd.read_csv(output_file, header=None, sep="\t")
        return annot

    if save:
        if output_dir is None:
            output_dir = os.path.join(os.path.abspath(os.path.curdir), "reference")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

    attrs = [
        "ensembl_gene_id",
        "chromosome_name",
        "start_position",
        "end_position",
        "exon_chrom_start",
        "exon_chrom_end",
        "5_utr_start",
        "5_utr_end",
        "3_utr_start",
        "3_utr_end",
        "strand",
        "external_gene_name",
        "gene_biotype",
    ]
    res = query_biomart(
        attributes=attrs,
        species=default_organisms[organism]["species"],
        ensembl_version=ensembl_genome_assembly,
    )

    if chr_prefix:
        res["chromosome_name"] = "chr" + res["chromosome_name"]

    # Keep only desired gene types
    res = res.loc[res["gene_biotype"].isin(gene_types), :]
    # Get strand
    res.loc[:, "strand"] = res["strand"].replace("1", "+").replace("-1", "-")
    # remove exons which are also UTRs

    bed_cols = ["chromosome_name", "start_position", "end_position", "strand"]
    # # Promoters = genebody:start +/- N
    promoter = res[bed_cols].drop_duplicates()
    for col in promoter.columns:
        if ("start" in col) or ("end" in col):
            promoter.loc[:, col] = promoter.loc[:, col].astype(int)

    promoter.loc[promoter["strand"] == "+", "start_position"] = promoter.loc[
        promoter["strand"] == "+", "start_position"
    ] - (promoter_width / 2)
    promoter.loc[promoter["strand"] == "-", "start_position"] = promoter.loc[
        promoter["strand"] == "-", "end_position"
    ] - (
        promoter_width / 2
    )  # + 1
    promoter.loc[promoter["strand"] == "+", "end_position"] = promoter.loc[:, "start_position"] + (
        promoter_width / 2
    )
    for col in bed_cols[1:3]:
        promoter.loc[:, col] = promoter.loc[:, col].astype(int)

    # # # fix where promoter start < 0
    promoter.loc[promoter["start_position"] < 0, "start_position"] = 0
    # # # remove end < start
    promoter = promoter.loc[~(promoter["start_position"] > promoter["end_position"])]

    # # genebody = start->end + promoter
    gb = res[bed_cols].drop_duplicates()
    for col in gb.columns:
        if ("start" in col) or ("end" in col):
            gb.loc[:, col] = gb.loc[:, col].astype(int)
    gb = gb.append(promoter)
    # for col in bed_cols[1:3]:
    #     gb.loc[:, col] = gb.loc[:, col].astype(int)
    gb = gb.sort_values(gb.columns.tolist())
    genebody_bed = BedTool.from_dataframe(gb)
    genebody_bed = genebody_bed.sort().merge()
    genebody = genebody_bed.to_dataframe()

    # # Exon
    exon = res[["chromosome_name", "exon_chrom_start", "exon_chrom_end"]].drop_duplicates().dropna()
    # # 5utr
    utr5 = res[["chromosome_name", "5_utr_start", "5_utr_end"]].drop_duplicates().dropna()
    # # 3utr
    utr3 = res[["chromosome_name", "3_utr_start", "3_utr_end"]].drop_duplicates().dropna()
    for d in [exon, utr5, utr3]:
        for col in d.columns:
            if ("start" in col) or ("end" in col):
                d.loc[:, col] = d.loc[:, col].astype(int)

    # # Introns = genebody - (promoter + exon + 5utr + 3utr)
    promoter = promoter.drop(["strand"], axis=1)
    exon.columns = utr3.columns = utr5.columns = promoter.columns = bed_cols[:3]
    others = promoter.append(exon).append(utr5).append(utr3)
    for col in others.columns[1:]:
        others.loc[:, col] = others.loc[:, col].astype(int)
    intron = genebody_bed.subtract(BedTool.from_dataframe(others))
    intron = intron.to_dataframe()

    # # Intergenic = ChromSizes - genebody
    cs = pybedtools.get_chromsizes_from_ucsc(genome_assembly)
    cs = pd.DataFrame(cs).T.reset_index()
    if not chr_prefix:
        cs["index"] = cs["index"].str.subtract("chr", "")
    cs = BedTool.from_dataframe(cs)
    intergenic = cs.subtract(genebody_bed)
    intergenic = intergenic.to_dataframe()

    # join all
    features = ["promoter", "genebody", "exon", "intron", "utr5", "utr3", "intergenic"]
    annot = list()
    for name in features:
        cur = locals()[name]
        cur = cur.iloc[:, list(range(3))]
        cur.columns = bed_cols[:3]
        cur.loc[:, "region_type"] = name
        if save:
            cur.sort_values(cur.columns.tolist()).drop("region_type", axis=1).to_csv(
                os.path.join(
                    output_dir,
                    "{}.{}.genomic_context.{}.bed".format(organism, ensembl_genome_assembly, name),
                ),
                index=False,
                header=False,
                sep="\t",
            )
        annot.append(cur)

    annot = pd.concat(annot)
    annot = annot.sort_values(annot.columns.tolist())

    # save
    if save:
        annot.to_csv(output_file, index=False, header=False, sep="\t")
    return annot


def get_chromosome_sizes(organism, genome_assembly=None, output_dir=None, overwrite=True):
    """
    Get a file with the sizes of chromosomes in a given organism/genome assembly.
    Saves results to disk and returns a path to a text file.

    Parameters
    ----------
    organism : :obj:`str`
        Organism to get chromosome sizes for.
        Currently supported: "human" and "mouse".

    genome_assembly : :obj:`str`, optional
        Ensembl assembly/version to use.
        Default for "human" is "hg38/grch38" and for "mouse" is "mm10/grcm38".

    output_dir : :obj:`str`, optional
        Directory to write output to.

        Defaults to "reference" in current directory.
    overwrite: :obj:`bool`, optional
        Whether existing files should be overwritten by new ones.
        Otherwise they will be kept and no action is made.

        Defaults to :obj:`True`.

    Returns
    -------
    str
        Path to text file with chromosome sizes.
    """
    import pybedtools

    if output_dir is None:
        output_dir = os.path.join(os.path.abspath(os.path.curdir), "reference")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    default_organisms = {"human": "hg38", "mouse": "mm10"}
    if genome_assembly is None:
        genome_assembly = default_organisms[organism]
        _LOGGER.warning(
            "Genome assembly not selected. Using assembly '{}' for '{}'.".format(
                genome_assembly, organism
            )
        )

    output_file = os.path.join(
        output_dir, "{}.{}.chromosome_sizes.txt".format(organism, genome_assembly)
    )
    if os.path.exists(output_file) and (not overwrite):
        msg = "Annotation file already exists and 'overwrite' is set to False."
        hint = " Returning existing annotation file: {}".format(output_file)
        _LOGGER.warning(msg + hint)
        return output_file

    sizes = pybedtools.get_chromsizes_from_ucsc(genome_assembly)
    with open(output_file, "w") as handle:
        for chrom, (_, size) in sizes.items():
            handle.write("{}\t{}\n".format(chrom, size))
    return output_file


def deseq_analysis(
    count_matrix,
    experiment_matrix,
    comparison_table,
    formula,
    output_dir,
    output_prefix,
    overwrite=True,
    alpha=0.05,
    independent_filtering=False,
    create_subdirectories=True,
    save_inputs=True,
    **kwargs
):
    """
    Perform differential comparison analysis with DESeq2.

    .. note::
        Do not include hyphens ("-") in any of the samples or groups names!
        R freaks out with this.

    # TODO: fix hyphens in names issue

    Parameters
    ----------
    count_matrix : :obj:`pandas.DataFrame`
        Data frame of shape (samples, variables) with raw read counts.

    experiment_matrix : :obj:`pandas.DataFrame`
        Data frame with columns "sample_name" and any
        other variables used in the `formula`.

    comparison_table : :obj:`pandas.DataFrame`
        Data frame with columns "comparison_name", "sample_group" and sample_name".

    formula : :obj:`str`
        Formula to test in R/patsy notation.
        Usually something like "~ batch + group".

    output_dir : :obj:`str`
        Output directory for produced files.

    output_prefix : :obj:`str`
        Prefix to add to produced files.

    overwrite: :obj:`bool`, optional
        Whether files existing should be overwritten. Defaults to :obj:`True`.

    alpha : number, optional
        Significance level to reject null hypothesis.
        This in practice has no effect as results for all features will be returned.
        Defaults to 0.05.

    create_subdirectories: :obj:`bool`
        Whether to create subdirectories for the result of each comparison.

    **kwargs: :obj:`dict`
        Additional keyword arguments to be passed to the DESeq function of DESeq2.

    Returns
    -------
    :obj:`pandas.DataFrame`
        Data frame with results, statistics for each feature.
    """
    from tqdm import tqdm
    from ngs_toolkit.utils import r2pandas_df, recarray2pandas_df

    from rpy2.robjects import numpy2ri, pandas2ri, r
    from rpy2.robjects.packages import importr

    numpy2ri.activate()
    pandas2ri.activate()

    importr("DESeq2")

    # order experiment and count matrices in same way
    experiment_matrix = experiment_matrix.set_index("sample_name").loc[count_matrix.columns, :]

    # save the matrices just in case
    if save_inputs:
        count_matrix.to_csv(os.path.join(output_dir, output_prefix + ".count_matrix.tsv"), sep="\t")
        experiment_matrix.to_csv(
            os.path.join(output_dir, output_prefix + ".experiment_matrix.tsv"), sep="\t"
        )
        comparison_table.to_csv(
            os.path.join(output_dir, output_prefix + ".comparison_table.tsv"), sep="\t"
        )

    # Rename samples to avoid R errors with sample names containing symbols
    count_matrix.columns = ["S{}".format(i) for i in range(len(count_matrix.columns))]
    experiment_matrix.index = ["S{}".format(i) for i in range(len(experiment_matrix.index))]
    # Replace hyphens with underscores
    experiment_matrix = experiment_matrix.replace("-", "_")
    comparison_table = comparison_table.replace("-", "_")

    # Run DESeq analysis
    dds = r.DESeqDataSetFromMatrix(
        countData=count_matrix.astype(int),
        colData=experiment_matrix,
        design=r("as.formula")(formula),
    )
    dds = r.DESeq(dds, parallel=True, **kwargs)
    # _save(dds, file=os.path.join(output_dir, output_prefix + ".deseq_dds_object.Rdata"))

    results = list()
    comps = comparison_table["comparison_name"].drop_duplicates().sort_values()
    for comp in tqdm(comps, total=len(comps), desc="Comparison"):
        if create_subdirectories:
            if not os.path.exists(os.path.join(output_dir, comp)):
                os.makedirs(os.path.join(output_dir, comp))
            out_file = os.path.join(
                output_dir, comp, output_prefix + ".deseq_result.{}.csv".format(comp)
            )
        else:
            out_file = os.path.join(output_dir, output_prefix + ".deseq_result.{}.csv".format(comp))

        if not overwrite and os.path.exists(get_this_file_or_timestamped(out_file)):
            continue
        _LOGGER.info("Doing comparison '{}'".format(comp))
        a = (
            comparison_table.loc[
                (comparison_table["comparison_name"] == comp)
                & (comparison_table["comparison_side"] >= 1),
                "sample_group",
            ]
            .drop_duplicates()
            .squeeze()
        )
        b = (
            comparison_table.loc[
                (comparison_table["comparison_name"] == comp)
                & (comparison_table["comparison_side"] <= 0),
                "sample_group",
            ]
            .drop_duplicates()
            .squeeze()
        )

        contrast = np.array(["sample_group", a, b])
        try:
            res = r("as.data.frame")(
                r.results(
                    dds,
                    contrast=contrast,
                    alpha=alpha,
                    independentFiltering=independent_filtering,
                    parallel=True,
                )
            )
        except Exception as e:
            _LOGGER.warning("DESeq2 group contrast '{}' didn't work!".format(contrast))
            _LOGGER.error(e)
            contrast = ["sample_group" + a, "sample_group" + b]
            try:
                res = r("as.data.frame")(
                    r.results(
                        dds,
                        contrast=contrast,
                        alpha=alpha,
                        independentFiltering=independent_filtering,
                        parallel=True,
                    )
                )
                _LOGGER.warning(
                    "DESeq2 group contrast '{}' did work now!".format(", ".join(contrast))
                )
            except Exception as e2:
                _LOGGER.warning(
                    "DESeq2 group contrast '{}' didn't work either!".format(", ".join(contrast))
                )
                _LOGGER.error(e2)
                raise e2

        if isinstance(res, np.recarray):
            res = recarray2pandas_df(res)
            res.index = count_matrix.index

        if not isinstance(res, pd.DataFrame):
            # convert to pandas dataframe
            res = r2pandas_df(res)

        res.loc[:, "comparison_name"] = comp
        res.index.name = "index"

        # save
        res.sort_values("pvalue").to_csv(out_file)
        # append
        results.append(res)

    # save all
    results = pd.concat(results)
    results.index.name = "index"
    results.to_csv(
        os.path.join(output_dir, output_prefix + ".deseq_result.all_comparisons.csv"), index=True,
    )

    # return
    return results


def least_squares_fit(
    matrix,
    design_matrix,
    test_model,
    null_model="~ 1",
    standardize_data=True,
    multiple_correction_method="fdr_bh",
):
    """
    Fit a least squares model with only categorical predictors.
    Computes p-values by comparing the log likelihood ratio of the chosen model to a `null_model`.

    Parameters
    ----------
    matrix : :obj:`pandas.DataFrame`
        A Data frame of shape (samples, variables).

    design_matrix : :obj:`pandas.DataFrame`
        A Data frame of shape (samples, variables) with all the variables in `test_model`.

    test_model : :obj:`str`
        Model design to test in R/patsy notation.

    null_model : :obj:`str`, optional
        Null model design in R/patsy notation. Defaults to "~ 1".

    standardize_data: :obj:`bool`, optional
        Whether data should be standardized prior to fitting. Defaults to :obj:`True`.

    multiple_correction_method : :obj:`str`, optional
        Method to use for multiple test correction.
        See statsmodels.sandbox.stats.multicomp.multipletests. Defaults to "fdr_bh".

    Returns
    -------
    :obj:`pandas.DataFrame`
        Statistics of model fitting and comparison between models for each feature.

    :Example:

    matrix = np.random.random(10000000).reshape(100, 100000)
    P = np.concatenate([[0] * 50, [1] * 50])  # dependent variable
    Q = np.concatenate([[0] * 25, [1] * 25] + [[0] * 25, [1] * 25])  # covariate
    design_matrix = pd.DataFrame([P, Q], index=["P", "Q"]).T
    matrix = matrix.T * ((1 + design_matrix.sum(axis=1)) * 4).values
    matrix = pd.DataFrame(matrix.T)
    test_model = "~ Q + P"
    null_model = "~ Q"
    res = least_squares_fit(matrix, design_matrix, test_model, null_model)
    res.head()

    """
    import patsy
    from scipy import stats
    from scipy.linalg import lstsq
    from sklearn.preprocessing import StandardScaler
    from statsmodels.sandbox.stats.multicomp import multipletests

    if standardize_data:
        norm = StandardScaler()
        matrix = pd.DataFrame(
            norm.fit_transform(matrix), index=matrix.index, columns=matrix.columns
        )

    a1 = patsy.dmatrix(test_model, design_matrix)
    betas1, residuals1, _, _ = lstsq(a1, matrix)

    a0 = patsy.dmatrix(null_model, design_matrix)
    betas0, residuals0, _, _ = lstsq(a0, matrix)

    results = pd.DataFrame(betas1.T, columns=a1.design_info.column_names, index=matrix.columns)

    # Calculate the log-likelihood ratios
    n = float(matrix.shape[0])
    results["model_residuals"] = residuals1
    results["null_residuals"] = residuals0
    results["model_log_likelihood"] = (
        (-n / 2.0) * np.log(2 * np.pi) - n / 2.0 * np.log(results["model_residuals"] / n) - n / 2.0
    )
    results["null_log_likelihood"] = (
        (-n / 2.0) * np.log(2 * np.pi) - n / 2.0 * np.log(results["null_residuals"] / n) - n / 2.0
    )

    results["log_likelihood_ratio"] = (
        results["model_log_likelihood"] - results["null_log_likelihood"]
    )
    results["D_statistic"] = 2 * results["log_likelihood_ratio"]
    results["p_value"] = stats.chi2.sf(
        results["log_likelihood_ratio"], df=betas1.shape[0] - betas0.shape[0]
    )
    results["q_value"] = multipletests(results["p_value"], method=multiple_correction_method)[1]

    if not standardize_data:
        results["mean"] = matrix.mean(axis=0)

    return results


def differential_from_bivariate_fit(
    comparison_table,
    matrix,
    output_dir,
    output_prefix,
    n_bins=250,
    multiple_correction_method="fdr_bh",
    plot=True,
    palette="colorblind",
    make_values_positive=False,
):
    """
    Perform differential analysis using a bivariate gaussian fit
    on the relationship between mean and fold-change for each comparison.

    Parameters
    ----------
    comparison_table : :obj:`pandas.DataFrame`
        Dataframe with 'comparison_name', 'comparison_side' and 'sample_name', 'sample_group' columns.

    matrix : :obj:`pandas.DataFrame`
        Matrix of `n_features, n_samples` with normalized, log-transformed values to perform analysis on.

    output_dir : :obj:`str`
        Output directory

    output_prefix : :obj:`str`
        Prefix for outputs.

    n_bins : :obj:`int`
        Number of bins of mean values along which to standardize fold-changes.

    multiple_correction_method : :obj:`str`
        Multiple correction method from `statsmodels.sandbox.stats.multicomp.multipletests`.

    plot: :obj:`bool`
        Whether to generate plots.

    palette : :obj:`str`
        Color palette to use. This can be any matplotlib palette and is passed to `sns.color_palette`.

    make_values_positive: :obj:`bool`
        Whether to transform `matrix` to have minimum value 0. Default False.

    Returns
    -------
    :obj:`pandas.DataFrame`
        Results of fitting and comparison between groups for each feature.
    """
    import matplotlib.pyplot as plt
    from ngs_toolkit.graphics import savefig
    from scipy.stats import gaussian_kde
    import seaborn as sns
    from statsmodels.sandbox.stats.multicomp import multipletests

    os.makedirs(output_dir, exist_ok=True)

    comparisons = comparison_table["comparison_name"].drop_duplicates().sort_values()

    if make_values_positive:
        matrix = matrix + abs(matrix.min().min())

    results = list()
    for i, comparison in enumerate(comparisons):
        _LOGGER.info("Doing comparison '{}'".format(comparison))
        out_file = os.path.join(output_dir, output_prefix + ".fit_result.{}.csv".format(comparison))

        comp = comparison_table.query("comparison_name == '{}'".format(comparison))
        sa = comp.query("comparison_side == 1")["sample_name"]
        ga = comp.query("comparison_side == 1")["sample_group"].unique()[0]
        sb = comp.query("comparison_side == 0")["sample_name"]
        gb = comp.query("comparison_side == 0")["sample_group"].unique()[0]
        a = matrix.loc[:, sa].mean(axis=1)
        a.name = ga
        b = matrix.loc[:, sb].mean(axis=1)
        b.name = gb

        # assemble stats
        res = a.to_frame()
        res = res.join(b)
        res["global_mean"] = matrix.mean(axis=1)
        res["global_std"] = matrix.std(axis=1)
        res["comparison_mean"] = res.mean(axis=1)
        res["comparison_std"] = res.mean(axis=1)
        res["log2FoldChange"] = np.log2(res[ga] / res[gb])
        res["comparison_name"] = comparison

        # standardize fold change
        bounds = np.linspace(0, res["comparison_mean"].max(), n_bins)
        for (start, end) in zip(bounds[:-2], bounds[1:-1]):
            r = res.loc[(res["comparison_mean"] > start) & (res["comparison_mean"] < end)].index
            v = res.loc[r, "log2FoldChange"]
            res.loc[r, "norm_log2FoldChange"] = (v - np.nanmean(v)) / np.nanstd(v)

        # let's try a bivariate gaussian kernel
        # separately for positive and negative to avoid biases in center of mass
        kernel = gaussian_kde(
            res.loc[
                res["norm_log2FoldChange"] > 0, ["comparison_mean", "norm_log2FoldChange"],
            ].T.values
        )
        res.loc[res["norm_log2FoldChange"] > 0, "density"] = kernel(
            res.loc[
                res["norm_log2FoldChange"] > 0, ["comparison_mean", "norm_log2FoldChange"],
            ].T.values
        )
        kernel = gaussian_kde(
            res.loc[
                res["norm_log2FoldChange"] <= 0, ["comparison_mean", "norm_log2FoldChange"],
            ].T.values
        )
        res.loc[res["norm_log2FoldChange"] <= 0, "density"] = kernel(
            res.loc[
                res["norm_log2FoldChange"] <= 0, ["comparison_mean", "norm_log2FoldChange"],
            ].T.values
        )

        # Let's calculate something like an empirical p-value on the density
        res["pvalue"] = (res["density"] - res["density"].min()) / (
            res["density"].max() - res["density"].min()
        )
        res["padj"] = multipletests(res["pvalue"].fillna(1), method=multiple_correction_method)[1]

        res["direction"] = (res["norm_log2FoldChange"] >= 0).astype(int).replace(0, -1)
        res.to_csv(out_file)
        results.append(res)

    # save all
    results = pd.concat(results)
    results.to_csv(
        os.path.join(output_dir, output_prefix + ".deseq_result.all_comparisons.csv"), index=True,
    )

    if not plot:
        return
    fig, axis = plt.subplots(
        2,
        len(comparisons),
        figsize=(4 * len(comparisons), 4 * 2),
        sharex=True,
        sharey="row",
        squeeze=False,
    )

    for i, comparison in enumerate(comparisons):
        axis[0, i].scatter(
            res["comparison_mean"],
            res["log2FoldChange"],
            alpha=0.2,
            s=5,
            color=sns.color_palette(palette)[0],
            rasterized=True,
        )
        axis[0, i].axhline(0, color="black", linestyle="--")
        axis[1, i].scatter(
            res["comparison_mean"],
            res["norm_log2FoldChange"],
            alpha=0.2,
            s=5,
            color=sns.color_palette(palette)[0],
            rasterized=True,
        )
        diff = res.loc[(res["pvalue"] < 0.05) & (res["norm_log2FoldChange"].abs() >= 2), :].index
        axis[0, i].scatter(
            res.loc[diff, "comparison_mean"],
            res.loc[diff, "log2FoldChange"],
            alpha=0.2,
            s=5,
            color=sns.color_palette(palette)[1],
            rasterized=True,
        )
        axis[1, i].scatter(
            res.loc[diff, "comparison_mean"],
            res.loc[diff, "norm_log2FoldChange"],
            alpha=0.2,
            s=5,
            color=sns.color_palette(palette)[1],
            rasterized=True,
        )
        axis[1, i].axhline(0, color="black", linestyle="--")
        axis[0, i].set_title(comparison + "\n" + ga + " vs " + gb)
        axis[1, i].set_xlabel("Comparison mean")
        if i == 0:
            axis[0, i].set_ylabel("log2 fold-change")
            axis[1, i].set_ylabel("Norm(log2 fold-change)")

    # save figure
    savefig(
        fig, os.path.join(output_dir, output_prefix + ".deseq_result.all_comparisons.scatter.svg"),
    )
    return results


# def independent_filtering(df, alpha=0.05, n_quantiles=100):
#     """
#     """
#     req_columns = ["pvalue", "baseMean"]
#     assert all([x in df.columns for x in req_columns])
#     # compute quantiles accross mean and pvalue distributions
#     stats = pd.DataFrame()
#     p = (np.arange(n_quantiles) / float(n_quantiles)) * 100
#     p = np.append(p, 100.)
#     for start, end in zip(p, p[1:]):
#         m = np.log2(1 + df['baseMean'])
#         i = df.loc[
#             (m >= np.percentile(m, start)) &
#             (m <= np.percentile(m, end)), :].index
#         stats.loc[start, "n"] = i.shape[0]
#         stats.loc[start, "mean"] = df.loc[i, "baseMean"].mean()
#         stats.loc[start, "mean_p"] = df.loc[i, "pvalue"].mean()
#         stats.loc[start, "n_sig_p"] = (df.loc[i, "pvalue"] < alpha).sum()
#     # plot
#     fig, axis = plt.subplots(1, 2, figsize=(4 * 2, 4 * 1))
#     axis[0].scatter(stats.index, stats.loc[:, "n_sig_p"])
#     axis[1].scatter(stats.index, -np.log10(stats.loc[:, "mean_p"]))
#     # choose inflection point
#     return


# def fit_curve():
#     def func(x, a, b, c):
#         return a * np.exp(-b * x) + c
#     plt.plot(stats.index, -np.log10(stats['mean_p']), 'b-', label='data')
#     popt, pcov = curve_fit(func, stats.index, -np.log10(stats['mean_p']))
#     plt.plot(stats.index, func(stats.index, *popt), 'r-', label='fit')
#     popt, pcov = curve_fit(func, stats.index, stats['mean_p'], bounds=(0, [3., 2., 1.]))
#     plt.plot(stats.index, func(stats.index, *popt), 'g--', label='fit-with-bounds')
#     plt.xlabel('x')
#     plt.ylabel('y')
#     plt.legend()
#     plt.show()


def lola(bed_files, universe_file, output_folder, genome, output_prefixes=None, cpus=8):
    """
    Perform location overlap analysis (LOLA).

    If bed_files is a list with more than one element, use ``output_prefixes``
    to pass a list of prefixes to label the output files for each input BED file.

    Files will be created in ``output_folder`` mimicking the output that the
    R function LOLA::writeCombinedEnrichment writes.

    Requires the R package "LOLA" to be installed:

    .. highlight:: R
    .. code-block:: R

        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("LOLA")

    Parameters
    ----------
    bed_files : {str, list}
        A string path to a BED file or a list of paths.

    universe_file : :obj:`str`
        A path to a BED file representing the universe from where the BED
        file(s) come from.

    output_folder : :obj:`str`
        Output folder for resulting files.

    genome : :obj:`str`, optional
        Genome assembly from which the BED files come from.
        This is used to get the LOLA databases from the
        :obj:`ngs_toolkit._CONFIG` parameters.

    output_prefixes : :obj:`list`, optional
        A list of strings with prefixes to be used in case
        ``bed_files`` is a list.

    cpus : :obj:`int`, optional
        Number of CPUs/threads to use.

        Defaults to 8.
    """
    from ngs_toolkit.utils import r2pandas_df
    from rpy2.robjects import numpy2ri, pandas2ri, r
    from rpy2.robjects.packages import importr

    numpy2ri.activate()
    pandas2ri.activate()

    importr("LOLA")

    # Get region databases from config
    _LOGGER.info("Getting LOLA databases for genome '%s' from configuration.", genome)

    msg = (
        "LOLA database values in configuration could not be found or understood. "
        "Please add a list of value(s) to this section "
        "'resources:lola:region_databases:%s'. "
        "For an example, see "
        "https://github.com/afrendeiro/toolkit/tree/master/ngs_toolkit/config/example.yaml"
    )
    try:
        databases = _CONFIG["resources"]["lola"]["region_databases"][genome]
    except KeyError:
        _LOGGER.error(msg, genome)
        raise

    if not isinstance(databases, list):
        if isinstance(databases, str):
            databases = list(databases)
        else:
            _LOGGER.error(msg, genome)
            raise KeyError(msg % genome)

    if len(databases) < 1:
        _LOGGER.error(msg)
        raise KeyError(msg)

    if isinstance(bed_files, str):
        bed_files = [bed_files]
    if output_prefixes is None:
        if len(bed_files) > 1:
            msg = (
                "Running more than one BED file at once while only specifying "
                "`output_folder` argument will cause output files to be named "
                "in the form "
                "'{output_folder}/{region_database}.{input_file}.tsv'."
                " To prevent this behaviour, pass a list of arguments to "
                "`output_prefixes`."
            )
            _LOGGER.warning(msg)
            output_prefixes = [
                rr.replace(os.path.sep, "__").replace(".bed", ".") for rr in bed_files
            ]
        else:
            output_prefixes = ["."]

    _LOGGER.info("Reading up universe file '{}'.".format(universe_file))
    universe = r("LOLA::readBed")(universe_file)
    _LOGGER.info("Loading region set databases.")
    _regionDB = r("LOLA::loadRegionDB")(np.array(databases))
    for suffix, bed_file in zip(output_prefixes, bed_files):
        _LOGGER.info("Reading up BED file '{}'.".format(bed_file))
        user_set = r("LOLA::readBed")(bed_file)
        _LOGGER.info("Running LOLA testing for file '{}'.".format(bed_file))
        _lola_results = r("LOLA::runLOLA")(user_set, universe, _regionDB, cores=cpus)
        _LOGGER.info("Converting results from R to Python")
        lola_results = r2pandas_df(_lola_results)
        _LOGGER.info("Saving all results for file '{}'.".format(bed_file))
        lola_results.to_csv(
            os.path.join(output_folder, "allEnrichments" + suffix + "tsv"), index=False, sep="\t",
        )
        for region_set in lola_results["collection"].drop_duplicates():
            _LOGGER.info("Saving results for collection '%s' only.", region_set)
            lola_results[lola_results["collection"] == region_set].to_csv(
                os.path.join(output_folder, "col_" + region_set + suffix + "tsv"),
                index=False,
                sep="\t",
            )


def meme_ame(
    input_fasta, output_dir, background_fasta=None, organism="human", motif_database_file=None,
):
    import subprocess

    if motif_database_file is None:
        # Get region databases from config
        _LOGGER.info("Getting 2bit reference genome for genome '%s' from configuration.", organism)

        msg = (
            "Reference genome in 2bit format value in configuration could not"
            "be found or understood. "
            "Please add a list of value(s) to this section "
            "'resources:meme:motif_databases:%s'. "
            "For an example, see "
            "https://github.com/afrendeiro/toolkit/tree/master/ngs_toolkit/config/example.yaml"
        )
        try:
            motif_database_file = _CONFIG["resources"]["meme"]["motif_databases"][organism]
        except KeyError:
            _LOGGER.error(msg, organism)
            return

        if not isinstance(motif_database_file, str):
            _LOGGER.error(msg, organism)
            return

    # shuffle input in no background is provided
    if background_fasta is None:
        shuffled = input_fasta + ".shuffled"
        cmd = """
        fasta-dinucleotide-shuffle -c 1 -f {0} > {1}
        """.format(
            input_fasta, shuffled
        )
        subprocess.call(cmd.split(" "))

    cmd = (
        "ame --bgformat 1 --scoring avg --method ranksum "
        "--pvalue-report-threshold 0.05 --control {0} -o {1} {2} {3}"
    ).format(
        background_fasta if background_fasta is not None else shuffled,
        output_dir,
        input_fasta,
        motif_database_file,
    )
    subprocess.call(cmd.split(" "))
    # subprocess.call("rm {}".format(shuffled).split(" "))


def homer_motifs(bed_file, output_dir, genome_assembly):
    import subprocess

    default_opts = "-size 1000 -h -p 2 -len 8,10,12,14 -noknown"
    cmd = "findMotifsGenome.pl {bed} {genome}r {out_dir} {opts}".format(
        bed=get_this_file_or_timestamped(bed_file),
        genome=genome_assembly,
        out_dir=output_dir,
        opts=default_opts,
    )
    subprocess.call(cmd.split(" "))


def homer_combine_motifs(
    comparison_dirs,
    output_dir,
    region_prefix="differential_analysis",
    reduce_threshold=0.6,
    match_threshold=10,
    info_value=0.6,
    p_value_threshold=1e-25,
    fold_enrichment=None,
    cpus=8,
    run=True,
    distributed=True,
    genome="hg38",
    motif_database=None,
):
    """
    Create consensus of de novo discovered motifs from HOMER

    Parameters
    ----------
    comparison_dirs : :obj:`list`
        Iterable of comparison directories where homer was run.
        Should contain a "homerMotifs.all.motifs" file.

    output_dir : :obj:`str`
        Output directory.

    p_value_threshold : number, optional
        Threshold for inclusion of a motif in the consensus set.
        Defaults to 1e-5

    cpus : number, optional
        Number of available CPUS/threads for multithread processing.
        Defaults to 8

    run: :obj:`bool`, optional
        Whether to run enrichment of each comparison in the consensus motifs.
        Default is True

    distributed: :obj:`bool`, optional
        Whether to run enrichment as a cluster job.
        Default is True

    genome : :obj:`str`
        Genome assembly of the data.
        Default is 'hg38'.

    motif_database : :obj:`str`
        Motif database to restrict motif matching too.

    Returns
    -------
    {str,None}
        If `run` is `False`, returns path to consensus motif file. Otherwise `None`.
    """
    from glob import glob
    import subprocess

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # concatenate files
    out_file = os.path.join(output_dir, "homerMotifs.combined.motifs")
    open(out_file, "w")  # make sure file is empty in the beginning
    with open(out_file, "a") as outfile:
        for dir_ in comparison_dirs:
            with open(os.path.join(dir_, "homerMotifs.all.motifs"), "r") as infile:
                for line in infile:
                    outfile.write(line)

    # Filter and get motif consensus
    extra = ""
    if motif_database:
        extra = " -known {}".format(motif_database)
    if fold_enrichment is None:
        fold_enrichment = ""
    else:
        fold_enrichment = " -F " + str(fold_enrichment)
    subprocess.call(
        "compareMotifs.pl {} {} -reduceThresh {}"
        " -matchThresh {}{} -pvalue {}"
        " -info {}{} -nofacts -cpu {}".format(
            out_file,
            output_dir,
            reduce_threshold,
            match_threshold,
            extra,
            p_value_threshold,
            info_value,
            fold_enrichment,
            cpus,
        ).split(" ")
    )

    # concatenate consensus motif files
    files = glob(os.path.join(output_dir, "homerResults/*motif"))
    combined_motifs = os.path.join(output_dir, "homerMotifs.filtered.motifs")
    with open(combined_motifs, "w") as outfile:
        for f in files:
            if ("similar" in f) or ("RV" in f):
                continue
            _LOGGER.info(f)
            with open(f, "r") as infile:
                for line in infile:
                    outfile.write(line)

    if run:
        for dir_ in comparison_dirs:
            # delete previous results if existing
            results = os.path.join(dir_, "knownResults.txt")
            if os.path.exists(results):
                _LOGGER.warning("Deleting previously existing results file: '{}'".format(results))
                os.remove(results)
            # prepare enrichment command with consensus set
            cmd = "findMotifsGenome.pl {bed} {genome}r {dir} -p {cpus} -nomotif -mknown {motif_file}".format(
                bed=get_this_file_or_timestamped(
                    os.path.join(dir_, region_prefix + "_regions.bed")
                ),
                genome=genome,
                cpus=cpus,
                dir=dir_,
                motif_file=combined_motifs,
            )
            # run
            # TODO: send this as a job through submit_job
            if distributed:
                subprocess.call(
                    "sbatch -J homer.{d} -o {dir}.homer.log -p shortq -c 8 --mem 20000".format(
                        d=os.path.basename(dir_), dir=dir_
                    ).split(" ")
                    + ["--wrap", cmd]
                )
            else:
                subprocess.call(cmd.split(" "))


@MEMORY.cache
def enrichr(gene_set, gene_set_libraries=None, kind="genes", max_attempts=5):
    """
    Use Enrichr on a list of genes (currently only genes supported through the API).
    If input contains <NaN> values, they will be ignored.

    The output of this function is cached to disk using joblib.
    To clear the cache of this function please remove the contents of
    `~/.ngs_toolkit/joblib/general/enrichr` or call
    :func:`ngs_toolkit.MEMORY.clear()` to remove all ngs_toolkit cached values.

    Parameters
    ----------
    gene_set : {:obj:`list`, :obj:`set`, :obj:`~pandas.Series`, :obj:`~pandas.DataFrame`}
        Set of genes to get enrichment for.
        Can be pandas DataFrame with single column or column "gene_name".

    gene_set_libraries : :obj:`list`, optional
        Gene set libraries to use.
        Defaults to values in initial configuration file.
        To see them, do: ``ngs_toolkit._CONFIG['resources']['enrichr']['gene_set_libraries']``

    kind : :obj:`str`, optional
        Type of input.
        Right now, only "genes" is supported.
        Defaults to "genes"

    max_attempts : :obj:`int`, optional
        Number of times to try a call to Enrichr API.
        Defaults to 5

    Returns
    -------
    :obj:`pandas.DataFrame`
        Results of enrichment analysis

    Raises
    -------
    Exception
        If `max_attempts` is exceeded
    """
    import json
    import requests

    from tqdm import tqdm

    ENRICHR_ADD = "http://amp.pharm.mssm.edu/Enrichr/addList"
    ENRICHR_RETRIEVE = "http://amp.pharm.mssm.edu/Enrichr/enrich"
    query_string = "?userListId={}&backgroundType={}"

    if gene_set_libraries is None:
        # Get region databases from config
        _LOGGER.debug("Getting Enrichr gene set libraries from configuration.")

        msg = "Enrichr gene set libraries value in configuration could not be found or understood. "
        msg += "Please add a list of value(s) to this section 'resources:mem:motif_databases'. "
        msg += "For an example, see https://github.com/afrendeiro/toolkit/tree/master/ngs_toolkit/config/example.yaml"
        try:
            gene_set_libraries = _CONFIG["resources"]["enrichr"]["gene_set_libraries"]
        except KeyError:
            _LOGGER.error(msg)
            return

        if not isinstance(gene_set_libraries, list):
            if isinstance(gene_set_libraries, str):
                gene_set_libraries = list(gene_set_libraries)
            else:
                _LOGGER.error(msg)
                return

        if len(gene_set_libraries) < 1:
            _LOGGER.error(msg)
            return

    # Handle various input types
    if isinstance(gene_set, (list, set)):
        genes = [x for x in gene_set if not pd.isnull(x)]
    elif isinstance(gene_set, pd.Series):
        genes = gene_set.dropna().tolist()
    elif isinstance(gene_set, pd.DataFrame):
        if gene_set.shape[1] == 1:
            genes = gene_set.squeeze().dropna().tolist()
        else:
            try:
                genes = gene_set["gene_name"].dropna().tolist()
            except KeyError:
                msg = "Provided input dataframe does not have a gene_name column!"
                _LOGGER.error(msg)
                raise ValueError(msg)

    if kind == "genes":
        # Build payload with bed file
        attr = "\n".join(genes)
    elif kind == "regions":
        raise NotImplementedError
        # Build payload with bed file
        # attr = "\n".join(
        #     dataframe[["chrom", "start", "end"]]
        #     .apply(lambda x: "\t".join([str(i) for i in x]), axis=1)
        #     .tolist()
        # )

    payload = {"list": (None, attr), "description": (None, "")}
    # Request adding gene set
    response = requests.post(ENRICHR_ADD, files=payload)
    if not response.ok:
        raise Exception("Error analyzing gene list")

    # Track gene set ID
    user_list_id = json.loads(response.text)["userListId"]

    results = pd.DataFrame()
    for gene_set_library in tqdm(
        gene_set_libraries, total=len(gene_set_libraries), desc="Gene set library"
    ):
        _LOGGER.debug("Using Enricher on {} gene set library.".format(gene_set_library))

        # Request enriched sets in gene set
        i = 0
        okay = False
        while not okay:
            if i == max_attempts:
                _LOGGER.error("API request failed {} times, exceeding `max_attempts`.".format(i))
                raise Exception("Fetching enrichment results maxed `max_attempts`.")
            response = requests.get(
                ENRICHR_RETRIEVE + query_string.format(user_list_id, gene_set_library)
            )
            okay = response.ok
            if not okay:
                _LOGGER.warning("API request failed. Retrying.")

        # Get enriched sets in gene set
        res = json.loads(response.text)
        # If there's no enrichemnt, continue
        if len(res) < 0:
            continue

        # Put in dataframe
        res = pd.DataFrame([pd.Series(s) for s in res[gene_set_library]])
        if res.shape[0] == 0:
            continue
        cols = [
            "rank",
            "description",
            "p_value",
            "z_score",
            "combined_score",
            "genes",
            "adjusted_p_value",
        ]
        if len(res.columns) == 7:
            res.columns = cols
        elif len(res.columns) == 9:
            res.columns = cols + ["old_p_value", "old_adjusted_p_value"]

        # Remember gene set library used
        res["gene_set_library"] = gene_set_library

        # Append to master dataframe
        results = results.append(res, ignore_index=True)

    return results


def run_enrichment_jobs(
    results_dir,
    genome,
    background_bed,
    steps=["lola", "meme", "homer", "enrichr"],
    overwrite=True,
    pep_config=None,
):
    """
    Submit parallel enrichment jobs for a specific analysis.

    Parameters
    ----------
    :param results_dir:
        Directory with files prepared by ngs_toolkit.general.run_enrichment_jobs

    :param genome:
        Genome assembly of the analysis.

    background_bed : :obj:`str`
        BED file to use as background for LOLA analysis.
        Typically the analysis' own consensus region set.

    steps : :obj:`list`, optional
        Steps of the analysis to perform.
        Defaults to ["region", lola", "meme", "homer", "enrichr"].

    :param overwrite: bool, optional
        Whether output should be overwritten.
        In this case no jobs will be submitted for jobs with existing output files.
        Defaults to True

    :param pep_config: :obj:`str`, optional
        Pickle file of the analysis.
        Only required for "region" enrichment.
    """
    # TODO: replace hardcoded paths with info from resources
    # TODO: remove pep_config requirement to "region_enrichment"
    import sys
    from glob import glob

    from ngs_toolkit.utils import submit_job

    dbs = {
        "human": "~/resources/motifs/motif_databases/HUMAN/HOCOMOCOv10.meme",
        "mouse": "~/resources/motifs/motif_databases/MOUSE/uniprobe_mouse.meme",
    }
    omap = {"hg38": "human", "hg19": "human", "mm10": "mouse"}

    jobs = list()
    # list of tuples with: job_name, log, exec, requirements (partition, cpu, mem, time), cmd

    # REGION
    if "region" in steps:
        files = glob(results_dir + "/*/*_regions*bed")
        for file in files:
            dir_ = os.path.dirname(file)
            name = os.path.basename(dir_)
            output_ = os.path.join(dir_, "region_type_enrichment*csv")
            if os.path.exists(output_) and (not overwrite):
                continue
            jobs.append(
                [
                    name + "_region",
                    os.path.join(dir_, name + ".region.log"),
                    os.path.join(dir_, name + ".region.sh"),
                    ("shortq", 1, 8000, "08:00:00"),
                    "{} -m ngs_toolkit.recipes.region_enrichment --output-file {} {} {}".format(
                        sys.executable, output_, file, pep_config
                    ),
                ]
            )

    # LOLA
    if "lola" in steps:
        # the star here is to support timestamped files
        files = glob(results_dir + "/*/*_regions*bed")
        for file in files:
            dir_ = os.path.dirname(file)
            name = os.path.basename(dir_)
            output_ = os.path.join(dir_, "allEnrichments*tsv")
            if os.path.exists(output_) and (not overwrite):
                continue
            jobs.append(
                [
                    name + "_lola",
                    os.path.join(dir_, name + ".lola.log"),
                    os.path.join(dir_, name + ".lola.sh"),
                    ("shortq", 2, 12000, "08:00:00"),
                    "{} -m ngs_toolkit.recipes.lola {} {} {} {} -c 2".format(
                        sys.executable, file, background_bed, dir_, genome
                    ),
                ]
            )

    # AME
    if "meme" in steps:
        files = glob(results_dir + "/*/*_regions.fa")
        for file in files:
            dir_ = os.path.dirname(file)
            name = os.path.basename(dir_)
            output_ = os.path.join(dir_, "ame.html")
            if os.path.exists(output_) and (not overwrite):
                continue
            jobs.append(
                [
                    name + "_meme",
                    os.path.join(dir_, name + ".meme_ame.log"),
                    os.path.join(dir_, name + ".meme_ame.sh"),
                    ("shortq", 1, 4000, "08:00:00"),
                    "fasta-dinucleotide-shuffle -c 1 -f {f} > {f}.shuffled.fa\n".format(f=file)
                    + " ame --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 0.05"
                    + " --control {f}.shuffled.fa -o {d} {f} {motifs}".format(
                        f=file, d=dir_, motifs=dbs[omap[genome]]
                    ),
                ]
            )

    # HOMER
    if "homer" in steps:
        files = glob(results_dir + "/*/*_regions*bed")
        for file in files:
            dir_ = os.path.dirname(file)
            name = os.path.basename(dir_)
            output_ = os.path.join(dir_, "homerResults.html")
            if os.path.exists(output_) and (not overwrite):
                continue
            jobs.append(
                [
                    name + "_homer",
                    os.path.join(dir_, name + ".homer.log"),
                    os.path.join(dir_, name + ".homer.sh"),
                    ("shortq", 8, 12000, "08:00:00"),
                    "findMotifsGenome.pl {f} {genome}r {d} -size 1000 -h -p 2 -len 8,10,12,14 -noknown".format(
                        f=file, d=dir_, genome=genome
                    ),
                ]
            )

    # Enrichr
    if "enrichr" in steps:
        files = glob(results_dir + "/*/*.gene_symbols.txt")
        for file in files:
            dir_ = os.path.dirname(file)
            name = os.path.basename(dir_)
            output_ = file.replace(".gene_symbols.txt", ".enrichr.csv")
            if os.path.exists(get_this_file_or_timestamped(output_)) and (not overwrite):
                continue
            jobs.append(
                [
                    name + "_enrichr",
                    os.path.join(dir_, name + ".enrichr.log"),
                    os.path.join(dir_, name + ".enrichr.sh"),
                    ("shortq", 1, 4000, "08:00:00"),
                    "{e} -m ngs_toolkit.recipes.enrichr {f} {o}".format(
                        e=sys.executable, f=file, o=output_
                    ),
                ]
            )

    for jobname, log_file, job_file, (partition, cores, mem, time), task in jobs:
        submit_job(
            task,
            job_file,
            log_file=log_file,
            jobname=jobname,
            partition=partition,
            cores=cores,
            mem=mem,
            time=time,
        )


def project_to_geo(
    project,
    output_dir="geo_submission",
    steps=["bam", "bigwig", "peaks"],
    samples=None,
    distributed=False,
    dry_run=False,
    **kwargs
):
    """
    Prepare raw sequencing files for submission to GEO.
    Files will be copied or generated in a new directory ``output_dir``.
    It will get the raw BAM file(s) of each sample, and in case of
    ATAC-seq/ChIP-seq samples, the bigWig and peak files. If multiple BAM
    files exist for each sample, all will be copied and sequencially named
    with the "fileN" suffix, where "N" is the file number.

    For each copied file a md5sum will be calculated.

    A pandas DataFrame with info on the sample's files and md5sums will be returned.

    Attributes
    ----------
    project : :obj:`peppy.Project`
        A :class:`peppy.Project` object to process.

    output_dir : :obj:`str`, optional
        Directory to create output. Will be created/overwriten if existing.

        Defaults to "geo_submission".
    samples : :obj:`list`, optional
        List of :class:`peppy.Sample` objects in project to restrict to.

        Defaults to all samples in project.
    distributed: :obj:`bool`, optional
        Whether processing should be distributed as jobs in a
        computing cluster for each sample.
        Currently available implementation supports a SLURM cluster only.

        Defaults is :obj:`False`.
    dry_run: :obj:`bool`, optional
        Whether copy/execution/submisison commands should be not be run, to test.

        Default is :obj:`False`.
    **kwargs : :obj:`dict`
        Additional keyword arguments will be passed to
        :meth:`ngs_toolkit.utils.submit_job` if ``distributed`` is :obj:`True`,
        and on to a :class:`divvy` submission template.
        Pass for example:

            * computing_configuration="slurm"
            * jobname="job"
            * cores=2
            * mem=8000
            * partition="longq"

    Returns
    ----------
    :class:`pandas.DataFrame`
        Annotation of samples and their BAM, BigWig,
        narrowPeak files and respective md5sums.
    """
    from ngs_toolkit.utils import submit_job

    output_dir = os.path.abspath(output_dir)
    if samples is None:
        samples = project.samples
    os.makedirs(output_dir, exist_ok=True)

    annot = pd.DataFrame(index=pd.Index([], name="sample_name"))
    for sample in samples:
        various = len(sample.data_source.split(" ")) > 1
        cmd = ""
        if "bam" in steps:
            for i, file in enumerate(sample.data_source.split(" ")):
                suffix = ".file{}".format(i) if various else ""
                # Copy raw file
                bam_file = os.path.join(output_dir, sample.name + "{}.bam".format(suffix))
                cmd += "cp {} {};\n".format(file, bam_file)
                cmd += "chmod 644 {};\n".format(bam_file)
                annot.loc[sample.name, "bam_file{}".format(i)] = bam_file

                # Copy or generate md5sum
                md5_file = bam_file + ".md5"
                if os.path.exists(file + ".md5"):
                    cmd += "cp {} {};\n".format(file + ".md5", md5_file)
                else:
                    cmd += "md5sum {} > {};\n".format(bam_file, md5_file)
                cmd += "chmod 644 {};\n".format(md5_file)
                annot.loc[sample.name, "bam_file{}_md5sum".format(i)] = md5_file

        # Copy bigWig files
        if sample.protocol in ["ATAC-seq", "ChIP-seq"]:
            if "bigwig" in steps:
                if hasattr(sample, "bigwig"):
                    bigwig_file = os.path.join(output_dir, sample.name + ".bigWig")
                    cmd += "cp {} {};\n".format(sample.bigwig, bigwig_file)
                    cmd += "chmod 644 {};\n".format(bigwig_file)
                    annot.loc[sample.name, "bigwig_file"] = bigwig_file

                    # Copy or generate md5sum
                    md5_file = bigwig_file + ".md5"
                    if os.path.exists(sample.bigwig + ".md5"):
                        cmd += "cp {} {};\n".format(sample.bigwig + ".md5", md5_file)
                    else:
                        cmd += "md5sum {} > {};\n".format(bigwig_file, md5_file)
                    cmd += "chmod 644 {};\n".format(md5_file)
                    annot.loc[sample.name, "bigwig_file_md5sum"] = md5_file
                else:
                    _LOGGER.warning(
                        "'{}' sample '{}' does not have a 'bigwig'".format(
                            sample.protocol, sample.name
                        )
                        + " attribute set. Skipping bigWig file."
                    )
        # Copy peaks
        if sample.protocol == "ATAC-seq":
            if "peaks" in steps:
                if hasattr(sample, "peaks"):
                    peaks_file = os.path.join(output_dir, sample.name + ".peaks.narrowPeak")
                    cmd += "cp {} {};\n".format(sample.peaks, peaks_file)
                    cmd += "chmod 644 {};\n".format(peaks_file)
                    annot.loc[sample.name, "peaks_file"] = peaks_file

                    # Copy or generate md5sum
                    md5_file = peaks_file + ".md5"
                    if os.path.exists(sample.peaks + ".md5"):
                        cmd += "cp {} {};\n".format(sample.peaks + ".md5", md5_file)
                    else:
                        cmd += "md5sum {} > {};\n".format(peaks_file, md5_file)
                    cmd += "chmod 644 {};\n".format(md5_file)
                    annot.loc[sample.name, "peaks_file_md5sum"] = md5_file
                else:
                    _LOGGER.warning(
                        "'{}' sample '{}' does not have a 'peaks' attribute set.".format(
                            sample.protocol, sample.name
                        )
                        + " Skipping peaks file."
                    )

        # Assemble job
        job_name = "project_to_geo.{}".format(sample.name)
        log_file = os.path.join(output_dir, job_name + ".log")
        job_file = os.path.join(output_dir, job_name + ".sh")
        submit_job(
            cmd,
            job_file,
            log_file=log_file,
            jobname=job_name,
            cores=1,
            mem=8000,
            dry_run=dry_run,
            **kwargs
        )

    return annot


def rename_sample_files(
    annotation_mapping,
    old_sample_name_column="old_sample_name",
    new_sample_name_column="new_sample_name",
    tmp_prefix="rename_sample_files",
    results_dir="results_pipeline",
    dry_run=False,
):
    """
    Rename existing directories with pipeline outputs for samples
    based on mapping of old/new sample names.

    All files within the directory with the old sample name will be renamed recursively.
    Old and new sample names can overlap - this procedure will handle these cases correctly
    by a 2-step process with temporary sample names with prefix ``tmp_prefix``.

    Attributes
    ----------
    annotation_mapping : :obj:`pandas.DataFrame`
        DataFrame with mapping of old (column "previous_sample_name")
        vs new ("new_sample_name") sample names.

    old_sample_name_column : :obj:`str`, optional
        Name of column with old sample names.

        Defaults to "old_sample_name".
    new_sample_name_column : :obj:`str`, optional
        Name of column with new sample names.

        Defaults to "new_sample_name".
    tmp_prefix : :obj:`str`, optional
        Prefix for temporary files to avoid overlap between old and new names.

        Defaults to "rename_sample_files".
    results_dir : :obj:`str`, optional
        Pipeline output directory containing sample output directories.
        This is usually the `data_dir` attribute of Analysis objects.

        Defaults to "results_pipeline".
    dry_run: :obj:`bool`, optional
        Whether to print commands instead of running them.

        Defaults to :obj:`False`.
    """

    def find_replace(from_pattern, to_pattern, root_dir=None, dry_run=False):
        import shutil

        root_dir = os.path.abspath(os.path.curdir if root_dir is None else root_dir)
        os.chdir(root_dir)
        for file in os.listdir(root_dir):  # parse through file list in the current directory
            if from_pattern in file:  # if pattern is found
                new_file = file.replace(from_pattern, to_pattern)
                if not dry_run:
                    shutil.move(file, new_file)  # rename
                else:
                    print(file, new_file)
                file = new_file
            if os.path.isdir(file):
                os.chdir(file)
                find_replace(from_pattern, to_pattern, ".")
                os.chdir("..")

    # 1) move to tmp name
    for i, series in annotation_mapping.iterrows():
        o = series[old_sample_name_column]
        t = "{}-{}".format(tmp_prefix, i)
        _LOGGER.debug("# Moving old sample '%s' to temporary name '%s'." % (o, t))
        find_replace(o, t, root_dir=results_dir, dry_run=dry_run)

    # 2) move to new name
    for i, series in annotation_mapping.iterrows():
        t = "{}-{}".format(tmp_prefix, i)
        n = series[new_sample_name_column]
        _LOGGER.debug("# Moving temporary named sample '%s' to '%s'." % (t, n))
        find_replace(t, n, root_dir=results_dir, dry_run=dry_run)


@MEMORY.cache
def query_biomart(attributes=None, species="hsapiens", ensembl_version="grch38", max_api_retries=5):
    """
    Query Biomart for gene attributes (https://www.ensembl.org/biomart/martview/).

    Available attributes can be see at Biomart. Each should be snakecase
    (lowecase with only underscore separators) as defined in the Biomart API.

    Example attributes:

        * "ensembl_gene_id"
        * "external_gene_name"
        * "hgnc_symbol"

    If a certain field contains multiple values separated by commas,
    it will attemp to return dataframe but it might fail.

    The output of this function is cached to disk using joblib.
    To clear the cache of this function please remove the contents of
    `~/.ngs_toolkit/joblib/general/query_biomart` or call
    :func:`ngs_toolkit.MEMORY.clear()` to remove all ngs_toolkit cached values.

    Parameters
    ----------
    attributes : :obj:`list`, optional
        List of ensembl atrributes to query.

        Defaults to ["ensembl_gene_id", "external_gene_name", "hgnc_id", "hgnc_symbol"].
    species : :obj:`str`, optional
        Ensembl string of species to query. Must be vertebrate.

        Defaults to "hsapiens".
    ensembl_version : :obj:`str`, optional
        Ensembl version to query. Currently "grch37", "grch38" and "grcm38" are tested.

        Defaults to "grch38".
    max_api_retries : :obj:`int`, optional
        How many times to try .

        Defaults to 5.

    Returns
    -------
    :obj:`pandas.DataFrame`
        Dataframe with required attributes for each entry.


    Raises
    ------
    :obj:`ngs_toolkit.exceptions.NetworkError:
        If API call to Enrichr is unsusscessfull for `max_api_retries`.
    """
    import requests
    import time

    supported = ["grch37", "grch38", "grcm38"]
    if ensembl_version not in supported:
        msg = "Ensembl version might not be supported."
        msg += " Tested versions are '{}'.".format("','".join(supported))
        msg += " Will try anyway."
        _LOGGER.warning(msg)

    if attributes is None:
        attributes = ["ensembl_gene_id", "external_gene_name", "hgnc_id", "hgnc_symbol"]

    # Build request XML
    ens_ver = "" if ensembl_version.endswith("38") else ensembl_version + "."
    url_query = "".join(
        [
            """http://{}ensembl.org/biomart/martservice?query=""".format(ens_ver),
            """<?xml version="1.0" encoding="UTF-8"?>""",
            """<!DOCTYPE Query>""",
            """<Query  virtualSchemaName="default" formatter="CSV" header="0" """
            """uniqueRows="0" count="" datasetConfigVersion="0.6" >""",
            """<Dataset name="{}_gene_ensembl" interface="default" >""".format(species),
        ]
        + ["""<Attribute name="{}" />""".format(attr) for attr in attributes]
        + ["""</Dataset>""", """</Query>"""]
    )

    n_fails = 0
    success = False
    while not success:
        req = requests.get(url_query, stream=True)
        if not req.ok:
            n_fails += 1
            msg = "Request to Biomart API was not successful."
            if n_fails == max_api_retries:
                _LOGGER.error(msg)
                raise NetworkError(msg)
            else:
                _LOGGER.warning(msg + " Retrying.")
                time.sleep(10)
        else:
            try:
                content = list(req.iter_lines())
                success = True
            except (
                requests.exceptions.ChunkedEncodingError,
                requests.urllib3.exceptions.ProtocolError,
            ):
                msg = "Retrieving Biomart request failed."
                n_fails += 1
                if n_fails == max_api_retries:
                    _LOGGER.error(msg)
                    raise NetworkError(msg)
                else:
                    _LOGGER.warning(msg + " Retrying.")
                    time.sleep(10)

    if isinstance(content[0], bytes):
        content = [x.decode("utf-8") for x in content]

    if (len(content) == 1) and (content[0].startswith("Query ERROR")):
        msg = (
            "Request to Biomart API was not successful. " "Please check your input: \n\t{}"
        ).format(content[0])
        # _LOGGER.error(msg)
        raise ValueError(msg)

    try:
        mapping = pd.DataFrame([x.strip().split(",") for x in content], columns=attributes)
    except AssertionError as e:
        msg = "Could not simply return dataframe with results."
        msg += (
            " It is likely this is because one of the requested attributes has commas as is quoted."
        )
        msg += " Or because querying an organism not present in vertebrate database."
        msg += " Will try to replace these fields and parse."
        _LOGGER.warning(msg)
        # input probably contains commas inside fields
        c = pd.Series([x.strip() for x in content])
        # well's assume that fields with commas have been quoted
        # get text inside double quotes
        cc = c.str.extract('"(.*)"')
        if cc.shape[1] > 1:
            _LOGGER.error("Attempt to fix failed.")
            raise e
        cc = cc.squeeze().str.replace(",", "")
        # now get text until first quote and concatenate with clean text
        f = c.str.extract('(.*),"').fillna("").squeeze() + "," + cc
        # now get back together with instances that didn't have quotes
        c.update(f)
        try:
            mapping = pd.DataFrame([x.strip().split(",") for x in c], columns=attributes)
        except AssertionError:
            _LOGGER.error("Attempt to fix failed.")
            raise e
        _LOGGER.info("Attempt to fix successful.")

    return mapping.replace("", np.nan)


def subtract_principal_component(
    x,
    pc=1,
    standardize=False,
    plot=True,
    plot_name="PCA_based_batch_correction.svg",
    max_pcs_to_plot=6,
):
    """
    Given a matrix (n_samples, n_variables), remove `pc` (1-based) from matrix.
    """
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from ngs_toolkit.graphics import savefig

    pc -= 1

    # All regions
    if standardize:
        x = StandardScaler().fit_transform(x)

    # PCA
    pca = PCA()
    x_hat = pca.fit_transform(x)

    # Remove PC
    x2 = x - np.outer(x_hat[:, pc], pca.components_[pc, :])

    max_pcs_to_plot = min(max_pcs_to_plot, x2.shape[0]) - 1

    # plot
    if plot:
        x2_hat = pca.fit_transform(x2)
        fig, axis = plt.subplots(max_pcs_to_plot, 2, figsize=(4 * 2, 4 * max_pcs_to_plot))
        axis[0, 0].set_title("Original")
        axis[0, 1].set_title("PC {} removed".format(pc + 1))
        for pc in range(max_pcs_to_plot):
            # before
            for j, _ in enumerate(x.index):
                axis[pc, 0].scatter(x_hat[j, pc], x_hat[j, pc + 1], s=50, rasterized=True)
            axis[pc, 0].set_xlabel("PC{}".format(pc + 1))
            axis[pc, 0].set_ylabel("PC{}".format(pc + 2))
            # after
            for j, _ in enumerate(x2.index):
                axis[pc, 1].scatter(
                    x2_hat[j, pc], x2_hat[j, pc + 1], s=35, alpha=0.8, rasterized=True
                )
            axis[pc, 1].set_xlabel("PC{}".format(pc + 1))
            axis[pc, 1].set_ylabel("PC{}".format(pc + 2))
        savefig(fig, plot_name)

    return x2


def fix_batch_effect_limma(matrix, batch_variable="batch", covariates=None):
    """
    Fix batch effect in matrix using limma.

    Requires the R package "limma" to be installed:

    .. highlight:: R
    .. code-block:: R

        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("limma")

    Parameters
    ----------
    matrix : :obj:`pandas.DataFrame`
        DataFrame with MultiIndex for potential covariate annotations

    formula : :obj:`str`, optional
        Model formula to regress out
        Defaults to "~batch"

    Returns
    -------
    :obj:`pandas.DataFrame`
        Regressed out matrix
    """
    import patsy
    from rpy2.robjects import numpy2ri, pandas2ri, r
    from rpy2.robjects.packages import importr

    numpy2ri.activate()
    pandas2ri.activate()

    importr("limma")

    if covariates is None:
        covariates = []

    if len(covariates) > 0:
        cov = patsy.dmatrix("~{} - 1".format(" + ".join(covariates)), matrix.columns.to_frame())
        fixed = r("removeBatchEffect")(
            x=matrix.values, batch=matrix.columns.get_level_values(batch_variable), design=cov,
        )
    else:
        fixed = r("removeBatchEffect")(
            x=matrix.values, batch=matrix.columns.get_level_values(batch_variable)
        )
    fixed = pd.DataFrame(np.asarray(fixed), index=matrix.index, columns=matrix.columns)
    return fixed
    # fixed.to_csv(os.path.join(analysis.results_dir, analysis.name + "_peaks.limma_fixed.csv"))
