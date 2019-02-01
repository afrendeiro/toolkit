#!/usr/bin/env python

import glob
import os
import re

import pandas as pd


def parse_ame(ame_dir):
    """
    Parse results of MEME-AME motif enrichment.

    :param ame_dir: Directory with MEME-AME results.
    :type ame_dir: str
    :returns: Data frame with enrichment statistics for each found TF motif.
    :rtype: pandas.DataFrame
    :raises: IOError
    """

    with open(os.path.join(ame_dir, "ame.txt"), 'r') as handle:
        lines = handle.readlines()

    output = list()
    for line in lines:
        # skip header lines
        if line[0] not in [str(i) for i in range(10)]:
            continue

        # get motif string and the first half of it (simple name)
        motif = line.strip().split(" ")[5].split("_")[0]
        # get corrected p-value
        q_value = float(line.strip().split(" ")[-2])
        # append
        output.append((motif, q_value))

    r = pd.Series(dict(output)).reset_index()
    r.columns = ["TF", "p_value"]
    return r


def parse_homer(homer_dir):
    """
    Parse results of HOMER findMotifs.pl de novo motif enrichment.

    :param homer_dir: Directory with HOMER results.
    :type homer_dir: str
    :returns: Data frame with enrichment statistics for each found TF motif.
    :rtype: pandas.DataFrame
    :raises: IOError
    """
    motif_htmls = sorted(glob.glob(os.path.join(homer_dir, "motif*.info.html")))

    if len(motif_htmls) < 1:
        raise IOError("Homer directory does not contain any discovered motifs.")

    output = pd.DataFrame()
    for motif_html in motif_htmls:

        motif = int(re.sub(".info.html", "", re.sub(os.path.join(homer_dir, "motif"), "", motif_html)))

        with open(motif_html, 'r') as handle:
            content = handle.read()

        # Parse table with motif info
        info_table = content[
            re.search("""<TABLE border="1" cellpading="0" cellspacing="0">""", content).end():
            re.search("</TABLE>", content).start()].strip()

        info_table = pd.DataFrame([x.split("</TD><TD>") for x in info_table.replace("<TR><TD>", "").split("</TD></TR>")])
        info_table.columns = ["description", "value"]
        info_table["description"] = info_table["description"].str.strip()
        info_table["motif"] = motif

        # Add most probable known motif name
        info_table["known_motif"] = content[
            re.search("<H4>", content).end():
            re.search("</H4>", content).start()]

        # append
        output = output.append(info_table, ignore_index=True)

    return output.sort_values("motif")


def parse_great_enrichment(input_tsv):
    """
    Parse output from GREAT enrichment (http://great.stanford.edu).

    :param input_tsv: TSV file exported from GREAT through the option "All data as .tsv" in "Global Controls".
    :type input_tsv: str
    :returns: Pandas dataframe with enrichment results.
    :rtype: pandas.DataFrame
    """
    df = pd.read_csv(input_tsv, sep="\t", skiprows=3)
    df.columns = df.columns.str.replace("# ", "")
    return df.loc[~df.iloc[:, 0].str.startswith("#")]
