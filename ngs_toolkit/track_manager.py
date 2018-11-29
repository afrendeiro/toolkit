import argparse
import os
import re
import string
import sys
import textwrap

import pandas as pd
from peppy import Project
import ngs_toolkit


_LOGGER = ngs_toolkit._LOGGER


def parse_arguments():
    """
    Argument Parsing.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="project_config_file")
    parser.add_argument(
        '-a', '--attrs', dest="attributes", nargs="?",
        help="Sample attributes (annotation sheet columns) to use to order tracks. " +
             "Add attributes comma-separated with no whitespace. " +
             "Default is list of 'group_attributes' specified in `project_config_file`, " +
             "if none are found, will default to 'sample_name'")
    parser.add_argument(
        '-c', '--color-attr', dest="color_attribute", default=None,
        help="Sample attribute to use to color tracks with. Default is first attribute passed.")
    parser.add_argument(
        '-r', '--overlay-replicates', dest="overlay_replicates", action="store_true",
        help="Whether replicate samples should be overlaied in same track. Default=False.")
    parser.add_argument(
        '-l', '--link', dest="link", action="store_true",
        help="Whether bigWig files should be soft-linked to the track database directory. Default=False.")
    parser.add_argument("-V", "--version", action="version",
                        version=ngs_toolkit.__version__)

    args = parser.parse_args()

    return args


def get_colors(level, pallete="gist_rainbow", nan_color=[0.5, 0.5, 0.5]):
    """
    Given a level (list or iterable) with length n,
    return a same-sized list of RGB colors  for each unique value in the level.
    """
    import matplotlib.pyplot as plt

    pallete = plt.get_cmap(pallete)
    level_unique = list(set(level))
    n = len(level_unique)
    # get n equidistant colors
    p = [pallete(1. * i / n) for i in range(n)]
    # convert to integer RGB strings with no alpha channel
    p = [",".join([str(int(y * 255)) for y in x[:-1]]) for x in p]

    color_dict = dict(zip(list(set(level)), p))
    # color for nan cases
    color_dict[pd.np.nan] = nan_color
    return [color_dict[x] for x in level]


def make_ucsc_trackhub(args, prj, track_hub, bigwig_dir, genomes_hub, proj_name, proj_desc, user_email):
    """
    Make UCSC trackHub for project
    """
    # Start building hub
    text = """hub {proj}
shortLabel {description}
longLabel {description}
genomesFile genomes.txt
email {email}
""".format(proj=proj_name, description=proj_desc, email=user_email)
    with open(track_hub, 'w') as handle:
        handle.write(text)
    os.chmod(track_hub, 0755)

    # Templates for various levels
    track_parent = """track {proj}
container multiWig
shortLabel {description}
longLabel {description}
container multiWig
aggregate none
showSubtrackColorOnUi on
type bigWig
autoScale on
visibility full
maxHeightPixels 32:32:8
{0}{1}{2}
"""

    track_middle = """
    track {track}
    shortLabel {desc}
    longLabel {desc}
    parent {parent}
    container multiWig
    aggregate transparentOverlay
    showSubtrackColorOnUi on
    type bigWig
    visibility full
    maxHeightPixels 32:32:8
    {subgroups}
"""

    track_final = """
        track {name}
        shortLabel {name}
        longLabel {name}
        parent {parent}
        type bigWig
        graphTypeDefault bar
        visibility full
        height 32
        maxHeightPixels 32:32:8
        windowingFunction mean
        autoScale on
        bigDataUrl {bigwig}
        color {color}
        {subgroups}

    """

    # Make dataframe for groupby
    df = pd.DataFrame([s.as_series() for s in prj.samples]).fillna("none").drop_duplicates(subset="sample_name")

    # Keep only samples that have appropriate types to be displayed
    if 'library' in df.columns:
        var_ = 'library'
    elif 'protocol' in df.columns:
        var_ = 'protocol'
    else:
        raise ValueError("Samples must contain either a 'library' or 'protocol' attribute.")
    df = df[df[var_].isin(["ATAC-seq", "ChIP-seq", "ChIPmentation"])]

    # Keep only attributes that have at least one value
    d = df[args.attributes].loc[:, df[args.attributes].isnull().all()].columns.tolist()
    attributes = [a for a in args.attributes if a not in d]

    # Create a trackHub for each genome
    for genome in df['genome'].unique():
        if not os.path.exists(os.path.join(bigwig_dir, genome)):
            os.makedirs(os.path.join(bigwig_dir, genome))
        os.chmod(os.path.join(bigwig_dir, genome), 0755)

        df_g = df[(df['genome'] == genome)]

        # Genomes
        text = """genome {g}
trackDb {g}/trackDb.txt
    """.format(g=genome)
        with open(genomes_hub, 'a') as handle:
            handle.write(text)
        os.chmod(genomes_hub, 0755)

        # TrackDB
        track_db = os.path.join(os.path.join(bigwig_dir, genome, "trackDb.txt"))
        open(track_db, 'w').write("")

        # Create subgroups, an attribute sort order and an experiment matrix
        subgroups = "\n".join([
            """subGroup{0} {1} {1} \\\n    {2}""".format(i + 1, attr, " \\\n    ".join(
                [x + "=" + x for x in sorted(df_g[attr].unique())])) for i, attr in enumerate(attributes)])
        subgroups = subgroups.replace("    = \\\n", "")  # remove values of subgroups with no value
        sort_order = "\n" + "sortOrder " + " ".join([x + "=+" for x in attributes])
        dimensions = "\n" + "dimensions " + " ".join(
            ["".join(x) for x in zip(["dim{}=".format(x) for x in ["X", "Y"] + [string.ascii_uppercase[:-3]]], attributes)])

        track = track_parent.format(subgroups,
                                    sort_order,
                                    dimensions,
                                    proj=proj_name, description=proj_desc)

        # Get unique colors for the given attribute
        df_g['track_color'] = get_colors(df[args.color_attribute])

        # Group by requested attributes, add tracks
        for labels, indices in df_g.groupby(attributes).groups.items():
            subgroups = "subGroups " + " ".join(["{}={}".format(k, v) for k, v in zip(attributes, labels)])

            if len(indices) == 1:
                sample_attrs = df_g.ix[indices].squeeze()

                track += textwrap.dedent(track_final.format(
                    name=sample_attrs["sample_name"],
                    color=sample_attrs['track_color'],
                    parent=proj_name,
                    bigwig=os.path.basename(sample_attrs['bigwig']),
                    subgroups=subgroups))

                # Make symbolic link to bigWig
                dest = os.path.join(os.path.join(bigwig_dir, genome, os.path.basename(sample_attrs['bigwig'])))
                if not os.path.exists(dest) and args.link:
                    try:
                        os.symlink(sample_attrs['bigwig'], dest)
                        os.chmod(dest, 0755)
                    except OSError:
                        print("Sample {} track file does not exist!".format(sample_attrs["sample_name"]))
                        continue

            else:
                name = "_".join([x for x in labels if x != ""])
                desc = " ".join([x for x in labels if x != ""])

                track += track_middle.format(
                    track=name,
                    desc=desc,
                    parent=proj_name,
                    subgroups=subgroups)

                for index in indices:
                    sample_attrs = df_g.ix[index].squeeze()
                    track += track_final.format(
                        name=sample_attrs["sample_name"],
                        color=sample_attrs['track_color'],
                        parent=name,
                        bigwig=os.path.basename(sample_attrs['bigwig']),
                        subgroups=subgroups)

                # Make symbolic link to bigWig
                dest = os.path.join(os.path.join(bigwig_dir, genome, os.path.basename(sample_attrs['bigwig'])))
                if not os.path.exists(dest) and args.link:
                    try:
                        os.symlink(sample_attrs['bigwig'], dest)
                        os.chmod(dest, 0755)
                    except OSError:
                        print("Sample {} track file does not exist!".format(sample_attrs["sample_name"]))
                        continue
                # Make directories readable and executable
                os.chmod(os.path.join(bigwig_dir, genome), 0755)
                os.chmod(bigwig_dir, 0755)

        track = re.sub("_none", "", track)

        # write trackDb to file
        with open(track_db, 'w') as handle:
            handle.write(textwrap.dedent(track))

    msg = "\n".join([
        "Finished producing trackHub!",
        "----------------------------",
        "Add the following URL to your UCSC trackHubs:",
        "{url}/hub.txt".format(url=prj['trackhubs']['url']),
        "or alternatively follow this URL: " +
        "http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db={genome}&hubUrl={url}/hub.txt".format(
            # the link can only point to one genome, so by default it will be the last one used.
            genome=genome, url=prj['trackhubs']['url'])
    ]) + "\n"

    if 'trackhubs' in prj:
        if "url" in prj['trackhubs']:
            print(msg)


def make_igv_tracklink(prj, track_file, track_url):
    """
    Make IGV track link for project
    """
    if not os.path.exists(os.path.dirname(track_file)):
        os.makedirs(os.path.dirname(track_file))

    # Start building hub
    link_header = "http://localhost:60151/load"

    # Make dataframe
    df = pd.DataFrame([s.as_series() for s in prj.samples]).fillna("none")

    # Keep only samples that have appropriate types to be displayed
    if 'library' in df.columns:
        var_ = 'library'
    elif 'protocol' in df.columns:
        var_ = 'protocol'
    else:
        raise ValueError("Samples must contain either a 'library' or 'protocol' attribute.")
    df = df[df[var_].isin(["ATAC-seq", "ChIP-seq", "ChIPmentation"])]

    text = "<html><head></head><body>"

    # Create a trackHub for each genome
    for genome in df['genome'].unique():
        df_g = df[(df['genome'] == genome)]

        text += "<a href="
        text += link_header
        text += "?file=" + ",".join(prj['trackhubs']['url'] + "/{}/".format(genome) + df_g['sample_name'] + ".bigWig")
        text += "?names=" + ",".join(df_g['sample_name'])
        text += "?genome={}".format(genome)
        text += ">Open {} tracks in IGV browser.</a>\n".format(genome)
    text += "</body></html>"

    # write to file
    with open(track_file, 'w') as handle:
        handle.write(text + "\n")
    os.chmod(track_file, 0655)
    os.chmod(os.path.dirname(track_file), 0755)

    msg = "\n".join([
        "Finished producing IGV track file!", "'{}'".format(track_file),
        "You can follow this URL to open tracks in a local IGV session: " +
        "{url}\n".format(url=track_url)
    ]) + "\n"
    print(msg)


def main():
    args = parse_arguments()

    # Start project object
    prj = Project(args.project_config_file)

    if "trackhubs" not in prj:
        raise ValueError("Project configuration does not have a trackhub section.")
    if "trackhub_dir" not in prj.trackhubs:
        raise ValueError("Project trackhub configuration does not have a trackhub_dir attribute.")

    if args.attributes is None:
        # If no attributes are set try to use group attributes specified in the project config
        if hasattr(prj, "group_attributes"):
            print("Using group attributes from project configuration file: '{}'.".format(",".join(prj.group_attributes)))
            args.attributes = prj.group_attributes
        # If none are set in project config, use just "sample_name"
        else:
            args.attributes = ['sample_name']
    else:
        args.attributes = args.attributes.split(',')

    if args.color_attribute is None:
        args.color_attribute = args.attributes[0]

    # Setup paths and hub files
    bigwig_dir = os.path.join(prj.trackhubs.trackhub_dir)
    track_hub = os.path.join(bigwig_dir, "hub.txt")
    genomes_hub = os.path.join(bigwig_dir, "genomes.txt")
    open(genomes_hub, 'w').write("")

    # Setup attributes
    proj_name = prj['project_name'] if "project_name" in prj else os.path.basename(prj['paths']['output_dir'])
    proj_desc = prj['project_description'] if "project_description" in prj else proj_name
    user_email = prj['email'] if "email" in prj else ""

    # In the future there will be more actions than this
    make_ucsc_trackhub(args, prj, track_hub, bigwig_dir, genomes_hub, proj_name, proj_desc, user_email)

    track_file = os.path.join(bigwig_dir, "igv", "index.html")
    track_url = os.path.join(prj['trackhubs']['url'], "igv")
    make_igv_tracklink(prj, track_file, track_url)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
