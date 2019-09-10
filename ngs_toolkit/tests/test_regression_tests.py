#!/usr/bin/env python

import os
from ngs_toolkit import Analysis
import pybedtools
import pandas as pd


def test_pybedtools_to_from_dataframe():

    class T(Analysis):

        def __init__(self, *args, **kwargs):
            super(T, self).__init__(
                *args,
                **kwargs
            )
            self.samples = list()

        def trigger(self):
            d = pd.DataFrame(
                [['chr1', 1, 10], ['chr2', 1, 10]],
                columns=['chrom', 'start', 'end'])
            b = pybedtools.BedTool.from_dataframe(d)
            return b.to_dataframe()

    t = T()
    assert isinstance(t.trigger(), pd.DataFrame)


def test_get_right_timestamped_file(tmpdir):
    from ngs_toolkit.utils import get_this_file_or_timestamped

    target = os.path.join(tmpdir, "human.grch38.genomic_context.bed")
    assert get_this_file_or_timestamped(target) == target

    outs = [
        "human.grch38.genomic_context.2019-09-03-11:46:42.bed",
        "human.grch38.genomic_context.exon.2019-09-03-11:46:36.bed",
        "human.grch38.genomic_context.genebody.2019-09-03-11:46:36.bed",
        "human.grch38.genomic_context.intergenic.2019-09-03-11:46:41.bed",
        "human.grch38.genomic_context.intron.2019-09-03-11:46:38.bed",
        "human.grch38.genomic_context.promoter.2019-09-03-11:46:36.bed",
        "human.grch38.genomic_context.utr3.2019-09-03-11:46:40.bed",
        "human.grch38.genomic_context.utr5.2019-09-03-11:46:39.bed"]
    outs = [os.path.join(tmpdir, f) for f in outs]

    # Now with several existing files that also match the regex
    for f in outs:
        with open(f, "w") as handle:
            handle.write(f)

    assert get_this_file_or_timestamped(target) == outs[0]
