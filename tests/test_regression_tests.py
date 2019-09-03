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
