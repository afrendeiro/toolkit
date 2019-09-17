#!/usr/bin/env python


import pytest


class Test_check_has_attributes:
    # here we use 'get_resources' as en example
    # decorated function that won't fail for some other
    # reason on a fairly empty analysis object
    def test_empty_analysis(self, empty_analysis):
        # Make sure it raises AttributeError
        with pytest.raises(AttributeError):
            empty_analysis.get_resources(steps=[])

    def test_null_analysis(self, null_analysis):
        # Make sure it raises AttributeError
        with pytest.raises(AttributeError):
            null_analysis.get_resources(steps=[])

    def test_full_analysis(self, full_analysis):
        full_analysis.get_resources(steps=[])

    # here we use 'calculate_peak_support' as en example
    # decorated function. It will however fail for another
    # reason due to the fairly empty analysis object (last test)
    def test_empty_analysis_2(self, empty_analysis):
        # Make sure it raises AttributeError
        with pytest.raises(AttributeError):
            empty_analysis.calculate_peak_support()

    def test_null_analysis_2(self, null_analysis):
        # Make sure it raises AttributeError
        with pytest.raises(AttributeError):
            null_analysis.calculate_peak_support()

    def test_full_analysis_2(self, atac_analysis):
        # This passes on the decorator
        # but raises IOError specific to the function
        with pytest.raises(IOError):
            atac_analysis.calculate_peak_support()

    def test_iterable_attributes(self, atac_analysis):
        from ngs_toolkit import Analysis
        from ngs_toolkit.decorators import check_has_attributes

        class TestAnalysis(Analysis):
            @check_has_attributes(['samples'], [list])
            def test_function(self):
                print(self.samples)
                return True

        a = TestAnalysis()

        # has not samples set
        del a.samples
        with pytest.raises(AttributeError):
            a.test_function()

        # samples is None
        a.samples = None
        with pytest.raises(AttributeError):
            a.test_function()

        # samples is empty list
        a.samples = list()
        with pytest.raises(AttributeError):
            a.test_function()

        # has samples
        a.samples = [1, 2, 3]
        assert a.test_function()
