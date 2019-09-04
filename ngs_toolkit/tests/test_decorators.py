#!/usr/bin/env python


import pytest
import os


travis = "TRAVIS" in os.environ


class Test_check_organism_genome:
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


class Test_check_has_sites:
    # here we use 'calculate_peak_support' as en example
    # decorated function. It will however fail for another
    # reason due to the fairly empty analysis object (last test)
    def test_empty_analysis(self, empty_analysis):
        # Make sure it raises AttributeError
        with pytest.raises(AttributeError):
            empty_analysis.calculate_peak_support()

    def test_null_analysis(self, null_analysis):
        # Make sure it raises AttributeError
        with pytest.raises(AttributeError):
            null_analysis.calculate_peak_support()

    def test_full_analysis(self, atac_analysis):
        # This passes on the decorator
        # but raises IOError specific to the function
        with pytest.raises(IOError):
            atac_analysis.calculate_peak_support()
