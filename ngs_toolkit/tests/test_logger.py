#!/usr/bin/env python

import os
import pytest


@pytest.fixture
def log():
    return os.path.join(os.path.expanduser("~"), ".ngs_toolkit.log.txt")


def test_config_has_all_required_fields(log):
    import logging
    from ngs_toolkit import _LOGGER

    assert isinstance(_LOGGER, logging.Logger)
    previous_size = os.stat(log).st_size
    _LOGGER.info("Testing logger")
    new_size = os.stat(log).st_size
    assert new_size > previous_size


def test_clear_log(log):
    from ngs_toolkit import clear_log

    previous_size = os.stat(log).st_size
    clear_log()
    new_size = os.stat(log).st_size
    assert new_size < previous_size
    assert new_size == 0
