#!/usr/bin/env python


def test_config_has_all_required_fields():
    import logging
    from ngs_toolkit import _LOGGER
    import os

    assert isinstance(_LOGGER, logging.Logger)
    previous_size = os.stat(os.path.join(
        os.path.expanduser("~"), ".ngs_toolkit.log.txt")).st_size
    _LOGGER.info("Testing logger")
    new_size = os.stat(os.path.join(
        os.path.expanduser("~"), ".ngs_toolkit.log.txt")).st_size
    assert new_size > previous_size
