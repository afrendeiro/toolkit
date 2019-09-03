#!/usr/bin/env python

from functools import wraps
from ngs_toolkit import _LOGGER


def check_has_samples(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        msg = "Analysis does not have a 'samples' attributes."
        if not hasattr(args[0], "samples"):
            _LOGGER.error(msg)
            raise AttributeError(msg)

        msg = "Analysis 'samples' attribute is not a list."
        if not isinstance(args[0].samples, list):
            _LOGGER.error(msg)
            raise AttributeError(msg)

        msg = "Analysis 'samples' attribute empty."
        if len(args[0].samples) == 0:
            _LOGGER.error(msg)
            raise AttributeError(msg)
        return f(*args, **kwargs)

    return wrapper


def check_organism_genome(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        attrs = ["organism", "genome"]
        msg = "Analysis does not have 'organism' and 'genome' attributes set."
        hint = " You can set them with analysis.set_organism_genome, for example."
        r1 = all([hasattr(args[0], attr) for attr in attrs])
        r2 = all([getattr(args[0], attr) is not None for attr in attrs])
        if not all([r1, r2]):
            _LOGGER.error(msg + hint)
            raise AttributeError(msg)
        return f(*args, **kwargs)

    return wrapper


def check_has_sites(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        attrs = ["sites"]
        msg = "Analysis object does not have a `sites` attribute."
        hint = " Produce one with analysis.get_consensus_sites for example."
        r1 = all([hasattr(args[0], attr) for attr in attrs])
        r2 = all([getattr(args[0], attr) is not None for attr in attrs])
        if not all([r1, r2]):
            _LOGGER.error(msg + hint)
            raise AttributeError(msg)
        return f(*args, **kwargs)

    return wrapper


def read_csv_timestamped(f):
    from ngs_toolkit.utils import get_this_file_or_timestamped
    @wraps(f)
    def wrapper(*args, **kwargs):
        for i, _ in enumerate(args):
            if isinstance(args[i], str):
                args = args[:i] + (
                    get_this_file_or_timestamped(args[i]),) + args[i + 1:]
        return f(*args, **kwargs)
    return wrapper


def to_csv_timestamped(f, exclude_functions=[]):
    from ngs_toolkit.utils import (
        record_analysis_output, get_timestamp,
        is_analysis_descendent)
    from ngs_toolkit import _CONFIG

    @wraps(f)
    def wrapper(*args, **kwargs):
        if is_analysis_descendent(exclude_functions=exclude_functions):
            # Add timestamp
            if _CONFIG["preferences"]["report"]["timestamp_tables"]:
                if len(args) > 1:
                    if isinstance(args[1], str):
                        s = args[1].split(".")
                        end = s[-1]
                        body = ".".join(s[:-1])
                        args = (args[0], ".".join([body, get_timestamp(), end])) + args[2:]
                        record_analysis_output(args[1], permissive=True)
                else:
                    if isinstance(args[0], str):
                        s = args[0].split(".")
                        end = s[-1]
                        body = ".".join(s[:-1])
                        args = (".".join([body, get_timestamp(), end])) + args[1:]
                        record_analysis_output(args[0], permissive=True)
        return f(*args, **kwargs)
    return wrapper


def timestamped_input(f):
    from ngs_toolkit.utils import get_this_file_or_timestamped
    @wraps(f)
    def wrapper(file):
        return f(get_this_file_or_timestamped(file))
    return wrapper
