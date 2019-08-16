#!/usr/bin/env python

from functools import wraps
from ngs_toolkit import _LOGGER


def check_has_samples(f):
    @wraps(f)
    def wrapper(*args, **kwds):
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
        return f(*args, **kwds)

    return wrapper


def check_organism_genome(f):
    @wraps(f)
    def wrapper(*args, **kwds):
        attrs = ["organism", "genome"]
        msg = "Analysis does not have 'organism' and 'genome' attributes set."
        hint = " You can set them with analysis.set_organism_genome, for example."
        r1 = all([hasattr(args[0], attr) for attr in attrs])
        r2 = all([getattr(args[0], attr) is not None for attr in attrs])
        if not all([r1, r2]):
            _LOGGER.error(msg + hint)
            raise AttributeError(msg)
        return f(*args, **kwds)

    return wrapper


def check_has_sites(f):
    @wraps(f)
    def wrapper(*args, **kwds):
        attrs = ["sites"]
        msg = "Analysis object does not have a `sites` attribute."
        hint = " Produce one with analysis.get_consensus_sites for example."
        r1 = all([hasattr(args[0], attr) for attr in attrs])
        r2 = all([getattr(args[0], attr) is not None for attr in attrs])
        if not all([r1, r2]):
            _LOGGER.error(msg + hint)
            raise AttributeError(msg)
        return f(*args, **kwds)

    return wrapper
