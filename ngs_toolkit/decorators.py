#!/usr/bin/env python

from functools import wraps


def check_organism_genome(f):
    from ngs_toolkit.general import Analysis

    @wraps(f)
    def wrapper(*args, **kwds):
        Analysis._check_organism_genome(args[0])
        return f(*args, **kwds)
    return wrapper


def check_has_sites(f):
    from ngs_toolkit.general import Analysis

    @wraps(f)
    def wrapper(*args, **kwds):
        Analysis._check_has_sites(args[0])
        return f(*args, **kwds)
    return wrapper
