#!/usr/bin/env python

from functools import wraps
from ngs_toolkit import _LOGGER
from ngs_toolkit.utils import warn_or_raise


def check_has_attributes(attributes=None, object_types=None, permissive=False):
    attributes = [] or attributes
    object_types = [None] * len(attributes) or object_types
    if len(attributes) != len(object_types):
        msg = "`attributes` and `object_types` arguments must be the same length."
        _LOGGER.error(msg)
        raise ValueError(msg)

    def decorator(f):
        @wraps(f)
        def wrapper(*args, **kwargs):
            import pandas as pd

            # check attributes are set
            msg = "Analysis '{}' attribute(s) are not set."
            has = pd.Series(
                [hasattr(args[0], attr) for attr in attributes],
                index=attributes)
            if not has.all():
                warn_or_raise(AttributeError(msg.format(",".join(has[~has].index))), permissive)

            # check attributes are not None
            msg = "Analysis '{}' attribute(s) are None."
            not_none = pd.Series(
                [getattr(args[0], attr) is not None for attr in attributes],
                index=attributes)
            if not not_none.all():
                warn_or_raise(AttributeError(msg.format(",".join(not_none[~not_none].index))), permissive)

            # check the type of attribute values matches requested
            msg = "Analysis '{}' attribute(s) are not of requested types '{}'."
            t_attributes = [a for a, t in zip(attributes, object_types) if t is not None]
            t_object_types = [t for a, t in zip(attributes, object_types) if t is not None]
            not_type = pd.Series(
                [isinstance(getattr(args[0], attr), t) is not None
                 for attr, t in zip(t_attributes, t_object_types)],
                index=t_attributes, dtype=object)
            if not not_type.all():
                warn_or_raise(
                    AttributeError(msg.format(
                        ",".join(not_type[~not_type].index),
                        ",".join([str(t) for t in t_object_types]))),
                    permissive)

            # for iterable types, check length > 0
            msg = "Analysis '{}' attribute(s) have 0 elements."
            i_attributes = [a for a, t in zip(attributes, object_types) if hasattr(a, "__iter__")]
            i_object_types = [t for a, t in zip(attributes, object_types) if hasattr(a, "__iter__")]
            not_empty = pd.Series(
                [len(getattr(args[0], attr)) > 0 for attr in i_attributes],
                index=i_attributes)
            if not not_empty.all():
                warn_or_raise(
                    AttributeError(msg.format(
                        ",".join(not_empty[~not_empty].index),
                        ",".join([str(t) for t in i_object_types]))),
                    permissive)
            return f(*args, **kwargs)
        return wrapper
    return decorator


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


def to_csv_timestamped(f, exclude_functions=None):

    # TODO: fix to files without "." (dot)
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
                        record_analysis_output(args[1])
                else:
                    if isinstance(args[0], str):
                        s = args[0].split(".")
                        end = s[-1]
                        body = ".".join(s[:-1])
                        args = (".".join([body, get_timestamp(), end])) + args[1:]
                        record_analysis_output(args[0])
        return f(*args, **kwargs)
    return wrapper


def timestamped_input(f):
    from ngs_toolkit.utils import get_this_file_or_timestamped

    @wraps(f)
    def wrapper(file):
        return f(get_this_file_or_timestamped(file))
    return wrapper
