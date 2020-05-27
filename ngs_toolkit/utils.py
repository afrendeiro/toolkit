#!/usr/bin/env python


import os

import numpy as np
import pandas as pd


def _format_string_with_environment_variables(string):
    """
    Given a string, containing curly braces with dollar sign,
    format it with the environment variables.

    Parameters
    ----------
    string : :obj:`str`
        String to format.

    Returns
    ----------
    :obj:`str`
        Formated string.

    Raises
    -------
    :obj:`ValueError`
        If not all patterns are set environment variables.
    """
    if string is None:
        return string

    not_to_format = pd.Series(string).str.extractall(r"[^$]{(.*?)}")[0].values
    to_format = pd.Series(string).str.extractall(r"\${(.*?)}")[0].values
    attrs = os.environ
    if not all([x in attrs for x in to_format]):
        msg = "Not all required patterns were found in the environment variables."
        _LOGGER.error(msg)
        raise ValueError(msg)
    # first, make sure there are no unintended strings being formatted
    for attr in set(not_to_format):
        string = string.replace("{" + attr + "}", "<>" + attr + "<>")

    # now format the intended ones
    for attr in set(to_format):
        string = string.replace("${" + attr + "}", r"{" + attr + r"}")

    # now format
    string = string.format(**attrs)

    # now return the not ones to original
    for attr in set(not_to_format):
        string = string.replace("<>" + attr + "<>", "{" + attr + "}")
    return string


def warn_or_raise(exception, permissive=False):
    from ngs_toolkit import _LOGGER

    msg = exception.args[0]
    if permissive:
        _LOGGER.warning(msg)
    else:
        _LOGGER.error(msg)
        raise exception


def have_unbuffered_output():
    """Set unbuffered output for current session."""
    import sys

    class Unbuffered:
        def __init__(self, stream):
            self.stream = stream

        def write(self, data):
            self.stream.write(data)
            self.stream.flush()

        def writelines(self, datas):
            self.stream.writelines(datas)
            self.stream.flush()

        def __getattr__(self, attr):
            return getattr(self.stream, attr)

    sys.stdout = Unbuffered(sys.stdout)


def is_running_inside_ipython():
    """Check whether code is running inside an IPython session."""
    try:
        cfg = get_ipython().config
        if cfg["IPKernelApp"]["parent_appname"] == "ipython-notebook":
            return True
        else:
            return False
    except NameError:
        return False


def filter_kwargs_by_callable(kwargs, callabl, exclude=None):
    from inspect import signature

    args = signature(callabl).parameters.keys()
    return {k: v for k, v in kwargs.items() if (k in args) and k not in (exclude or [])}


def get_timestamp(fmt="%Y-%m-%d-%H:%M:%S"):
    """Get current timestamp in ``fmt`` format."""
    from datetime import datetime

    return datetime.today().strftime(fmt)


def remove_timestamp_if_existing(file):
    """Remove timestamp from path if matching timestamp pattern exists."""
    import re

    return re.sub(r"\d{4}-\d{2}-\d{2}-\d{2}:\d{2}:\d{2}\.", "", file)


def get_this_file_or_timestamped(file, permissive=True):
    """
    Get a path to an existing timestamped file based on an non-timestamped path.

    Parameters
    ----------
    file_name : :obj:`str`
        File name of analysis output to record.
    permissive : :obj:`bool`
        Whether failure to find timestamped file should return the original file
        or raise a IndexError.

    Raises
    ----------
    IndexError
        If not `permissive` and can't find timestamped file.
    """
    from glob import glob
    import re

    from ngs_toolkit.utils import sorted_nicely
    from ngs_toolkit import _LOGGER

    split = file.split(".")
    body = ".".join(split[:-1])
    end = split[-1]

    res = sorted_nicely(glob(body + "*" + end))
    res = [x for x in res if re.search(body + r"\.\d{4}-\d{2}-\d{2}-\d{2}:\d{2}:\d{2}\.", x)]
    if len(res) > 1:
        _LOGGER.warning(
            "Could not get unequivocal timestamped file for '{}'.".format(file)
            + " Returning latest: '{}'.".format(res[-1])
        )
    try:
        # get newest file
        return res[-1]
    except IndexError:
        if permissive:
            return file
        else:
            msg = "Could not remove timestamp from file path."
            msg += " Probabably it does not exist."
            _LOGGER.error(msg)
            raise IndexError(msg)


def is_analysis_descendent(exclude_functions=None):
    """
    Check whether any call in the traceback comes from a function part of a
    :meth:`ngs_toolkit.Analysis` object.

    Parameters
    ----------
    exclude_functions : :obj:`list`
        List of function names to exclude from.

    Returns
    ----------
    tuple
        If is descentent, returns tuple of (Analysis instance, function name),
        othewise returns :obj:`False`.
    """
    import inspect
    from ngs_toolkit import Analysis

    if exclude_functions is None:
        exclude_functions = list()

    # Let's get the object that called the function previous to this one
    # # the use case is often:
    # # Analysis().do_work() <- do work will produce a plot and save it using savefig.
    # # If savefig(track=True), this function will be called and we can trace which Analysis object did so

    # Go up the stack until an Analysis object is found:
    for s in inspect.stack():
        if s.function in exclude_functions:
            return False
        # # Function does not have "self" in locals and therefore
        # # unlikely part of Analysis
        if "self" not in s.frame.f_locals:
            continue
        # # Get Analysis object
        if isinstance(s.frame.f_locals["self"], Analysis):
            return (s.frame.f_locals["self"], s.function)
    return False


def record_analysis_output(file_name, **kwargs):
    """
    Register a file that is an output of an Analysis.
    The file will be associated with the function that produced it and
    saved in the attribute ``output_files``.

    Parameters
    ----------
    file_name : :obj:`str`
        File name of analysis output to record.
    **kwargs : :obj:`bool`
        Keyword arguments passed to :meth:`ngs_toolkit.utils.is_analysis_descendent`.
    """
    out = is_analysis_descendent(**kwargs)
    if out:
        out[0].record_output_file(file_name, out[1])


def submit_job(
    code,
    job_file,
    log_file=None,
    computing_configuration=None,
    dry_run=False,
    limited_number=False,
    total_job_lim=500,
    refresh_time=10,
    in_between_time=5,
    **kwargs
):
    """
    Submit a job to be run.
    Uses divvy to allow running on a local computer or distributed computing resources.

    Parameters
    ----------
    code : :obj:`str`
        String of command(s) to be run.
    job_file : :obj:`str`
        File to write job ``code`` to.
    log_file : :obj:`str`
        Log file to write job output to.

        Defaults to ``job_file`` with ".log" ending.
    computing_configuration : :obj:`str`
        Name of :class:`divvy` computing configuration to use.

        Defaults to 'default' which is to run job in localhost.
    dry_run: :obj:`bool`
        Whether not to actually run job.

        Defaults to :obj:`False`.
    limited_number: :obj:`bool`
        Whether to restrict jobs to a maximum number.
        Currently only possible if using "slurm".

        Defaults to :obj:`False`.
    total_job_lim : :obj:`int`
        Maximum number of jobs to restrict to.

        Defaults to 500.
    refresh_time : :obj:`int`
        Time in between checking number of jobs in seconds.

        Defaults to 10.
    in_between_time : :obj:`int`
        Time in between job submission in seconds.

        Defaults to 5.
    **kwargs : :obj:`dict`
        Additional keyword arguments will be passed to the chosen submission template according to `computing_configuration`.
        Pass for example: jobname="job", cores=2, mem=8000, partition="longq".
    """
    import time
    import subprocess

    import divvy
    from ngs_toolkit import _CONFIG, _LOGGER

    # reduce level of logging from divvy
    # only for divvy <=0.
    if "logging" in divvy.__dict__.keys():
        divvy.logging.getLogger("divvy").setLevel("ERROR")

    def count_jobs_running(check_cmd="squeue", sep="\n"):
        """
        Count running jobs on a cluster by invoquing a command that lists the jobs.
        """
        return subprocess.check_output(check_cmd).split(sep).__len__()

    def submit_job_if_possible(
        cmd, check_cmd="squeue", total_job_lim=800, refresh_time=10, in_between_time=5
    ):
        submit = count_jobs_running(check_cmd) < total_job_lim
        while not submit:
            time.sleep(refresh_time)
            submit = count_jobs_running(check_cmd) < total_job_lim
        subprocess.call(cmd)
        time.sleep(in_between_time)

    if log_file is None:
        log_file = ".".join(job_file.split(".")[:-1]) + ".log"

    # Get computing configuration from config
    if computing_configuration is None:
        try:
            computing_configuration = _CONFIG["preferences"]["computing_configuration"]
        except KeyError:
            msg = "'computing_configuration' was not given"
            msg += " and default could not be get from config."
            hint = " Pass a value or add one to the section"
            hint += " preferences:computing_configuration'"
            hint += " in the ngs_toolkit config file."
            _LOGGER.error(msg + hint)
            raise

    dcc = divvy.ComputingConfiguration()
    if computing_configuration is not None:
        dcc.activate_package(computing_configuration)

    # Generate job script
    d = {"code": code, "logfile": log_file}
    d.update(kwargs)
    dcc.write_script(job_file, d)

    # Submit job
    if not dry_run:
        scmd = dcc["compute"]["submission_command"]
        cmd = scmd.split(" ") + [job_file]

        # simply submit if not limiting submission to the number of already running jobs
        if not limited_number:
            subprocess.call(cmd)
        else:
            # otherwise, submit only after `total_job_lim` is less than number of runnning jobs
            # this is only possible for slurm now though
            if scmd != "slurm":
                subprocess.call(cmd)
            else:
                submit_job_if_possible(cmd, check_cmd="slurm")


def chunks(l, n):
    """
    Partition iterable `l` in chunks of size `n`.

    Parameters
    ----------
    l : iterable
        Iterable (e.g. list or numpy array).
    n : :obj:`int`
        Size of chunks to generate.
    """
    n = max(1, n)
    return list(l[i : i + n] for i in range(0, len(l), n))


def sorted_nicely(l):
    """
    Sort an iterable in the way that humans expect.

    Parameters
    ----------
    l : iterable
        Iterable to be sorted

    Returns
    -------
    iterable
        Sorted iterable
    """
    import re

    def convert(text):
        return int(text) if text.isdigit() else text

    def alphanum_key(key):
        return [convert(c) for c in re.split("([0-9]+)", key)]

    return sorted(l, key=alphanum_key)


def standard_score(x):
    """
    Compute a standard score, defined as (x - min(x)) / (max(x) - min(x)).

    Parameters
    ----------
    x : :class:`numpy.ndarray`
        Numeric array.

    Returns
    -------
    :class:`numpy.ndarray`
        Transformed array.
    """
    return (x - x.min()) / (x.max() - x.min())


def z_score(x):
    """
    Compute a Z-score, defined as (x - mean(x)) / std(x).

    Parameters
    ----------
    x : :class:`numpy.ndarray`
        Numeric array.

    Returns
    -------
    :class:`numpy.ndarray`
        Transformed array.
    """
    return (x - x.mean()) / x.std()


def logit(x):
    """
    Compute the logit of x, defined as log(x / (1 - x)).

    Parameters
    ----------
    x : :class:`numpy.ndarray`
        Numeric array.

    Returns
    -------
    :class:`numpy.ndarray`
        Transformed array.
    """
    return np.log(x / (1 - x))


def count_dataframe_values(x):
    """
    Count number of non-null values in a dataframe.

    Parameters
    ----------
    x : :class:`pandas.DataFrame`
        Pandas DataFrame

    Returns
    -------
    int
        Number of non-null values.
    """
    return np.multiply(*x.shape) - x.isnull().sum().sum()


def location_index_to_bed(index):
    """
    Get a pandas DataFrame with columns "chrom", "start", "end"
    from an pandas Index of strings in form "chrom:start-end".

    Parameters
    ----------
    index : {:obj:`list`, :class:`pandas.Index`, :class:`pandas.Series)`, :class:`pandas.DataFrame)`}
        Index strings of the form "chrom:start-end".

    Returns
    -------
    :class:`pandas.DataFrame`
        Pandas dataframe.
    """
    bed = pd.DataFrame(index=index)
    if isinstance(index, list):
        index = pd.Series(index=index, name="region").index
    elif isinstance(index, pd.DataFrame):
        index = index.to_series(name="region").index
    elif isinstance(index, pd.Index):
        pass
    else:
        msg = "index is not list, Series or Index"
        TypeError(msg)
    bed.loc[:, "chrom"] = index.str.split(":").str[0]
    index2 = index.str.split(":").str[1]
    bed.loc[:, "start"] = index2.str.split("-").str[0].astype(int)
    bed.loc[:, "end"] = index2.str.split("-").str[1].astype(int)
    return bed


def bed_to_index(df):
    """
    Get an index of the form chrom:start-end
    from a a dataframe with such columns.

    Parameters
    ----------
    df : {:class:`pandas.DataFrame`, :class:`pybedtools.bedtool.BedTool`, :obj:`str`}
        DataFrame with columns "chrom", "start" and "end".

    Returns
    -------
    :class:`pandas.Index`
        Pandas index.
    """
    import pybedtools

    if isinstance(df, pybedtools.BedTool):
        df = df.to_dataframe()
    elif isinstance(df, str):
        df = pybedtools.BedTool(df).to_dataframe()
    cols = ["chrom", "start", "end"]
    if not all([x in df.columns for x in cols]):
        raise AttributeError("DataFrame does not have '{}' columns.".format("', '".join(cols)))
    index = (
        df["chrom"].astype(str)
        + ":"
        + df["start"].astype(int).astype(str)
        + "-"
        + df["end"].astype(int).astype(str)
    )
    return pd.Index(index, name="region")


def bedtool_to_index(bedtool):
    """Convert bedtool or path to a bed file to list of region identifiers"""
    import pybedtools

    if isinstance(bedtool, str):
        bedtool = pybedtools.BedTool(bedtool)
    elif isinstance(bedtool, pybedtools.BedTool):
        pass
    else:
        msg = "Input not pybedtools.BedTool or string to BED file."
        raise ValueError(msg)
    return [str(i.chrom) + ":" + str(i.start) + "-" + str(i.stop) for i in bedtool]


def to_bed_index(sites):
    """Convert bedtool, BED file or dataframe to list of region identifiers"""
    import pybedtools

    msg = "Input not pybedtools.BedTool, pandas.DataFrame or path to BED file."
    if isinstance(sites, str):
        return bedtool_to_index(pybedtools.BedTool(sites))
    elif isinstance(sites, pybedtools.BedTool):
        return bedtool_to_index(sites)
    elif isinstance(sites, pd.core.frame.DataFrame):
        return bed_to_index(sites)
    elif isinstance(sites, str):
        return bedtool_to_index(sites)
    else:
        raise ValueError(msg)


def sort_bed_nicely(bed_file):
    """Sorts BED file but in sorted_nicely order"""
    import pybedtools
    import tempfile
    from ngs_toolkit.utils import sorted_nicely

    if isinstance(bed_file, str):
        sites = pybedtools.BedTool(bed_file)
    elif isinstance(bed_file, pybedtools.BedTool):
        pass
    chrom_names = tempfile.NamedTemporaryFile()
    with open(chrom_names.name, "w") as handle:
        for chrom in sorted_nicely(sites.to_dataframe()["chrom"].unique()):
            handle.write(chrom + "\n")
    return sites.sort(g=chrom_names.name)


def filter_bed_file(input_bed, filter_bed, output_bed):
    """
    Filter BED file for entries that overlap another BED file.

    Parameters
    ----------
    input_bed : :obj:`str`
        BED file to filter.

    filter_bed : :obj:`str`
        BED file with entries to filter from input_bed.

    output_bed : :obj:`str`
        Output BED file.
    """
    import pybedtools

    input_bed = get_this_file_or_timestamped(input_bed)

    (
        pybedtools.BedTool(input_bed)
        .intersect(pybedtools.BedTool(filter_bed), v=True)
        .saveas(output_bed)
    )


def timedelta_to_years(x):
    """
    Convert a timedelta to years.

    Parameters
    ----------
    x : :class:`pandas.Timedelta`
        A Timedelta object.

    Returns
    -------
    float
        Years.
    """
    return x / np.timedelta64(60 * 60 * 24 * 365, "s")


def signed_max(x, f=0.66, axis=0):
    """
    Return maximum or minimum of array ``x`` depending on the sign of the
    majority of values.
    If there isn't a clear majority (at least ``f`` fraction in one side),
    return mean of values.
    If given a pandas DataFrame or 2D numpy array, will apply this across rows
    (columns-wise, axis=0) or across columns (row-wise, axis=1).
    Will return :obj:`numpy.nan` for non-numeric values.

    Parameters
    ----------
    x : {:class:`numpy.ndarray`, :class:`pandas.DataFrame`, :class:`pandas.Series`}
        Input values to reduce
    f : :obj:`float`
        Threshold fraction of majority agreement.

        Default is 0.66.
    axis : :obj:`int`
        Whether to apply across rows (0, column-wise) or
        across columns (1, row-wise).

        Default is 0.

    Returns
    -------
    :class:`pandas.Series`
        Pandas Series with values reduced to the signed maximum.
    """
    if axis not in [0, 1]:
        raise ValueError("Axis must be one of 0 (columns) or 1 (rows).")

    if len(x.shape) == 1:
        # Return nan if not numeric
        if x.dtype not in [np.float_, np.int_]:
            return np.nan

        types = [type(i) for i in x]
        if not all(types):
            return np.nan
        else:
            if types[0] not in [np.float_, float, np.int_, int]:
                return np.nan
        ll = float(len(x))
        neg = sum(x < 0)
        pos = sum(x > 0)
        obs_f = max(neg / ll, pos / ll)
        if obs_f >= f:
            if neg > pos:
                return min(x)
            else:
                return max(x)
        else:
            return np.mean(x)
    else:
        if not isinstance(x, pd.DataFrame):
            x = pd.DataFrame(x)

        if axis == 1:
            x = x.T
        res = pd.Series(np.empty(x.shape[1]), index=x.columns)
        for v in x.columns:
            res[v] = signed_max(x.loc[:, v], f=f)

        return res


def log_pvalues(x, f=0.1):
    """
    Calculate -log10(p-value) of array.

    Replaces infinite values with:

    .. highlight:: python
    .. code-block:: python

        max(x) + max(x) * f

    that is, fraction ``f`` more than the maximum non-infinite -log10(p-value).

    Parameters
    ----------
    x : :class:`pandas.Series`
        Series with numeric values
    f : :obj:`float`
        Fraction to augment the maximum value by if ``x`` contains infinite values.

        Defaults to 0.1.

    Returns
    -------
    :class:`pandas.Series`
        Transformed values.
    """
    ll = -np.log10(x)
    rmax = ll[ll != np.inf].max()
    return ll.replace(np.inf, rmax + rmax * f)


def fix_dataframe_header(df, force_dtypes=float):
    # First try to see whether it is not a MultiIndex after all
    if df.index.isna().any().sum() == 1:
        df.columns = df.loc[df.index.isna(), :].squeeze()
        df.columns.name = None
        df.index.name = None
        return df.loc[~df.index.isna()]
    # If so, get potential categorical columns and prepare MultiIndex with them
    cols = list()
    for i in df.index[:300]:
        try:
            df.loc[i, :].astype(float)
        except ValueError:
            cols.append(i)
    if len(cols) == 0:
        pass
    elif len(cols) == 1:
        df.columns = df.loc[cols[0]]
        df = df.loc[~df.index.isin(cols)]
    else:
        df.columns = pd.MultiIndex.from_arrays(df.loc[cols].values, names=cols)
        df = df.loc[~df.index.isin(cols)]
    df.index.name = None
    if force_dtypes is not None:
        return df.astype(force_dtypes)
    else:
        return df


def r2pandas_df(r_df):
    """Make :class:`pandas.DataFrame` from a ``R`` dataframe given by :class:`rpy`."""
    from rpy2.robjects import pandas2ri

    if isinstance(r_df, pd.DataFrame):
        return r_df
    if isinstance(r_df, np.recarray):
        return pd.DataFrame.from_records(r_df)
    df = pd.DataFrame(np.asarray(r_df)).T
    df.columns = [str(x) for x in r_df.colnames]
    df.index = [str(x) for x in r_df.rownames]
    return df


def recarray2pandas_df(recarray):
    """Make :class:`pandas.DataFrame` from :class:`numpy.recarray`."""
    df = pd.DataFrame.from_records(recarray, index=list(range(recarray.shape[0])))
    return df


# def detect_peaks(x, mph=None, mpd=1, threshold=0, edge='rising',
#                  kpsh=False, valley=False):
#     """Detect peaks in data based on their amplitude and other features.

#     Parameters
#     ----------
#     x : 1D array_like
#         data.
#     mph : {None, number},optional (default = None)
#         detect peaks that are greater than minimum peak height.
#     mpd : positive integer,optional (default = 1)
#         detect peaks that are at least separated by minimum peak distance (in
#         number of data).
#     threshold : positive number,optional (default = 0)
#         detect peaks (valleys) that are greater (smaller) than `threshold`
#         in relation to their immediate neighbors.
#     edge : {None, 'rising', 'falling', 'both'},optional (default = 'rising')
#         for a flat peak, keep only the rising edge ('rising'), only the
#         falling edge ('falling'), both edges ('both'), or don't detect a
#         flat peak (None).
#     kpsh: :obj:`bool`,optional (default = False)
#         keep peaks with same height even if they are closer than `mpd`.
#     valley: :obj:`bool`,optional (default = False)
#         if True (1), detect valleys (local minima) instead of peaks.
#     show: :obj:`bool`,optional (default = False)
#         if True (1), plot data in matplotlib figure.
#     ax : a matplotlib.axes.Axes instance,optional (default = None).

#     Returns
#     -------
#     ind : 1D array_like
#         indeces of the peaks in `x`.

#     Notes
#     -----
#     The detection of valleys instead of peaks is performed internally by simply
#     negating the data: `ind_valleys = detect_peaks(-x)`

#     The function can handle NaN's

#     See this IPython Notebook [1]_.

#     References
#     ----------
#     .. [1] http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb

#     Examples
#     --------
#     >>> from detect_peaks import detect_peaks
#     >>> x = np.random.randn(100)
#     >>> x[60:81] = np.nan
#     >>> # detect all peaks and plot data
#     >>> ind = detect_peaks(x, show=True)
#     >>> print(ind)

#     >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
#     >>> # set minimum peak height = 0 and minimum peak distance = 20
#     >>> detect_peaks(x, mph=0, mpd=20, show=True)

#     >>> x = [0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0]
#     >>> # set minimum peak distance = 2
#     >>> detect_peaks(x, mpd=2, show=True)

#     >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
#     >>> # detection of valleys instead of peaks
#     >>> detect_peaks(x, mph=0, mpd=20, valley=True, show=True)

#     >>> x = [0, 1, 1, 0, 1, 1, 0]
#     >>> # detect both edges
#     >>> detect_peaks(x, edge='both', show=True)

#     >>> x = [-2, 1, -2, 2, 1, 1, 3, 0]
#     >>> # set threshold = 2
#     >>> detect_peaks(x, threshold = 2, show=True)

#     x = np.atleast_1d(x).astype('float64')
#     if x.size < 3:
#         return np.array([], dtype=int)
#     if valley:
#         x = -x
#     # find indices of all peaks
#     dx = x[1:] - x[:-1]
#     # handle NaN's
#     indnan = np.where(np.isnan(x))[0]
#     if indnan.size:
#         x[indnan] = np.inf
#         dx[np.where(np.isnan(dx))[0]] = np.inf
#     ine, ire, ife = np.array([[], [], []], dtype=int)
#     if not edge:
#         ine = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) > 0))[0]
#     else:
#         if edge.lower() in ['rising', 'both']:
#             ire = np.where((np.hstack((dx, 0)) <= 0) & (np.hstack((0, dx)) > 0))[0]
#         if edge.lower() in ['falling', 'both']:
#             ife = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) >= 0))[0]
#     ind = np.unique(np.hstack((ine, ire, ife)))
#     # handle NaN's
#     if ind.size and indnan.size:
#         # NaN's and values close to NaN's cannot be peaks
#         ind = ind[np.in1d(ind, np.unique(np.hstack((indnan, indnan - 1, indnan + 1))), invert=True)]
#     # first and last values of x cannot be peaks
#     if ind.size and ind[0] == 0:
#         ind = ind[1:]
#     if ind.size and ind[-1] == x.size - 1:
#         ind = ind[:-1]
#     # remove peaks < minimum peak height
#     if ind.size and mph is not None:
#         ind = ind[x[ind] >= mph]
#     # remove peaks - neighbors < threshold
#     if ind.size and threshold > 0:
#         dx = np.min(np.vstack([x[ind] - x[ind - 1], x[ind] - x[ind + 1]]), axis=0)
#         ind = np.delete(ind, np.where(dx < threshold)[0])
#     # detect small peaks closer than minimum peak distance
#     if ind.size and mpd > 1:
#         ind = ind[np.argsort(x[ind])][::-1]  # sort ind by peak height
#         idel = np.zeros(ind.size, dtype=bool)
#         for i in range(ind.size):
#             if not idel[i]:
#                 # keep peaks with the same height if kpsh is True
#                 idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
#                     & (x[ind[i]] > x[ind] if kpsh else True)
#                 idel[i] = 0  # Keep current peak
#         # remove the small peaks and sort back the indices by their occurrence
#         ind = np.sort(ind[~idel])

#     return ind


def collect_md5_sums(df):
    """
    Given a dataframe with columns with paths to md5sum files ending in '_md5sum',
    replace the paths to the md5sum files with the actual checksums.

    Useful to use in combination with :func:`~ngs_toolkit.general.project_to_geo`.

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        A dataframe with columns ending in '_md5sum'.

    Returns
    -------
    :class:`pandas.DataFrame`
        DataFrame with md5sum columns replaced with the actual md5sums.
    """
    cols = df.columns[df.columns.str.endswith("_md5sum")]
    for col in cols:
        for i, path in df.loc[:, col].iteritems():
            if not pd.isnull(path):
                cont = open(path, "r").read().strip()
                if any([x.isspace() for x in cont]):
                    cont = cont.split(" ")[0]
                df.loc[i, col] = cont
    return df


def decompress_file(file, output_file=None):
    """
    Decompress a gzip-compressed file out-of-memory.
    Output default is same as ``file`` without ".gz" ending.
    """
    import shutil
    import gzip

    from ngs_toolkit import _LOGGER

    if output_file is None:
        if not file.endswith(".gz"):
            msg = "`output_file` not given and input_file does not end in '.gz'."
            _LOGGER.error(msg)
            raise ValueError(msg)
        output_file = file.replace(".gz", "")
    # decompress
    with gzip.open(file, "rb") as _in:
        with open(output_file, "wb") as _out:
            shutil.copyfileobj(_in, _out)
            # for line in _in.readlines():
            # _out.write(line.decode('utf-8'))


def compress_file(file, output_file=None):
    """
    Compress a gzip-compressed file out-of-memory.
    Output default is same as ``file`` but with ".gz" ending.
    """
    import shutil
    import gzip

    if output_file is None:
        output_file = file + ".gz"
    # compress
    with open(file, "rb") as _in:
        with gzip.open(output_file, "wb") as _out:
            shutil.copyfileobj(_in, _out)
            # for line in _in.readlines():
            # _out.write(line.decode('utf-8'))


def download_file(url, output_file, chunk_size=1024):
    """
    Download a file and write to disk in chunks (not in memory).

    Parameters
    ----------
    url : :obj:`str`
        URL to download from.
    output_file : :obj:`str`
        Path to file as output.
    chunk_size : :obj:`int`
        Size in bytes of chunk to write to disk at a time.
    """
    if url.startswith("ftp://"):
        import shutil
        import urllib.request as request
        from contextlib import closing

        with closing(request.urlopen(url)) as r:
            with open(output_file, "wb") as f:
                shutil.copyfileobj(r, f)
    else:
        import requests

        response = requests.get(url, stream=True)
        with open(output_file, "wb") as outfile:
            outfile.writelines(response.iter_content(chunk_size=chunk_size))


def download_gzip_file(url, output_file, **kwargs):
    """
    Download a gzip compressed file and write uncompressed
    file to disk in chunks (not in memory).

    Parameters
    ----------
    url : :obj:`str`
        URL to download from.
    output_file : :obj:`str`
        Path to file as output.
    **kwargs : :obj:`dict`
        Additional keyword arguments are passed to :meth:`ngs_toolkit.utils.download_file`.
    """
    if not output_file.endswith(".gz"):
        output_file += ".gz"
    download_file(url, output_file, **kwargs)
    decompress_file(output_file)
    if os.path.exists(output_file) and os.path.exists(output_file.replace(".gz", "")):
        os.remove(output_file)


def deseq_results_to_bed_file(
    deseq_result_file,
    bed_file,
    sort=True,
    ascending=False,
    normalize=False,
    significant_only=False,
    alpha=0.05,
    abs_fold_change=1.0,
):
    """
    Write BED file with fold changes from DESeq2 as score value.
    """
    from ngs_toolkit import _LOGGER
    from sklearn.preprocessing import MinMaxScaler

    df = pd.read_csv(deseq_result_file, index_col=0)

    msg = "DESeq2 results do not have a 'log2FoldChange' column."
    if not ("log2FoldChange" in df.columns.tolist()):
        _LOGGER.error(msg)
        raise AssertionError(msg)

    if sort is True:
        df = df.sort_values("log2FoldChange", ascending=ascending)

    if significant_only is True:
        df = df.loc[(df["padj"] < alpha) & (df["log2FoldChange"].abs() > abs_fold_change), :]

    # decompose index string (chrom:start-end) into columns
    df["chrom"] = map(lambda x: x[0], df.index.str.split(":"))
    r = pd.Series(map(lambda x: x[1], df.index.str.split(":")))
    df["start"] = map(lambda x: x[0], r.str.split("-"))
    df["end"] = map(lambda x: x[1], r.str.split("-"))
    df["name"] = df.index
    if normalize:
        MinMaxScaler(feature_range=(0, 1000)).fit_transform(df["log2FoldChange"])
    df["score"] = df["log2FoldChange"]

    df[["chrom", "start", "end", "name", "score"]].to_csv(
        bed_file, sep="\t", header=False, index=False
    )


def homer_peaks_to_bed(homer_peaks, output_bed):
    """
    Convert HOMER peak calls to BED format.
    The fifth column (score) is the -log10(p-value) of the peak.

    Parameters
    ----------
    homer_peaks : :obj:`str`
        HOMER output with peak calls.
    output_bed : :obj:`str`
        Output BED file.
    """
    df = pd.read_csv(homer_peaks, sep="\t", comment="#", header=None)
    df["-log_pvalue"] = (-np.log10(df.iloc[:, -2])).replace(pd.np.inf, 1000)
    df["name"] = df[1] + ":" + df[2].astype(str) + "-" + df[3].astype(str)

    (
        df[[1, 2, 3, "name", "-log_pvalue"]]
        .sort_values([1, 2, 3], ascending=True)
        .to_csv(output_bed, header=False, index=False, sep="\t")
    )


def macs2_call_chipseq_peak(signal_samples, control_samples, output_dir, name, distributed=True):
    """
    Call ChIP-seq peaks with MACS2 in a slurm job.

    Parameters
    ----------
    signal_samples : :obj:`list`
        Signal Sample objects.

    control_samples : :obj:`list`
        Background Sample objects.

    output_dir : :obj:`list`
        Parent directory where MACS2 outputs will be stored.

    name : :obj:`str`
        Name of the MACS2 comparison being performed.

    distributed: :obj:`bool`
        Whether to submit a SLURM job or to return a string with the runnable.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    runnable = """macs2 callpeak -t {0} -c {1} -n {2} --outdir {3}""".format(
        " ".join([s.aligned_filtered_bam for s in signal_samples]),
        " ".join([s.aligned_filtered_bam for s in control_samples]),
        name,
        output_dir,
    )

    if distributed:
        job_name = "macs2_{}".format(name)
        job_file = os.path.join(output_dir, name + ".macs2.sh")
        submit_job(runnable, job_file, cores=4, jobname=job_name)
    else:
        return runnable


def homer_call_chipseq_peak_job(
    signal_samples, control_samples, output_dir, name, distributed=True
):
    """
    Call ChIP-seq peaks with MACS2 in a slurm job.

    Parameters
    ----------
    signal_samples : :obj:`list`
        Signal Sample objects.

    control_samples : :obj:`list`
        Background Sample objects.

    output_dir : :obj:`list`
        Parent directory where MACS2 outputs will be stored.

    name : :obj:`str`
        Name of the MACS2 comparison being performed.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # make tag directory for the signal and background samples separately
    signal_tag_directory = os.path.join(output_dir, "homer_tag_dir_" + name + "_signal")
    fs = " ".join([s.aligned_filtered_bam for s in signal_samples])
    runnable = """makeTagDirectory {0} {1}\n""".format(signal_tag_directory, fs)
    background_tag_directory = os.path.join(output_dir, "homer_tag_dir_" + name + "_background")
    fs = " ".join([s.aligned_filtered_bam for s in control_samples])
    runnable += """makeTagDirectory {0} {1}\n""".format(background_tag_directory, fs)

    # call peaks
    output_file = os.path.join(output_dir, name, name + "_homer_peaks.factor.narrowPeak")
    runnable += """findPeaks {signal} -style factor -o {output_file} -i {background}\n""".format(
        output_file=output_file, background=background_tag_directory, signal=signal_tag_directory,
    )
    output_file = os.path.join(output_dir, name, name + "_homer_peaks.histone.narrowPeak")
    runnable += """findPeaks {signal} -style histone -o {output_file} -i {background}""".format(
        output_file=output_file, background=background_tag_directory, signal=signal_tag_directory,
    )

    if distributed:
        job_name = "homer_findPeaks_{}".format(name)
        job_file = os.path.join(output_dir, name + ".homer.sh")
        submit_job(runnable, job_file, cores=4, jobname=job_name)
    else:
        return runnable


def bed_to_fasta(input_bed, output_fasta, genome_file):
    """
    Retrieves DNA sequence underlying specific region.
    Names of FASTA entries will be of form ``chr:start-end``.

    Parameters
    ----------
    input_bed : :obj:`str`
        Path to input BED file.

    output_fasta : :obj:`str`
        Path to resulting FASTA file.

    genome_file : :obj:`str`
        Path to genome file in either 2bit or FASTA format.
        Will be guessed based on file ending.

    Raises
    ----------
    ValueError
        If `genome_file` format cannot be guessed or is not supported.
    """

    if genome_file.endswith(".2bit"):
        bed_to_fasta_through_2bit(input_bed, output_fasta, genome_file)
    elif genome_file.endswith(".fa") or genome_file.endswith(".fasta"):
        bed_to_fasta_through_fasta(input_bed, output_fasta, genome_file)
    else:
        msg = "Format of `genome_file` must be one of FASTA or 2bit, "
        msg += "with file ending in either '.fa', '.fasta' or '.2bit'"
        raise ValueError(msg)


def read_bed_file_three_columns(input_bed: str) -> pd.DataFrame:
    """Read BED file into dataframe, make 'name' field from location."""
    bed = pd.read_csv(
        input_bed, sep="\t", header=None, usecols=[0, 1, 2], names=["chrom", "start", "end"]
    )
    bed["name"] = bed_to_index(bed)
    return bed


def bed_to_fasta_through_2bit(input_bed, output_fasta, genome_2bit):
    """
    Retrieves DNA sequence underlying specific region.
    Requires the `twoBitToFa` command from UCSC kent tools.
    Names of FASTA entries will be of form "chr:start-end".

    Parameters
    ----------
    input_bed : :obj:`str`
        Path to input BED file.

    output_fasta : :obj:`str`
        Path to resulting FASTA file.

    genome_2bit : :obj:`str`
        Path to genome 2bit file.
    """
    import subprocess

    tmp_bed = input_bed + ".tmp.bed"
    bed = read_bed_file_three_columns(input_bed)
    bed.to_csv(tmp_bed, sep="\t", header=None, index=False)

    cmd = "twoBitToFa {0} -bed={1} {2}".format(genome_2bit, tmp_bed, output_fasta)
    subprocess.call(cmd.split(" "))
    os.remove(tmp_bed)


def bed_to_fasta_through_fasta(input_bed, output_fasta, genome_fasta):
    """
    Retrieves DNA sequence underlying specific region.
    Uses ``bedtools getfasta`` (internally through pybedtools.BedTool.sequence).
    Names of FASTA entries will be of form "chr:start-end".

    Parameters
    ----------
    input_bed : :obj:`str`
        Path to input BED file.

    output_fasta : :obj:`str`
        Path to resulting FASTA file.

    genome_fasta : :obj:`str`
        Path to genome FASTA file.
    """
    import pybedtools

    bed = read_bed_file_three_columns(input_bed)
    bed = pybedtools.BedTool.from_dataframe(bed)
    bed.sequence(fi=genome_fasta, fo=output_fasta, name=True)


def count_reads_in_intervals(bam, intervals, permissive=True):
    """
    Count total number of reads in a iterable holding strings
    representing genomic intervals of the form ``"chrom:start-end"``.

    Please make sure both ``intervals`` and ``bam`` file are
    zero- or one-indexed.

    Parameters
    ----------
    bam : :obj:`str`
        Path to BAM file.

    intervals : :obj:`list`
        List of strings with genomic coordinates in format
        ``"chrom:start-end"``.

    Returns
    -------
    :obj:`dict`
        Dict of read counts for each interval.
    """
    import pysam
    from ngs_toolkit import _LOGGER

    counts = dict()

    bam = pysam.AlignmentFile(bam, mode="rb")

    errors: int = 0
    for interval in intervals:
        try:
            counts[interval] = bam.count(region=interval)
        except ValueError:
            if permissive:
                errors += 1
            else:
                raise
            # if fix_off_by_one:
            #     i = interval.split(":")[1]
            #     s = (
            #         interval.split(":")[0] +
            #         ":" + str(int(i.split("-")[0]) + 1) +
            #         "-" + str(int(i.split("-")[1]) + 1))
            #     counts[interval] = bam.count(region=s)
    bam.close()
    if errors > 0:
        _LOGGER.warning("There have been %i errors. Beware.", errors)

    return counts


def normalize_quantiles_r(array):
    """
    Quantile normalization with a R implementation.

    Requires the R package "preprocessCore" to be installed:

    .. highlight:: R
    .. code-block:: R

        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("preprocessCore")

    Parameters
    ----------
    array : :class:`numpy.ndarray`
        Numeric array to normalize.

    Returns
    -------
    :class:`numpy.ndarray`
        Normalized numeric array.
    """
    from rpy2.robjects import numpy2ri, pandas2ri, r
    from rpy2.robjects.packages import importr

    numpy2ri.activate()
    pandas2ri.activate()

    importr("preprocessCore")

    return np.array(r("normalize.quantiles")(array))


def normalize_quantiles_p(df_input):
    """
    Quantile normalization with a ure Python implementation.
    Code from https://github.com/ShawnLYU/Quantile_Normalize.

    Parameters
    ----------
    df_input : :class:`pandas.DataFrame`
        Dataframe to normalize of shape (features, samples).

    Returns
    -------
    :class:`numpy.ndarray`
        Normalized numeric array.
    """
    df = df_input.copy()
    # compute rank
    dic = {}
    for col in df:
        dic.update({col: sorted(df[col])})
    sorted_df = pd.DataFrame(dic)
    rank = sorted_df.mean(axis=1).tolist()
    # sort
    for col in df:
        t = np.searchsorted(np.sort(df[col]), df[col])
        df[col] = [rank[i] for i in t]
    return df


def cqn(matrix, gc_content, lengths):
    """
    Conditional quantile normalization (CQN) with the ``cqn`` R library.
    It uses GC content and length of regulatory elements as covariates.

    Requires the R package "cqn" to be installed:

    .. highlight:: R
    .. code-block:: R

        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install("cqn")

    Parameters
    ----------
    matrix : :class:`pandas.DataFrame`
        DataFrame to normalize.
    gc_content : :class:`pandas.Series`
        Series with GC content of each feature in ``matrix``.
    lengths : :class:`pandas.Series`
        Series with length of each feature in ``matrix``.

    Returns
    ----------
    :class:`pandas.DataFrame`
        Normalized DataFrame
    """
    from rpy2.robjects import numpy2ri, pandas2ri, r
    from rpy2.robjects.packages import importr

    numpy2ri.activate()
    pandas2ri.activate()

    importr("cqn")

    cqn_out = r.cqn(matrix, x=gc_content, lengths=lengths)

    y_r = cqn_out[list(cqn_out.names).index("y")]
    y = pd.DataFrame(np.array(y_r), index=matrix.index, columns=matrix.columns)
    offset_r = cqn_out[list(cqn_out.names).index("offset")]
    offset = pd.DataFrame(np.array(offset_r), index=matrix.index, columns=matrix.columns)

    return y + offset


def count_bam_file_length(bam_file: str) -> int:
    """Get length of BAM indexed file"""
    import pysam

    return pysam.AlignmentFile(bam_file).count()


def count_lines(file: str) -> int:
    """Count lines of plain text file"""
    i = -1
    with open(file, "r") as f:
        for i, _ in enumerate(f):
            pass
    return i + 1


def get_total_region_area(bed_file: str) -> int:
    """Get sum of BED records"""
    peaks = pd.read_csv(bed_file, sep="\t", header=None)
    return int((peaks.iloc[:, 2] - peaks.iloc[:, 1]).sum())


def get_region_lengths(bed_file: str) -> pd.Series:
    """Get length of each record in BED file"""
    peaks = pd.read_csv(bed_file, sep="\t", header=None)
    return peaks.iloc[:, 2] - peaks.iloc[:, 1]


def get_regions_per_chromosomes(bed_file: str) -> pd.Series:
    """Count record per chromosome in BED file"""
    peaks = pd.read_csv(bed_file, sep="\t", header=None)
    return peaks.iloc[:, 0].value_counts()
