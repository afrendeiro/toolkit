#!/usr/bin/env python

import os

import numpy as np
import pandas as pd

from ngs_toolkit import _CONFIG, _LOGGER
from ngs_toolkit.decorators import check_has_attributes

# Bugs:
# TODO: unsupervised_analysis plotting fails if only one sample in one group
# TODO: plot_differential_enrichment fails if only one comparison
# TODO: make plotting functions that take analysis as arguments track too

# Improvements:
# TODO: move plot_features to be under Analysis, abstract. Call in plot_differential
# TODO: Make from_pep intialization a static method returning an instance
# TODO: Actually use ngs_toolkit.constants
# TODO: Add PAGE as enrichment method
# TODO: Analysis.annotate_samples: reimplement to support CNV dict of resolutions
# TODO: Analysis.annotate_samples: implement connection to analysis' numeric_attributes or another way of preserving dtypes
# TODO: Add function to complete comparison_table information such as "comparison_genome" and "data_type" and call it when setting automatically
# TODO: Add function to create comparison_table from samples' group_attributes
# TODO: Make recipe to get all (or subset through CLI) resources
# TODO: add matrix_features to RNASeqAnalysis class
# TODO: make pandas dataframes indexable by Sample objects

# Testing:
# TODO: test having no config set
# TODO: test differential analysis with many factors
# TODO: test subproject initialization

# Ideas:
# TODO: idea: make Analysis.annotate() call both annotate_features(), annotate_samples() and their ancestors with `steps`
# TODO: Idea: merging analysis. If same type, merge matrices, otherwise use dicts?
# TODO: Idea: if genome of analysis is set, get required static files for that genome assembly automatically
# TODO: Idea: Analysis.load_data: get default output_map by having functions declare what they output perhaps also with a dict of kwargs to pass to pandas.read_csv
# TODO: for recipes.ngs_analysis implement steps, build lock file system to measure progress, allow continuing

# Code:
# TODO: add type hinting (this implies adding all imports up in the file)
# TODO: replace _LOGGER.message("Message: {}".format("value")) with _LOGGER.message("Message: %s", "value")


class Analysis(object):
    """
    Generic class holding functions and data from a typical NGS analysis.

    Other modules implement classes inheriting from this that in general contain
    data type-specific functions (e.g. :class:`~ngs_toolkit.atacseq.ATACSeqAnalysis`
    has a :func:`~ngs_toolkit.atacseq.ATACSeqAnalysis.get_consensus_sites` function to generate a peak consensus map).

    Objects of this type can be used to store data (e.g. dataframes), variables
    (e.g. paths to files or configurations) and can easily be filled with
    existing data using :func:`~ngs_toolkit.analysis.Analysis.load_data`
    for cross-environment portability, or serialized (saved to a file as a
    pickle object) for rapid loading in the same environment.
    See the :func:`~ngs_toolkit.analysis.Analysis.to_pickle`,
    :func:`~ngs_toolkit.analysis.Analysis.from_pickle` and
    :func:`~ngs_toolkit.analysis.Analysis.update` functions for this.

    Parameters
    ----------
    name : :obj:`str`, optional
        Name of the analysis.

        Defaults to "analysis".
    from_pep : :obj:`str`, optional
        PEP configuration file to initialize analysis from.
        The analysis will adopt as much attributes from the PEP as possible
        but keyword arguments passed at initialization will still have priority.
        Keyword arguments will be passed to :obj:`peppy.Project` if matching its
        initialization signature.

        Defaults to :obj:`None` (no PEP used).
    from_pickle : :obj:`str`, optional
        Pickle file of an existing serialized analysis object
        from which the analysis should be loaded.

        Defaults to :obj:`None` (will not load from pickle).
    root_dir : :obj:`str`, optional
        Base directory for the project.

        Defaults to current directory or to what is specified in PEP if :attr:`~ngs_toolkit.analysis.Analysis.from_pep`.
    data_dir : :obj:`str`, optional
        Directory containing processed data (e.g. by looper) that will
        be input to the analysis. This is in principle not required.

        Defaults to "data".
    results_dir : :obj:`str`, optional
        Directory to contain outputs produced by the analysis.

        Defaults to "results".
    prj : :class:`peppy.Project`, optional
        A :class:`peppy.Project` object that this analysis is tied to.

        Defaults to :obj:`None`.
    samples : :obj:`list`, optional
        List of :class:`peppy.Sample` objects that this analysis is tied to.

        Defaults to :obj:`None`.
    subset_to_data_type : :obj:`bool`, optional
        Whether to keep only samples that match the data type of the analysis.

        Defaults to :obj:`True`.
    kwargs : :obj:`dict`, optional
        Additional keyword arguments will simply be stored as object attributes.
    """

    _data_type = None

    def __init__(
        self,
        name=None,
        from_pep=False,
        from_pickle=False,
        root_dir=None,
        data_dir="data",
        results_dir="results",
        prj=None,
        samples=None,
        subset_to_data_type=True,
        **kwargs
    ):
        # Add given args and kwargs to object
        _LOGGER.debug("Adding given arguments to analysis object.")
        self.name = name
        self.root_dir = root_dir
        self.pep = None if not from_pep else from_pep
        self.prj = prj
        self.samples = samples

        # Add default values for matrices
        self.raw_matrix_name = "matrix_raw"
        self.norm_matrix_name = "matrix_norm"
        self.feature_matrix_name = "matrix_features"
        _LOGGER.debug("Adding additional kwargs to analysis object.")
        self.__dict__.update(kwargs)

        # If from_pickle, load and return
        if from_pickle is not False:
            _LOGGER.info("Updating analysis object from pickle file: '{}'.".format(from_pickle))
            self.update(pickle_file=from_pickle)
            return None

        _LOGGER.debug("Setting data type-specific attributes to None.")
        attrs = [
            "data_type",
            "__data_type__",
            "var_unit_name",
            "quantity",
            "norm_units",
            "raw_matrix_name",
            "norm_matrix_name",
            "annot_matrix_name",
            "norm_method",
        ]
        for attr in attrs:
            if not hasattr(self, attr):
                setattr(self, attr, None)
        if not hasattr(self, "thresholds"):
            self.thresholds = {"alpha": 0.05, "log2_fold_change": 0}
        if not hasattr(self, "output_files"):
            self.output_files = list()

        # Generate from PEP configuration file
        if from_pep is not False:
            from peppy import Project
            from ngs_toolkit.utils import filter_kwargs_by_callable

            self.from_pep(pep_config=from_pep, **filter_kwargs_by_callable(kwargs, Project))

        # Store projects attributes in self
        _LOGGER.debug("Trying to set analysis attributes.")
        self.set_project_attributes(overwrite=False, subset_to_data_type=subset_to_data_type)

        # Get name
        if self.name is None:
            self.name = "analysis"

        # Try to set genome if not set
        self.organism, self.genome = (None, None)
        _LOGGER.debug("Trying to get analysis genome.")
        self.set_organism_genome()

        # Set root_dir
        if self.root_dir is None:
            self.root_dir = os.curdir
        self.root_dir = os.path.abspath(self.root_dir)

        # Set data and results dir
        self.data_dir = os.path.join(self.root_dir, data_dir)
        self.results_dir = os.path.join(self.root_dir, results_dir)

        # # if given absolute paths, keep them, otherwise append to root directory
        for _dir, attr in [(data_dir, "data_dir"), (results_dir, "results_dir")]:
            if not os.path.isabs(_dir):
                _dir = os.path.join(self.root_dir, _dir)
            setattr(self, attr, _dir)

        # Try to make directories
        for _dir in [self.data_dir, self.results_dir]:
            if not os.path.exists(_dir):
                try:
                    os.makedirs(_dir)
                except OSError:
                    _LOGGER.debug("Could not make directory for Analysis: '%s'", _dir)

        # Add sample input file locations
        _LOGGER.debug("Trying to set sample input file attributes.")
        self.set_samples_input_files()

        # Set pickle file
        if not hasattr(self, "pickle_file"):
            self.pickle_file = os.path.join(self.results_dir, self.name + ".pickle")

    def __repr__(self):
        t = "'{}' analysis".format(self.data_type) if self.data_type is not None else "Analysis"
        samples = " with {} samples".format(len(self.samples)) if self.samples is not None else ""
        organism = " of organism '{}'".format(self.organism) if self.organism is not None else ""
        genome = " ({})".format(self.genome) if self.genome is not None else ""
        suffix = "."
        return t + " '{}'".format(self.name) + samples + organism + genome + suffix

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass

    @staticmethod
    def _overwride_sample_representation():
        """
        Make :class:`peppy.Sample` objects have a more pretty representation.
        """
        from peppy import Sample

        def r(self):
            return self.name

        Sample.__repr__ = r

    @staticmethod
    def _check_data_type_is_supported(data_type):
        """
        Check if data type is allowed in ngs_toolkit configuration.

        Parameters
        ----------
        data_type : :obj:`str`
            Data type to check.

        Returns
        ----------
        :obj:`bool`
            Whether data_type is supported.
        """
        supported = _CONFIG["supported_data_types"]
        return data_type in supported

    def _get_data_type(self, data_type=None):
        """
        Get data type of current object, throwing an error if
        data_type is not allowed in ngs_toolkit configuration.

        Parameters
        ----------
        data_type : :obj:`str`, optional
            Data type to check.

        Returns
        ----------
        :obj:`str`
            A supported data_type.
        """
        if data_type is None:
            msg = "Data type not defined and Analysis object does not have"
            msg += " a `data_type` attribute."
            try:
                data_type = self.data_type
            except AttributeError as e:
                _LOGGER.error(msg)
                raise e
            if data_type is None:
                _LOGGER.error(msg)
                raise ValueError(msg)
        else:
            msg = "Data type is not supported."
            hint = " Check which data types are supported in the 'supported_data_types'"
            hint += " section of the configuration file."
            if not self._check_data_type_is_supported(data_type):
                raise ValueError(msg + hint)
        return data_type

    def _check_samples_have_file(self, attr, f=all, samples=None):
        """
        Checks that there are existing files for an attribute of the analysis samples.
        This requires a reducing function such as 'all' or 'any' to evaluate how to
        return a value across all samples.

        Parameters
        ----------
        attr : :obj:`str`
            An attribute of the analysis' samples to check existence of files.
        f : function, optional
            Function to reduce output across samples.

            Defaults to :obj:`all`.
        samples : :obj:`list`, optional
            Samples to consider.

            Defaults to all in analysis.

        Returns
        ----------
        :obj:`bool`
            Whether samples have file.

        Raises
        -------
        :obj:`AttributeError`
            If attribute does not exist in samples
        """
        if samples is None:
            samples = self.samples

        try:
            files = [str(getattr(sample, attr)) for sample in samples]
        except AttributeError:
            msg = "Sample did not have attribute '{}'".format(attr)
            _LOGGER.error(msg)
            raise

        return f([os.path.exists(file) for file in files])

    def _get_samples_have_file(self, attr, samples=None):
        """
        Get samples with an existing file under `attr`.

        Parameters
        ----------
        attr : :obj:`str`
            Attribute to check
        samples : :obj:`list`, optional
            Samples to consider.

            Defaults to all in analysis.

        Returns
        -------
        list
            List of :class:`peppy.Sample` objects.

        Raises
        -------
        :obj:`AttributeError`
            If attribute does not exist in samples
        """
        if samples is None:
            samples = self.samples
        return [sample for sample in samples if os.path.exists(str(getattr(sample, attr)))]

    def _get_samples_missing_file(self, attr, samples=None):
        """
        Get samples without an existing file under `attr`.

        Parameters
        ----------
        attr : :obj:`str`
            Attribute to check

        samples : :obj:`list`, optional
            Samples to consider.

            Defaults to all in analysis.

        Returns
        -------
        list
            List of :class:`peppy.Sample` objects.

        Raises
        -------
        :obj:`AttributeError`
            If attribute does not exist in samples
        """
        if samples is None:
            samples = self.samples
        return [
            sample
            for sample in samples
            if sample not in self._get_samples_have_file(attr, samples=samples)
        ]

    def _get_samples_with_input_file(self, input_file, permissive=False, samples=None):
        """
        Get samples with existing files of attribute `input_file`.

        If none has, raise error. Else, if permissive, return samples with existing file.
        Otherwise return only with all samples have it, otherwise throw IOError.

        Parameters
        ----------
        input_file : :obj:`str`
            Attribute to check.
        permissive : :obj:`bool`, optional
            Whether to allow returning a subset of samples if not all have file.

            Defaults to :obj:`False`.
        samples : :obj:`list`, optional
            Samples to consider.

            Defaults to all in analysis.

        Returns
        -------
        list
            List of :class:`peppy.Sample` objects.

        Raises
        -------
        :obj:`IOError`
            If not permissive and not all sample input files are found.
        """
        if samples is None:
            samples = self.samples
        check = self._check_samples_have_file(attr=input_file, f=all, samples=samples)
        if check:
            return samples

        missing = self._get_samples_missing_file(attr=input_file, samples=samples)

        msg = "None of the samples have '{}' files.".format(input_file)
        if all([s in missing for s in samples]):
            if permissive:
                _LOGGER.warning(msg)
                return []
            else:
                _LOGGER.error(msg)
                raise IOError(msg)

        msg = "Not all samples have '{}' files.".format(input_file)
        hint = " Samples missing files: '{}'".format(", ".join([s.name for s in missing]))
        if permissive:
            _LOGGER.warning(msg + hint)
            return [s for s in samples if s not in missing]
        else:
            _LOGGER.error(msg + hint)
            raise IOError(msg)

    @staticmethod
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
        from ngs_toolkit.utils import _format_string_with_environment_variables

        return _format_string_with_environment_variables(string)

    def _format_string_with_attributes(self, string):
        """
        Given a string, containing curly braces, format it with the attributes from self.

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
            If not all patterns are analysis variables.
        """
        if string is None:
            return string
        to_format = pd.Series(string).str.extractall(r"{(.*?)}")[0].values
        attrs = self.__dict__.keys()
        if not all([x in attrs for x in to_format]):
            msg = "Not all required patterns were found as attributes of object '{}'.".format(self)
            _LOGGER.error(msg)
            raise ValueError(msg)
        return string.format(**self.__dict__)

    def from_pep(self, pep_config, **kwargs):
        """
        Create a peppy.Project from a PEP configuration file
        and associate is with the analysis.

        Parameters
        ----------
        pep_config : :obj:`str`
            PEP configuration file.

        Attributes
        ----------
        prj : :obj:`peppy.Project`
            Project object from given PEP configuration file.
        """
        import peppy

        # peppy.project.logging.disable()
        self.prj = peppy.Project(cfg=pep_config, **kwargs)

    def update(self, pickle_file=None):
        """
        Update all of the object"s attributes with the attributes
        from a serialized object (ie object stored in a file) object.

        Parameters
        ----------
        pickle_file : :obj:`str`, optional
            Pickle file to load.

            Defaults to the analysis' ``pickle_file``.
        """
        self.__dict__.update(self.from_pickle(pickle_file=pickle_file).__dict__)

    def set_organism_genome(self):
        """
        Attempt to derive the analysis' organism and genome assembly
        by inspecting the same attributes of its samples.

        Attributes
        ----------
        organism : :obj:`str`
            Organism of the analysis if all samples agree in these attributes.
        genome : :obj:`str`
            Genome assembly of the analysis if all samples agree in these attributes.
        """
        if self.samples is None or self.samples == []:
            _LOGGER.warning(
                "Genome assembly for analysis was not set and cannot be derived from samples."
            )
        else:
            hint = "Will not set an organism for analysis."
            organisms = [
                x for x in {getattr(s, "organism", None) for s in self.samples} if x is not None
            ]
            if len(organisms) == 1:
                _LOGGER.info("Setting analysis organism as '{}'.".format(organisms[0]))
                self.organism = organisms[0]
            elif organisms:
                msg = "Did not found any organism in the analysis samples. "
                _LOGGER.warning(msg + hint)
            else:
                msg = "Found several organism for the various analysis samples. "
                _LOGGER.warning(msg + hint)

            hint = "Will not set a genome for analysis."
            genomes = [
                x for x in {getattr(s, "genome", None) for s in self.samples} if x is not None
            ]
            if len(genomes) == 1:
                _LOGGER.info("Setting analysis genome as '{}'.".format(genomes[0]))
                self.genome = genomes[0]
            elif genomes:
                msg = "Did not found any genome assembly in the analysis samples. "
                _LOGGER.warning(msg + hint)
            else:
                msg = "Found several genome assemblies for the various analysis samples. "
                _LOGGER.warning(msg + hint)

    def set_project_attributes(self, overwrite=True, subset_to_data_type=True):
        """
        Set Analysis object attributes ``samples``, ``sample_attributes``
        and ``group_atrributes`` to the values in the associated Project
        object if existing.

        Parameters
        ----------
        overwrite: :obj:`bool`, optional
            Whether to overwrite attribute values if existing.

            Defaults to :obj:`True`.
        subset_to_data_type: :obj:`bool`, optional
            Whether to subset samples and comparison_table to entries
            of same ``data_type`` as analysis.

            Defaults to :obj:`True`.

        Attributes
        ----------
        samples : :obj:`list`
            List of peppy.Samples if contained in the PEP configuration.
        sample_attributes : :obj:`list`
            Sample attributes if specified in the PEP configuration.
        group_attributes : :obj:`list`
            Groups attributes if specified in the PEP configuration.
        comparison_table : :obj:`pandas.DataFrame`
            Comparison table if specified in the PEP configuration.
        """
        hint = " Adding a '{}' section to your project configuration file allows the analysis"
        hint += " object to use those attributes during the analysis."
        if self.prj is None:
            _LOGGER.warning(
                "Analysis object does not have an attached Project. "
                + "Will not add special attributes to analysis such as "
                + "samples, their attributes and comparison table."
            )
            return

        try:
            self.prj.root_dir = self.prj._config.root_dir
        except AttributeError:
            tmp = os.path.dirname(self.prj._config.sample_table)
            if os.path.basename(tmp) == "metadata":
                tmp = os.path.abspath(os.path.join(tmp, ".."))
            self.prj.root_dir = tmp

        for attr, parent in [
            ("name", self.prj),
            ("root_dir", self.prj._config),
            ("samples", self.prj),
            ("sample_attributes", self.prj._config),
            ("group_attributes", self.prj._config),
            ("comparison_table", self.prj._config),
        ]:
            if not hasattr(parent, attr):
                _LOGGER.warning(
                    "Associated project does not have any '{}'.".format(attr) + hint.format(attr)
                    if attr != "samples"
                    else ""
                )
            else:
                msg = "Setting project's '{0}' as the analysis '{0}'.".format(attr)
                if overwrite:
                    _LOGGER.info(msg)
                    setattr(self, attr, getattr(parent, attr))
                else:
                    if not hasattr(self, attr):
                        _LOGGER.info(msg)
                        setattr(self, attr, getattr(parent, attr))
                    else:
                        if getattr(self, attr) is None:
                            setattr(self, attr, getattr(parent, attr))
                        else:
                            _LOGGER.debug(
                                "{} already exist for analysis, not overwriting.".format(
                                    attr.replace("_", " ").capitalize()
                                )
                            )

        if hasattr(self, "comparison_table"):
            if isinstance(getattr(self, "comparison_table"), str):
                _LOGGER.debug("Reading up comparison table.")
                self.comparison_table = pd.read_csv(self.comparison_table)

        if subset_to_data_type:
            if hasattr(self, "samples"):
                if self.data_type is not None:
                    _LOGGER.info(
                        "Subsetting samples for samples of type '{}'.".format(self.data_type)
                    )
                    self.samples = [s for s in self.samples if s.protocol == self.data_type]
            if hasattr(self, "comparison_table"):
                if (
                    (self.comparison_table is not None)
                    and (self.data_type is not None)
                    and ("data_type" in self.comparison_table.columns)
                ):
                    _LOGGER.info(
                        "Subsetting comparison_table for comparisons of type '{}'.".format(
                            self.data_type
                        )
                    )
                    self.comparison_table.query("data_type == @self.data_type", inplace=True)

    def set_samples_input_files(self, overwrite=True):
        """
        Add input file values to sample objects dependent on data type.
        These are specified in the ``ngs_toolkit`` configuration file
        under ``sample_input_files:<data type>:<attribute>``.

        Parameters
        ----------
        overwrite: :obj:`bool`, optional
            Whether to overwrite attribute values if existing.

            Defaults to :obj:`True`.
        """

        def _format_string_with_sample_attributes(sample, string):
            """
            Given a string, containing curly braces, format it with the
            attributes from the sample object or the analysis.
            """
            if string is None:
                return string

            # fix around sample.__dict__ being overwritten
            sample_dict = {x: sample[x] for x in sample}

            to_format = pd.Series(string).str.extractall(r"{(.*?)}")[0].values
            attrs = [x for x in sample_dict.keys()]

            if not all([x in attrs for x in to_format]):
                d = sample_dict.copy()
                d.update(self.__dict__)
                return string.format(**d)
            else:
                return string.format(**sample_dict)

        def _set(attr, value, obj, sample):
            # This first bit here is to handle recursively
            # cases where value is a dict or (PathExAttMap in fact)
            if isinstance(value, dict):
                setattr(obj, attr, value)
                for k, v in getattr(obj, attr).items():
                    _set(k, v, getattr(obj, attr), sample)
                return
            try:
                value = _format_string_with_sample_attributes(sample, value)
            except KeyError:
                _LOGGER.error(
                    "Failed formatting for sample '{}', value '{}'.".format(sample.name, value)
                )
            if overwrite:
                _LOGGER.debug(msg.format(attr, sample.name, value))
                setattr(obj, attr, value)
            else:
                if not hasattr(sample, attr):
                    _LOGGER.debug(msg.format(attr, sample.name, value))
                    setattr(obj, attr, value)
                else:
                    if getattr(sample, attr) is None:
                        _LOGGER.debug(msg.format(attr, sample.name, value))
                        setattr(obj, attr, value)
                    else:
                        _LOGGER.debug(
                            "{} already exists in sample, not overwriting.".format(
                                attr.replace("_", " ").capitalize()
                            )
                        )

        if self.samples is None or self.samples == []:
            _LOGGER.warning(
                "Analysis object does not have attached Samples. "
                + "Will not add special attributes to samples such as "
                + "input file locations."
            )
            return

        for sample in self.samples:
            sample.sample_root = os.path.join(
                sample.project.root_dir, sample.project._config.results_subdir, sample.sample_name
            )

        msg = "Setting '{}' in sample {} as '{}'."
        for data_type in _CONFIG["sample_input_files"]:
            for sample in [s for s in self.samples if s.protocol == data_type]:
                if ("name" in sample) and ("sample_name" not in sample):
                    sample.sample_name = sample.name
                if ("sample_name" in sample) and ("name" not in sample):
                    sample.name = sample.sample_name
                for attr, value in _CONFIG["sample_input_files"][data_type].items():
                    _set(attr, value, sample, sample)

    def to_pickle(self, timestamp=False):
        """
        Serialize object (ie save to disk) to pickle format.

        Parameters
        ----------
        timestamp: :obj:`bool`, optional
            Whether to timestamp the file.

            Defaults to :obj:`False`.
        """
        import datetime
        import time
        import pickle

        if timestamp:
            ts = datetime.datetime.fromtimestamp(time.time()).strftime("%Y%m%d-%H%M%S")
            p = self.pickle_file.replace(".pickle", ".{}.pickle".format(ts))
        else:
            p = self.pickle_file
        pickle.dump(self, open(p, "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    def from_pickle(self, pickle_file=None):
        """
        Load object from pickle file.

        Parameters
        ----------
        pickle_file : :obj:`str`, optional
            Pickle file to load.

            Default is the object's attribute ``pickle_file``.

        Returns
        -------
        :class:`~ngs_toolkit.Analysis`
            The analysis serialized in the pickle file.
        """
        import pickle

        if pickle_file is None:
            pickle_file = self.pickle_file
        return pickle.load(open(pickle_file, "rb"))

    def get_sample_annotation(self, attributes=None, samples=None):
        """
        Get dataframe annotation of sample attributes.

        Attributes
        -------
        attributes : :obj:`None`, optional
            Attributes to include.

            Defaults to the union of sample_attributes and group_attributes in Analysis.
        samples : :obj:`None`, optional
            Samples to subset.

            Defaults to all samples in Analysis.

        Returns
        -------
        pandas.DataFrame
            Dataframe with requested attributes (columns) for each sample (rows).
        """
        if attributes is None:
            attributes = list(set(self.sample_attributes).union(set(self.group_attributes)))

        if samples is None:
            samples = self.samples

        try:
            v = [getattr(s, p, np.nan) for s in samples for p in attributes]
        except AttributeError:
            msg = "All samples must have all attributes specified!"
            _LOGGER.error(msg)
            raise AttributeError(msg)
        df = pd.DataFrame(
            np.array(v).reshape(len(samples), len(attributes)),
            index=[s.name for s in samples],
            columns=attributes,
        )
        return df

    def load_data(
        self, output_map=None, only_these_keys=None, prefix="{results_dir}/{name}", permissive=True,
    ):
        """
        Load the output files of the major functions of the Analysis.

        Parameters
        ----------
        output_map : :obj:`dict`
            Dictionary with {attribute_name: (file_path, kwargs)} to load the files.
            The kwargs in the tuple will be passed to :func:`pandas.read_csv`.

            Default is the required to read the keys in ``only_these_keys``.
        only_these_keys : :obj:`list`, optional
            Iterable of analysis attributes to load up.
            Possible attributes:

                * "matrix_raw"
                * "matrix_norm"
                * "matrix_features"
                * "differential_results"
                * "differential_enrichment"

            Default is all of the above.
        prefix : :obj:`str`, optional
            String prefix of files to load.
            Variables in curly braces will be formated with attributes of analysis.

            Default is "{results_dir}/{name}".
        permissive : :obj:`bool`, optional
            Whether an error should be ignored if reading a file causes IOError.

            Default is :obj:`True`.

        Attributes
        ----------
        <various> : :class:`pandas.DataFrame`
            Dataframes holding the respective data, available as attributes
            described in the ``only_these_keys`` parameter.

        Raises
        ----------
        :obj:`IOError`
            If not permissive and a file is not found.
        """
        from ngs_toolkit.utils import get_this_file_or_timestamped

        prefix = self._format_string_with_attributes(prefix)

        if output_map is None:
            kwargs = {"index_col": 0}
            output_map = {
                "matrix_raw": (prefix + ".matrix_raw.csv", kwargs),
                "matrix_norm": (prefix + ".matrix_norm.csv", kwargs),
                "matrix_features": (prefix + ".matrix_features.csv", kwargs),
                "differential_results": (
                    os.path.join(
                        self.results_dir,
                        "differential_analysis_{}".format(self.data_type),
                        "differential_analysis.deseq_result.all_comparisons.csv",
                    ),
                    kwargs,
                ),
            }

        if only_these_keys is None:
            only_these_keys = output_map.keys()

        output_map = {k: v for k, v in output_map.items() if k in only_these_keys}

        for name, (file, kwargs) in output_map.items():
            _LOGGER.info("Loading '{}' analysis attribute.".format(name))
            try:
                setattr(self, name, pd.read_csv(get_this_file_or_timestamped(file), **kwargs))

                # Fix possible multiindex for matrix_norm
                if name == "matrix_norm":
                    if not getattr(self, name).dtypes.all().name == "float64":
                        msg = (
                            "`matrix_norm` value has non-float values."
                            " This is likely because it is a pandas.DataFrame"
                            " with MultiIndex columns.\n"
                            " Since the length and type of MultiIndex values"
                            " cannot be safely inferred from the CSV format"
                            " please determine the number of rows part of the "
                            " MultiIndex and read in the file explicitely:\n"
                            "``analysis.matrix_norm = pandas.read_csv('{file}'"
                            ", index_col=0, header=list(range(6)))``"
                        ).format(file=output_map["matrix_norm"][0])
                        _LOGGER.error(msg)
            except IOError as e:
                if not permissive:
                    raise e
                _LOGGER.warning(e)

        # if "differential_enrichment" in output_map:
        #     self.enrichment_results = dict()
        #     self.enrichment_results["enrichr"] = pd.read_csv(
        #         os.path.join(
        #             "results",
        #             "differential_analysis_RNA-seq",
        #             "enrichments",
        #             "differential_analysis.enrichr.csv",
        #         )
        #     )

    def record_output_file(
        self,
        file_name,
        name="analysis",
        dump_yaml=True,
        output_yaml="{root_dir}/{name}.analysis_record.yaml",
    ):
        """
        Record an analysis output.

        Will also write all records to a YAML file and call `Analysis.generate_report`
        if specified in general configuration.

        Parameters
        ----------
        file_name : :obj:`str`
            Filename of output to report.
        name : :obj:`str`, optional
            Name of the output to report.

            Defaults to "analysis".
        dump_yaml : :obj:`bool`, optional
            Whether to dump records to yaml file.

            Defaults to :obj:`True`.
        output_yaml : :obj:`str`, optional
            YAML file to dump records to.
            Will be formated with Analysis variables.

            Defaults to "{root_dir}/{name}.analysis_record.yaml".

        Attributes
        ----------
        output_files : :obj:`list`
            Appends a tuple of (``name``, ``file_name``) to ``output_files``.
        """
        import yaml

        self.output_files.append((name, file_name))
        if dump_yaml:
            yaml.safe_dump(
                self.output_files, open(self._format_string_with_attributes(output_yaml), "w")
            )

        if _CONFIG["preferences"]["report"]["continuous_generation"]:
            self.generate_report(pip_versions=False)

    def generate_report(
        self, output_html="{root_dir}/{name}.analysis_report.html", template=None, pip_versions=True
    ):
        """
        Record an analysis output.

        Parameters
        ----------
        output_html : :obj:`str`
            Filename of output to report.

            Defaults to "{root_dir}/{name}.analysis_report.html".
        template : :obj:`None`, optional
            Name of the output to report.

            Default is the HTML template distributed with ngs-toolkit.
        pip_versions: :obj:`bool`, optional
            Whether the versions of Python packages should be included
            in the report by using pip freeze.

            Default is :obj:`True`.
        """
        import time
        import sys
        from collections import OrderedDict

        from jinja2 import Template
        import pkg_resources
        from ngs_toolkit import __version__

        try:
            from pip._internal.operations import freeze
        except ImportError:  # pip < 10.0
            from pip.operations import freeze

        def fix_name(x, name):
            return (
                " - ".join(os.path.basename(x).replace(name, "").split(".")[:-1])
                .replace("_", " ")
                .capitalize()
            )

        output_html = self._format_string_with_attributes(output_html)

        # Lets reorganize the output_files
        # into a dict of {name: list(file_names)}
        keys = OrderedDict()
        for x in self.output_files:
            keys[x[0]] = x[1]
        outputs = {k: list() for k in keys}

        # get relative paths:
        # this enables viewing linked files in html independent of the machine
        for key, file in self.output_files:
            outputs[key].append(os.path.relpath(file, self.root_dir))

        # Select image outputs (non-CSV files)
        images = {k: [x for x in v if x.endswith(".svg")] for k, v in outputs.items()}
        # Generate dict = {"section": [(caption, file), ...]}
        images = {
            k.capitalize().replace("_", " "): [(fix_name(x, self.name), x,) for x in v]
            for k, v in images.items()
        }
        csvs = {k: [x for x in v if x.endswith(".csv")] for k, v in outputs.items()}
        csvs = {
            k.capitalize().replace("_", " "): [(fix_name(x, self.name), x,) for x in v]
            for k, v in csvs.items()
        }

        # Get template
        if template is None:
            resource_package = "ngs_toolkit"
            resource_path = "/".join(("templates", "report.html"))
            template = Template(
                pkg_resources.resource_string(resource_package, resource_path).decode()
            )
        else:
            template = Template(open(template, "r").read())

        # Format
        output = template.render(
            analysis=self,
            project_repr={k: v for k, v in self.__dict__.items() if isinstance(v, str)},
            samples=[s.to_dict() for s in self.samples] if self.samples is not None else [],
            time=time.asctime(),
            images=images,
            csvs=csvs,
            python_version=sys.version,
            library_version=__version__,
            freeze=[] if not pip_versions else list(freeze.freeze()),
        )

        # Write
        with open(output_html, "w") as handle:
            handle.write(output)

    def set_matrix(self, matrix_name, csv_file, prefix="{results_dir}/{name}", **kwargs):
        """
        Set an existing CSV file as the value of the analysis' matrix.

        Parameters
        ----------
        matrix_name : :obj:`str`
            The attribute name of the matrix.

            Options are "matrix_raw" and "matrix_norm".
        csv_file : :obj:`str`
            Path to valid CSV file to be used as matrix.
            Assumes header and index column.
            Customize additional overwriding options to read CSV by passing kwargs.
        prefix : :obj:`str`, optional
            String prefix of paths to save files.
            Variables in curly braces will be formated with attributes of analysis.

            Defaults to "{results_dir}/{name}".
        **kwargs : :obj:`dict`
            Additional keyword arguments will be passed to :obj:`pandas.read_csv`.

        Attributes
        ----------
        matrix_name : :obj:`pandas.DataFrame`
            An attribute named `matrix_name` holding the respecive matrix.
        """
        from ngs_toolkit.utils import fix_dataframe_header

        prefix = self._format_string_with_attributes(prefix)
        output_file = prefix + "." + matrix_name + ".csv"

        options = {"index_col": 0}
        options.update(kwargs)
        df = pd.read_csv(csv_file, **options)
        if isinstance(df.columns, pd.MultiIndex):
            df = fix_dataframe_header(df)
        _LOGGER.info("Set {} as '{}' attribute.".format(csv_file, matrix_name))
        setattr(self, matrix_name, df)
        _LOGGER.info("Saving '{}' attribute to {}.".format(matrix_name, output_file))
        df.to_csv(output_file, index=True)

    @check_has_attributes(["organism", "genome"])
    def get_resources(
        self,
        steps=["blacklist", "tss", "genomic_context"],
        organism=None,
        genome_assembly=None,
        output_dir=None,
        overwrite=False,
    ):
        """
        Get genome-centric resources used by several ``ngs_toolkit`` analysis
        functions.

        Parameters
        ----------
        steps : :obj:`list`, optional
            What kind of annotations to get. Options are:

                 * "genome": Genome sequence (2bit format)
                 * "blacklist": Locations of blacklisted regions for genome
                 * "tss": Locations of gene"s TSSs
                 * "genomic_context": Genomic context of genome
                 * "chromosome_sizes": Sizes of chromosomes

            Defaults to ["blacklist", "tss", "genomic_context"].
        organism : :obj:`str`, optional
            Organism to get for. Currently supported are "human" and "mouse".

            Defaults to analysis' own organism.
        genome_assembly : :obj:`str`, optional
            Genome assembly to get resources for.
            Currently supported are "hg19", "hg38" and "mm10".

            Defaults to the genome assembly of the analysis.
        output_dir : :obj:`str`, optional
            Directory to save results to.

            Defaults to the value of ``preferences:root_reference_dir`` in the
            configuration, if that is not set, to a directory called "reference"
            in the analysis root directory.
        overwrite: :obj:`bool`, optional
            Whether existing files should be overwritten by new ones.
            Otherwise they will be kept and no action is made.

            Defaults to :obj:`False`.

        Returns
        -------
        : dict
            Dictionary with keys same as the options as steps, containing paths
            to the requested files.

            The values of the 'genome' step are also a dictionary with keys
            "2bit" and "fasta" for each file type respectively.
        """
        from ngs_toolkit.general import (
            get_genome_reference,
            get_blacklist_annotations,
            get_tss_annotations,
            get_genomic_context,
            get_chromosome_sizes,
        )
        from ngs_toolkit.utils import get_this_file_or_timestamped

        if organism is None:
            organism = self.organism
        if genome_assembly is None:
            genome_assembly = self.genome
        if output_dir is None:
            output_dir = self._format_string_with_environment_variables(
                _CONFIG["preferences"]["root_reference_dir"]
            )
            if output_dir is None:
                output_dir = os.path.join(self.root_dir, "reference")

        kwargs = {
            "organism": organism,
            "genome_assembly": genome_assembly,
            "output_dir": output_dir,
            "overwrite": overwrite,
        }
        mapping = {"hg19": "grch37", "hg38": "grch38", "mm10": "grcm38"}
        output = dict()

        if "genome" in steps:
            output["genome_file"] = dict()
            output["genome_file"]["2bit"] = get_genome_reference(file_format="2bit", **kwargs)
            fasta = output["genome_file"]["2bit"].replace(".2bit", ".fa")
            if os.path.exists(fasta):
                output["genome_file"]["fasta"] = fasta
            else:
                output["genome_file"]["fasta"] = get_genome_reference(file_format="fasta", **kwargs)
        if "blacklist" in steps:
            output["blacklist_file"] = get_blacklist_annotations(**kwargs)
        if "tss" in steps:
            get_tss_annotations(**kwargs)
            output["tss_file"] = get_this_file_or_timestamped(
                os.path.join(
                    output_dir,
                    "{}.{}.gene_annotation.protein_coding.tss.bed".format(
                        self.organism, mapping[self.genome]
                    ),
                )
            )
        if "genomic_context" in steps:
            get_genomic_context(**kwargs)
            output["genomic_context_file"] = os.path.join(
                output_dir, "{}.{}.genomic_context.bed".format(self.organism, mapping[self.genome]),
            )

        if "chromosome_sizes" in steps:
            output["chromosome_sizes_file"] = get_chromosome_sizes(**kwargs)
        return output

    def normalize_rpm(
        self,
        matrix="matrix_raw",
        samples=None,
        mult_factor=1e6,
        log_transform=True,
        pseudocount=1,
        save=True,
        assign=True,
    ):
        """
        Normalization of matrix of (n_features, n_samples) by total in each sample.

        Parameters
        ----------
        matrix : :obj:`str`, optional
            Attribute name of matrix to normalize.

            Defaults to "matrix_raw".
        samples : :obj:`list`, optional
            Iterable of :class:`peppy.Sample` objects to restrict matrix to.

            Defaults to all samples in matrix.
        mult_factor : :obj:`float`, optional
            A constant to multiply values for.

            Defaults to 1e6.
        log_transform: :obj:`bool`, optional
            Whether to log transform values or not.

            Defaults to :obj:`True`.
        pseudocount : {:obj:`int`, :obj:`float`}, optional
            A constant to add to values.

            Defaults to 1.
        save : :obj:`bool`, optional
            Whether to write normalized DataFrame to disk.

            Defaults to :obj:`True`.
        assign : :obj:`bool`, optional
            Whether to assign the normalized DataFrame to ``matrix_norm``.

            Defaults to :obj:`True`.

        Attributes
        ----------
        matrix_norm : :class:`pandas.DataFrame`
            If ``assign``, a :class:`pandas.DataFrame` normalized with respective method.
        norm_method : :obj:`str`
            If ``assign``, it is the name of method used to normalize: "rpm".

        Returns
        -------
        :class:`pandas.DataFrame`
            Normalized dataframe.
        """
        to_norm = self.get_matrix(matrix=matrix, samples=samples)
        # apply normalization over total
        matrix_norm = (to_norm / to_norm.sum()) * mult_factor

        # Make non-negative
        if matrix_norm.min().min() <= 0:
            matrix_norm += np.absolute(matrix_norm.min().min())

        # Log2 transform
        if log_transform:
            matrix_norm = np.log2(pseudocount + matrix_norm)

        if save:
            matrix_norm.to_csv(
                os.path.join(self.results_dir, self.name + ".matrix_norm.csv"), index=True,
            )
        if assign:
            self.matrix_norm = matrix_norm
            self.norm_method = "rpm"

        return matrix_norm

    def normalize_quantiles(
        self,
        matrix="matrix_raw",
        samples=None,
        implementation="Python",
        log_transform=True,
        pseudocount=1,
        save=True,
        assign=True,
    ):
        """
        Quantile normalization of matrix of (n_features, n_samples).

        Parameters
        ----------
        matrix : :obj:`str`
            Attribute name of matrix to normalize.

            Defaults to "matrix_raw".
        samples : :obj:`list`, optional
            Iterable of :class:`peppy.Sample` objects to restrict matrix to.

            Defaults to all in matrix.
        implementation : :obj:`str`, optional
            One of ``Python`` or ``R``.
            Dictates which implementation is to be used.
            The R implementation comes from the `preprocessCore` package,
            and the Python one is from https://github.com/ShawnLYU/Quantile_Normalize.
            They give very similar results.

            Default is "Python".
        log_transform: :obj:`bool`, optional
            Whether to log transform values or not.

            Default is :obj:`True`.
        pseudocount : float, optional
            A constant to add before log transformation.

            Default is 1.
        save : :obj:`bool`, optional
            Whether to write normalized DataFrame to disk.

            Default is :obj:`True`.
        assign : :obj:`bool`, optional
            Whether to assign the normalized DataFrame to an attribute `matrix_norm`.

            Default is :obj:`True`.

        Attributes
        ----------
        matrix_norm : :class:`pandas.DataFrame`
            If ``assign``, a pandas DataFrame normalized with respective method.
        norm_method : :obj:`str`
            If ``assign``, it is the name of method used to normalize: "quantile".

        Returns
        -------
        :class:`pandas.DataFrame`
            Normalized dataframe.
        """
        from ngs_toolkit.utils import normalize_quantiles_r, normalize_quantiles_p

        to_norm = self.get_matrix(matrix=matrix, samples=samples)

        if implementation == "R":
            matrix_norm = pd.DataFrame(
                normalize_quantiles_r(to_norm.values), index=to_norm.index, columns=to_norm.columns,
            )
        elif implementation == "Python":
            matrix_norm = normalize_quantiles_p(to_norm)
        else:
            msg = "Implementation of quantile normalization must be one of 'R' of 'Python'"
            _LOGGER.error(msg)
            raise ValueError(msg)

        # Make non-negative
        if matrix_norm.min().min() <= 0:
            matrix_norm += np.absolute(matrix_norm.min().min())

        # Log2 transform
        if log_transform:
            matrix_norm = np.log2(pseudocount + matrix_norm)

        if save:
            matrix_norm.to_csv(
                os.path.join(self.results_dir, self.name + ".matrix_norm.csv"), index=True,
            )
        if assign:
            self.matrix_norm = matrix_norm
            self.norm_method = "quantile"

        return matrix_norm

    def normalize_median(
        self,
        matrix="matrix_raw",
        samples=None,
        function=np.nanmedian,
        fillna=True,
        save=True,
        assign=True,
    ):
        """
        Normalization of matrices of (n_features, n_samples)
        by subtracting the median from each sample/feature.
        Most appopriate for CNV data.

        Parameters
        ----------
        matrix : :obj:`str`, optional
            Attribute name of dictionary of matrices to normalize.

            Defaults to "matrix_raw".
        samples : :obj:`list`, optional
            Samples to restrict analysis to.

            Defaults to all samples in Analysis object.
        function : function, optional
            An alternative function to calculate across samples. Data will be subtracted by this.

            Defaults to ``numpy.nanmedian``.
        fillna: :obj:`bool`, optional
            Whether to fill NaN with zero.

            Defaults to :obj:`True`.
        save : :obj:`bool`, optional
            Whether results should be saved to disc.

            Defaults to :obj:`True`.
        assign : :obj:`bool`, optional
            Whether to assign the normalized DataFrame to an attribute `matrix_norm`.

            Default is :obj:`True`.

        Attributes
        ----------
        matrix_norm : :class:`pandas.DataFrame`
            If ``assign``, a pandas DataFrame normalized with respective method.
        norm_method : :obj:`str`
            If ``assign``, it is the name of method used to normalize: "median".

        Returns
        -------
        :class:`pandas.DataFrame`
            Normalized dataframe.
        """
        matrix = self.get_matrix(matrix, samples=samples)

        matrix_norm = dict()

        to_norm = self.get_matrix(matrix=matrix, samples=samples)
        matrix_norm = (to_norm.T - function(to_norm, axis=1)).T
        if fillna:
            matrix_norm = matrix_norm.fillna(0)
        if save:
            matrix_norm.to_csv(
                os.path.join(self.results_dir, self.name + ".matrix_norm.csv"), index=True,
            )
        if assign:
            self.matrix_norm = matrix_norm
            self.norm_method = "median"

        return matrix_norm

    def normalize_pca(
        self, pc, matrix="matrix_raw", samples=None, save=True, assign=True, **kwargs
    ):
        """
        Normalization of a matrix by subtracting the
        contribution of Principal Component `pc` from each sample/feature.

        Parameters
        ----------
        pc : :obj:`int`
            Principal Component to remove. 1-based.

        matrix : :obj:`str`, optional
            Attribute name of dictionary of matrices to normalize.

            Defaults to "matrix_raw".
        samples : :obj:`list`, optional
            Samples to restrict analysis to.

            Defaults to all samples.
        save : :obj:`bool`, optional
            Whether results should be saved to disc.

            Defaults to :obj:`True`.
        assign : :obj:`bool`, optional
            Whether to assign the normalized DataFrame to an attribute `matrix_norm`.

            Default is :obj:`True`.

        **kwargs : :obj:`dict`, optional
            Additional keyword arguments will be passed to
            :class:`ngs_toolkit.general.subtract_principal_component`.

        Attributes
        ----------
        matrix_norm : :class:`pandas.DataFrame`
            If ``assign``, a pandas DataFrame normalized with respective method.
        norm_method : :obj:`str`
            If ``assign``, it is the name of method used to normalize: "pca".

        Returns
        -------
        :class:`pandas.DataFrame`
            Normalized dataframe.
        """
        from ngs_toolkit.general import subtract_principal_component

        if pc is None:
            raise ValueError("Principal Component to remove must be specified!")

        matrix = self.get_matrix(matrix, samples=samples)

        # first make sure data is centered
        to_norm = self.normalize_median(matrix, samples=samples, save=False, assign=False)
        # then remove the PC
        default_kwargs = {
            "plot_name": os.path.join(self.results_dir, "PCA_based_batch_correction.svg")
        }
        default_kwargs.update(kwargs)
        matrix_norm = subtract_principal_component(to_norm.T.fillna(0), pc=pc, **default_kwargs).T

        if save:
            matrix_norm.to_csv(
                os.path.join(self.results_dir, self.name + ".matrix_norm.csv"), index=True,
            )
        if assign:
            self.matrix_norm = matrix_norm
            self.norm_method = "pca"

        return matrix_norm

    def normalize_vst(self, matrix="matrix_raw", samples=None, save=True, assign=True, **kwargs):
        """
        Normalization of a matrix using
        Variance Stabilization Transformation (VST) method from DESeq2.

        Parameters
        ----------
        matrix : :obj:`str`, optional
            Attribute name of dictionary of matrices to normalize.

            Defaults to "matrix_raw".
        samples : :obj:`list`, optional
            Samples to restrict analysis to.

            Defaults to all samples.
        save : :obj:`bool`, optional
            Whether results should be saved to disc.

            Defaults to :obj:`True`.
        assign : :obj:`bool`, optional
            Whether to assign the normalized DataFrame to an attribute `matrix_norm`.

            Default is :obj:`True`.

        **kwargs : :obj:`dict`
            Additional keywork arguments will be passed
            to `DESeq2::varianceStabilizingTransformation`.

        Attributes
        ----------
        matrix_norm : :class:`pandas.DataFrame`
            If ``assign``, a DataFrame normalized with VST method.
        norm_method : :obj:`str`
            If ``assign``, it is the name of method used to normalize: "vst".

        Returns
        -------
        :class:`pandas.DataFrame`
            Normalized dataframe.
        """
        from rpy2.robjects import numpy2ri, pandas2ri, r
        from rpy2.robjects.packages import importr

        numpy2ri.activate()
        pandas2ri.activate()

        importr("DESeq2")

        matrix = self.get_matrix(matrix, samples=samples)

        # Apply VST
        matrix_norm = pd.DataFrame(
            r.varianceStabilizingTransformation(matrix.values, **kwargs),
            index=matrix.index,
            columns=matrix.columns,
        )

        if save:
            matrix_norm.to_csv(
                os.path.join(self.results_dir, self.name + ".matrix_norm.csv"), index=True,
            )
        if assign:
            self.matrix_norm = matrix_norm
            self.norm_method = "vst"

        return matrix_norm

    def normalize(
        self, method="quantile", matrix="matrix_raw", samples=None, save=True, assign=True, **kwargs
    ):
        """
        Normalization of matrix of (n_features, n_samples).

        Parameters
        ----------
        method : :obj:`str`, optional
            Normalization method to apply. One of:
             - ``rpm``: Reads per million normalization (RPM).
             - ``vst``: Variance stabilization transformation (uses ``DESeq2`` R package).
             - ``quantile``: Quantile normalization and log2 transformation.
             - ``cqn``: Conditional quantile normalization (uses ``cqn`` R package).
                      Only available for ATAC-seq.
             - ``median``: Substraction of median per feature.
                      Only useful for CNV.
             - ``pca``: Subtraction of Principal Component from matrix.
                      Requires which PC to subtract. ``pc`` must be passed as kwarg.

            Defaults to "quantile".
        matrix : :obj:`str`, optional
            Attribute name of matrix to normalize.

            Defaults to "matrix_raw".
        samples : :obj:`list`, optional
            Iterable of :class:`peppy.Sample` objects to restrict matrix to.

            Default is all samples in matrix.
        save : :obj:`bool`, optional
            Whether to write normalized DataFrame to disk.

            Defaults to :obj:`True`.
        assign : :obj:`bool`, optional
            Whether to assign the normalized DataFrame to an attribute `matrix_norm`.

            Default is :obj:`True`.
        **kwargs : :obj:`dict`
            Additional keyword arguments will be passed to the respective
            normalization function.

        Attributes
        ----------
        matrix_norm : :class:`pandas.DataFrame`
            If ``assign``, a pandas DataFrame normalized with respective method.
        norm_method : :obj:`str`
            If ``assign``, it is the ``method`` used to normalize.

        Returns
        -------
        :class:`pandas.DataFrame`
            Normalized dataframe.
        """
        if method == "rpm":
            return self.normalize_rpm(matrix=matrix, samples=samples, save=save, assign=assign)
        elif method == "quantile":
            return self.normalize_quantiles(
                matrix=matrix, samples=samples, save=save, assign=assign
            )
        elif method == "cqn":
            if self.data_type == "RNA-seq":
                raise ValueError(
                    "Cannot use `cqn` normalization with this data_type: {}".format(self.data_type)
                )
            return self.normalize_cqn(matrix=matrix, samples=samples, save=save, assign=assign)
        elif method == "median":
            return self.normalize_median(matrix=matrix, samples=samples, save=save, assign=assign)
        elif method == "pca":
            if "pc" not in kwargs:
                raise ValueError("`pca` normalization requires `pc` as kwarg")
            return self.normalize_pca(
                matrix=matrix, samples=samples, save=save, assign=assign, pc=kwargs["pc"],
            )
        elif method == "vst":
            return self.normalize_vst(
                matrix=matrix, samples=samples, save=save, assign=assign, **kwargs
            )
        else:
            msg = "Requested normalization method is not available!"
            _LOGGER.error(msg)
            raise ValueError(msg)

    def remove_factor_from_matrix(
        self,
        factor,
        method="combat",
        covariates=None,
        matrix="matrix_norm",
        samples=None,
        save=True,
        assign=True,
        make_positive=True,
    ):
        """
        Remove an annotated factor from a matrix using Combat.

        Requires the Python port of "Combat" to be installed. Install for example the following fork:

        .. highlight:: shell
        .. code-block:: shell

            pip install git+https://github.com/afrendeiro/combat.git

        Parameters
        ----------
        factor : :obj:`str`
            The name of the factor to remove from matrix.

        method : :obj:`str`
            The method to use to remove the factor.

            Default is "combat".
        covariates : :obj:`list`
            Covariates to consider when removing factor.
            These will be kept in the data.

        matrix : {str, pandas.DataFrame}
            The name of the attribute with the matrix or a DataFrame.

            Defaults to "matrix_norm".
        samples : :obj:`list`
            Iterable of :class:`peppy.Sample` objects to restrict matrix to.

            Default is not to subset matrix.
        save : :obj:`bool`, optional
            Whether to write normalized DataFrame to disk.

            Defaults to :obj:`True`.
        assign : :obj:`bool`
            Whether to assign the result to "matrix_norm".

            Defaults to :obj:`True`.
        make_positive  : :obj:`bool`
            Whether to make resulting matrix non-negative.
            Not implemented yet.

            Defaults to :obj:`True`.

        Returns
        -------
        :class:`pandas.DataFrame`
            Requested matrix (dataframe).
        """
        from combat import combat
        from patsy import dmatrix

        if method != "combat":
            msg = "Only implemented method is 'combat'."
            _LOGGER.error(msg)
            raise NotImplementedError(msg)

        matrix = self.get_matrix(matrix, samples=samples)

        # Include only variables with variance
        std = matrix.std(axis=1)
        if any(std == 0):
            msg = "Matrix contains features with zero variance. Removing those."
            _LOGGER.warning(msg)
            matrix = matrix.loc[std > 0, :]

        # make vector of factor to remove
        if samples is None:
            samples = [s for s in self.samples if s.name in matrix.columns]
        batch = pd.Series([getattr(s, factor) for s in samples], index=[s.name for s in samples])

        # make design model of covariates
        if covariates is not None:
            _LOGGER.debug("Generating design matrix for covariates.")
            if not isinstance(matrix.columns, pd.MultiIndex):
                _LOGGER.debug(
                    "Matrix was not MultiIndex, annotating matrix with sample attributes."
                )
                matrix = self.annotate_samples(matrix=matrix, save=False, assign=False)
            d = matrix.columns.to_frame().set_index("sample_name")
            covariates = dmatrix("~ " + " + ".join(covariates), d)

        matrix_norm = combat(matrix, batch=batch, model=covariates)

        if save:
            matrix_norm.to_csv(
                os.path.join(self.results_dir, self.name + ".matrix_norm.csv"), index=True,
            )
        if assign:
            self.matrix_norm = matrix_norm
            self.norm_method = "{} + {}".format(self.norm_method, method).replace("None + ", "")

        return matrix_norm

    def get_matrix(self, matrix, samples=None):
        """
        Get a matrix that is an attribute of self subsetted for the requested samples.

        Parameters
        ----------
        matrix : {str, pandas.DataFrame}
            The name of the attribute with the matrix or a DataFrame already.
        samples : :obj:`list`
            Iterable of :class:`peppy.Sample` objects to restrict matrix to.

            Default is not to subset matrix.

        Returns
        -------
        :class:`pandas.DataFrame`
            Requested matrix (dataframe).
        """
        if isinstance(matrix, str):
            matrix = getattr(self, matrix)
        if samples is None:
            # all samples in matrix
            return matrix
        else:
            # subset to requested samples
            return matrix.loc[:, [s.name for s in samples]]

    def get_matrix_stats(
        self,
        matrix="matrix_raw",
        samples=None,
        save=True,
        output_prefix="stats_per_feature",
        assign=True,
    ):
        """
        Gets a matrix of feature-wise (ie for every gene or region) statistics such
        across samples such as mean, variance, deviation, dispersion and amplitude.

        Parameters
        ----------
        matrix : :obj:`str`
            Attribute name of matrix to normalize.

            Defaults to "matrix_raw".
        samples : :obj:`list` [:class:`peppy.Sample`]
            Subset of samples to use.

            Defaults to all in analysis.
        save : :obj:`bool`, optional
            Whether to write the annotated DataFrame to disk.

            Default is :obj:`True`.
        output_prefix : :obj:`str`, optional
            Prefix to add to output file when save is True.

            Default is "matrix_features".
        assign : :obj:`bool`, optional
            Whether to assign the annoatated DataFrame to "matrix_features".

            Default is :obj:`True`.

        Returns
        -------
        :class:`pandas.DataFrame`
            Statistics for each feature.

        Attributes
        ----------
        stats : :class:`pandas.DataFrame`
            A DataFrame with statistics for each feature.
        """
        matrix = self.get_matrix(matrix=matrix, samples=samples)

        metrics = pd.DataFrame(index=pd.Index(matrix.index, name=self.var_unit_name))
        # calculate mean coverage
        metrics.loc[:, "mean"] = matrix.mean(axis=1)
        # calculate coverage variance
        metrics.loc[:, "variance"] = matrix.var(axis=1)
        # calculate std deviation (sqrt(variance))
        metrics.loc[:, "std_deviation"] = np.sqrt(metrics.loc[:, "variance"])
        # calculate dispersion (variance / mean)
        metrics.loc[:, "dispersion"] = metrics.loc[:, "variance"] / metrics.loc[:, "mean"]
        # calculate qv2 (std / mean) ** 2
        metrics.loc[:, "qv2"] = (metrics.loc[:, "std_deviation"] / metrics.loc[:, "mean"]) ** 2
        # calculate "amplitude" (max - min)
        metrics.loc[:, "amplitude"] = metrics.max(axis=1) - metrics.min(axis=1)
        # calculate interquantile range
        metrics.loc[:, "iqr"] = metrics.quantile(0.75, axis=1) - metrics.quantile(0.25, axis=1)
        metrics.index.name = "index"
        if save:
            metrics.to_csv(
                os.path.join(self.results_dir, self.name + ".{}.csv".format(output_prefix)),
                index=True,
            )

        if assign:
            self.stats = metrics
        return metrics

    def annotate_features(
        self,
        samples=None,
        matrix="matrix_norm",
        feature_tables=None,
        permissive=True,
        save=True,
        assign=True,
        output_prefix="matrix_features",
    ):
        """
        Annotates analysis features (regions/genes) by aggregating annotations
        per feature (genomic context, chromatin state, gene annotations and statistics)
        if present and relevant depending on the data type of the Analysis.

        The numeric matrix to be used is specified in `matrix`.
        If any two two annotation dataframes have equally named columns (e.g. chrom, start, end),
        the value of the first is kept.

        Parameters
        ----------
        samples : :obj:`list`
            Iterable of :class:`peppy.Sample` objects to restrict matrix to.
            Calculated metrics will be restricted to these samples.

            Defaults to all in analysis (the matrix will not be subsetted).
        matrix : :obj:`str`
            Attribute name of matrix to annotate.

            Defaults to "matrix_norm".
        feature_tables : :obj:`list`
            Attribute names of dataframes used to annotate the numeric dataframe.

            Default is ["gene_annotation", "region_annotation", "chrom_state_annotation", "support", "stats"]
            for ATAC-seq and ChIP-seq and ["stats"] for all others.

        permissive : :obj:`bool`
            Whether DataFrames that do not exist should be simply skipped or an error will be thrown.

            Defaults to :obj:`True`.
        save : :obj:`bool`, optional
            Whether to write the annotated DataFrame to disk.

            Default is :obj:`True`.
        assign : :obj:`bool`, optional
            Whether to assign the annoatated DataFrame to "matrix_features".

            Default is :obj:`True`.
        output_prefix : :obj:`str`, optional
            Prefix to add to output file when ``save`` is :obj:`True`.

            Default is "matrix_features".

        Raises
        ----------
        AttributeError
            If not `permissive` a required DataFrame does not exist as an object attribute.

        Attributes
        ----------
        matrix_features : :class:`pandas.DataFrame`
            A pandas DataFrame containing annotations of the region features.
        """
        if samples is None:
            samples = self.samples
        matrix = getattr(self, matrix)

        next_matrix = matrix
        # add closest gene
        msg = "`{}` attribute does not exist."

        if feature_tables is None:
            if self.data_type in ["ATAC-seq", "ChIP-seq"]:
                feature_tables = [
                    "gene_annotation",
                    "region_annotation",
                    "chrom_state_annotation",
                    "support",
                    "stats",
                ]
            else:
                feature_tables = ["stats"]

        for matrix_name in feature_tables:
            if hasattr(self, matrix_name):
                cur_matrix = getattr(self, matrix_name)
                matrix_features = pd.merge(
                    next_matrix,
                    cur_matrix[cur_matrix.columns.difference(next_matrix.columns)],
                    left_index=True,
                    right_index=True,
                    how="left",
                )
                next_matrix = matrix_features
            else:
                if not permissive:
                    _LOGGER.error(msg.format(matrix_name))
                    raise AttributeError(msg.format(matrix_name))
                else:
                    _LOGGER.warning(msg.format(matrix_name) + " Proceeding anyway.")

        if "matrix_features" not in locals():
            matrix_features = next_matrix

        # Pair indexes
        msg = "Annotated matrix does not have same feature length as matrix_raw matrix."
        if not matrix.shape[0] == matrix_features.shape[0]:
            _LOGGER.error(msg)
            raise AssertionError(msg)
        matrix_features.index = matrix.index
        matrix_features.index.name = "index"

        # Save
        if save:
            matrix_features.to_csv(
                os.path.join(self.results_dir, self.name + ".{}.csv".format(output_prefix)),
                index=True,
            )
        if assign:
            self.matrix_features = matrix_features

        return matrix_features

    def annotate_samples(
        self,
        matrix="matrix_norm",
        attributes=None,
        numerical_attributes=None,
        save=False,
        assign=False,
    ):
        """
        Annotate matrix ``(n_features, n_samples)`` with sample metadata
        (creates MultiIndex on columns). Numerical attributes can be pass as a iterable
        to ``numerical_attributes`` to be converted.

        Parameters
        ----------
        matrix : :obj:`str`, optional
            Attribute name of matrix to annotate.

            Defaults to "matrix_norm".
        attributes : :obj:`list`, optional
            Desired attributes to be annotated.

            Defaults to all attributes in the original sample annotation sheet of the analysis' Project.
        numerical_attributes : :obj:`list`, optional
            Attributes which are numeric even though they
            might be so in the samples" attributes. Will attempt
            to convert values to numeric.
        save : :obj:`bool`, optional
            Whether to write normalized DataFrame to disk.

            Default is :obj:`True`.
        assign : :obj:`bool`, optional
            Whether to assign the normalized DataFrame to "matrix_norm".

            Default is :obj:`True`.

        Returns
        -------
        :class:`pandas.DataFrame`
            Annotated dataframe with requested sample attributes.

        Attributes
        ----------
        matrix_norm : :obj:`pandas.DataFrame`
            A pandas DataFrame with  MultiIndex column index containing the sample's attributes specified.
        """
        if (attributes is None) and hasattr(self, "sample_attributes"):
            _LOGGER.info(
                "Using 'sample_attributes' from analysis to annotate matrix: {}".format(
                    ",".join(self.sample_attributes)
                )
            )
            attributes = self.sample_attributes
        if (attributes is None) and hasattr(self, "prj"):
            _LOGGER.warning(
                "Analysis has no 'sample_attributes' set. "
                + "Will use all columns from project annotation sheet: {}".format(
                    ",".join(self.prj.sheet.columns)
                )
            )
            attributes = self.prj.sheet.columns
        if attributes is None:
            msg = "Attributes not given and could not be set from Analysis or Project."
            raise ValueError(msg)

        if self.data_type == "CNV":
            if matrix is None:
                if not isinstance(getattr(self, matrix), pd.DataFrame):
                    _LOGGER.error(
                        "For CNV data type, the matrix to be annotated must be"
                        + " directly passed to the function throught the `matrix` argument!"
                    )
                    raise ValueError

        matrix = self.get_matrix(matrix)

        if isinstance(matrix.columns, pd.core.indexes.multi.MultiIndex):
            matrix.columns = matrix.columns.get_level_values("sample_name")

        samples = [s for s in self.samples if s.name in matrix.columns.tolist()]

        attrs = list()
        for attr in attributes:
            _LOGGER.debug("Attribute: '{}'".format(attr))
            ll = list()
            for sample in samples:  # keep order of samples in matrix
                try:
                    ll.append(getattr(sample, attr))
                except AttributeError:
                    ll.append(np.nan)
            if numerical_attributes is not None:
                if attr in numerical_attributes:
                    ll = pd.Series(ll).replace("", np.nan).astype(float).tolist()
            attrs.append(ll)

        # Generate multiindex columns
        index = pd.MultiIndex.from_arrays(attrs, names=attributes)
        df = matrix[[s.name for s in samples]]
        df.columns = index

        # Save
        if save:
            df.to_csv(
                os.path.join(self.results_dir, self.name + ".matrix_norm.csv"), index=True,
            )
        if assign:
            self.matrix_norm = df
        return df

    def annotate_matrix(self, **kwargs):
        """
        Convinience function to create dataframes annotated with feature and samples attributes.

        Simply calls :func:`Analysis.annotate_features` and :func:`Analysis.annotate_samples`.

        Parameters
        ----------
        kwargs : :obj:`dict`
            Additional keyword arguments are passed to the above mentioned functions.
        """
        self.annotate_features(**kwargs)
        self.annotate_samples(**kwargs)

    def get_level_colors(
        self,
        index=None,
        matrix="matrix_norm",
        levels=None,
        pallete="tab20",
        uniform_cmap="plasma",
        diverging_cmap="RdYlBu_r",
        nan_color=(0.662745, 0.662745, 0.662745, 1.0),
        as_dataframe=False,
    ):
        """
        Get tuples of floats representing a colour for a sample in a given variable in a
        dataframe"s index (particularly useful with MultiIndex dataframes).

        If given, will use the provieded ``index`` argument, otherwise, the columns
        and its levels of an attribute of self named ``matrix``.
        ``levels`` can be passed to subset the levels of the index.

        Will try to guess if each variable is categorical or numerical and return either colours
        from a colour ``pallete`` or a ``cmap``, respectively with null values set to ``nan_color``
        (a 4-value tuple of floats).

        Parameters
        ----------
        index : pandas.Index, optional
            Pandas Index to use.

            Default is to use the column Index of the provided ``matrix``.
        matrix : :obj:`str`, optional
            Name of analysis attribute containing a dataframe with pandas.MultiIndex columns to use.

            Default is to use the provided ``index``.
        levels : :obj:`list`, optional
            Levels of multiindex to restrict to.

            Defaults to all in index under use.
        pallete : :obj:`str`, optional
            Name of matplotlib color palete to use with categorical levels.
            See matplotlib.org/examples/color/colormaps_reference.html.

            Defaults to "tab20".
        {uniform_cmap, diverging_cmap} : :obj:`str`, optional
            Name of matplotlib color paletes to use with numerical levels.
            Uniform will be used if values in level are non-negative, while diverging if including negative.
            See matplotlib.org/examples/color/colormaps_reference.html.

            Defaults to "plasma" and "RdYlBu_r", respectively.
        nan_color : tuple, optional
            Color for missing (i.e. NA) values.

            Defaults to ``(0.662745, 0.662745, 0.662745, 0.5)`` == ``grey``.
        as_dataframe: :obj:`bool`, optional
            Whether a dataframe should be returned.

            Defaults to :obj:`False`.

        Returns
        -------
        {list, pandas.DataFrame}
            Matrix of shape (level, sample) with rgb values of each of the variable.
            If as_dataframe, this will be a pandas.DataFrame otherwise, list of lists.
        """
        import matplotlib
        import matplotlib.pyplot as plt

        if index is None:
            index = self.get_matrix(matrix).columns

        if levels is not None:
            drop = [l.name for l in index.levels if l.name not in levels]
            index = index.droplevel(drop)

        # Handle special case of single level
        if not isinstance(index, pd.core.indexes.multi.MultiIndex):
            index = pd.MultiIndex.from_arrays([index.values], names=[index.name])

        _pallete = plt.get_cmap(pallete)
        _uniform_cmap = plt.get_cmap(uniform_cmap)
        _diverging_cmap = plt.get_cmap(diverging_cmap)

        colors = list()
        for l, level in enumerate(index.levels):
            # For empty levels (all values nan), return nan colour
            if level.empty:
                colors.append([nan_color] * len(index))
                _LOGGER.warning("Level {} has only NaN values.".format(level.name))
                continue
            # determine the type of data in each level
            # TODO: check this works in all cases
            values = index.get_level_values(level.name)
            if sum([isinstance(x, str) for x in values]) >= 1:
                dtype = "categorical"
            else:
                dtype = "numerical or boolean"

            # Add either colors based on categories or numerical scale
            if dtype == "categorical":
                _LOGGER.debug("Level '{}' has a categorical type.".format(level.name, dtype))
                n = len(set(values))
                # get n equidistant colors
                p = [_pallete(1.0 * i / n) for i in range(n)]
                color_dict = dict(zip(list(set(index.get_level_values(level.name))), p))
                # color for nan cases
                color_dict[np.nan] = nan_color
                col = [color_dict[x] for x in index.get_level_values(level.name)]
            else:
                # Create a range of either 0-max if only positive values are found
                # or symmetrically from the maximum absolute value found
                if all((values.dropna() == True) | (values.dropna() == False)):
                    # boolean colormap
                    _LOGGER.debug("Level '{}' has a boolean type.".format(level.name, dtype))
                    norm = matplotlib.colors.Normalize(vmin=-0.2, vmax=1.2)
                    col = _diverging_cmap(norm(values.astype(float)))
                elif not any(values.dropna() < 0):
                    # uniform colormap with no negative values
                    # works with True/False and nan (the later are replaced after anyway)
                    _LOGGER.debug(
                        "Level '{}' has a numeric type with no negative values.".format(
                            level.name, dtype
                        )
                    )
                    norm = matplotlib.colors.Normalize(vmin=values.min(), vmax=values.max())
                    col = _uniform_cmap(norm(values.astype(float)))
                else:
                    # numeric diverging centered on zero
                    _LOGGER.debug(
                        "Level '{}' has a numeric type with negative values.".format(
                            level.name, dtype
                        )
                    )
                    v = np.nanmax(np.absolute(values))
                    norm = matplotlib.colors.Normalize(vmin=-v, vmax=v)
                    col = _diverging_cmap(norm(values.astype(float)))

                # replace color for nan cases
                col[
                    np.where(index.get_level_values(level.name).to_series().isnull().tolist())
                ] = nan_color
            # append vector (list) of sample values to list of levels
            colors.append([tuple(x) for x in col])

        if as_dataframe:
            return pd.DataFrame(colors, index=index.names, columns=index).T
        return colors

    def unsupervised_analysis(
        self,
        steps=["correlation", "manifold", "pca", "pca_association"],
        matrix="matrix_norm",
        samples=None,
        attributes_to_plot=None,
        output_dir="{results_dir}/unsupervised_analysis_{data_type}",
        output_prefix="all_{var_unit_name}s",
        standardize_matrix=True,
        manifold_algorithms=[
            "MDS",
            "Isomap",
            "LocallyLinearEmbedding",
            "SpectralEmbedding",
            "TSNE",
        ],
        maniford_kwargs={},
        display_corr_values=False,
        plot_max_pcs=4,
        save_additional=False,
        prettier_sample_names=True,
        rasterized=False,
        dpi=300,
        **kwargs
    ):
        """
        General unsupervised analysis of a matrix.

        Apply unsupervised clustering, manifold learning and dimensionality reduction
        methods on numeric matrix.
        Colours and labels samples by their attributes as given in `attributes_to_plot`.

        This analysis has 4 possible steps:
         - "correlation":
                Pairwise sample correlation with 2 distance metrics plotted as heatmap.
         - "manifold":
                Manifold learning of latent spaces for projection of samples.
                See here available algorithms:
                http://scikit-learn.org/stable/modules/classes.html#module-sklearn.manifold
         - "pca":
                For PCA analysis, if `test_pc_association` is `True`, will compute association of PCs
                with sample attributes given in `attributes_to_plot`. For numeric attributes,
                the Pearson correlation will be computed and for categoriacal, a pairwise
                Kruskal-Wallis H-test (ANOVA).
         - "pca_association":
                For PCA analysis, if `test_pc_association` is `True`, will compute association of PCs
                with sample attributes given in `attributes_to_plot`. For numeric attributes,
                the Pearson correlation will be computed and for categoriacal, a pairwise
                Kruskal-Wallis H-test (ANOVA).

        Parameters
        ----------
        steps : :obj:`list`, optional
            List of step keywords to be performed as described above.

            Defaults to all available.
        matrix : :obj:`str`, optional
            Name of analysis attribute contatining the numeric dataframe to perform analysis on.
            Must have a pandas.MultiIndex as column index.

            Defaults to "matrix_norm".
        samples : :obj:`list`, optional
            List of sample objects to restrict analysis to.

            Defaults to all in analysis.
        attributes_to_plot : :obj:`list`, optional
            List of attributes shared between sample groups should be plotted.

            Defaults to attributes in analysis.group_attributes.
        output_dir : :obj:`str`, optional
            Directory for generated files and plots.

            Defaults to "{results_dir}/unsupervised_analysis_{data_type}".
        output_prefix : :obj:`str`, optional
            Prefix for output files.

            Defaults to "all_regions" if data_type is ATAC-seq and "all_genes" if data_type is RNA-seq.
        standardize_matrix: :obj:`bool`, optional
            Whether to standardize variables in `matrix` by removing the mean and scaling to unit variance.
            It is not applied to the "correlation" step.

            Default is :obj:`True`.
        manifold_algorithms : :obj:`list`, optional
            List of manifold algorithms to use. See available algorithms here:
            http://scikit-learn.org/stable/modules/classes.html#module-sklearn.manifold

            Defaults to ['MDS', 'Isomap', 'LocallyLinearEmbedding', 'SpectralEmbedding', 'TSNE'],
        maniford_kwargs : :obj:`dict`, optional
            Dictionary of keyword arguments to pass to the algorithms in ``manifold_algorithms``.
            Should be of the form {"algorithm_name": {"key": value}}
        display_corr_values: :obj:`bool`, optional
            Whether values in heatmap of sample correlations should be displayed overlaid on top of colours.

            Defaults to :obj:`False`.
        save_additional: :obj:`bool`, optional
            Whether additional results such as PCA projection, loadings should be saved.

            Defaults to :obj:`False`.
        prettier_sample_names: :obj:`bool`, optional
            Whether it should attempt to prettify sample names by removing the data type from plots.

            Defaults to :obj:`True`.
        rasterized: :obj:`bool`, optional
            Whether elements with many objects should be rasterized.

            Defaults to :obj:`False`.
        dpi : :obj:`int`, optional
            Definition of rasterized image in dots per inch (dpi).

            Defaults to 300.
        **kwargs: optional
            kwargs are passed to :func:`~ngs_toolkit.Analysis.get_level_colors` and
            :func:`~ngs_toolkit.graphics.plot_projection`.
        """
        # TODO: Treat PCA as just one of several possible decomposition methods
        # by either:
        #  - creating a new, independent decomposition section
        #  - lumping all under common decomposition/manifold and handle it internally,
        #    add association testing under that for any arbitrary
        from collections import defaultdict
        import itertools
        import matplotlib.pyplot as plt
        from ngs_toolkit.graphics import savefig, plot_projection
        import seaborn as sns
        from sklearn import manifold
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
        from statsmodels.sandbox.stats.multicomp import multipletests
        from scipy.stats import kruskal, pearsonr

        output_dir = self._format_string_with_attributes(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        matrix = self.get_matrix(matrix)

        if ("pca_association" in steps) and ("pca" not in steps):
            steps.append("pca")

        output_prefix = self._format_string_with_attributes(output_prefix)

        if not isinstance(matrix.columns, pd.core.indexes.multi.MultiIndex):
            msg = "Provided quantification matrix must have columns with MultiIndex."
            hint = " Will try to use `analysis.annotate_samples` to do that."
            _LOGGER.info(msg + hint)
            matrix = self.annotate_samples(matrix=matrix, save=False, assign=False)

        if samples is None:
            samples = [
                s for s in self.samples if s.name in matrix.columns.get_level_values("sample_name")
            ]
        else:
            samples = [
                s for s in samples if s.name in matrix.columns.get_level_values("sample_name")
            ]
        if len(samples) == 0:
            msg = "None of the samples could be found in the quantification matrix."
            _LOGGER.error(msg)
            raise ValueError(msg)
        if len(samples) == 1:
            msg = "Only one sample could be found in the quantification matrix."
            hint = " Function needs more than one."
            _LOGGER.error(msg + hint)
            raise ValueError(msg)

        msg = "`attributes_to_plot` were not specified and the analysis does not have a "
        msg += " 'group_attributes' variable."
        if attributes_to_plot is None:
            try:
                attributes_to_plot = self.group_attributes
            except AttributeError:
                _LOGGER.error(msg)
                raise
        # Raise error when requested factor is not known
        miss = [attr for attr in attributes_to_plot if attr not in matrix.columns.names]
        if len(miss) > 0:
            msg = "Requested '{}' value is not present as column level of matrix.".format(
                ", ".join(miss)
            )
            _LOGGER.error(msg)
            raise ValueError(msg)
        # remove attributes not in matrix
        attributes_to_plot = [attr for attr in attributes_to_plot if attr in matrix.columns.names]
        # remove attributes with all NaNs
        attributes_to_plot = [
            attr
            for attr in attributes_to_plot
            if not pd.isnull(matrix.columns.get_level_values(attr)).all()
        ]
        if len(attributes_to_plot) == 0:
            msg = (
                "None of the factors in `attributes_to_plot` could be found in the "
                + "quantification matrix index or they are all NaN."
            )
            _LOGGER.error(msg)
            raise ValueError(msg)

        # All regions, matching samples (provided samples in matrix)
        x = matrix.loc[
            :, matrix.columns.get_level_values("sample_name").isin([s.name for s in samples]),
        ]

        # Get matrix of samples vs levels with colors as values
        cd_kwargs = {
            k: v for k, v in kwargs.items() if k in ["pallete", "uniform_cmap", "diverging_cmap"]
        }

        # # it will always be a matrix for all samples
        color_dataframe = self.get_level_colors(
            index=matrix.columns,
            levels=["sample_name"] + attributes_to_plot
            if "sample_name" not in attributes_to_plot
            else attributes_to_plot,
            as_dataframe=True,
            **cd_kwargs
        )
        if "sample_name" not in attributes_to_plot:
            color_dataframe = color_dataframe.drop("sample_name", axis=1)
        # will be filtered now by the requested samples if needed
        color_dataframe = color_dataframe.loc[x.columns.get_level_values("sample_name"), :]

        projection_kwargs = {
            k: v
            for k, v in kwargs.items()
            if k
            in [
                "plot_max_dims",
                "rasterized",
                "plot_group_centroids",
                "axis_ticklabels",
                "axis_ticklabels_name",
                "axis_lines",
                "legends",
                "always_legend",
            ]
        }

        if isinstance(x.columns, pd.MultiIndex):
            sample_display_names = x.columns.get_level_values("sample_name")
        else:
            sample_display_names = x.columns

        if "correlation" in steps:
            # Pairwise correlations
            for method in ["pearson", "spearman"]:
                _LOGGER.info("Plotting pairwise correlation with '{}' metric.".format(method))
                xp = x.copy()
                xp.columns = xp.columns.get_level_values("sample_name")
                xp = xp.astype(float).corr(method)

                cd = color_dataframe
                cd.index = cd.index.get_level_values("sample_name")
                g = sns.clustermap(
                    xp,
                    xticklabels=False,
                    yticklabels=sample_display_names,
                    annot=display_corr_values,
                    cmap="Spectral_r",
                    figsize=(0.4 * x.shape[1], 0.4 * x.shape[1]),
                    cbar_kws={"label": "{} correlation".format(method.capitalize())},
                    row_colors=cd,
                    col_colors=cd,
                )
                g.ax_heatmap.set_yticklabels(
                    g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small"
                )
                g.ax_heatmap.set_xlabel(None, visible=False)
                g.ax_heatmap.set_ylabel(None, visible=False)
                savefig(
                    g.fig,
                    os.path.join(
                        output_dir,
                        "{}.{}.{}_correlation.clustermap.svg".format(
                            self.name, output_prefix, method
                        ),
                    ),
                )

        if standardize_matrix:
            std = StandardScaler()
            x = pd.DataFrame(std.fit_transform(x.T).T, index=x.index, columns=x.columns)

        if "manifold" in steps:
            # Manifolds
            # TODO: test usage of non default manifolds from sklearn
            params = defaultdict(dict)
            params.update(
                {
                    "MDS": {"n_jobs": -1},
                    "Isomap": {"n_jobs": -1},
                    "LocallyLinearEmbedding": {},
                    "SpectralEmbedding": {"n_jobs": -1},
                    "TSNE": {"init": "pca"},
                }
            )
            params.update(maniford_kwargs)
            for algo in manifold_algorithms:
                msg = "Learning manifold with '{}' algorithm".format(algo)
                _LOGGER.info(msg + ".")

                manif = getattr(manifold, algo)(**params[algo])
                try:
                    x_new = manif.fit_transform(x.T)
                except (TypeError, ValueError):
                    hint = " Number of samples might be too small to perform '{}'".format(algo)
                    _LOGGER.error(msg + " failed!" + hint)
                    continue

                x_new = pd.DataFrame(x_new, index=x.columns, columns=list(range(x_new.shape[1])))
                if save_additional:
                    for d, label in [(x_new, "embedding")]:
                        _LOGGER.debug("Saving {} {} matrix to disk.".format(algo, label))
                        d.to_csv(
                            os.path.join(
                                output_dir,
                                "{}.{}.{}.{}.csv".format(
                                    self.name, output_prefix, algo.lower(), label
                                ),
                            )
                        )

                _LOGGER.info("Plotting projection of manifold with '{}' algorithm.".format(algo))
                plot_projection(
                    df=x_new,
                    color_dataframe=color_dataframe,
                    dims=1,
                    output_file=os.path.join(
                        output_dir, "{}.{}.{}.svg".format(self.name, output_prefix, algo.lower()),
                    ),
                    attributes_to_plot=attributes_to_plot,
                    axis_ticklabels_name=algo,
                    **projection_kwargs
                )

        if "pca" in steps:
            # PCA
            pcs = min(*x.shape) - 1
            _LOGGER.info("Decomposing data with 'PCA' algorithm for {} dimensions.".format(pcs))
            pca = PCA(n_components=pcs, svd_solver="arpack")
            x_new = pca.fit_transform(x.T)

            pcs_order = range(pca.n_components_)
            x_new = pd.DataFrame(x_new, index=x.columns, columns=pcs_order)
            comps = pd.DataFrame(pca.components_.T, index=x.index, columns=pcs_order)

            if save_additional:
                for d, label in [(x_new, "embedding"), (comps, "loading")]:
                    _LOGGER.debug("Saving PCA {} matrix to disk.".format(label))
                    d.to_csv(
                        os.path.join(
                            output_dir, "{}.{}.pca.{}.csv".format(self.name, output_prefix, label)
                        )
                    )

            # Write % variance expained to disk
            variance = pd.Series(
                pca.explained_variance_ratio_ * 100, name="percent_variance", index=pcs_order,
            ).to_frame()
            variance["log_variance"] = np.log10(pca.explained_variance_)
            variance.index = variance.index.values
            variance.index.name = "PC"
            variance.to_csv(
                os.path.join(
                    output_dir, "{}.{}.pca.explained_variance.csv".format(self.name, output_prefix),
                )
            )

            # plot % explained variance per PC
            _LOGGER.info("Plotting variance explained with PCA.")
            fig, axis = plt.subplots(1, 3, figsize=(4 * 3, 4))
            axis[0].plot(variance.index.values, variance["percent_variance"], "o-")
            axis[0].set_ylim(
                (0, variance["percent_variance"].max() + variance["percent_variance"].max() * 0.1,)
            )
            axis[1].plot(variance.index.values, variance["log_variance"], "o-")
            axis[2].plot(variance.index.values, variance["percent_variance"].cumsum(), "o-")
            axis[2].set_ylim((0, 100))
            for ax in axis:
                ax.axvline(len(attributes_to_plot), linestyle="--")
                ax.set_xlabel("PC")
            axis[0].set_ylabel("% variance")
            axis[1].set_ylabel("log variance")
            axis[2].set_ylabel("Cumulative % variance")
            sns.despine(fig)
            savefig(
                fig,
                os.path.join(
                    output_dir, "{}.{}.pca.explained_variance.svg".format(self.name, output_prefix),
                ),
            )

            # plot pca
            pcs = min(x_new.shape[1] - 1, plot_max_pcs)

            if pcs > 0:
                _LOGGER.info("Plotting PCA up to ({} + 1) dimensions.".format(pcs))
                plot_projection(
                    df=x_new,
                    color_dataframe=color_dataframe,
                    dims=pcs,
                    output_file=os.path.join(
                        output_dir, "{}.{}.pca.svg".format(self.name, output_prefix)
                    ),
                    attributes_to_plot=attributes_to_plot,
                    axis_ticklabels_name="PCA",
                    **projection_kwargs
                )
            else:
                _LOGGER.warning("Only one PC selected, cannot plot projection.")

        if "pca_association" in steps:
            # Test association of PCs with attributes
            _LOGGER.info("Computing association of given attributes with principal components.")
            associations = list()
            for attr in attributes_to_plot:
                # Get all values of samples for this attr
                groups = x_new.index.get_level_values(attr).unique().dropna()

                if groups.nunique() == 1:
                    _LOGGER.warning(
                        "Attribute '{}' cannot be tested because all values are equal or null.".format(
                            attr
                        )
                    )
                    continue

                # Determine if attr is categorical or continuous
                if all([isinstance(i, (str, bool)) for i in groups]):
                    variable_type = "categorical"
                elif all([isinstance(i, (int, float, np.int, np.float)) for i in groups]):
                    variable_type = "numerical"
                else:
                    _LOGGER.warning(
                        "Attribute '{}' cannot be tested because data type is not undestood.".format(
                            attr
                        )
                    )
                    variable_type = "not-detected"

                _LOGGER.debug("Attribute '{}' is of type {}.".format(attr, variable_type))

                for pc in pcs_order:
                    _LOGGER.debug("Attribute '{}'; PC {}.".format(attr, pc + 1))
                    if variable_type == "categorical":
                        # It categorical, test pairwise combinations of attributes
                        for group1, group2 in itertools.combinations(groups, 2):
                            _LOGGER.debug("Testing group '{}' with '{}'.".format(group1, group2))
                            g1_mask = x_new.index.get_level_values(attr) == group1
                            g2_mask = x_new.index.get_level_values(attr) == group2

                            g1_values = x_new.loc[g1_mask, pc]
                            g2_values = x_new.loc[g2_mask, pc]

                            # Test ANOVA (or Kruskal-Wallis H-test)
                            p = kruskal(g1_values, g2_values)[1]

                            # Append
                            associations.append([pc + 1, attr, variable_type, group1, group2, p])

                    elif variable_type == "numerical":
                        # It numerical, calculate pearson correlation
                        pc_values = x_new.loc[:, pc]
                        trait_values = x_new.index.get_level_values(attr)
                        p = pearsonr(pc_values, trait_values)[1]

                        associations.append([pc + 1, attr, variable_type, np.nan, np.nan, p])
                    else:
                        associations.append([pc + 1, attr, variable_type, np.nan, np.nan, np.nan])

            associations = pd.DataFrame(
                associations,
                columns=["pc", "attribute", "variable_type", "group_1", "group_2", "p_value",],
            )

            if associations.empty:
                msg = "Couldn't test any associations between PCs and factors."
                hint = " Perhaps PCA produced only 1 PC or the type of the factors could not be detected?"
                _LOGGER.warning(msg + hint)
                return

            # correct p-values
            associations.loc[:, "adj_pvalue"] = multipletests(
                associations["p_value"], method="fdr_bh"
            )[1]

            # write
            _LOGGER.info("Saving associations.")
            associations.to_csv(
                os.path.join(
                    output_dir,
                    "{}.{}.pca.variable_principle_components_association.csv".format(
                        self.name, output_prefix
                    ),
                ),
                index=False,
            )

            if associations["attribute"].nunique() < 2:
                _LOGGER.info("Few attributes tested, can't plot associations.")
                return

            # Plot
            for var in ["p_value", "adj_pvalue"]:
                pivot = (
                    associations.groupby(["pc", "attribute"])[var]
                    .min()
                    .reset_index()
                    .pivot(index="pc", columns="attribute", values=var)
                    .dropna(axis=1)
                )

                # skip if dataframe has no variation
                if (pivot.nunique() == 1).all().all():
                    continue

                # heatmap of -log p-values
                g = sns.clustermap(
                    -np.log10(pivot),
                    row_cluster=False,
                    annot=True,
                    cbar_kws={"label": "-log10(p_value) of association"},
                    rasterized=rasterized,
                    vmin=0,
                )
                g.ax_heatmap.set_xticklabels(
                    g.ax_heatmap.get_xticklabels(), rotation=45, ha="right"
                )
                savefig(
                    g.fig,
                    os.path.join(
                        output_dir,
                        "{}.{}.pca.variable_principle_components_association.{}.svg".format(
                            self.name, output_prefix, var
                        ),
                    ),
                )

                # heatmap of masked significant
                g = sns.clustermap(
                    (pivot < 0.05).astype(int),
                    row_cluster=False,
                    cbar_kws={"label": "significant association"},
                    rasterized=rasterized,
                    vmin=0,
                    vmax=1,
                    cmap="Paired",
                )
                g.ax_heatmap.set_xticklabels(
                    g.ax_heatmap.get_xticklabels(), rotation=45, ha="right"
                )
                savefig(
                    g.fig,
                    os.path.join(
                        output_dir,
                        "{}.{}.pca.variable_principle_components_association.{}.masked.svg".format(
                            self.name, output_prefix, var
                        ),
                    ),
                )

    def differential_analysis(
        self,
        comparison_table=None,
        samples=None,
        covariates=None,
        filter_support=False,
        output_dir="{results_dir}/differential_analysis_{data_type}",
        output_prefix="differential_analysis",
        overwrite=True,
        distributed=False,
        deseq_kwargs=None,
        **kwargs
    ):
        """
        Perform differential regions/genes across samples that are
        associated with a certain trait.
        Currently the only implementation is with DESeq2.
        This implies the rpy2 library and the respective R library
        are installed.

        Requires the R package "DESeq2" to be installed:

        .. highlight:: R
        .. code-block:: R

            if (!requireNamespace("BiocManager", quietly = TRUE))
                install.packages("BiocManager")
            BiocManager::install("DESeq2")

        For other implementations of differential analysis see
        `ngs_toolkit.general.least_squares_fit`
        and `ngs_toolkit.general.differential_from_bivariate_fit`.

        Parameters
        ----------
        comparison_table : :obj:`pandas.DataFrame`
            A dataframe with "comparison_name", "comparison_side" and
            "sample_name", "sample_group" columns.

            Defaults to the analysis' own "comparison_table" attribute.
        samples : :obj:`list`, optional
            Samples to limit analysis to.

            Defaults to all samples in analysis object.
        covariates : :obj:`list`, optional
            Additional variables to take into account in the model fitting.

            Defaults to None.
        filter_support: :obj:`bool`, optional
            Whether features not supported in a given comparison should
            be removed (i.e. regions with no peaks in any sample in a
            comparison are not tested).
            Applies only to ATAC-/ChIP-seq data.

            Default is :obj:`True`.
        output_dir : :obj:`str`, optional
            Output directory for analysis.
            Variables in curly braces will be formated
            with attributes from analysis.

            Defaults to "{results_dir}/differential_analysis_{data_type}".
        output_prefix : :obj:`str`, optional
            Prefix for output files.

            Defaults to "differential_analysis".
        overwrite: :obj:`bool`, optional
            Whether results should be overwritten in case they already exist.

            Defaults to :obj:`True`.
        distributed: :obj:`bool`, optional
            Whether analysis should be distributed in a computing cluster
            for each comparison.
            Additional configuration can be passed in ``kwargs``.

            Defaults to :obj:`False`.
        deseq_kwargs : :obj:`dict`, optional
            Additional keyword arguments to be passed to the
            `DESeq` function of DESeq2.

        kwargs : :obj:`dict`, optional
            Additional keyword arguments are passed to
            :func:`~ngs_toolkit.utils.submit_job` and then to the
            chosen `divvy` submission template according to
            `computing_configuration`.
            Pass for example `cores=4, mem=8000, partition="longq",
            time="08:00:00"`.

        Returns
        -------
        :class:`pandas.DataFrame`
            Results for all comparisons.
            Will be :obj:`None` if `distributed` is `True`.

        Attributes
        ----------
        differential_results : :obj:`pandas.DataFrame`
            Pandas dataframe with results.
        """
        import sys

        from ngs_toolkit.general import deseq_analysis
        from ngs_toolkit.utils import submit_job

        if comparison_table is None:
            msg = "`comparison_table` was not given and is not set in analysis object."
            hint = "Add a `comparison_table` attribute to the analysis object."
            try:
                comparison_table = self.comparison_table
            except AttributeError as e:
                _LOGGER.error(msg)
                _LOGGER.info(hint)
                raise e

        if deseq_kwargs is None:
            deseq_kwargs = {}

        # Check comparisons
        # check comparison table has required columns
        req_attrs = [
            "comparison_name",
            "comparison_side",
            "sample_name",
            "sample_group",
        ]
        if not all([x in comparison_table.columns for x in req_attrs]):
            raise AssertionError(
                "Given comparison table does not have all of '{}' columns.".format(
                    "', '".join(req_attrs)
                )
            )
        # check all comparisons have samples in two sides
        if not all(comparison_table.groupby("comparison_name")["comparison_side"].nunique() == 2):
            msg = "All comparisons must have samples in each side of the comparison."
            raise AssertionError(msg)
        # check if any comparison and sample group has samples disagreeing in side
        if not all(
            comparison_table.groupby(["comparison_name", "sample_group"])[
                "comparison_side"
            ].nunique()
            == 1
        ):
            msg = "Samples in same comparison and group must agree"
            msg += " on their side of the comparison."
            raise AssertionError(msg)

        # Handle samples under self
        if samples is None:
            samples = self.samples
        samples = [s for s in samples if s.name in comparison_table["sample_name"].tolist()]

        # Make output dir
        output_dir = self._format_string_with_attributes(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Get matrix and samples
        count_matrix = self.matrix_raw
        if samples is None:
            samples = self.samples
        samples = [
            s
            for s in samples
            if (s.name in comparison_table["sample_name"].tolist())
            & (s.name in count_matrix.columns)
        ]
        count_matrix = count_matrix[[s.name for s in samples]]

        # Get experiment matrix
        # by getting any other relevant covariates as required
        if covariates is not None:
            sample_table = pd.DataFrame([s.as_series() for s in samples])
            # check all covariates are in the samples and none is null
            if not all([x in sample_table.columns for x in covariates]):
                msg = "Not all of the specified covariate variables are in the selected samples."
                raise AssertionError(msg)
            if sample_table[covariates].isnull().any().any():
                msg = "None of the selected samples can have a Null value in the specified covariate variables."
                _LOGGER.error(msg)
                raise AssertionError(msg)

            # add covariates to comparison table
            comparison_table = (
                comparison_table.set_index("sample_name")
                .join(sample_table.set_index("sample_name")[covariates])
                .reset_index()
            )

        # Make table for DESeq2
        experiment_matrix = comparison_table.loc[
            :, ["sample_name", "sample_group"] + (covariates if covariates is not None else []),
        ].drop_duplicates()

        # Check whether the is a complex design
        complx = (comparison_table.groupby("sample_name")["sample_group"].nunique() > 1).any()

        if complx and (not distributed):
            distributed = True
            msg = "Detected complex design with samples in various groups."
            msg += " Will run analysis in distributed mode."
            _LOGGER.info(msg)

        # Make formula for DESeq2
        formula = "~ {}sample_group".format(
            " + ".join(covariates) + " + " if covariates is not None else ""
        )

        # Run DESeq2 analysis
        if not distributed:
            # filter features without support for given comparison
            if (self.data_type in ["ATAC-seq", "ChIP-seq"]) and filter_support:
                if not hasattr(self, "support"):
                    msg = "`filter_support` enabled but analysis has no `support` attribute!"
                    _LOGGER.error(msg)
                    raise ValueError(msg)
                _LOGGER.debug("Filtering out unsupported regions.")
                sup = self.support.loc[
                    :,
                    [
                        s.name
                        for s in self.samples
                        if s.name in comparison_table["sample_name"].tolist()
                    ],
                ]
                count_matrix = count_matrix.loc[(sup > 0).any(axis=1), :]
            results = deseq_analysis(
                count_matrix,
                experiment_matrix,
                comparison_table,
                formula,
                output_dir,
                output_prefix,
                overwrite=overwrite,
                **deseq_kwargs
            )
            try:
                results = results.set_index("index")
            except KeyError:
                pass
            _LOGGER.info(
                "Setting results of differential analysis to a variable 'differential_results'."
            )
            self.differential_results = results
            return results

        else:
            for comparison_name in comparison_table["comparison_name"].drop_duplicates():
                # make directory for comparison input/output
                out = os.path.join(os.path.abspath(output_dir), comparison_name)
                if not os.path.exists(out):
                    os.makedirs(out)

                comp = comparison_table.loc[
                    comparison_table["comparison_name"] == comparison_name, :
                ]
                comp.to_csv(os.path.join(out, "comparison_table.csv"), index=False)

                exp = experiment_matrix.loc[
                    experiment_matrix["sample_name"].isin(comp["sample_name"].tolist())
                    & experiment_matrix["sample_group"].isin(comp["sample_group"].tolist()),
                    :,
                ]
                exp.to_csv(os.path.join(out, "experiment_matrix.csv"), index=False)

                count = count_matrix.loc[:, comp["sample_name"].drop_duplicates()]
                # filter features without support for given comparison
                if self.data_type in ["ATAC-seq", "ChIP-seq", "ChIPmentation"] and filter_support:
                    if not hasattr(self, "support"):
                        msg = "`filter_support` enabled by analysis has no `support` attribute!"
                        _LOGGER.error(msg)
                        raise ValueError(msg)
                    _LOGGER.debug(
                        "Filtering out unsupported regions for comparison '{}'.".format(
                            comparison_name
                        )
                    )
                    sup = self.support.loc[
                        :, [s.name for s in self.samples if s.name in comp["sample_name"].tolist()],
                    ]
                    sup = sup.reindex(count.index).dropna()
                    count = count.reindex(sup[(sup > 0).any(axis=1)].index)
                count.to_csv(os.path.join(out, "count_matrix.csv"), index=True)

                # Assemble job and submit
                job_name = "deseq_job.{}".format(comparison_name)
                log_file = os.path.join(out, job_name + ".log")
                job_file = os.path.join(out, job_name + ".sh")
                cmd = (
                    "{executable} -m ngs_toolkit.recipes.deseq2 "
                    "--no-save-inputs --output-prefix {output_prefix} "
                    "--formula '{formula}' "
                    "{overwrite} {out}"
                ).format(
                    executable=sys.executable,
                    output_prefix=output_prefix,
                    formula=formula,
                    overwrite=" --overwrite" if overwrite else "",
                    out=out,
                )
                submit_job(cmd, job_file, log_file=log_file, jobname=job_name, **kwargs)
            if "computing_configuration" in kwargs:
                if kwargs["computing_configuration"] not in ["localhost", "default"]:
                    return
            return self.collect_differential_analysis(
                comparison_table=comparison_table, overwrite=overwrite
            )

    def collect_differential_analysis(
        self,
        comparison_table=None,
        input_dir="{results_dir}/differential_analysis_{data_type}",
        input_prefix="differential_analysis",
        output_dir="{results_dir}/differential_analysis_{data_type}",
        output_prefix="differential_analysis",
        permissive=True,
        save=True,
        assign=True,
        overwrite=False,
    ):
        """
        Collect results from DESeq2 differential analysis.
        Particularly useful when running ``differential_analysis`` in distributed mode.

        Parameters
        ----------
        comparison_table : :obj:`pandas.DataFrame`
            A dataframe with "comparison_name", "comparison_side" and "sample_name", "sample_group" columns.

            Defaults to the analysis's own "comparison_table" attribute.
        input_dir, output_dir : :obj:`str`, optional
            In-/Output directory of files.
            Values within curly brackets "{data_type}", will be formated with attributes from analysis.

            Defaults to "{results_dir}/differential_analysis_{data_type}".
        input_prefix, output_prefix : :obj:`str`, optional
            Prefix of the in-/output files.

            Defaults for both is "differential_analysis".
        permissive : :obj:`bool`, optional
            Whether non-existing files should be skipped or an error be thrown.

            Defaults to :obj:`True`.
        save : :obj:`bool`, optional
            Whether to save results to disk.

            Defaults to :obj:`True`.
        assign : :obj:`bool`, optional
            Whether to add results to a `differential_results` attribute.

            Defaults to :obj:`True`.
        overwrite: :obj:`bool`, optional
            Whether results should be overwritten in case they already exist.

            Defaults to :obj:`False`.

        Returns
        -------
        :class:`pandas.DataFrame`
            Results for all comparisons.
            Will be :obj:`None` if ``overwrite`` is :obj:`False` and a results file already exists.

        Attributes
        ----------
        differential_results : :obj:`pandas.DataFrame`
            Pandas dataframe with results.
        """
        from tqdm import tqdm

        if comparison_table is None:
            msg = "`comparison_table` was not given and is not set in analysis object."
            hint = "Add a `comparison_table` attribute to the analysis object."
            try:
                comparison_table = self.comparison_table
            except AttributeError as e:
                _LOGGER.error(msg)
                _LOGGER.info(hint)
                raise e

        input_dir = self._format_string_with_attributes(input_dir)
        # Make output dir
        output_dir = self._format_string_with_attributes(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        results_file = os.path.join(output_dir, output_prefix + ".deseq_result.all_comparisons.csv")
        if not overwrite and os.path.exists(results_file):
            msg = "Differential analysis results '{}' already exist and argument `overwrite` is False."
            hint = " Will not do ``anything``."
            _LOGGER.warning(msg.format(results_file) + hint)
            return

        # Subset comparisons by data type
        if ("data_type" in comparison_table.columns) and (self.data_type is not None):
            _LOGGER.debug("Subsetting comparisons for data type '{}'".format(self.data_type))
            comps = (
                comparison_table.loc[
                    comparison_table["data_type"] == self.data_type, "comparison_name"
                ]
                .drop_duplicates()
                .sort_values()
            )
            if comps.empty:
                msg = "After subsetting comparison data_types for '{}', table was empty."
                msg += " Continuing with all comparisons in table."
                _LOGGER.warning(msg)
                comps = comparison_table["comparison_name"].drop_duplicates().sort_values()
        else:
            msg = "Comparison table does not have a 'data_type' column."
            msg += " Collecting all comparisons in table."
            _LOGGER.warning(msg)
            comps = comparison_table["comparison_name"].drop_duplicates().sort_values()

        results = list()
        for comp in tqdm(comps, total=len(comps), desc="Comparison"):
            res_file = os.path.join(
                input_dir, comp, input_prefix + ".deseq_result.{}.csv".format(comp)
            )
            _LOGGER.debug("Collecting comparison '{}'".format(comp))
            try:
                res2 = pd.read_csv(res_file, index_col=0)
            except IOError as e:
                if permissive:
                    _LOGGER.warning(
                        "Results file for comparison '{}' do not exist. Skipping.".format(comp)
                    )
                    continue
                else:
                    raise e
            results.append(res2.reset_index())

        if not results:
            msg = "No comparison had a valid results file!"
            if permissive:
                _LOGGER.warning(msg)
                return
            else:
                _LOGGER.error(msg)
                raise IOError(msg)

        results = pd.concat(results)

        # set index
        if "index" in results.columns:
            results = results.set_index("index")

        if save:
            results.to_csv(results_file, index=True)

        if assign:
            self.differential_results = results
        return results

    def plot_differential(
        self,
        steps=[
            "distributions",
            "counts",
            "scatter",
            "volcano",
            "ma",
            "stats_heatmap",
            "correlation",
            "heatmap",
        ],
        results=None,
        comparison_table=None,
        samples=None,
        matrix="matrix_norm",
        only_comparison_samples=False,
        alpha=0.05,
        corrected_p_value=True,
        fold_change=None,
        diff_based_on_rank=False,
        max_rank=1000,
        ranking_variable="pvalue",
        respect_stat_thresholds=True,
        output_dir="{results_dir}/differential_analysis_{data_type}",
        output_prefix="differential_analysis",
        plot_each_comparison=True,
        mean_column="baseMean",
        log_fold_change_column="log2FoldChange",
        p_value_column="pvalue",
        adjusted_p_value_column="padj",
        comparison_column="comparison_name",
        rasterized=True,
        robust=False,
        feature_labels=False,
        group_colours=True,
        group_attributes=None,
        **kwargs
    ):
        """
        Plot differential features (eg chromatin region, genes) discovered with supervised
        group comparisons by ``ngs_toolkit.general.differential_analysis``.
        This will plot number and direction of discovered features, scatter, MA and volcano
        plots for each comparison and joint heatmaps of log fold changes, normalized values
        or Z-scores of individual samples or groups in the differential features.

        Parameters
        ----------
        steps : :obj:`list`, optional
            Types of plots to make:
                - "distributions": Distribution of p-values and fold-changes
                - "counts" - Count of differential features per comparison given certain thresholds.
                - "scatter" - Scatter plots (group 1 vs group 2).
                - "volcano" - Volcano plots (log fold change vs -log p-value)
                - "ma" - MA plots (log mean vs log fold change)
                - "stats_heatmap" - Heatmap of p-values and fold-changes for comparisons.
                - "correlation" - Correlation of samples or sample groups in differential features.
                - "heatmap" - Heatmaps of samples or sample groups in differential features.

            Defaults to all of the above.
        results : :obj:`pandas.DataFrame`, optional
            Data frame with differential analysis results.
            See ``ngs_toolkit.general.differential_analysis`` for more information.
        comparison_table : :obj:`pandas.DataFrame`, optional
            Comparison table. If provided, group-wise plots will be produced.

            Defaults to the analysis' "comparison_table" attribute.
        samples : :obj:`list`, optional
            List of sample objects to restrict analysis to.

            Defaults to all samples in analysis.
        matrix : :obj:`str`, optional
            Matrix of quantification to use for plotting feature values across samples/groups.

            Defaults to "matrix_norm".
        only_comparison_samples: :obj:`bool`, optional
            Whether to use only samples present in the `comparison_table` and `results` table.

            Defaults to :obj:`False`.
        alpha : float, optional
            Significance level to consider a feature differential.

            Defaults to 0.05.
        corrected_p_value: :obj:`bool`, optional
            Whether to use a corrected p-valueto consider a feature differential.

            Defaults to :obj:`True`.
        fold_change : float, optional
            Effect size (log2 fold change) to consider a feature differential. Considers absolute values.

            Default is no log2 fold change threshold.
        diff_based_on_rank: :obj:`bool`, optional
            Whether a feature should be considered differential based on its rank.
            Use in combination with `max_rank`, `ranking_variable` and `respect_stat_thresholds`.

            Defaults to :obj:`False`.
        max_rank : :obj:`int`, optional
            Rank to use when using `diff_based_on_rank`.

            Defaults to 1000.
        ranking_variable : :obj:`str`, optional
            Which variable to use for ranking when using `diff_based_on_rank`.

            Defaults to "pvalue".
        respect_stat_thresholds: :obj:`bool`, optional
            Whether the statistical thresholds from `alpha` and `fold_change` should still be
            respected when using `diff_based_on_rank`.

            Defaults to :obj:`True`.
        output_dir : :obj:`str`, optional
            Directory to create output files.

            Defaults to "{results_dir}/differential_analysis_{data_type}"
        output_prefix : :obj:`str`, optional
            Prefix to use when creating output files.

            Defaults to "differential_analysis".
        plot_each_comparison: :obj:`bool`, optional
            Whether each comparison should be plotted in scatter, MA and volcano plots.
            Useful to turn off with many comparisons.

            Defaults to :obj:`True`.
        mean_column : :obj:`str`, optional
            Column  in `results` data frame containing values for mean values across samples.

            Defaults to "baseMean".
        log_fold_change_column : :obj:`str`, optional
            Column in `results` data frame containing values for log2FoldChange values across samples.

            Defaults to "log2FoldChange".
        p_value_column : :obj:`str`, optional
            Column  in `results` data frame containing values for p-values across samples.

            Defaults to "pvalue".
        adjusted_p_value_column : :obj:`str`, optional
            Column  in `results` data frame containing values for adjusted p-values across samples.

            Defaults to "padj".
        comparison_column : :obj:`str`, optional
            Column  in `results` data frame containing the name of the comparison.

            Defaults to "comparison_name".
        rasterized: :obj:`bool`, optional
            Whether plots with many objects should be rasterized.

            Defaults to :obj:`True`.
        robust: :obj:`bool`, optional
            Whether heatmap color scale ranges should be robust (using quantiles) rather than extreme values.
            Useful for noisy/extreme data.

            Defaults to :obj:`False`.
        feature_labels: :obj:`bool`, optional
            Whether features (regions/genes) should be labeled in heatmaps.

            Defaults to :obj:`False`.
        group_colours: :obj:`bool`, optional
            Whether groups of samples should be coloured in heatmaps.

            Defaults to :obj:`True`.
        group_attributes : :obj:`list`, optional
            Which variables to colour if `group_colours` if :obj:`True`.

            Defaults to all of analysis.group_attributes.
        **kwargs: :obj:`dict`, optional
            Additional keyword arguments will be passed to `Analysis.get_level_colors`.
        """
        # TODO: split plotting into smaller parts
        import matplotlib.pyplot as plt
        from ngs_toolkit.graphics import add_colorbar_to_axis, savefig
        import seaborn as sns

        if results is None:
            msg = "Differential results dataframe not given and Analysis object does not"
            msg += " have a `differential_results` attribute."
            hint = " Run differential_analysis to produce differential results."
            try:
                results = self.differential_results
            except AttributeError as e:
                _LOGGER.error(msg + hint)
                raise e
            if (results is None) or (not isinstance(results, pd.DataFrame)):
                hint = " Run differential_self to produce differential results."
                _LOGGER.error(msg)
                raise ValueError
        results = results.copy()

        req_attrs = [
            mean_column,
            log_fold_change_column,
            p_value_column,
            adjusted_p_value_column,
            comparison_column,
        ]
        if not all([x in results.columns for x in req_attrs]):
            raise AssertionError(
                "Results dataframe must have '{}' columns.".format(", ".join(req_attrs))
            )

        # Make output dir
        output_dir = self._format_string_with_attributes(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Get matrix and samples
        if samples is None:
            samples = self.samples

        matrix = self.get_matrix(matrix, samples)

        if not isinstance(matrix.columns, pd.core.indexes.multi.MultiIndex):
            msg = "Provided quantification matrix must have columns with MultiIndex."
            hint = " Will try to use `analysis.annotate_samples` to do that."
            _LOGGER.info(msg + hint)
            matrix = self.annotate_samples(matrix=matrix, save=False, assign=False)

        samples = [s for s in samples if s.name in matrix.columns]

        # Get labels
        var_name = self.var_unit_name
        quantity = self.quantity

        if comparison_table is None:
            msg = "`comparison_table` was not given and is not set in analysis object."
            hint = "Will skip certain plots done at comparison level. "
            hint += "Add a `comparison_table` attribute to the analysis object."
            try:
                comparison_table = self.comparison_table
            except AttributeError:
                _LOGGER.warning(msg)
                _LOGGER.info(hint)

        if only_comparison_samples and comparison_table is not None:
            # select only comparisons from results table
            comparison_table = comparison_table.loc[
                comparison_table["comparison_name"].isin(
                    results["comparison_name"].unique().tolist()
                )
            ]
            samples = [s for s in samples if s.name in comparison_table["sample_name"].tolist()]
        matrix = matrix[[s.name for s in samples]]

        # Handle group colouring
        if group_colours:
            if group_attributes is None:

                try:
                    group_attributes = self.group_attributes
                except AttributeError:
                    msg = "`group_colours` is True, and `group_attributes` was not given and cannot be get from analysis!"
                    raise AssertionError(msg)

            # This will always be a matrix for all samples

            # Get matrix of samples vs levels with colors as values
            cd_kwargs = {
                k: v
                for k, v in kwargs.items()
                if k in ["pallete", "uniform_cmap", "diverging_cmap"]
            }

            color_dataframe = pd.DataFrame(
                self.get_level_colors(index=matrix.columns, levels=group_attributes, **cd_kwargs),
                index=group_attributes,
                columns=matrix.columns,
            ).T
            # will be filtered/ordered now by the requested samples if needed
            color_dataframe = color_dataframe.loc[matrix.columns, :]
            color_dataframe = color_dataframe.loc[[s.name for s in samples], :]
            color_dataframe.index = color_dataframe.index.get_level_values("sample_name")

        # Extract significant based on p-value and fold-change
        if fold_change is not None:
            fc = results[log_fold_change_column].abs() > fold_change
        else:
            fc = [True] * results.shape[0]
        if corrected_p_value:
            p_var = adjusted_p_value_column
        else:
            p_var = p_value_column
        results.loc[(results[p_var] < alpha) & fc, "diff"] = True
        results.loc[:, "diff"] = results.loc[:, "diff"].fillna(False)
        # Declare significant based on top ranked features
        if diff_based_on_rank:
            for comparison in results[comparison_column].unique():
                if ranking_variable == log_fold_change_column:
                    i = (
                        results.loc[results[comparison_column] == comparison, ranking_variable]
                        .abs()
                        .sort_values()
                        .tail(max_rank)
                        .index
                    )
                else:
                    i = (
                        results.loc[results[comparison_column] == comparison, ranking_variable]
                        .sort_values()
                        .head(max_rank)
                        .index
                    )
                results.loc[
                    (results[comparison_column] == comparison) & results.index.isin(i), "diff_rank",
                ] = True
            results.loc[:, "diff_rank"] = results.loc[:, "diff_rank"].fillna(False)
            if respect_stat_thresholds:
                results.loc[:, "diff"] = (results.loc[:, "diff"].isin([True])) & (
                    results.loc[:, "diff_rank"].isin([True])
                )
            else:
                results.loc[:, "diff"] = results.loc[:, "diff_rank"]

        # Annotate direction of change
        results.loc[:, "direction"] = results.loc[:, log_fold_change_column].apply(
            lambda x: "up" if x >= 0 else "down"
        )

        # PLOTS
        _LOGGER.info("Starting to generate plots for differential comparisons.")
        comparisons = sorted(results[comparison_column].drop_duplicates())
        n_side = int(np.ceil(np.sqrt(len(comparisons))))

        # P-value and Fold-change distributions
        if "distributions" in steps:
            for variable, label, axvline in [
                (p_value_column, "P-value", False),
                (adjusted_p_value_column, "Adjusted p-value", False),
                (log_fold_change_column, "log2(fold-change)", True),
            ]:
                _LOGGER.info("Plotting distribution of {}.".format(label))

                # log fold-changes distributions
                fig, axis = plt.subplots(1, 1, figsize=(4, 4))
                sns.distplot(results[variable].dropna(), kde=False, hist=True, ax=axis)
                if axvline:
                    axis.axvline(0, color="black", alpha=0.5)
                axis.set_xlabel(label)
                axis.set_ylabel(var_name.capitalize() + "s (frequency)")
                sns.despine(fig)
                savefig(
                    fig,
                    os.path.join(output_dir, output_prefix + "." + variable + ".distribution.svg"),
                )

                if plot_each_comparison:
                    # per comparison
                    g = sns.FacetGrid(
                        data=results, col=comparison_column, col_wrap=n_side, height=3, aspect=1,
                    )
                    g.map(sns.distplot, variable, kde=False, hist=True)
                    for ax in g.axes:
                        ax.set_yscale("log")
                        if axvline:
                            ax.axvline(0, color="black", alpha=0.5)
                        ax.set_xlabel(label)
                        ax.set_ylabel(var_name.capitalize() + "s (frequency)")
                    sns.despine(g.fig)
                    savefig(
                        g.fig,
                        os.path.join(
                            output_dir,
                            output_prefix + "." + variable + ".distribution.per_comparison.svg",
                        ),
                    )

        # Number of differential vars
        if "counts" in steps:
            _LOGGER.info("Calculating number of differential {}s per comparison.".format(var_name))
            n_vars = float(matrix.shape[0])
            total_diff = (
                results.groupby([comparison_column])["diff"]
                .sum()
                .sort_values(ascending=False)
                .reset_index()
            )
            split_diff = (
                results.groupby([comparison_column, "direction"])["diff"]
                .sum()
                .sort_values(ascending=False)
                .reset_index()
            )
            split_diff.loc[split_diff["direction"] == "down", "diff"] *= -1
            split_diff["label"] = (
                split_diff[comparison_column].astype(str) + ", " + split_diff["direction"]
            )
            total_diff["diff_perc"] = (total_diff["diff"] / n_vars) * 100
            split_diff["diff_perc"] = (split_diff["diff"] / n_vars) * 100

            _LOGGER.info("Plotting number of differential {}s per comparison.".format(var_name))
            fig, axis = plt.subplots(2, 2, figsize=(4 * 2, 4 * 2))
            sns.barplot(
                data=total_diff, x="diff", y=comparison_column, orient="h", ax=axis[0, 0],
            )
            sns.barplot(
                data=total_diff, x="diff_perc", y=comparison_column, orient="h", ax=axis[0, 1],
            )
            sns.barplot(
                data=split_diff,
                x="diff",
                y=comparison_column,
                hue="direction",
                dodge=False,
                orient="h",
                ax=axis[1, 0],
            )
            sns.barplot(
                data=split_diff,
                x="diff_perc",
                y=comparison_column,
                hue="direction",
                dodge=False,
                orient="h",
                ax=axis[1, 1],
            )
            for ax in axis[0, :]:
                ax.set_xlabel("", visible=False)
            for ax in axis[:, 1]:
                ax.set_yticklabels(ax.get_yticklabels(), visible=False)
            axis[-1, 0].set_xlabel("Frequency of differential {}s".format(var_name))
            axis[-1, 1].set_xlabel("Frequency of differential {}s (% of total)".format(var_name))
            for ax in axis[1, :].flatten():
                ax.axvline(0, linestyle="--", color="black", alpha=0.6)
            m = split_diff["diff"].abs().max()
            axis[1, 0].set_xlim((-m, m))
            m = split_diff["diff_perc"].abs().max()
            axis[1, 1].set_xlim((-m, m))
            sns.despine(fig)
            savefig(
                fig,
                os.path.join(output_dir, output_prefix + ".number_differential.directional.svg"),
            )

        if plot_each_comparison:
            _LOGGER.debug("Doing detailed plotting per comparison:")

            # Add same colour scale to all plots/comparisons
            smallest_p_value = -np.log10(np.nanpercentile(results[p_value_column], 1e-5))
            if smallest_p_value in [np.inf, np.nan]:
                smallest_p_value = 300
            _LOGGER.debug(
                "Maximum -log10(p-value) across comparisons is {}".format(smallest_p_value)
            )
            pval_cmap = "Reds"

            # Pairwise scatter plots
            if (comparison_table is not None) and ("scatter" in steps):
                _LOGGER.info(
                    "Plotting scatter of {} distribution for each group in each comparison.".format(
                        var_name
                    )
                )
                fig, axes = plt.subplots(
                    n_side, n_side, figsize=(n_side * 4, n_side * 4), sharex=True, sharey=True,
                )
                if n_side > 1 or n_side > 1:
                    axes = iter(axes.flatten())
                else:
                    axes = iter([axes])
                for comparison in comparisons:
                    _LOGGER.debug("Comparison '{}'...".format(comparison))
                    c = comparison_table.loc[comparison_table[comparison_column] == comparison, :]
                    a = c.loc[c["comparison_side"] >= 1, "sample_name"]
                    b = c.loc[c["comparison_side"] <= 0, "sample_name"]

                    a = matrix.loc[:, [s.name for s in samples if s.name in a.tolist()],].mean(
                        axis=1
                    )
                    b = matrix.loc[:, [s.name for s in samples if s.name in b.tolist()],].mean(
                        axis=1
                    )

                    # Hexbin plot
                    ax = next(axes)
                    try:
                        ax.hexbin(
                            b,
                            a,
                            alpha=0.85,
                            cmap="Greys",
                            color="black",
                            edgecolors="white",
                            linewidths=0,
                            bins="log",
                            mincnt=1,
                            rasterized=True,
                        )
                    except ValueError:
                        _LOGGER.warning(
                            "Couldn't plot scatter for comparison '{}'.".format(comparison)
                        )
                        continue

                    # Scatter for significant features
                    diff_vars = results.loc[
                        (results[comparison_column] == comparison) & (results["diff"].isin([True])),
                        :,
                    ]
                    if diff_vars.shape[0] > 0:
                        # get color vector based on p-value
                        col = -np.log10(
                            results.loc[
                                (results[comparison_column] == comparison)
                                & (results["diff"].isin([True])),
                                p_value_column,
                            ].squeeze()
                        )
                        _LOGGER.debug("Shapes: {} {} {}".format(a.shape, b.shape, diff_vars.shape))
                        # in case there's just one significant feature:
                        if isinstance(col, np.float_):
                            col = np.array([col])
                        collection = ax.scatter(
                            b.loc[diff_vars.index],
                            a.loc[diff_vars.index],
                            alpha=0.5,
                            s=2,
                            c=col,
                            cmap=pval_cmap,
                            vmin=0,
                            vmax=smallest_p_value,
                        )
                        add_colorbar_to_axis(collection, label="-log10(p-value)")
                    ax.set_title(comparison)
                    # Name groups
                    xl = (
                        c.loc[c["comparison_side"] <= 0, "sample_group"].drop_duplicates().squeeze()
                    )
                    yl = (
                        c.loc[c["comparison_side"] >= 1, "sample_group"].drop_duplicates().squeeze()
                    )
                    if not (isinstance(xl, str) and isinstance(yl, str)):
                        xl = "Down-regulated"
                        yl = "Up-regulated"
                    ax.set_xlabel(xl)
                    ax.set_ylabel(yl)

                    # x = y square
                    lims = [
                        np.min([ax.get_xlim(), ax.get_ylim()]),
                        np.max([ax.get_xlim(), ax.get_ylim()]),
                    ]
                    ax.plot(lims, lims, linestyle="--", alpha=0.5, zorder=0, color="black")
                    ax.set_aspect("equal")
                    ax.set_xlim(lims)
                    ax.set_ylim(lims)
                for ax in axes:
                    ax.set_visible(False)
                sns.despine(fig)
                savefig(fig, os.path.join(output_dir, output_prefix + ".scatter_plots.svg"))

            # Volcano plots
            if "volcano" in steps:
                _LOGGER.info("Plotting volcano plots for each comparison.")
                fig, axes = plt.subplots(
                    n_side, n_side, figsize=(n_side * 4, n_side * 4), sharex=False, sharey=False,
                )
                if n_side > 1 or n_side > 1:
                    axes = iter(axes.flatten())
                else:
                    axes = iter([axes])
                for comparison in comparisons:
                    _LOGGER.debug("Comparison '{}'...".format(comparison))
                    t = results.loc[results[comparison_column] == comparison, :]

                    # Hexbin plot
                    ax = next(axes)
                    ax.hexbin(
                        t[log_fold_change_column],
                        -np.log10(t[p_value_column]),
                        alpha=0.85,
                        cmap="Greys",
                        color="black",
                        edgecolors="white",
                        linewidths=0,
                        bins="log",
                        mincnt=1,
                        rasterized=True,
                    )

                    # Scatter for significant
                    diff_vars = t.loc[t["diff"].isin([True]), :]
                    if diff_vars.shape[0] > 0:
                        collection = ax.scatter(
                            t.loc[diff_vars.index, log_fold_change_column],
                            -np.log10(t.loc[diff_vars.index, p_value_column]),
                            alpha=0.5,
                            s=2,
                            c=-np.log10(t.loc[diff_vars.index, p_value_column]),
                            cmap=pval_cmap,
                            vmin=0,
                            vmax=smallest_p_value,
                        )
                        add_colorbar_to_axis(collection, label="-log10(p-value)")
                    ax.set_title(comparison)
                    ax.set_xlabel("log2(fold-change)")
                    ax.set_ylabel("-log10(p-value)")
                    ax.axvline(0, linestyle="--", alpha=0.5, zorder=0, color="black")
                    ll = np.max([abs(ii) for ii in ax.get_xlim()])
                    ax.set_xlim(-ll, ll)

                    # Add lines of significance
                    ax.axhline(
                        -np.log10(t.loc[t["diff"].isin([True]), p_value_column].max()),
                        linestyle="--",
                        alpha=0.5,
                        zorder=0,
                        color="black",
                    )
                    if fold_change is not None:
                        ax.axvline(
                            -fold_change, linestyle="--", alpha=0.5, zorder=0, color="black",
                        )
                        ax.axvline(
                            fold_change, linestyle="--", alpha=0.5, zorder=0, color="black",
                        )
                for ax in axes:
                    ax.set_visible(False)
                sns.despine(fig)
                savefig(fig, os.path.join(output_dir, output_prefix + ".volcano_plots.svg"))

            # MA plots
            if "ma" in steps:
                _LOGGER.info("Plotting MA plots for each comparison.")
                fig, axes = plt.subplots(
                    n_side, n_side, figsize=(n_side * 4, n_side * 4), sharex=False, sharey=False,
                )
                if n_side > 1 or n_side > 1:
                    axes = iter(axes.flatten())
                else:
                    axes = iter([axes])
                for comparison in comparisons:
                    _LOGGER.debug("Comparison '{}'...".format(comparison))
                    t = results.loc[results[comparison_column] == comparison, :]

                    # Hexbin plot
                    ax = next(axes)
                    ax.hexbin(
                        np.log10(t[mean_column]),
                        t[log_fold_change_column],
                        alpha=0.85,
                        cmap="Greys",
                        color="black",
                        edgecolors="white",
                        linewidths=0,
                        bins="log",
                        mincnt=1,
                        rasterized=True,
                    )

                    # Scatter for significant
                    diff_vars = t.loc[t["diff"].isin([True]), :]
                    if diff_vars.shape[0] > 0:
                        collection = ax.scatter(
                            np.log10(t.loc[diff_vars.index, mean_column]),
                            t.loc[diff_vars.index, log_fold_change_column],
                            alpha=0.5,
                            s=2,
                            c=-np.log10(t.loc[diff_vars.index, p_value_column]),
                            cmap=pval_cmap,
                            vmin=0,
                            vmax=smallest_p_value,
                        )
                        add_colorbar_to_axis(collection, label="-log10(p-value)")
                    ax.set_title(comparison)
                    ax.set_xlabel("Mean {}".format(quantity.lower()))
                    ax.set_ylabel("log2(fold-change)")
                    ax.axhline(0, linestyle="--", alpha=0.5, zorder=0, color="black")
                    ll = np.max([abs(ii) for ii in ax.get_ylim()])
                    ax.set_ylim(-ll, ll)

                    # Add lines of significance
                    if fold_change is not None:
                        ax.axhline(
                            -fold_change, linestyle="--", alpha=0.5, zorder=0, color="black",
                        )
                        ax.axhline(
                            fold_change, linestyle="--", alpha=0.5, zorder=0, color="black",
                        )
                for ax in axes:
                    ax.set_visible(False)
                sns.despine(fig)
                savefig(fig, os.path.join(output_dir, output_prefix + ".ma_plots.svg"))

        if results.loc[:, "diff"].sum() < 1:
            msg = "No significantly different regions found in any comparison."
            msg += " Skipping heatmap plots on differential {}s.".format(var_name)
            _LOGGER.warning(msg)
            return

        # Observe values of variables across all comparisons
        all_diff = results[results["diff"].isin([True])].index.drop_duplicates()
        if isinstance(matrix.columns, pd.MultiIndex):
            sample_cols = matrix.columns.get_level_values("sample_name").tolist()
        else:
            sample_cols = matrix.columns.tolist()

        if (comparison_table is not None) and ("heatmap" in steps):
            _LOGGER.info("A comparison table was given, will try to plot values per sample group.")
            if results[comparison_column].drop_duplicates().shape[0] > 1:
                _LOGGER.info("Getting per-group values for each comparison.")
                groups = pd.DataFrame()
                for sample_group in comparison_table["sample_group"].drop_duplicates():
                    c = comparison_table.loc[
                        comparison_table["sample_group"] == sample_group, "sample_name"
                    ].drop_duplicates()
                    if c.shape[0] > 0:
                        groups.loc[:, sample_group] = matrix[
                            [d for d in c if d in sample_cols]
                        ].mean(axis=1)

                if groups.empty:
                    # It seems comparisons were not done in a all-versus-all fashion
                    for group in comparison_table["sample_group"].drop_duplicates():
                        c = comparison_table.loc[
                            comparison_table["sample_group"] == group, "sample_name"
                        ].drop_duplicates()
                        if c.shape[0] > 0:
                            groups.loc[:, group] = matrix[c].mean(axis=1)

                # Select only differential regions from groups
                groups = groups.loc[all_diff, :].sort_index(axis=1)

                n = groups.isnull().sum().sum()
                if n > 0:
                    _LOGGER.warning(
                        "{} {}s (across all comparisons) were not found in quantification matrix!".format(
                            n, var_name
                        )
                    )
                    m = groups.columns[groups.isnull().sum() == groups.shape[0]]
                    if len(m) > 0:
                        _LOGGER.warning(
                            "{} comparison groups were not found in quantification matrix: '{}'!".format(
                                len(m), ", ".join(m)
                            )
                            + " Proceeding without those."
                        )
                        groups = groups.loc[:, ~groups.columns.isin(m)]
                f = groups.index[groups.isnull().sum(1) == groups.shape[1]]
                if len(f) > 0:
                    _LOGGER.warning(
                        "{} {}s were not found in quantification matrix!".format(len(m), var_name)
                        + " Proceeding without those."
                    )
                    groups = groups.dropna()
                n = groups.isnull().sum().sum()
                if n != 0:
                    _LOGGER.error(
                        "{} {}s (across all comparisons) still have NaNs. Cannot proceed!".format(
                            n, var_name
                        )
                    )
                else:
                    _LOGGER.info(
                        "Plotting clustered heatmaps of sample groups in all differential {}s found.".format(
                            var_name
                        )
                    )
                    figsize = (max(5, 0.12 * groups.shape[1]), 5)
                    # Heatmaps
                    # Comparison level
                    g = sns.clustermap(
                        groups.corr(),
                        xticklabels=False,
                        yticklabels=True,
                        cbar_kws={
                            "label": "Pearson correlation\non differential {}s".format(var_name)
                        },
                        metric="correlation",
                        rasterized=True,
                        figsize=(figsize[0], figsize[0]),
                    )
                    g.ax_heatmap.set_yticklabels(
                        g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small",
                    )
                    savefig(
                        g.fig,
                        os.path.join(
                            output_dir,
                            output_prefix + ".diff_{}.groups.clustermap.corr.svg".format(var_name),
                        ),
                    )

                    g = sns.clustermap(
                        groups,
                        xticklabels=True,
                        yticklabels=feature_labels,
                        cbar_kws={"label": "{} of\ndifferential {}s".format(quantity, var_name)},
                        robust=robust,
                        metric="correlation",
                        rasterized=True,
                        figsize=figsize,
                    )
                    g.ax_heatmap.set_ylabel(
                        "Differential {}s (n = {})".format(var_name, groups.shape[0])
                    )
                    g.ax_heatmap.set_xticklabels(
                        g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small",
                    )
                    g.ax_heatmap.set_yticklabels(
                        g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small",
                    )
                    savefig(
                        g.fig,
                        os.path.join(
                            output_dir,
                            output_prefix + ".diff_{}.groups.clustermap.svg".format(var_name),
                        ),
                    )

                    g = sns.clustermap(
                        groups,
                        xticklabels=True,
                        yticklabels=feature_labels,
                        z_score=0,
                        cbar_kws={
                            "label": "Z-score of {}\non differential {}s".format(quantity, var_name)
                        },
                        cmap="RdBu_r",
                        center=0,
                        robust=robust,
                        metric="correlation",
                        rasterized=True,
                        figsize=figsize,
                    )
                    g.ax_heatmap.set_ylabel(
                        "Differential {}s (n = {})".format(var_name, groups.shape[0])
                    )
                    g.ax_heatmap.set_xticklabels(
                        g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small",
                    )
                    g.ax_heatmap.set_yticklabels(
                        g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small",
                    )
                    savefig(
                        g.fig,
                        os.path.join(
                            output_dir,
                            output_prefix + ".diff_{}.groups.clustermap.z0.svg".format(var_name),
                        ),
                    )

                    # same without clustering
                    g = sns.clustermap(
                        groups,
                        col_cluster=False,
                        xticklabels=True,
                        yticklabels=feature_labels,
                        cbar_kws={"label": "{} of\ndifferential {}s".format(quantity, var_name)},
                        robust=robust,
                        metric="correlation",
                        rasterized=True,
                        figsize=figsize,
                    )
                    g.ax_heatmap.set_ylabel(
                        "Differential {}s (n = {})".format(var_name, groups.shape[0])
                    )
                    g.ax_heatmap.set_xticklabels(
                        g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small",
                    )
                    g.ax_heatmap.set_yticklabels(
                        g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small",
                    )
                    savefig(
                        g.fig,
                        os.path.join(
                            output_dir,
                            output_prefix
                            + ".diff_{}.groups.sorted.clustermap.svg".format(var_name),
                        ),
                    )

                    g = sns.clustermap(
                        groups,
                        col_cluster=False,
                        xticklabels=True,
                        yticklabels=feature_labels,
                        z_score=0,
                        cbar_kws={
                            "label": "Z-score of {}\non differential {}s".format(quantity, var_name)
                        },
                        cmap="RdBu_r",
                        center=0,
                        robust=robust,
                        metric="correlation",
                        rasterized=True,
                        figsize=figsize,
                    )
                    g.ax_heatmap.set_ylabel(
                        "Differential {}s (n = {})".format(var_name, groups.shape[0])
                    )
                    g.ax_heatmap.set_xticklabels(
                        g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small",
                    )
                    g.ax_heatmap.set_yticklabels(
                        g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small",
                    )
                    savefig(
                        g.fig,
                        os.path.join(
                            output_dir,
                            output_prefix
                            + ".diff_{}.groups.sorted.clustermap.z0.svg".format(var_name),
                        ),
                    )

        # Fold-changes and P-values
        # pivot table of genes vs comparisons
        if "stats_heatmap" in steps:
            _LOGGER.info("Getting fold-change and p-value values per comparison.")
            if results.index.name is None:
                results.index.name = "index"
            fold_changes = pd.pivot_table(
                results.loc[all_diff, :].reset_index(),
                index=results.index.name,
                columns=comparison_column,
                values=log_fold_change_column,
            ).fillna(0)
            p_values = (
                -np.log10(
                    pd.pivot_table(
                        results.loc[all_diff, :].reset_index(),
                        index=results.index.name,
                        columns=comparison_column,
                        values=adjusted_p_value_column,
                    )
                )
            ).fillna(0)

            # get a signed p-value
            if fold_changes.shape == p_values.shape:
                p_values *= (fold_changes > 0).astype(int).replace(0, -1)

            for matrix_, label, desc in [
                (fold_changes, "log fold change", "log fold change"),
                (p_values, "p value", "-log10(signed p-value)"),
            ]:

                if matrix_.isnull().sum().sum() > 0:
                    _LOGGER.warning(
                        "Some sample groups or regions in {} matrix contain NaNs. Removing those.".format(
                            label
                        )
                    )
                    matrix_ = matrix_.dropna().T.dropna().T

                if matrix_.shape[1] > 1:
                    _LOGGER.info(
                        "Plotting group-wise correlation of {}s per sample groups in all differential {}s found.".format(
                            var_name, label
                        )
                    )
                    figsize = (max(5, 0.12 * matrix_.shape[1]), 5)
                    grid = sns.clustermap(
                        matrix_.corr(),
                        xticklabels=False,
                        yticklabels=True,
                        cbar_kws={"label": "Pearson correlation\non {}s".format(desc)},
                        vmin=0,
                        vmax=1,
                        metric="correlation",
                        rasterized=True,
                        figsize=(figsize[0], figsize[0]),
                    )
                    grid.ax_heatmap.set_yticklabels(
                        grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small",
                    )
                    grid.ax_heatmap.set_xlabel("Comparison groups")
                    grid.ax_heatmap.set_ylabel("Comparison groups")
                    savefig(
                        grid.fig,
                        os.path.join(
                            output_dir,
                            output_prefix
                            + ".diff_{}.groups.{}.clustermap.corr.svg".format(
                                var_name, label.replace(" ", "_")
                            ),
                        ),
                    )
                    _LOGGER.info(
                        "Plotting clustered heatmaps of {}s per sample groups in all differential {}s found.".format(
                            var_name, label
                        )
                    )
                    try:
                        grid = sns.clustermap(
                            matrix_.loc[all_diff, :],
                            xticklabels=True,
                            yticklabels=feature_labels,
                            cbar_kws={"label": "{} of\ndifferential {}s".format(desc, var_name)},
                            cmap="RdBu_r",
                            center=0,
                            robust=robust,
                            metric="correlation",
                            rasterized=True,
                            figsize=figsize,
                        )
                        grid.ax_heatmap.set_ylabel(
                            "Differential {}s (n = {})".format(
                                var_name, matrix_.loc[all_diff, :].shape[0]
                            )
                        )
                        grid.ax_heatmap.set_xticklabels(
                            grid.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small",
                        )
                        grid.ax_heatmap.set_yticklabels(
                            grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small",
                        )
                        grid.ax_heatmap.set_xlabel("Comparison groups")
                        savefig(
                            grid.fig,
                            os.path.join(
                                output_dir,
                                output_prefix
                                + ".diff_{}.groups.{}.clustermap.svg".format(
                                    var_name, label.replace(" ", "_")
                                ),
                            ),
                        )
                    except FloatingPointError:
                        _LOGGER.error(
                            "{} likely contains null or infinite values. Cannot plot.".format(label)
                        )

        # Sample level
        if "heatmap" in steps:
            _LOGGER.info(
                "Getting per sample values of {} in all differential {}s found.".format(
                    quantity, var_name
                )
            )
            if isinstance(matrix.columns, pd.core.indexes.multi.MultiIndex):
                matrix.columns = matrix.columns.get_level_values("sample_name")

            matrix2 = matrix.loc[all_diff, :].sort_index(axis=1)

            n = matrix2.isnull().sum().sum()
            if n > 0:
                _LOGGER.warning(
                    "WARNING! {} {} (across all comparisons) were not found in quantification matrix!".format(
                        n, var_name
                    )
                    + " Proceeding without those."
                )
                matrix2 = matrix2.dropna()
            figsize = (max(5, 0.12 * matrix2.shape[1]), 5)
            if group_colours:
                extra = {"col_colors": color_dataframe}
            else:
                extra = {}

            _LOGGER.info(
                "Plotting sample-wise correlation heatmaps of {} values per sample in all differential {}s found.".format(
                    quantity, var_name
                )
            )
            grid = sns.clustermap(
                matrix2.corr(),
                yticklabels=True,
                xticklabels=False,
                cbar_kws={"label": "Pearson correlation\non differential {}s".format(var_name)},
                metric="correlation",
                figsize=(figsize[0], figsize[0]),
                rasterized=rasterized,
                robust=robust,
                **extra
            )
            grid.ax_heatmap.set_yticklabels(
                grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small"
            )
            savefig(
                grid.fig,
                os.path.join(
                    output_dir,
                    output_prefix + ".diff_{}.samples.clustermap.corr.svg".format(var_name),
                ),
            )

            _LOGGER.info(
                "Plotting clustered heatmaps of {} values per sample in all differential {}s found.".format(
                    quantity, var_name
                )
            )
            grid = sns.clustermap(
                matrix2,
                yticklabels=feature_labels,
                cbar_kws={"label": "{} of\ndifferential {}s".format(quantity, var_name)},
                xticklabels=True,
                vmin=0,
                metric="correlation",
                figsize=figsize,
                rasterized=rasterized,
                robust=robust,
                **extra
            )
            grid.ax_heatmap.set_ylabel(
                "Differential {}s (n = {})".format(var_name, matrix2.shape[0])
            )
            grid.ax_heatmap.set_xticklabels(
                grid.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small"
            )
            grid.ax_heatmap.set_yticklabels(
                grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small"
            )
            savefig(
                grid.fig,
                os.path.join(
                    output_dir, output_prefix + ".diff_{}.samples.clustermap.svg".format(var_name),
                ),
            )

            grid = sns.clustermap(
                matrix2,
                yticklabels=feature_labels,
                z_score=0,
                cbar_kws={"label": "Z-score of {}\non differential {}s".format(quantity, var_name)},
                xticklabels=True,
                cmap="RdBu_r",
                center=0,
                metric="correlation",
                figsize=figsize,
                rasterized=rasterized,
                robust=robust,
                **extra
            )
            grid.ax_heatmap.set_ylabel(
                "Differential {}s (n = {})".format(var_name, matrix2.shape[0])
            )
            grid.ax_heatmap.set_xticklabels(
                grid.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small"
            )
            grid.ax_heatmap.set_yticklabels(
                grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small"
            )
            savefig(
                grid.fig,
                os.path.join(
                    output_dir,
                    output_prefix + ".diff_{}.samples.clustermap.z0.svg".format(var_name),
                ),
            )

            grid = sns.clustermap(
                matrix2,
                col_cluster=False,
                yticklabels=feature_labels,
                cbar_kws={"label": "{} of\ndifferential {}s".format(quantity, var_name)},
                xticklabels=True,
                vmin=0,
                metric="correlation",
                figsize=figsize,
                rasterized=rasterized,
                robust=robust,
                **extra
            )
            grid.ax_heatmap.set_ylabel(
                "Differential {}s (n = {})".format(var_name, matrix2.shape[0])
            )
            grid.ax_heatmap.set_xticklabels(
                grid.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small"
            )
            grid.ax_heatmap.set_yticklabels(
                grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small"
            )
            savefig(
                grid.fig,
                os.path.join(
                    output_dir,
                    output_prefix + ".diff_{}.samples.sorted.clustermap.svg".format(var_name),
                ),
            )

            grid = sns.clustermap(
                matrix2,
                col_cluster=False,
                yticklabels=feature_labels,
                z_score=0,
                cbar_kws={"label": "Z-score of {}\non differential {}s".format(quantity, var_name)},
                xticklabels=True,
                cmap="RdBu_r",
                center=0,
                metric="correlation",
                figsize=figsize,
                rasterized=rasterized,
                robust=robust,
                **extra
            )
            grid.ax_heatmap.set_ylabel(
                "Differential {}s (n = {})".format(var_name, matrix2.shape[0])
            )
            grid.ax_heatmap.set_xticklabels(
                grid.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small"
            )
            grid.ax_heatmap.set_yticklabels(
                grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small"
            )
            savefig(
                grid.fig,
                os.path.join(
                    output_dir,
                    output_prefix + ".diff_{}.samples.sorted.clustermap.z0.svg".format(var_name),
                ),
            )

    def differential_overlap(
        self,
        differential=None,
        output_dir="{results_dir}/differential_analysis_{data_type}",
        output_prefix="differential_analysis",
    ):
        """
        Visualize intersection of sets of differential regions/genes.

        Parameters
        ----------
        differential : :obj:`pandas.DataFrame`, optional
            DataFrame containing result of comparisons filtered for features considered as differential.

            Defaults to the ``differential_results`` attribute, subset by the object's ``thresholds``.
        output_dir : :obj:`str`, optional
            Directory to create output files.

            Defaults to "{results_dir}/differential_analysis_{data_type}".
        output_prefix : :obj:`str`, optional
            Prefix to use when creating output files.

            Defaults to "differential_analysis".
        """
        # Make output dir
        import itertools
        import matplotlib
        import matplotlib.pyplot as plt
        from ngs_toolkit.graphics import savefig
        from ngs_toolkit.utils import log_pvalues
        from scipy.stats import fisher_exact
        import seaborn as sns
        from statsmodels.sandbox.stats.multicomp import multipletests
        from tqdm import tqdm

        output_dir = self._format_string_with_attributes(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        total = self.matrix_raw.shape[0]
        unit = self.var_unit_name

        if differential is None:
            differential = self.differential_results.loc[
                (self.differential_results["padj"] < self.thresholds["alpha"])
            ]

        if "direction" not in differential.columns:
            differential.loc[:, "direction"] = differential["log2FoldChange"].apply(
                lambda x: "up" if x > 0 else "down"
            )

        differential.index.name = "index"
        differential.loc[:, "intersect"] = 1
        piv = pd.pivot_table(
            differential.reset_index(),
            index="index",
            columns=["comparison_name", "direction"],
            values="intersect",
            fill_value=0,
        )

        intersections = pd.DataFrame(
            columns=["group1", "group2", "dir1", "dir2", "size1", "size2", "intersection", "union",]
        )
        perms = list(
            itertools.permutations(
                piv.T.groupby(level=["comparison_name", "direction"]).groups.items(), 2
            )
        )
        for ((k1, dir1), i1), ((k2, dir2), i2) in tqdm(
            perms, total=len(perms), desc="Permutations"
        ):
            i1 = set(piv[i1][piv[i1] == 1].dropna().index)
            i2 = set(piv[i2][piv[i2] == 1].dropna().index)
            intersections = intersections.append(
                pd.Series(
                    [
                        k1,
                        k2,
                        dir1,
                        dir2,
                        len(i1),
                        len(i2),
                        len(i1.intersection(i2)),
                        len(i1.union(i2)),
                    ],
                    index=[
                        "group1",
                        "group2",
                        "dir1",
                        "dir2",
                        "size1",
                        "size2",
                        "intersection",
                        "union",
                    ],
                ),
                ignore_index=True,
            )
        # convert to %
        intersections.loc[:, "intersection"] = intersections["intersection"].astype(float)
        intersections.loc[:, "perc_1"] = (
            intersections["intersection"] / intersections["size1"] * 100.0
        )
        intersections.loc[:, "perc_2"] = (
            intersections["intersection"] / intersections["size2"] * 100.0
        )
        intersections.loc[:, "intersection_max_perc"] = intersections[["perc_1", "perc_2"]].max(
            axis=1
        )

        # calculate p-value from Fisher"s exact test
        intersections.loc[:, "a"] = intersections["intersection"]
        intersections.loc[:, "b"] = intersections["size1"] - intersections["intersection"]
        intersections.loc[:, "c"] = intersections["size2"] - intersections["intersection"]
        intersections.loc[:, "d"] = total - intersections[["b", "c", "intersection"]].sum(axis=1)

        for i, row in intersections.loc[:, ["a", "b", "c", "d"]].astype(int).iterrows():
            odds, p = fisher_exact(row.values.reshape((2, 2)), alternative="greater")
            intersections.loc[i, "odds_ratio"] = odds
            intersections.loc[i, "p_value"] = p
        # intersections["q_value"] = intersections["p_value"] * intersections.shape[0]
        intersections.loc[:, "q_value"] = multipletests(intersections["p_value"])[1]
        intersections.loc[:, "log_p_value"] = log_pvalues(intersections["p_value"])
        intersections.loc[:, "log_q_value"] = log_pvalues(intersections["q_value"])

        # save
        intersections.to_csv(
            os.path.join(output_dir, output_prefix + ".differential_overlap.csv"), index=False,
        )
        intersections = pd.read_csv(
            os.path.join(output_dir, output_prefix + ".differential_overlap.csv")
        )

        for metric, label, description, fill_value in [
            ("intersection", "intersection", "total in intersection", 0),
            ("intersection_max_perc", "percentage_overlap", "max of intersection %", 0),
            ("log_p_value", "significance", "p-value", 0),
        ]:
            _LOGGER.debug(metric)
            # make pivot tables
            piv_up = pd.pivot_table(
                intersections[(intersections["dir1"] == "up") & (intersections["dir2"] == "up")],
                index="group1",
                columns="group2",
                values=metric,
            ).fillna(fill_value)
            piv_down = pd.pivot_table(
                intersections[
                    (intersections["dir1"] == "down") & (intersections["dir2"] == "down")
                ],
                index="group1",
                columns="group2",
                values=metric,
            ).fillna(fill_value)
            if metric == "intersection":
                piv_up = np.log10(1 + piv_up)
                piv_down = np.log10(1 + piv_down)
            np.fill_diagonal(piv_up.values, np.nan)
            np.fill_diagonal(piv_down.values, np.nan)

            # heatmaps
            if metric == "intersection_max_perc":
                extra = {"vmin": 0, "vmax": 100}
            else:
                extra = {}
            fig, axis = plt.subplots(1, 2, figsize=(8 * 2, 8), subplot_kw={"aspect": "equal"})
            sns.heatmap(
                piv_down,
                square=True,
                cmap="Blues",
                cbar_kws={"label": "Concordant {}s ({})".format(unit, description)},
                ax=axis[0],
                **extra
            )
            sns.heatmap(
                piv_up,
                square=True,
                cmap="Reds",
                cbar_kws={"label": "Concordant {}s ({})".format(unit, description)},
                ax=axis[1],
                **extra
            )
            axis[0].set_title("Downregulated {}s".format(unit))
            axis[1].set_title("Upregulated {}s".format(unit))
            axis[0].set_xticklabels(axis[0].get_xticklabels(), rotation=90, ha="center")
            axis[1].set_xticklabels(axis[1].get_xticklabels(), rotation=90, ha="center")
            axis[0].set_yticklabels(axis[0].get_yticklabels(), rotation=0, ha="right")
            axis[1].set_yticklabels(axis[1].get_yticklabels(), rotation=0, ha="right")
            savefig(
                fig,
                os.path.join(
                    output_dir,
                    output_prefix + ".differential_overlap.{}.up_down_split.svg".format(label),
                ),
            )

            # combined heatmap
            # with upregulated {}s in upper square matrix and downredulated in down square
            piv_combined = pd.DataFrame(
                np.triu(piv_up), index=piv_up.index, columns=piv_up.columns
            ).replace(0, np.nan)
            piv_combined.update(
                pd.DataFrame(
                    np.tril(-piv_down), index=piv_down.index, columns=piv_down.columns
                ).replace(0, np.nan)
            )
            piv_combined = piv_combined.fillna(fill_value)
            if metric == "intersection":
                piv_combined = np.log10(1 + piv_combined)
            np.fill_diagonal(piv_combined.values, np.nan)

            if metric == "intersection_max_perc":
                extra = {"vmin": -150, "vmax": 150}
            else:
                extra = {}
            fig, axis = plt.subplots(1, figsize=(8, 8))
            sns.heatmap(
                piv_combined,
                square=True,
                cmap="RdBu_r",
                center=0,
                cbar_kws={"label": "Concordant {}s ({})".format(unit, description)},
                ax=axis,
                **extra
            )
            axis.set_xticklabels(axis.get_xticklabels(), rotation=90, ha="center")
            axis.set_yticklabels(axis.get_yticklabels(), rotation=0, ha="right")
            savefig(
                fig,
                os.path.join(
                    output_dir,
                    output_prefix + ".differential_overlap.{}.up_down_together.svg".format(label),
                ),
            )

            # Rank plots
            if metric == "log_pvalue":
                r = pd.melt(
                    piv_combined.reset_index(),
                    id_vars=["group1"],
                    var_name="group2",
                    value_name="agreement",
                )
                r = r.dropna().sort_values("agreement")
                r = r.iloc[range(0, r.shape[0], 2)]
                r["rank"] = r["agreement"].rank(ascending=False)

                fig, axis = plt.subplots(1, 3, figsize=(3 * 4, 4))
                axis[0].scatter(r["rank"], r["agreement"])
                axis[0].axhline(0, linestyle="--", color="black")
                axis[1].scatter(r["rank"].tail(10), r["agreement"].tail(10))
                axis[2].scatter(r["rank"].head(10), r["agreement"].head(10))
                for i, row in r.tail(10).iterrows():
                    axis[1].text(
                        row["rank"],
                        row["agreement"],
                        s=row["group1"] + "\n" + row["group2"],
                        fontsize=5,
                    )
                for i, row in r.head(10).iterrows():
                    axis[2].text(
                        row["rank"],
                        row["agreement"],
                        s=row["group1"] + "\n" + row["group2"],
                        fontsize=5,
                    )
                for ax in axis:
                    ax.set_ylabel("Agreement (-log(p-value))")
                    ax.set_xlabel("Rank")
                sns.despine(fig)
                savefig(
                    fig,
                    os.path.join(
                        output_dir,
                        output_prefix + ".differential_overlap.{}.agreement.rank.svg".format(label),
                    ),
                )

            # Observe disagreement
            # (overlap of down-regulated with up-regulated and vice-versa)
            piv_up = pd.pivot_table(
                intersections[(intersections["dir1"] == "up") & (intersections["dir2"] == "down")],
                index="group1",
                columns="group2",
                values=metric,
            )
            piv_down = pd.pivot_table(
                intersections[(intersections["dir1"] == "down") & (intersections["dir2"] == "up")],
                index="group1",
                columns="group2",
                values=metric,
            )

            piv_disagree = pd.concat([piv_up, piv_down]).groupby(level=0).max()
            if metric == "intersection":
                piv_disagree = np.log10(1 + piv_disagree)
            np.fill_diagonal(piv_disagree.values, np.nan)

            fig, axis = plt.subplots(1, 2, figsize=(16, 8), subplot_kw={"aspect": "equal"})
            sns.heatmap(
                piv_disagree,
                square=True,
                cmap="Greens",
                cbar_kws={"label": "Discordant {}s ({})".format(unit, description)},
                ax=axis[0],
            )

            norm = matplotlib.colors.Normalize(vmin=0, vmax=piv_disagree.max().max())
            cmap = plt.get_cmap("Greens")
            for j, g2 in enumerate(piv_disagree.index):
                for i, g1 in enumerate(piv_disagree.columns):
                    axis[1].scatter(
                        len(piv_disagree.index) - (j + 0.5),
                        len(piv_disagree.index) - (i + 0.5),
                        s=(100 ** (norm(piv_disagree.loc[g1, g2]))) - 1,
                        color=cmap(norm(piv_disagree.loc[g1, g2])),
                        marker="o",
                    )
                    axis[1].set_title("Rotate plot -90 degrees")
            axis[0].set_xticklabels(axis[0].get_xticklabels(), rotation=90, ha="center")
            axis[0].set_yticklabels(axis[0].get_yticklabels(), rotation=0, ha="right")
            axis[1].set_xlim((0, len(piv_disagree.index)))
            axis[1].set_ylim((0, len(piv_disagree.columns)))
            savefig(
                fig,
                os.path.join(
                    output_dir,
                    output_prefix + ".differential_overlap.{}.disagreement.svg".format(label),
                ),
            )

            # Rank plots
            if metric == "log_pvalue":
                r = pd.melt(
                    piv_disagree.reset_index(),
                    id_vars=["group1"],
                    var_name="group2",
                    value_name="disagreement",
                )
                r = r.dropna().sort_values("disagreement")
                r = r.iloc[range(0, r.shape[0], 2)]
                r["rank"] = r["disagreement"].rank(ascending=False)

                fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 4), subplot_kw={"aspect": "equal"})
                axis[0].scatter(r["rank"], r["disagreement"])
                axis[1].scatter(r["rank"].tail(10), r["disagreement"].tail(10))
                for i, row in r.tail(10).iterrows():
                    axis[1].text(
                        row["rank"],
                        row["disagreement"],
                        s=row["group1"] + "\n" + row["group2"],
                        fontsize=5,
                    )
                for ax in axis:
                    ax.set_ylabel("Disagreement (-log(p-value))")
                    ax.set_xlabel("Rank")
                sns.despine(fig)
                savefig(
                    fig,
                    os.path.join(
                        output_dir,
                        output_prefix
                        + ".differential_overlap.{}.disagreement.rank.svg".format(label),
                    ),
                )

    @check_has_attributes(["organism", "genome"])
    def differential_enrichment(
        self,
        differential=None,
        output_dir="{results_dir}/differential_analysis_{data_type}/enrichments",
        output_prefix="differential_analysis",
        genome=None,
        steps=["region", "lola", "meme", "homer", "enrichr"],
        directional=True,
        max_diff=1000,
        sort_var="pvalue",
        distributed=False,
        overwrite=False,
    ):
        """
        Perform various types of enrichment analysis given a dataframe of the results from differential analysis.
        Performs enrichment of gene sets (RNA-seq and ATAC-seq), genomic regions, chromatin states
        Location Overlap Analysis (LOLA) and TF motif enrichment (over-representation and de-novo search)
        (ATAC-seq only).

        Parameters
        ----------
        differential : :obj:`pandas.DataFrame`
            Data frame with differential results as produced by ``differential_analysis``,
            but filtered by some threshold for the relevant (significant regions).
            Must contain a "comparison_name" column.

            Defaults to ``analysis.differential_results``.
        output_dir : :obj:`str`, optional
            Directory to create output files.

            Defaults to "{results_dir}/differential_analysis_{data_type}".
        output_prefix : :obj:`str`, optional
            Prefix to use when creating output files.

            Defaults to "differential_analysis".
        genome : :obj:`str`, optional
            Genome assembly of the analysis.

            Defaults to Analysis's ``genome`` attribute.
        steps : :obj:`list`, optional
            Steps of the analysis to perform.

            Defaults to all possible: ["region", lola", "meme", "homer", "enrichr"].
        directional: :obj:`bool`, optional
            Whether enrichments should be performed in a direction-dependent way
            (up-regulated and down-regulated features separately).
            This requires a column named "log2FoldChange" to exist.

            Defaults to :obj:`True`.
        max_diff : :obj:`int`, optional
            Number of maximum features to perform enrichment for ranked by variable in `max_diff`.

            Defaults to 1000.
        sort_var : :obj:`str`, optional
            Variable to sort for when setting `max_diff`.

            Defaults to "pvalue".
        distributed: :obj:`bool`, optional
            Whether work should be submitted as jobs in a computing cluster.

            Defaults to :obj:`False`.
        overwrite: :obj:`bool`, optional
            Whether output files should be overwritten when `distributed` is :obj:`True`.

            Defaults to :obj:`False`.

        Attributes
        ----------
        enrichment_results : :obj:`dict`
            Dictionary with keys as in `steps` and values with pandas.DataFrame
            of enrichment results.
        """
        # TODO: inspect given matrix and warn if matrix hasn't been subset for 'significant' features
        # TODO: separate and fix mouse TF ids
        # TODO: separate homer_consensus output processing
        # TODO: add overwrite function when distributed==False
        from ngs_toolkit.parsers import parse_ame, parse_homer
        from ngs_toolkit.general import enrichr, run_enrichment_jobs
        from ngs_toolkit.utils import get_this_file_or_timestamped
        from tqdm import tqdm

        serial = not distributed

        if differential is None:
            msg = "Data frame with differential comparison results `differential` "
            msg += "was not passed and is not defined in the Analysis object as a "
            msg += "`differential_results` attribute."
            hint = " Run analysis.differential_analysis to get differential results."
            try:
                differential = self.differential_results
            except AttributeError as e:
                _LOGGER.error(msg)
                _LOGGER.info(hint)
                raise e
            if differential is None:
                _LOGGER.error(msg)
                _LOGGER.info(hint)
                raise ValueError

        if genome is None:
            genome = self.genome

        matrix = self.matrix_features

        known = ["region", "lola", "meme", "homer", "enrichr"]
        if not all([x in known for x in steps]):
            _LOGGER.warning("Not all provided steps for enrichment are known! Proceeding anyway.")

        output_dir = self._format_string_with_attributes(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        region_enr = list()
        lola_enr = list()
        meme_enr = list()
        homer_enr = list()
        pathway_enr = list()
        possible_steps = [
            # name, df, function, kwargs, suffix, joint_output_suffix
            (
                "region",
                region_enr,
                pd.read_csv,
                {},
                "region_type_enrichment.csv",
                ".region_type_enrichment.csv",
            ),
            ("meme", meme_enr, parse_ame, {}, "ame.txt", ".meme_ame.csv"),
            ("homer", homer_enr, parse_homer, {}, "homerResults", ".homer_motifs.csv"),
            ("lola", lola_enr, pd.read_csv, {"sep": "\t"}, "allEnrichments.tsv", ".lola.csv",),
            (
                "enrichr",
                pathway_enr,
                pd.read_csv,
                {"encoding": "utf-8"},
                output_prefix + ".enrichr.csv",
                ".enrichr.csv",
            ),
        ]

        gene_level_data_types = ["RNA-seq", "CRISPR"]
        region_level_data_types = ["ATAC-seq", "ChIP-seq", "CNV"]

        if self.data_type in gene_level_data_types:
            possible_steps = [x for x in possible_steps if x[0] == "enrichr"]

        # Examine each region cluster
        comps = differential["comparison_name"].drop_duplicates()
        for comp in tqdm(comps, total=len(comps), desc="Comparison"):
            if directional:
                # Separate in up/down-regulated genes
                params = [(np.less, 0, "down", "head"), (np.greater, 0, "up", "tail")]
            else:
                # get all genes
                params = [(np.less, np.inf, "all", "head")]

            for f, arg, direction, top in params:
                if directional:
                    diff = differential.loc[
                        (differential["comparison_name"] == comp)
                        & (f(differential["log2FoldChange"], arg)),
                        :,
                    ].index
                else:
                    diff = differential.loc[(differential["comparison_name"] == comp), :].index

                # Handle extremes of regions
                if diff.shape[0] < 1:
                    continue
                if diff.shape[0] > max_diff:
                    if directional:
                        diff = getattr(
                            differential[
                                (differential["comparison_name"] == comp)
                                & (f(differential["log2FoldChange"], arg))
                            ][sort_var].sort_values(),
                            top,
                        )(max_diff).index
                    else:
                        diff = getattr(
                            differential[(differential["comparison_name"] == comp)][
                                sort_var
                            ].sort_values(),
                            top,
                        )(max_diff).index

                # Add data_type specific info
                comparison_df = matrix.loc[diff, :]
                if comparison_df.shape != comparison_df.dropna().shape:
                    _LOGGER.warning(
                        "There are differential regions which are not in the set"
                        + " of annotated regions for comparison '{}'!".format(comp)
                        + " Continuing enrichment without those."
                    )
                    comparison_df = comparison_df.dropna()

                # Prepare output dir
                comparison_dir = os.path.join(output_dir, "{}.{}".format(comp, direction))
                if not os.path.exists(comparison_dir):
                    os.makedirs(comparison_dir)

                # Prepare files and run (if not distributed)
                if self.data_type in gene_level_data_types:
                    _LOGGER.debug(
                        "Doing genes of comparison '{}', direction '{}'.".format(comp, direction)
                    )
                    comparison_df.index.name = "gene_name"
                    # write gene names to file
                    clean = comparison_df.reset_index()["gene_name"].drop_duplicates().sort_values()
                    clean.to_csv(
                        os.path.join(comparison_dir, output_prefix + ".gene_symbols.txt"),
                        header=False,
                        index=False,
                    )

                    if "enrichr" in steps:
                        if serial:
                            if not os.path.exists(
                                os.path.join(comparison_dir, output_prefix + ".enrichr.csv")
                            ):
                                enr = enrichr(comparison_df["gene_name"])
                                enr.to_csv(
                                    os.path.join(comparison_dir, output_prefix + ".enrichr.csv"),
                                    index=False,
                                    encoding="utf-8",
                                )
                elif self.data_type in region_level_data_types:
                    _LOGGER.debug(
                        "Doing regions of comparison '{}', direction '{}'.".format(comp, direction)
                    )
                    # do the suite of region enrichment analysis
                    self.characterize_regions_function(
                        comparison_df,
                        output_dir=comparison_dir,
                        prefix=output_prefix,
                        run=serial,
                        genome=genome,
                        steps=steps,
                    )

                # collect enrichments
                if serial:
                    # read/parse, label and append
                    for name, df, function, kwargs, suffix, _ in possible_steps:
                        if name in steps:
                            enr = function(
                                get_this_file_or_timestamped(os.path.join(comparison_dir, suffix)),
                                **kwargs
                            )
                            enr.loc[:, "comparison_name"] = comp
                            enr.loc[:, "direction"] = direction
                            enr.loc[:, "label"] = "{}.{}".format(comp, direction)
                            df.append(enr)

        if serial:
            # write combined enrichments
            _LOGGER.info("Saving combined enrichments for all comparisons.")
            self.enrichment_results = dict()

            for name, df, function, kwargs, suffix, output_suffix in possible_steps:
                if name in steps:
                    self.enrichment_results[name] = pd.concat(df, axis=0)
                    self.enrichment_results[name].to_csv(
                        os.path.join(output_dir, output_prefix + output_suffix),
                        index=False,
                        encoding="utf-8",
                    )
        else:
            background = ""
            if self.data_type in region_level_data_types:
                try:
                    _LOGGER.info("Using background region set from analysis.sites")
                    background = getattr(self, "sites").fn
                except AttributeError:
                    _LOGGER.warning(
                        "Using no background region set because 'analysis.sites' is not set!"
                    )
            _LOGGER.info("Submitting enrichment jobs.")
            run_enrichment_jobs(
                results_dir=output_dir,
                genome=genome,
                background_bed=background,
                steps=steps,
                overwrite=overwrite,
                pep_config=self.prj.config_file,
            )

    def collect_differential_enrichment(
        self,
        steps=["region", "lola", "meme", "homer", "homer_consensus", "enrichr"],
        directional=True,
        permissive=True,
        output_dir="{results_dir}/differential_analysis_{data_type}/enrichments",
        input_prefix="differential_analysis",
        output_prefix="differential_analysis",
        differential=None,
    ):
        """
        Collect the results of enrichment analysis ran after a differential analysis.

        Parameters
        ----------
        steps : :obj:`list`, optional
            Steps of the enrichment analysis to collect results for.

            Defaults to ["region", "lola", "meme", "homer", "enrichr"].
        directional: :obj:`bool`, optional
            Whether enrichments were made in a direction-dependent way
            (up-regulated and down-regulated features separately).
            This implies a column named "direction" exists".

            Defaults to :obj:`True`.
        differential : :obj:`pandas.DataFrame`, optional
            Data frame with differential results to select which comparisons to collect
            enrichments for. Usually produced by ``ngs_toolkit.general.differential_analysis``.

            Defaults to analysis's ``differential_results`` attributes.
        output_dir : :obj:`str`, optional
            Directory to create output files.

            Defaults to "{results_dir}/differential_analysis_{data_type}".
        input_prefix, output_prefix : :obj:`str`, optional
            File prefix of input/output files.

            Defaults to "differential_analysis".
        permissive : :obj:`bool`, optional
            Whether to skip non-existing files, giving a warning.

            Defaults to :obj:`True`.

        Attributes
        ----------
        enrichment_results : :obj:`dict`
            Dictionary with keys as in ``steps`` and values with pandas.DataFrame
            of enrichment results.
        """
        # TODO: separate and fix mouse TF ids
        # TODO: separate homer_consensus output processing
        from ngs_toolkit.parsers import parse_ame, parse_homer
        from tqdm import tqdm

        output_dir = self._format_string_with_attributes(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        data_type_steps = {
            "ATAC-seq": ["region", "lola", "meme", "homer", "homer_consensus", "enrichr",],
            "ChIP-seq": ["region", "lola", "meme", "homer", "homer_consensus", "enrichr",],
            "RNA-seq": ["enrichr"],
        }
        if steps is None:
            steps = [s for s in steps if s in data_type_steps[self.data_type]]
        for step in steps:
            if step not in data_type_steps[self.data_type]:
                steps.pop(steps.index(step))
        if len(steps) == 0:
            msg = "No valid steps for the respective data type selected."
            _LOGGER.error(msg)
            raise ValueError(msg)

        if differential is None:
            differential = self.differential_results

        msg = "Collecting enrichments of comparison '{}', direction '{}'."
        error_msg = "{} results for comparison '{}', direction '{}' were not found!"

        region_enr = list()
        lola_enr = list()
        meme_enr = list()
        homer_enr = list()
        homer_consensus = list()
        pathway_enr = list()
        possible_steps = [
            # name, df, function, kwargs, suffix, joint_output_suffix
            (
                "region",
                region_enr,
                pd.read_csv,
                {},
                "region_type_enrichment.csv",
                ".region_type_enrichment.csv",
            ),
            ("meme", meme_enr, parse_ame, {}, "ame.txt", ".meme_ame.csv"),
            ("homer", homer_enr, parse_homer, {}, "homerResults", ".homer_motifs.csv"),
            (
                "homer_consensus",
                homer_consensus,
                pd.read_csv,
                {"sep": "\t"},
                "knownResults.txt",
                ".homer_consensus.csv",
            ),
            ("lola", lola_enr, pd.read_csv, {"sep": "\t"}, "allEnrichments.tsv", ".lola.csv",),
            (
                "enrichr",
                pathway_enr,
                pd.read_csv,
                {"encoding": "utf-8"},
                input_prefix + ".enrichr.csv",
                ".enrichr.csv",
            ),
        ]

        if self.data_type == "RNA-seq":
            possible_steps = [x for x in possible_steps if x[0] == "enrichr"]

        # Examine each region cluster
        comps = differential["comparison_name"].drop_duplicates()
        for comp in tqdm(comps, total=len(comps), desc="Comparison"):
            if directional:
                # Separate in up/down-regulated genes
                params = list()
                if (
                    differential[
                        (differential["comparison_name"] == comp)
                        & (differential["log2FoldChange"] > 0)
                    ].shape[0]
                    > 0
                ):
                    params.append("up")
                if (
                    differential[
                        (differential["comparison_name"] == comp)
                        & (differential["log2FoldChange"] < 0)
                    ].shape[0]
                    > 0
                ):
                    params.append("down")
                if len(params) == 0:
                    continue
            else:
                params = ["all"]

            for direction in params:
                comparison_dir = os.path.join(output_dir, "{}.{}".format(comp, direction))
                _LOGGER.debug(msg.format(comp, direction))

                # read/parse, label and append
                for name, df, function, kwargs, suffix, _ in possible_steps:
                    if name in steps:
                        try:
                            enr = function(os.path.join(comparison_dir, suffix), **kwargs)
                        except (IOError, AttributeError) as e:
                            if permissive:
                                _LOGGER.warning(error_msg.format(name, comp, direction))
                            else:
                                raise e
                        else:
                            if not enr.empty:
                                enr.loc[:, "comparison_name"] = comp
                                enr.loc[:, "direction"] = direction
                                enr.loc[:, "label"] = "{}.{}".format(comp, direction)
                                df.append(enr)
                            else:
                                _LOGGER.warning(
                                    "Comparison '{}' {} results are empty!".format(comp, name)
                                )

        # write combined enrichments
        _LOGGER.info("Saving combined enrichments for all comparisons.")
        self.enrichment_results = dict()

        for name, df, function, kwargs, suffix, output_suffix in possible_steps:
            if name in steps:
                if len(df) == 0:
                    msg = "No comparison has '{}' results!".format(name)
                    _LOGGER.warning(msg)
                else:
                    self.enrichment_results[name] = pd.concat(df, axis=0)
                    self.enrichment_results[name].to_csv(
                        os.path.join(output_dir, output_prefix + output_suffix),
                        index=False,
                        encoding="utf-8",
                    )

    def plot_differential_enrichment(
        self,
        steps=["region", "lola", "meme", "homer_consensus", "great", "enrichr"],
        plot_types=["barplots", "scatter", "correlation", "heatmap"],
        enrichment_type=None,
        enrichment_table=None,
        direction_dependent=True,
        output_dir="{results_dir}/differential_analysis_{data_type}/enrichments",
        comp_variable="comparison_name",
        output_prefix="differential_analysis",
        rasterized=True,
        clustermap_metric="correlation",
        top_n=5,
        z_score=0,
        cmap=None,
    ):
        """
        Make plots illustrating enrichment of features for various comparisons.

        Input can be the dictionary under `analysis.enrichment_results` or
        a single dataframe of enrichment terms across several comparisons for a given type of enrichment.
        In the later case both `enrichment_table` and `enrichment_type` must be given.

        Parameters
        ----------
        steps : :obj:`list`, optional
            Types of the enrichment analysis to plot.
            Options are ["region", "lola", "meme", "homer_consensus", great", "enrichr"].

            Defaults to all keys present in analysis.enrichment_results.
        plot_types : :obj:`list`, optional
            Types of plots to do for each enrichment type.
            One of ["barplot", "scatter", "correlation", "heatmap"].

            Defaults to all of the above.
        enrichment_type : :obj:`str`, optional
            Type of enrichment if run for a single type of enrichment.
            In this case `enrichment_table` must be given.
            One of {"region", "lola", "meme", "great", "enrichr"}.

            Default (``None``) is to run all keys present in analysis.enrichment_results.
        enrichment_table : :obj:`pandas.DataFrame`, optional
            Data frame with enrichment results as produced by
            ``differential_enrichment`` or ``collect_differential_enrichment``.
            If given, `enrichment_type` must be given too.

            Default (``None``) is the dataframes in all values present in analysis.enrichment_results.
        direction_dependent: :obj:`bool`, optional
            Whether enrichments were made in a direction-dependent way (up-regulated and down-regulated features separately).
            This implies a column named "direction" exists".

            Defaults to :obj:`True`.
        output_dir : :obj:`str`, optional
            Directory to create output files.

            Defaults to "{results_dir}/differential_analysis_{data_type}/enrichments".
        comp_variable : :obj:`str`, optional
            Column defining which comparison enrichment terms belong to.

            Defaults to "comparison_name".
        output_prefix : :obj:`str`, optional
            Prefix to use when creating output files.

            Defaults to "differential_analysis".
        rasterized: :obj:`bool`, optional
            Whether or not to rasterize heatmaps for efficient plotting.

            Defaults to :obj:`True`.
        clustermap_metric : :obj:`str`, optional
            Distance metric to use for clustermap clustering,
            See https://docs.scipy.org/doc/scipy/reference/spatial.distance.html for valid values.

            Default to "correlation" (Pearson's).
        top_n : :obj:`int`, optional
            Top terms to use to display in plots.

            Defaults to 5.
        z_score : {bool, int}, optional
            Which dimention/axis to perform Z-score transformation for.
            Pass :obj:`False` to skip plotting Z-score heatmaps.
            Numpy/Pandas conventions are used:
            `0` is row-wise (in this case across comparisons) and `1` is column-wise (across terms).

            Defaults to 0.
        cmap : :obj:`str`, optional
            Colormap to use in heatmaps.

            Defaults to :obj:`None`.
        """
        # TODO: split function in its smaller parts and call them appropriately.
        import matplotlib
        import matplotlib.pyplot as plt
        from ngs_toolkit.graphics import savefig
        from ngs_toolkit.utils import log_pvalues
        from scipy.stats import zscore
        import seaborn as sns

        def enrichment_barplot(input_df, x, y, group_variable, top_n, output_file):
            n = len(input_df[group_variable].drop_duplicates())
            n_side = int(np.ceil(np.sqrt(n)))

            top_data = (
                input_df.set_index(x).groupby(group_variable)[y].nlargest(top_n).reset_index()
            )

            fig, axis = plt.subplots(
                n_side,
                n_side,
                figsize=(4 * n_side, n_side * max(5, 0.12 * top_n)),
                sharex=False,
                sharey=False,
            )
            if isinstance(axis, np.ndarray):
                axis = iter(axis.flatten())
            else:
                axis = iter(np.array([axis]))
            for comp in top_data[group_variable].drop_duplicates().sort_values():
                df2 = top_data.loc[top_data[group_variable] == comp, :]
                ax = next(axis)
                sns.barplot(
                    df2[y],
                    df2[x],
                    estimator=max,
                    orient="horizontal",
                    ax=ax,
                    color=sns.color_palette("colorblind")[0],
                )
                ax.set_title(comp)
            for ax in axis:
                ax.set_visible(False)
            sns.despine(fig)
            savefig(fig, output_file)

        def enrichment_correlation_plot(
            input_df, output_file, label="Pearson correlation of enrichment"
        ):
            try:
                g = sns.clustermap(
                    input_df.T.corr(),
                    rasterized=rasterized,
                    xticklabels=True,
                    yticklabels=True,
                    cbar_kws={"label": label},
                )
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
                savefig(g.fig, output_file)
            except FloatingPointError:
                msg = "Plotting of correlation matrix failed: {}".format(output_file)
                _LOGGER.warning(msg)

        def enrichment_clustermap(
            input_df,
            output_file,
            label="Enrichment\nof differential regions",
            z_score=None,
            params=None,
        ):
            if params is None:
                params = dict()
            # plot clustered heatmap
            shape = input_df.shape

            # # fix some labels
            input_df.index = (
                input_df.index.str.replace(r"_Homo.*", "")
                .str.replace(r"_Mus.*", "")
                .str.replace(r" \(GO:.*", "")
                .str.replace("_", " ")
            )
            if z_score is not None:
                params.update({"cmap": "RdBu_r", "center": 0, "z_score": z_score})
            try:
                g = sns.clustermap(
                    input_df,
                    figsize=(0.25 * shape[1], 0.20 * shape[0]),
                    metric=clustermap_metric,
                    xticklabels=True,
                    yticklabels=True,
                    rasterized=rasterized,
                    cbar_kws={"label": label},
                    **params
                )
                g.ax_heatmap.set_xticklabels(
                    g.ax_heatmap.get_xticklabels(),
                    rotation=90,
                    va="top",
                    ha="center",
                    fontsize="xx-small",
                )
                g.ax_heatmap.set_yticklabels(
                    g.ax_heatmap.get_yticklabels(),
                    rotation=0,
                    ha="left",
                    va="center",
                    fontsize="xx-small",
                )
                g.ax_heatmap.set_xlabel("Comparison")
                g.ax_heatmap.set_ylabel("Term")
                savefig(g.fig, output_file)
            except FloatingPointError:
                msg = "Plotting of correlation matrix failed: {}".format(output_file)
                _LOGGER.warning(msg)

        if steps is None:
            steps = ["region", "lola", "meme", "homer_consensus", "great", "enrichr"]

        if (enrichment_table is None) and (enrichment_type is None):
            if not hasattr(self, "enrichment_results"):
                msg = "'enrichment_table' and 'enrichment_type' were not given "
                msg += "but analysis also does not have a 'enrichment_results' attribute."
                _LOGGER.error(msg)
                raise ValueError(msg)
            else:
                for step in steps:

                    if step not in self.enrichment_results:
                        msg = "'{}' in steps but it is not a value in analysis.enrichment_results!"
                        continue
                        _LOGGER.warning(msg)
                    enrichment_table = self.enrichment_results[step]
                    self.plot_differential_enrichment(
                        steps=[step],
                        plot_types=plot_types,
                        enrichment_table=enrichment_table,
                        enrichment_type=step,
                        direction_dependent=direction_dependent,
                        output_dir=output_dir,
                        comp_variable=comp_variable,
                        output_prefix=output_prefix,
                        rasterized=rasterized,
                        clustermap_metric=clustermap_metric,
                        top_n=top_n if step != "meme" else 300,
                        z_score=z_score,
                        cmap=cmap,
                    )
                return

        if enrichment_type not in [
            "region",
            "lola",
            "enrichr",
            "meme",
            "homer_consensus",
            "great",
        ]:
            raise AssertionError(
                "`enrichment_type` must be one of 'lola', 'enrichr', 'meme', 'homer_consensus', 'great'."
            )

        output_dir = self._format_string_with_attributes(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        if z_score == 0:
            z_score_label = "Row"
        elif z_score == 1:
            z_score_label = "Column"
        elif z_score in [None, False]:
            pass
        else:
            raise ValueError("Argument 'z_score' must be on of 0, 1 or None.")

        enrichment_table = enrichment_table.copy()
        if ("direction" in enrichment_table.columns) and direction_dependent:
            enrichment_table.loc[:, comp_variable] = (
                enrichment_table[comp_variable].astype(str).str.replace("_", " ")
                + " - "
                + enrichment_table["direction"].astype(str)
            )

        if enrichment_type == "region":
            _LOGGER.info("Plotting enrichments for 'region'")
            from ngs_toolkit.graphics import plot_region_context_enrichment

            enrichment_table["-log10(p-value)"] = log_pvalues(enrichment_table["p_value"])

            if not enrichment_table.index.name == "region":
                enrichment_table = enrichment_table.set_index("region")
            # Plot terms of each comparison in barplots and volcano plots
            plot_region_context_enrichment(
                enrichment_table,
                output_dir=output_dir,
                output_prefix=output_prefix,
                across_attribute=comp_variable,
                pvalue=0.05,
                top_n=top_n,
            )

            # pivot table
            if any([x in plot_types for x in ("correlation", "heatmaps")]):
                region_pivot = pd.pivot_table(
                    enrichment_table,
                    values="log2_odds_ratio",
                    columns=comp_variable,
                    index="region",
                ).fillna(0)

                # plot correlation
                if "correlation" in plot_types:
                    enrichment_correlation_plot(
                        input_df=region_pivot,
                        label="Correlation of enrichment\nof differential regions",
                        output_file=os.path.join(
                            output_dir, output_prefix + ".region_type_enrichment.correlation.svg",
                        ),
                    )

                # plot clustered heatmaps
                if "heatmap" in plot_types:
                    enrichment_clustermap(
                        region_pivot,
                        output_file=os.path.join(
                            output_dir,
                            output_prefix + ".region_type_enrichment.cluster_specific.svg",
                        ),
                        label="log2(odd ratio) of enrichment\nof differential regions",
                        params={"cmap": "RdBu_r", "center": 0},
                    )

        if enrichment_type == "lola":
            _LOGGER.info("Plotting enrichments for 'lola'")
            cols = _CONFIG["resources"]["lola"]["region_set_labeling_columns"]
            odds_col = _CONFIG["resources"]["lola"]["output_column_names"]["odds_ratio"]
            pval_col = _CONFIG["resources"]["lola"]["output_column_names"]["log_p_value"]
            # get a unique label for each lola region set
            cols = [col for col in cols if col in enrichment_table.columns]
            if not cols:
                raise ValueError(
                    "None of the columns present in were found in {} the LOLA results.".format(
                        "CONFIG:resources:lola:region_set_labeling_columns"
                    )
                )
            enrichment_table.loc[:, "label"] = (
                enrichment_table[cols].astype(str).apply(", ".join, axis=1)
            )
            enrichment_table.loc[:, "label"] = (
                enrichment_table["label"]
                .str.replace("nan", "")
                .str.replace("None", "")
                .str.replace(", $", "")
                .str.replace(", , ", ", ")
                .str.decode("unicode_escape")
            ).astype(str)

            # Replace inf values with maximum non-inf p-value observed
            r = enrichment_table.loc[enrichment_table[pval_col] != np.inf, pval_col].max()
            r += r * 0.1
            enrichment_table.loc[:, pval_col] = enrichment_table[pval_col].replace(np.inf, r)

            # Plot top_n terms of each comparison in barplots
            if "barplots" in plot_types:
                enrichment_barplot(
                    enrichment_table,
                    x="label",
                    y=pval_col,
                    group_variable=comp_variable,
                    top_n=top_n,
                    output_file=os.path.join(
                        output_dir, output_prefix + ".lola.barplot.top_{}.svg".format(top_n),
                    ),
                )

            # Significance vs fold enrichment over background
            if "scatter" in plot_types:
                n = len(enrichment_table[comp_variable].drop_duplicates())
                n_side = int(np.ceil(np.sqrt(n)))
                fig, axis = plt.subplots(
                    n_side,
                    n_side,
                    figsize=(3 * n_side, 3 * n_side),
                    sharex=False,
                    sharey=False,
                    squeeze=False,
                )
                axis = axis.flatten()
                for i, comp in enumerate(
                    enrichment_table[comp_variable].drop_duplicates().sort_values()
                ):
                    enr = enrichment_table[enrichment_table[comp_variable] == comp].reset_index(
                        drop=True
                    )
                    enr.loc[:, "combined"] = enr[[odds_col, pval_col]].apply(zscore).mean(axis=1)
                    axis[i].scatter(
                        enr[odds_col], enr[pval_col], c=enr["combined"], s=8, alpha=0.75,
                    )

                    # label top points
                    for j in enr[pval_col].sort_values().tail(5).index:
                        axis[i].text(
                            enr.loc[j, odds_col],
                            enr.loc[j, pval_col],
                            s=enr.loc[j, "label"],
                            ha="right",
                            fontsize=5,
                        )
                    axis[i].set_title(comp)

                for ax in axis.reshape((n_side, n_side))[:, 0]:
                    ax.set_ylabel("-log10(p-value)")
                for ax in axis.reshape((n_side, n_side))[-1, :]:
                    ax.set_xlabel("log odds ratio")
                sns.despine(fig)
                savefig(
                    fig, os.path.join(output_dir, output_prefix + ".lola.scatterplot.svg"),
                )

            # Plot heatmaps of terms for each comparison
            if len(enrichment_table[comp_variable].drop_duplicates()) < 2:
                return

            # pivot table
            if ("correlation" not in plot_types) and ("heatmap" not in plot_types):
                return
            lola_pivot = pd.pivot_table(
                enrichment_table, values=pval_col, columns=comp_variable, index="label",
            ).fillna(0)
            lola_pivot = lola_pivot.replace(np.inf, lola_pivot[lola_pivot != np.inf].max().max())

            # plot correlation
            if "correlation" in plot_types:
                enrichment_correlation_plot(
                    input_df=lola_pivot,
                    label="Correlation of enrichment\nof differential regions",
                    output_file=os.path.join(output_dir, output_prefix + ".lola.correlation.svg"),
                )

            # plot clustered heatmaps of top terms
            if "heatmap" in plot_types:
                top = (
                    enrichment_table.set_index("label")
                    .groupby(comp_variable)[pval_col]
                    .nlargest(top_n)
                )
                top_terms = top.index.get_level_values("label").unique()
                enrichment_clustermap(
                    lola_pivot.loc[top_terms, :],
                    output_file=os.path.join(
                        output_dir, output_prefix + ".lola.cluster_specific.svg"
                    ),
                    label="-log10(p-value) of enrichment\nof differential regions",
                )
                if z_score is not False:
                    enrichment_clustermap(
                        lola_pivot.loc[top_terms, :],
                        output_file=os.path.join(
                            output_dir,
                            output_prefix
                            + ".lola.cluster_specific.{}_z_score.svg".format(z_score_label),
                        ),
                        label="{} Z-score of enrichment\nof differential regions".format(
                            z_score_label
                        ),
                        z_score=z_score,
                    )

        if enrichment_type == "meme":
            _LOGGER.info("Plotting enrichments for 'meme'")
            enrichment_table.loc[:, "log_p_value"] = log_pvalues(enrichment_table["p_value"])

            # Plot top_n terms of each comparison in barplots
            if "barplots" in plot_types:
                enrichment_barplot(
                    enrichment_table,
                    x="TF",
                    y="log_p_value",
                    group_variable=comp_variable,
                    top_n=top_n,
                    output_file=os.path.join(
                        output_dir, output_prefix + ".motifs.barplot.top_{}.svg".format(top_n),
                    ),
                )

            if len(enrichment_table[comp_variable].drop_duplicates()) < 2:
                return
            # Plot heatmaps of terms for each comparison
            if ("correlation" not in plot_types) and ("heatmap" not in plot_types):
                return
            motifs_pivot = pd.pivot_table(
                enrichment_table, values="log_p_value", columns="TF", index=comp_variable,
            ).fillna(0)
            # plot correlation
            if "correlation" in plot_types:
                enrichment_correlation_plot(
                    input_df=motifs_pivot,
                    label="Correlation of enrichment\nof differential regions",
                    output_file=os.path.join(output_dir, output_prefix + ".motifs.correlation.svg"),
                )

            # plot clustered heatmaps of top terms
            if "heatmap" in plot_types:
                top = (
                    enrichment_table.set_index("TF")
                    .groupby(comp_variable)["log_p_value"]
                    .nlargest(top_n)
                )
                top_terms = top.index.get_level_values("TF").unique()
                enrichment_clustermap(
                    motifs_pivot.loc[:, top_terms],
                    output_file=os.path.join(
                        output_dir, output_prefix + ".motifs.cluster_specific.svg"
                    ),
                    label="-log10(p-value) of enrichment\nof differential regions",
                )
                if z_score is not False:
                    enrichment_clustermap(
                        motifs_pivot.loc[:, top_terms],
                        output_file=os.path.join(
                            output_dir,
                            output_prefix
                            + ".motifs.cluster_specific.{}_z_score.svg".format(z_score_label),
                        ),
                        label="{} Z-score of enrichment\nof differential regions".format(
                            z_score_label
                        ),
                        z_score=z_score,
                    )

        if enrichment_type == "homer_consensus":
            _LOGGER.info("Plotting enrichments for 'homer_consensus'")
            enrichment_table.loc[:, "enrichment_over_background"] = (
                enrichment_table["% of Target Sequences with Motif"]
                / enrichment_table["% of Background Sequences with Motif"]
            )
            enrichment_table.loc[:, "log_p_value"] = log_pvalues(enrichment_table["P-value"])

            # Plot top_n terms of each comparison in barplots
            top_n = min(
                top_n,
                enrichment_table.set_index("Motif Name")
                .groupby(comp_variable)["log_p_value"]
                .count()
                .min()
                - 1,
            )
            if "barplots" in plot_types:
                enrichment_barplot(
                    enrichment_table,
                    x="Motif Name",
                    y="log_p_value",
                    group_variable=comp_variable,
                    top_n=top_n,
                    output_file=os.path.join(
                        output_dir,
                        output_prefix + ".homer_consensus.barplot.top_{}.svg".format(top_n),
                    ),
                )

            # Significance vs fold enrichment over background
            if "scatter" in plot_types:
                n = len(enrichment_table[comp_variable].drop_duplicates())
                n_side = int(np.ceil(np.sqrt(n)))
                fig, axis = plt.subplots(
                    n_side, n_side, figsize=(3 * n_side, 3 * n_side), sharex=False, sharey=False,
                )
                axis = axis.flatten()
                for i, comp in enumerate(
                    enrichment_table[comp_variable].drop_duplicates().sort_values()
                ):
                    enr = enrichment_table[enrichment_table[comp_variable] == comp]
                    enr.loc[:, "Motif Name"] = (
                        enr["Motif Name"]
                        .str.replace(".*BestGuess:", "")
                        .str.replace(r"-ChIP-Seq.*", "")
                    )

                    enr.loc[:, "combined"] = (
                        enr[["enrichment_over_background", "log_p_value"]]
                        .apply(zscore)
                        .mean(axis=1)
                    )
                    axis[i].scatter(
                        enr["enrichment_over_background"],
                        enr["log_p_value"],
                        c=enr["combined"],
                        s=8,
                        alpha=0.75,
                    )

                    # label top points
                    for j in enr["combined"].sort_values().tail(5).index:
                        axis[i].text(
                            enr.loc[j, "enrichment_over_background"],
                            enr.loc[j, "log_p_value"],
                            s=enr.loc[j, "Motif Name"],
                            ha="right",
                            fontsize=5,
                        )
                    axis[i].set_title(comp)

                for ax in axis.reshape((n_side, n_side))[:, 0]:
                    ax.set_ylabel("-log10(p-value)")
                for ax in axis.reshape((n_side, n_side))[-1, :]:
                    ax.set_xlabel("Enrichment over background")
                sns.despine(fig)
                savefig(
                    fig,
                    os.path.join(output_dir, output_prefix + ".homer_consensus.scatterplot.svg"),
                )

            # Plot heatmaps of terms for each comparison
            if len(enrichment_table[comp_variable].drop_duplicates()) < 2:
                return

            for label, metric in [
                ("-log10(p-value) of enrichment\nin", "log_p_value"),
                ("Enrichment over background of\n", "enrichment_over_background"),
            ]:
                # pivot table
                if ("correlation" not in plot_types) and ("heatmap" not in plot_types):
                    return
                motifs_pivot = pd.pivot_table(
                    enrichment_table, values=metric, columns="Motif Name", index=comp_variable,
                ).fillna(0)

                # plot correlation
                if "correlation" in plot_types:
                    enrichment_correlation_plot(
                        input_df=motifs_pivot,
                        label="Correlation of enrichment\nof differential regions",
                        output_file=os.path.join(
                            output_dir, output_prefix + ".homer_consensus.correlation.svg",
                        ),
                    )

                top = (
                    enrichment_table.set_index("Motif Name")
                    .groupby(comp_variable)[metric]
                    .nlargest(top_n)
                )
                top_terms = top.index.get_level_values("Motif Name").unique()

                # plot clustered heatmaps of top terms
                if "heatmap" in plot_types:
                    enrichment_clustermap(
                        motifs_pivot.loc[:, top_terms].T,
                        output_file=os.path.join(
                            output_dir, output_prefix + ".homer_consensus.cluster_specific.svg",
                        ),
                        label=label + " differential regions",
                    )
                    if z_score is not False:
                        enrichment_clustermap(
                            motifs_pivot.loc[:, top_terms].T,
                            output_file=os.path.join(
                                output_dir,
                                output_prefix
                                + ".homer_consensus.cluster_specific.{}_z_score.svg".format(
                                    z_score_label
                                ),
                            ),
                            label="{} Z-score of {} differential regions".format(
                                z_score_label, label
                            ),
                            z_score=z_score,
                        )

        if enrichment_type == "enrichr":
            _LOGGER.info("Plotting enrichments for 'enrichr'")
            enrichment_table["log_p_value"] = log_pvalues(enrichment_table["p_value"])

            for gene_set_library in enrichment_table["gene_set_library"].unique():
                _LOGGER.debug(gene_set_library)

                # Plot top_n terms of each comparison in barplots
                n = len(enrichment_table[comp_variable].drop_duplicates())
                n_side = int(np.ceil(np.sqrt(n)))

                if "barplots" in plot_types:
                    enrichment_barplot(
                        enrichment_table.loc[
                            enrichment_table["gene_set_library"] == gene_set_library
                        ],
                        x="description",
                        y="log_p_value",
                        group_variable=comp_variable,
                        top_n=top_n,
                        output_file=os.path.join(
                            output_dir,
                            output_prefix
                            + ".enrichr.{}.barplot.top_{}.svg".format(gene_set_library, top_n),
                        ),
                    )

                    # # ^^ possible replacement
                    # grid = sns.catplot(
                    #     data=top_data, x="log_p_value", y="description",
                    #     order=top_data.groupby("description")["log_p_value"].mean().sort_values(ascending=False).index,
                    #     kind="bar", orient="horiz", col=comp_variable, col_wrap=n_side, palette="magma_r")
                    # grid.savefig(os.path.join(output_dir, output_prefix + ".enrichr.{}.barplot.top_{}.joint_comparisons.svg".format(
                    #         gene_set_library, top_n)), bbox_inches="tight", dpi=300)

                # Scatter plots of Z-score vs p-value vs combined score
                if "scatter" in plot_types:
                    fig, axis = plt.subplots(
                        n_side, n_side, figsize=(4 * n_side, 4 * n_side), sharex=True, sharey=True,
                    )
                    axis = axis.flatten()
                    # normalize color across comparisons
                    d = enrichment_table.loc[
                        (enrichment_table["gene_set_library"] == gene_set_library),
                        "combined_score",
                    ].describe()
                    norm = matplotlib.colors.Normalize(vmin=d["min"], vmax=d["max"])
                    for i, comparison in enumerate(enrichment_table[comp_variable].unique()):
                        enr = enrichment_table[
                            (enrichment_table["gene_set_library"] == gene_set_library)
                            & (enrichment_table[comp_variable] == comparison)
                        ]
                        sns.scatterplot(
                            data=enr,
                            x="z_score",
                            y="log_p_value",
                            size="combined_score",
                            hue="combined_score",
                            hue_norm=norm,
                            ax=axis[i],
                            rasterized=rasterized,
                            palette="magma",
                        )
                        axis[i].set_title(comparison)

                        for metric in ["log_p_value", "z_score", "combined_score"]:
                            f = (
                                pd.DataFrame.nsmallest
                                if metric == "z_score"
                                else pd.DataFrame.nlargest
                            )
                            for _, s in f(enr, 5, metric).iterrows():
                                axis[i].text(s["z_score"], s["log_p_value"], s=s["description"])
                    sns.despine(fig)
                    savefig(
                        fig,
                        os.path.join(
                            output_dir,
                            output_prefix
                            + ".enrichr.{}.zscore_vs_pvalue.scatterplot.svg".format(
                                gene_set_library
                            ),
                        ),
                    )

                # Plot heatmaps of terms for each comparison
                if len(enrichment_table[comp_variable].drop_duplicates()) < 2:
                    continue

                # pivot table
                if ("correlation" not in plot_types) and ("heatmap" not in plot_types):
                    return
                enrichr_pivot = pd.pivot_table(
                    enrichment_table[enrichment_table["gene_set_library"] == gene_set_library],
                    values="log_p_value",
                    columns="description",
                    index=comp_variable,
                ).fillna(0)

                # plot correlation
                if "correlation" in plot_types:
                    enrichment_correlation_plot(
                        input_df=enrichr_pivot,
                        label="Correlation of enrichment\nof differential gene sets",
                        output_file=os.path.join(
                            output_dir,
                            output_prefix + ".enrichr.{}.correlation.svg".format(gene_set_library),
                        ),
                    )

                top = (
                    enrichment_table[enrichment_table["gene_set_library"] == gene_set_library]
                    .set_index("description")
                    .groupby(comp_variable)["p_value"]
                    .nsmallest(top_n)
                )
                top_terms = top.index.get_level_values("description").unique()
                # top_terms = top_terms[top_terms.isin(lola_pivot.columns[lola_pivot.sum() > 5])]

                # plot clustered heatmap
                if "heatmap" in plot_types:
                    enrichment_clustermap(
                        enrichr_pivot[list(set(top_terms))].T,
                        output_file=os.path.join(
                            output_dir,
                            output_prefix
                            + ".enrichr.{}.cluster_specific.svg".format(gene_set_library),
                        ),
                        label="-log10(p-value) of enrichment\nof differential genes",
                    )
                    if z_score is not False:
                        enrichment_clustermap(
                            enrichr_pivot[list(set(top_terms))].T,
                            output_file=os.path.join(
                                output_dir,
                                output_prefix
                                + ".enrichr.{}.cluster_specific.{}_z_score.svg".format(
                                    gene_set_library, z_score_label
                                ),
                            ),
                            label="{} Z-score of enrichment\nof differential regions".format(
                                z_score_label
                            ),
                            z_score=z_score,
                        )

        if enrichment_type == "great":
            _LOGGER.info("Plotting enrichments for 'great'")
            enrichment_table["log_q_value"] = log_pvalues(enrichment_table["HyperFdrQ"])

            for gene_set_library in enrichment_table["Ontology"].unique():
                _LOGGER.info(gene_set_library)
                if "barplots" in plot_types:
                    # Plot top_n terms of each comparison in barplots
                    enrichment_barplot(
                        enrichment_table.loc[enrichment_table["Ontology"] == gene_set_library],
                        x="description",
                        y="log_p_value",
                        group_variable=comp_variable,
                        top_n=top_n,
                        output_file=os.path.join(
                            output_dir,
                            output_prefix
                            + ".great.{}.barplot.top_{}.svg".format(gene_set_library, top_n),
                        ),
                    )

                # Plot heatmaps of terms for each comparison
                if len(enrichment_table[comp_variable].drop_duplicates()) < 2:
                    return

                # pivot table
                if ("correlation" not in plot_types) and ("heatmap" not in plot_types):
                    return
                great_pivot = pd.pivot_table(
                    enrichment_table[enrichment_table["Ontology"] == gene_set_library],
                    values="log_q_value",
                    columns="Desc",
                    index=comp_variable,
                ).fillna(0)

                # plot correlation
                if "correlation" in plot_types:
                    enrichment_correlation_plot(
                        input_df=great_pivot,
                        label="Correlation of enrichment\nof differential gene sets",
                        output_file=os.path.join(
                            output_dir,
                            output_prefix + ".great.{}.correlation.svg".format(gene_set_library),
                        ),
                    )

                top = (
                    enrichment_table[enrichment_table["Ontology"] == gene_set_library]
                    .set_index("Desc")
                    .groupby(comp_variable)["HyperP"]
                    .nsmallest(top_n)
                )
                top_terms = top.index.get_level_values("Desc").unique()

                # plot clustered heatmaps
                if "heatmap" in plot_types:
                    enrichment_clustermap(
                        great_pivot[list(set(top_terms))].T,
                        output_file=os.path.join(
                            output_dir,
                            output_prefix
                            + ".great.{}.cluster_specific.svg".format(gene_set_library),
                        ),
                        label="-log10(p-value) of enrichment\nof differential genes",
                    )
                    if z_score is not False:
                        enrichment_clustermap(
                            great_pivot[list(set(top_terms))].T,
                            output_file=os.path.join(
                                output_dir,
                                output_prefix
                                + ".great.{}.cluster_specific.{}_z_score.svg".format(
                                    gene_set_library, z_score_label
                                ),
                            ),
                            label="{} Z-score of enrichment\nof differential regions".format(
                                z_score_label
                            ),
                            z_score=z_score,
                        )

    def run_full_analysis_recipe(self, **kwargs):
        """
        Run the :class:`ngs_toolkit.recipes.ngs_analysis` recipe on the current
        Analysis object.

        Parameters
        ----------
        **kwargs : :obj:`dict`
            Additional keyword arguments are passed to
            :func:`ngs_toolkit.recipes.ngs_analysis.main_analysis_pipeline`.
        """
        from ngs_toolkit.recipes.ngs_analysis import main_analysis_pipeline

        main_analysis_pipeline(self)
