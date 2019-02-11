#!/usr/bin/env python


from collections import Counter, defaultdict
import datetime
import itertools
import os
import pickle
import textwrap
import time
import matplotlib
import matplotlib.pyplot as plt

from ngs_toolkit import _CONFIG, _LOGGER
from ngs_toolkit.decorators import check_organism_genome
from ngs_toolkit.graphics import (
    add_colorbar_to_axis,
    savefig,
    plot_projection)
from ngs_toolkit.parsers import parse_ame, parse_homer
from ngs_toolkit.general import (
    deseq_analysis,
    enrichr,
    run_enrichment_jobs,
    get_genome_reference,
    get_blacklist_annotations,
    get_tss_annotations,
    get_genomic_context)
from ngs_toolkit.utils import log_pvalues

import numpy as np
import pandas as pd
from peppy import Project, Sample
from pypiper.ngstk import NGSTk
from scipy.stats import fisher_exact, kruskal, pearsonr, zscore
import seaborn as sns
from sklearn import manifold
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from statsmodels.sandbox.stats.multicomp import multipletests
from tqdm import tqdm


class Analysis(object):
    """
    Generic class holding functions and data from a typical NGS analysis.

    Other modules implement classes inheriting from this that in general contain
    data type-specific functions (e.g. ``ngs_toolkit.atacseq.ATACSeqAnalysis``
    has a ``get_consensus_sites`` function to generate a peak consensus map).

    Objects of this type can be used to store data (e.g. dataframes), variables
    (e.g. paths to files or configurations) and are easily serializable (saved
    to a file as an object) for rapid loading and cross-environment portability.
    See the ``ngs_toolkit.general.Analysis.to_pickle``,
    ``ngs_toolkit.general.Analysis.from_pickle`` and
    ``ngs_toolkit.general.Analysis.update`` functions for this.

    Parameters
    ----------
    name : str, optional
        Name of the analysis.
        Defaults to ``analysis``.

    samples : list, optional
        List of ``peppy.Sample`` objects that this analysis is tied to.
        Defaults to ``None``.

    prj : peppy.Project, optional
        A ``peppy.Project`` object that this analysis is tied to.
        Defaults to ``None``.

    data_dir : str, optional
        Directory containing processed data (e.g. by looper) that will
        be input to the analysis. This is in principle not required.
        Defaults to ``data``.

    results_dir : str, optional
        Directory to contain outputs produced by the analysis.
        Defaults to ``results``.

    pickle_file : str, optional
        A pickle file to serialize the object.
        Defaults to "`name`.pickle".

    from_pickle : bool, optional
        Whether the analysis should be loaded from an existing
        serialized analysis object in ``pickle_file``.
        Defaults to False.

    kwargs : dict, optional
        Additional keyword arguments will simply be stored as object attributes.
    """
    def __init__(
            self,
            name=None,
            samples=None,
            prj=None,
            root_dir=None,
            data_dir="data",
            results_dir="results",
            pickle_file=None,
            from_pickle=False,
            from_pep=False,
            **kwargs):
        # parse kwargs with default
        if root_dir is None:
            self.root_dir = os.curdir
        self.root_dir = os.path.abspath(self.root_dir)

        # # if given absolute paths, keep them, otherwise append to root directory
        for dir_, attr in [(data_dir, "data_dir"), (results_dir, "results_dir")]:
            if not os.path.isabs(dir_):
                setattr(self, attr, os.path.join(self.root_dir, dir_))
            else:
                setattr(self, attr, dir_)

        self.samples = samples
        self.prj = prj
        self.pickle_file = pickle_file

        # parse remaining kwargs
        _LOGGER.debug("Adding additional kwargs to analysis object.")
        self.__dict__.update(kwargs)

        _LOGGER.debug("Setting data type-specific attributes to None.")
        attrs = [
            "data_type", "var_names",
            "quantity", "norm_units", "raw_matrix_name",
            "norm_matrix_name", "annot_matrix_name"]
        for attr in attrs:
            if not hasattr(self, attr):
                setattr(self, attr, None)

        # Generate from PEP configuration file
        if from_pep is not False:
            self.from_pep(from_pep)
            if name is None:
                if hasattr(getattr(self, "prj"), "name"):
                    if getattr(self, "prj").name is not None:
                        name = getattr(getattr(self, "prj"), "name")
                        _LOGGER.info("Setting Project's name as the analysis's name: '{}'.".format(name))
                        self.name = name
        self.name = "analysis" if name is None else name

        # Set default location for the pickle
        if self.pickle_file is None:
            self.pickle_file = os.path.join(results_dir, "{}.pickle".format(self.name))
        # reload itself if required
        if from_pickle:
            _LOGGER.info("Updating analysis object from pickle file: '{}'.".format(self.pickle_file))
            self.update(pickle_file=self.pickle_file)

        for directory in [self.data_dir, self.results_dir]:
            if not os.path.exists(directory):
                os.makedirs(directory)

        # Store projects attributes in self
        _LOGGER.debug("Trying to set analysis attributes.")
        self.set_project_attributes(overwrite=False)

        # Try to set genome if not set
        self.organism, self.genome = (None, None)
        _LOGGER.debug("Trying to get analysis genome.")
        self.set_organism_genome()
        # TODO: if genome is set, get required static files for that genome assembly

        # Add sample input file locations
        _LOGGER.debug("Trying to set sample input file attributes.")
        self.set_samples_input_files()

    def __repr__(self):
        t = "'{}' analysis".format(self.data_type) if self.data_type is not None else "Analysis"
        samples = " with {} samples".format(len(self.samples)) if self.samples is not None else ""
        organism = " of organism '{}'".format(self.organism) if self.organism is not None else ""
        genome = " ({})".format(self.genome) if self.genome is not None else ""
        suffix = "."
        return t + " '{}'".format(self.name) + samples + organism + genome + suffix

    @staticmethod
    def _overwride_sample_representation():
        """
        Make peppy.Sample objects have a more pretty representation.
        """
        def r(self): return self.name
        Sample.__repr__ = r

    @staticmethod
    def _check_data_type_is_supported(data_type):
        """
        Check if data type is allowed in ngs_toolkit configuration.

        Parameters
        ----------
        data_type : str
            Data type to check.

        Parameters
        ----------
        bool
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
        data_type : str, optional
            Data type to check.
        """
        # TODO: test
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
        attr : str
            An attribute of the analysis' samples to check existence of files.

        f : function
            Function to reduce output across samples.
            Defaults to `all`.

        samples : list, optional
            Samples to consider.
            Defaults to all in analysis

        Returns
        ----------
        bool
            Whether samples have file
        """
        if samples is None:
            samples = self.samples
        return f([os.path.exists(str(getattr(sample, attr))) for sample in samples])

    def _get_samples_have_file(self, attr, samples=None):
        """
        Get samples with an existing file under `attr`.

        Parameters
        ----------
        attr : str
            Attribute to check

        samples : list, optional
            Samples to consider.
            Defaults to all in analysis

        Returns
        -------
        list
            List of peppy.Sample objects
        """
        if samples is None:
            samples = self.samples
        return [sample for sample in samples if os.path.exists(str(getattr(sample, attr)))]

    def _get_samples_missing_file(self, attr, samples=None):
        """
        Get samples without an existing file under `attr`.

        Parameters
        ----------
        attr : str
            Attribute to check

        samples : list, optional
            Samples to consider.
            Defaults to all in analysis

        Returns
        -------
        list
            List of peppy.Sample objects.
        """
        if samples is None:
            samples = self.samples
        return [sample for sample in samples if sample not in self._get_samples_have_file(attr, samples=samples)]

    def _get_samples_with_input_file(self, input_file, permissive=False, samples=None):
        """
        Get samples with existing files of attribute `input_file`.

        If none has, raise error. Else, if permissive, return samples with existing file.
        Otherwise return only with all samples have it, otherwise throw IOError.

        Parameters
        ----------
        input_file : str
            Attribute to check

        permissive : bool, optional
            Whether to allow returning a subset of samples if not all have file.
            Defaults to False

        samples : list, optional
            Samples to consider.
            Defaults to all in analysis

        Returns
        -------
        list
            List of peppy.Sample objects

        Raises
        -------
        IOError
            If not permissive and not all sample input files are found.
        """
        if samples is None:
            samples = self.samples
        check = self._check_samples_have_file(attr=input_file, f=any if permissive else all, samples=samples)
        missing = self._get_samples_missing_file(attr=input_file, samples=samples)

        msg = "None of the samples have '{}' files.".format(input_file)
        if all([s in missing for s in samples]):
            _LOGGER.error(msg)
            raise IOError(msg)

        msg = "Not all samples have '{}' files.".format(input_file)
        hint = " Samples missing files: {}".format(", ".join([s.name for s in missing]))
        if not check:
            if permissive:
                _LOGGER.warning(msg + hint)
                return [s for s in samples if s not in missing]
            else:
                _LOGGER.error(msg)
                raise IOError(msg)
        else:
            return samples

    @staticmethod
    def _format_string_with_environment_variables(string):
        """
        Given a string, containing curly braces with dollar sign,
        format it with the environment variables.

        Parameters
        ----------
        string : str
            String to format.

        Returns
        ----------
        str
            Formated string.
        """
        if string is None:
            return string
        to_format = pd.Series(string).str.extractall(r"\${(.*?)}")[0].values
        attrs = os.environ
        if not all([x in attrs for x in to_format]):
            msg = "Not all required patterns were found in the environment variables."
            _LOGGER.error(msg)
            raise ValueError(msg)
        for attr in to_format:
            string = string.replace("${" + attr + "}", "{" + attr + "}")
        return string.format(**attrs)

    def _format_string_with_attributes(self, string):
        """
        Given a string, containing curly braces, format it with the attributes from self.

        Parameters
        ----------
        string : str
            String to format.

        Returns
        ----------
        str
            Formated string.
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

    def from_pep(self, pep_config):
        """
        Create a peppy.Project from a PEP configuration file
        and associate is with the analysis.

        Parameters
        ----------
        pep_config : str
            PEP configuration file.

        Attributes:
        ----------
        prj : peppy.Project
            peppy.Project from given PEP configuration file.
        """
        self.prj = Project(pep_config)
        # msg = "Provided PEP configuration file could not be read."
        # try:
        #     self.prj = Project(pep_config)
        # except (KeyError, yaml.scanner.ScannerError):  # This is for a malformed yaml
        #     # Does not cover a bad path: IsADirectoryError (Python3) and IOError (Python2)
        #     _LOGGER.error(msg)
        #     raise

    def update(self, pickle_file=None):
        """
        Update all of the object"s attributes with the attributes from a serialized
        object (i.e. object stored in a file) object.

        Parameters
        ----------
        pickle_file : str, optional
            Pickle file to load. By default this is the object"s attribute `pickle_file`.
        """
        self.__dict__.update(self.from_pickle(pickle_file=pickle_file).__dict__)

    def set_organism_genome(self):
        """
        Attempt to derive the analysis' organism and genome assembly
        by inspecting the same attributes of its samples.


        Attributes
        ----------
        organism, genome : str
            Organism and genome assembly of the analysis
            if all samples agree in these attributes.

        """
        if self.samples is None:
            _LOGGER.warning("Genome assembly for analysis was not set and cannot be derived from samples.")
        else:
            organisms = list(set([s.organism for s in self.samples]))
            if len(organisms) == 1:
                _LOGGER.info("Setting analysis organism as '{}'.".format(organisms[0]))
                self.organism = organisms[0]
            else:
                _LOGGER.warning("Found several organism for the various analysis samples. " +
                                "Will not set a organism for analysis.")
            genomes = list(set([s.genome for s in self.samples]))
            if len(genomes) == 1:
                _LOGGER.info("Setting analysis genome as '{}'.".format(genomes[0]))
                self.genome = genomes[0]
            else:
                _LOGGER.warning("Found several genome assemblies for the various analysis samples. " +
                                "Will not set a genome for analysis.")

    def set_project_attributes(self, overwrite=True):
        """
        Set Analysis object attributes ``samples``, ``sample_attributes`` and ``group_atrributes``
        to the values in the associated Project object if existing.

        Parameters
        ----------
        overwrite : bool, optional
            Whether to overwrite attribute values if existing.
            Defaults to True

        Attributes
        ----------
        samples : list
            List of peppy.Samples if contained in the PEP configuration

        sample_attributes, group_attributes : list
            Sample attributes if specified in the PEP configuration

        comparison_table : pandas.DataFrame
            Comparison table if specified in the PEP configuration
        """
        hint = " Adding a '{}' section to your project configuration file allows the analysis"
        hint += " object to use those attributes during the analysis."
        if self.prj is not None:
            for attr, parent in [
                    ("samples", self.prj), ("sample_attributes", self.prj),
                    ("group_attributes", self.prj), ("comparison_table", self.prj.metadata)]:
                if not hasattr(parent, attr):
                    _LOGGER.warning("Associated project does not have any '{}'.".format(attr) +
                                    hint.format(attr) if attr != "samples" else "")
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
                                _LOGGER.debug("{} already exist for analysis, not overwriting."
                                              .format(attr.replace("_", " ").capitalize()))
            if hasattr(self, "comparison_table"):
                if isinstance(getattr(self, "comparison_table"), str):
                    _LOGGER.debug("Reading up comparison table.")
                    self.comparison_table = pd.read_csv(self.comparison_table)
        else:
            _LOGGER.warning("Analysis object does not have an attached Project. " +
                            "Will not add special attributes to analysis such as " +
                            "samples, their attributes and comparison table.")

    def set_samples_input_files(self, overwrite=True):
        """
        Add input file values to sample objects dependent on data type.
        These are specified in the `ngs_toolkit` configuration file
        under "sample_input_files:<data type>:<attribute>".

        Parameters
        ----------
        overwrite : bool, optional
            Whether to overwrite attribute values if existing.
            Defaults to True
        """
        if self.samples is None:
            _LOGGER.error("Analysis object does not have attached Samples. " +
                          "Will not add special attributes to samples such as " +
                          "input file locations.")
            return

        msg = "Setting '{}' in sample {} as '{}'."
        for data_type in _CONFIG["sample_input_files"]:
            for sample in [s for s in self.samples if s.protocol == data_type]:
                for attr, value in _CONFIG["sample_input_files"][data_type].items():
                    if value is None:
                        pass
                    elif ("{data_dir}" in value) and ("{sample_name}" in value):
                        value = value.format(data_dir=self.data_dir, sample_name=sample.name)
                    elif "{data_dir}" in value:
                        value = value.format(data_dir=self.data_dir)
                    elif "{sample_name}" in value:
                        value = value.format(sample_name=sample.name)
                    if overwrite:
                        _LOGGER.debug(msg.format(attr, sample.name, value))
                        setattr(sample, attr, value)
                    else:
                        if not hasattr(sample, attr):
                            _LOGGER.debug(msg.format(attr, sample.name, value))
                            setattr(sample, attr, value)
                        else:
                            if getattr(sample, attr) is None:
                                _LOGGER.debug(msg.format(attr, sample.name, value))
                                setattr(sample, attr, value)
                            else:
                                _LOGGER.debug("{} already exists in sample, not overwriting."
                                              .format(attr.replace("_", " ").capitalize()))

    def to_pickle(self, timestamp=False):
        """
        Serialize object (i.e. save to disk) to pickle format.

        Parameters
        ----------
        timestamp : bool
            Whether to timestamp the file.
        """
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
        pickle_file : str, optional
            Pickle file to load.
            By default this is the object"s attribute `pickle_file`.
        """
        if pickle_file is None:
            pickle_file = self.pickle_file
        return pickle.load(open(pickle_file, "rb"))

    @check_organism_genome
    def get_annotations(
            self,
            organism=None, genome_assembly=None, output_dir=None,
            steps=["blacklist", "tss", "genomic_context"], overwrite=False):
        """
        Get genome annotations and other resources for several ngs_toolkit analysis.

        Parameters
        ----------
        organism : str, optional
            Organism to get for. Currently supported are "human" and "mouse".
            Defaults to analysis' own organism.

        genome_assembly : str, optional
            Genome assembly to get for.
            Currently supported are "hg19", "hg38" and "mm10".
            Defaults to analysis' own genome assembly.

        output_dir : str, optional
            Directory to save results to.
            Defaults to the value of "preferences:root_reference_dir" in the configuration,
            if that is not set, to a directory called "reference" in the analysis root directory.

        steps : list, optional
            What kind of annotations to get.
            Options are:
                 - "genome": Genome sequence (2bit format)
                 - "blacklist": Locations of blacklisted regions for genome
                 - "tss": Locations of gene"s TSSs
                 - "genomic_context": Genomic context of genome
            Defaults to ["blacklist", "tss", "genomic_context"]

        overwrite : bool, optional
            Whether existing files should be overwritten by new ones.
            Otherwise they will be kept and no action is made.
            Defaults to False.

        Returns
        -------
        dict
            Dictionary with keys same as the options as steps, containing paths to the requested files.
            The values of the 'genome' step are also a dictionary with keys "2bit" and "fasta" for
            each file type respectively.
        """
        if organism is None:
            organism = self.organism
        if genome_assembly is None:
            genome_assembly = self.genome
        if output_dir is None:
            output_dir = self._format_string_with_environment_variables(
                _CONFIG['preferences']['root_reference_dir'])
            if output_dir is None:
                output_dir = os.path.join(self.root_dir, "reference")

        args = {"organism": organism,
                "genome_assembly": genome_assembly,
                "output_dir": output_dir,
                "overwrite": overwrite}
        mapping = {"hg19": "grch37", "hg38": "grch38", "mm10": "grcm38"}
        output = dict()

        if "genome" in steps:
            output["genome_file"] = dict()
            output["genome_file"]["2bit"] = get_genome_reference(**args)
            fasta = output["genome_file"]["2bit"].replace(".2bit", ".fa")
            if os.path.exists(fasta):
                output["genome_file"]["fasta"] = fasta
        if "blacklist" in steps:
            output["blacklist_file"] = get_blacklist_annotations(**args)
        if "tss" in steps:
            get_tss_annotations(**args)
            output["tss_file"] = os.path.join(
                output_dir, "{}.{}.gene_annotation.protein_coding.tss.bed"
                .format(self.organism, mapping[self.genome]))
        if "genomic_context" in steps:
            get_genomic_context(**args)
            output["genomic_context_file"] = os.path.join(
                output_dir, "{}.{}.genomic_context.bed"
                .format(self.organism, mapping[self.genome]))

        return output

    def get_matrix(self, matrix=None, matrix_name=None, samples=None):
        """
        Return a matrix that is an attribute of self subsetted for the requested samples.

        Parameters
        ----------
        matrix : pandas.DataFrame
            Pandas DataFrame.

        samples : list
            Iterable of peppy.Sample objects to restrict matrix to.
            If not provided (`None` is passed) the matrix will not be subsetted.

        matrix_name : str
            Name of the matrix that is an attribute of the object
            with values for samples in `samples`.

        Returns
        -------
        pandas.DataFrame
            Requested DataFrame.
        """
        if (matrix is None) & (matrix_name is None):
            msg = "Either arguments `matrix` or `matrix_name` must be provided."
            _LOGGER.error(msg)
            raise ValueError(msg)
        # default to matrix to be normalized
        if matrix is None:
            r_matrix = getattr(self, matrix_name)
        else:
            r_matrix = matrix
        # default to all samples in self with matching names in matrix
        if samples is None:
            r_matrix = r_matrix.loc[:, [s.name for s in self.samples]]
        else:
            r_matrix = r_matrix.loc[:, [s.name for s in samples]]

        return r_matrix

    def annotate_with_sample_metadata(
            self,
            quant_matrix=None,
            attributes=None,
            numerical_attributes=None,
            save=True,
            assign=True):
        """
        Annotate matrix ``(n_regions, n_samples)`` with sample metadata
        (creates MultiIndex on columns). Numerical attributes can be pass as a iterable
        to ``numerical_attributes`` to be converted.

        Parameters
        ----------
        quant_matrix : str
            Attribute name of matrix to annotate.
            Default is  infered from the analysis data_type in the following way:
                - ATAC-seq or ChIP-seq: ``coverage_annotated``;
                - RNA-seq: ``expression_annotated``.

        attributes : list
            Desired attributes to be annotated.
            Defaults to all attributes in the original sample annotation sheet of the analysis" Project.

        numerical_attributes : list
            Attributes which are numeric even though they
            might be so in the samples" attributes. Will attempt
            to convert values to numeric.

        save : bool
            Whether to write normalized DataFrame to disk.
            Default is True.

        assign : bool
            Whether to assign the normalized DataFrame to an attribute
             - ``accessibility`` if ``data_type`` is "ATAC-seq,
             - ``binding`` if ``data_type`` is "ChIP-seq, or
             - ``expression`` if ``data_type`` is "RNA-seq.

        Returns
        -------
        pandas.DataFrame
            Annotated dataframe with requested sample attributes.

        Attributes
        ----------
        {accessibility, binding, expression} : pandas.DataFrame
            A pandas DataFrame with  MultiIndex column index containing the sample's attributes specified.
        """
        if (attributes is None) and hasattr(self, "sample_attributes"):
            _LOGGER.info("Using 'sample_attributes' from analysis to annotate matrix: {}"
                         .format(",".join(self.sample_attributes)))
            attributes = self.sample_attributes
        if (attributes is None) and hasattr(self, "prj"):
            _LOGGER.warn(
                "Analysis has no 'sample_attributes' set. " +
                "Will use all columns from project annotation sheet: {}"
                .format(",".join(self.prj.sheet.columns)))
            attributes = self.prj.sheet.columns
        if attributes is None:
            msg = "Attributes not given and could not be set from Analysis or Project."
            raise ValueError(msg)

        if self.data_type == "ATAC-seq":
            output_matrix = "accessibility"
            if quant_matrix is None:
                quant_matrix = "coverage_annotated"
        elif self.data_type == "ChIP-seq":
            output_matrix = "binding"
            if quant_matrix is None:
                quant_matrix = "coverage_annotated"
        elif self.data_type == "CNV":
            output_matrix = "cnv"
            if quant_matrix is None:
                if not isinstance(getattr(self, quant_matrix), pd.DataFrame):
                    _LOGGER.error("For CNV data type, the matrix to be annotated must be" +
                                  " directly passed to the function throught the `quant_matrix` argument!")
                    raise ValueError
            if quant_matrix is None:
                quant_matrix = "coverage_norm"
        elif self.data_type == "RNA-seq":
            output_matrix = "expression"
            if quant_matrix is None:
                quant_matrix = "expression_annotated"
        else:
            _LOGGER.warning("Data type of object not known, will not set as attribute.")
            assign = False
            output_matrix = ""
            if quant_matrix is None:
                msg = "Data type of object not known, must specify `quant_matrix` to annotate!"
                raise ValueError(msg)

        matrix = getattr(self, quant_matrix)

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
                os.path.join(
                    self.results_dir,
                    self.name + "{}.annotated_metadata.csv"
                        .format("." + output_matrix if output_matrix != "" else output_matrix)),
                index=True)
        if assign:
            setattr(self, output_matrix, df)
        return df

    def get_level_colors(
            self, index=None, matrix=None, levels=None,
            pallete="tab20", cmap="RdBu_r", nan_color=(0.662745, 0.662745, 0.662745, 1.0),
            # TODO: test dataframe return
            as_dataframe=False):
        """
        Get tuples of floats representing a colour for a sample in a given variable in a
        dataframe"s index (particularly useful with MultiIndex dataframes).

        If given, will use the provieded ``index`` argument, otherwise, the the columns
        and its levels of an attribute of self named ``matrix``.
        ``levels`` can be passed to subset the levels of the index.

        Will try to guess if each variable is categorical or numerical and return either colours
        from a colour ``pallete`` or a ``cmap``, respectively with null values set to ``nan_color``
        (a 4-value tuple of floats).

        Parameters
        ----------

        index : pandas.Index, optional
            Pandas Index to use.
            If not provided will use the column Index of the provided ``matrix``.

        matrix : str, optional
            Name of analysis attribute containing a dataframe with pandas.MultiIndex columns to use.
            If not provided will use the provided ``index``.

        levels : list, optional
            Levels of multiindex to restrict to.
            Defaults to all in index under use.

        pallete : str, optional
            Name of matplotlib color palete to use with categorical levels.
            See matplotlib.org/examples/color/colormaps_reference.html.
            Defaults to "tab20".

        cmap : str, optional
            Name of matplotlib color palete to use with numerical levels.
            See matplotlib.org/examples/color/colormaps_reference.html.
            Defaults to "RdBu_r".

        nan_color : tuple, optional
            Color for missing (i.e. NA) values.
            Defaults to ``(0.662745, 0.662745, 0.662745, 0.5)`` == ``grey``.

        as_dataframe : bool, optional
            Whether a dataframe should be return.
            Defaults to False.
            Not implemented yet.

        Returns
        -------
        {list, pandas.DataFrame}
            Matrix of shape (level, sample) with rgb values of each of the variable.
            If as_dataframe, this will be a pandas.DataFrame otherwise, list of lists.
        """
        if index is None:
            if matrix is None:
                msg = "One of `index` or `matrix` must be provided."
                _LOGGER.error(msg)
                raise ValueError(msg)
            index = getattr(self, matrix).columns

        if levels is not None:
            drop = [l.name for l in index.levels if l.name not in levels]
            index = index.droplevel(drop)

        # Handle special case of single level
        if not isinstance(index, pd.core.indexes.multi.MultiIndex):
            index = pd.MultiIndex.from_arrays([index.values], names=[index.name])

        _cmap = plt.get_cmap(cmap)
        _pallete = plt.get_cmap(pallete)

        colors = list()
        for level in index.levels:
            # For empty levels (all values nan), return nan colour
            if len(level) == 0:
                colors.append([nan_color] * len(index))
                _LOGGER.warning("Level {} has only NaN values.".format(level.name))
                continue
            # determine the type of data in each level
            # TODO: check this works in all cases
            most_common = Counter([type(x) for x in level]).most_common()[0][0]
            _LOGGER.debug("Getting colours for level '{}', which has type '{}'."
                          .format(level.name, most_common))

            # Add either colors based on categories or numerical scale
            if most_common in [int, float, np.float32, np.float64, np.int32, np.int64]:
                values = index.get_level_values(level.name)
                # Create a range of either 0-100 if only positive values are found
                # or symmetrically from the maximum absolute value found
                if not any(values.dropna() < 0):
                    norm = matplotlib.colors.Normalize(vmin=values.min(), vmax=values.max())
                else:
                    r = max(abs(values.min()), abs(values.max()))
                    norm = matplotlib.colors.Normalize(vmin=-r, vmax=r)

                col = _cmap(norm(values))
                # replace color for nan cases
                col[np.where(index.get_level_values(level.name).to_series().isnull().tolist())] = nan_color
                colors.append(col.tolist())
            else:
                n = len(set(index.get_level_values(level.name)))
                # get n equidistant colors
                p = [_pallete(1. * i / n) for i in range(n)]
                color_dict = dict(zip(list(set(index.get_level_values(level.name))), p))
                # color for nan cases
                color_dict[np.nan] = nan_color
                col = [color_dict[x] for x in index.get_level_values(level.name)]
                colors.append(col)

        if as_dataframe:
            colors = pd.DataFrame(
                colors, index=index.levels, columns=index).T

        return colors

    def unsupervised_analysis(
            self,
            steps=["correlation", "manifold", "pca", "pca_association"],
            data_type=None,
            quant_matrix=None,
            samples=None,
            attributes_to_plot=None,
            plot_prefix=None,
            standardize_matrix=False,
            manifold_algorithms=["MDS", "Isomap", "LocallyLinearEmbedding", "SpectralEmbedding", "TSNE"],
            display_corr_values=False,
            plot_max_pcs=8,
            prettier_sample_names=True,
            pallete="tab20",
            cmap="RdBu_r",
            rasterized=False,
            dpi=300,
            output_dir="{results_dir}/unsupervised_analysis_{data_type}",
            **kwargs):
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
        analysis : ngs_toolkit.general.Analysis
            Analysis object to perform analysis for.

        steps : list, optional
            List of step keywords to be performed as described above.
            Defaults to all available.

        data_type : str, optional
            Data type. One of "ATAC-seq" or "RNA-seq".
            Defaults to "ATAC-seq".

        quant_matrix : str, optional
            Name of analysis attribute contatining the numeric dataframe to perform analysis on.
            Defaults to the value of "norm_matrix_name" which is data-type specific.
            Must have a pandas.MultiIndex as column index.

        samples : list, optional
            List of sample objects to restrict analysis to.
            Defaults to all in analysis.

        attributes_to_plot : list, optional
            List of attributes shared between sample groups should be plotted.
            Defaults to attributes in analysis.group_attributes.

        plot_prefix : str, optional
            Prefix for output files.
            Defaults to "all_sites" if data_type is ATAC-seq and "all_genes" if data_type is RNA-seq.

        standardize_matrix : bool, optional
            Whether to standardize variables in `quant_matrix` by removing the mean and scaling to unit variance.

        manifold_algorithms : list, optional
            List of manifold algorithms to use. See available algorithms here:
            http://scikit-learn.org/stable/modules/classes.html#module-sklearn.manifold

        display_corr_values : bool, optional
            Whether values in heatmap of sample correlations should be
            displayed overlaid on top of colours. Defaults to False.

        prettier_sample_names : bool, optional
            Whether it should attempt to prettify sample names by removing the data type from plots.
            Defaults to True.

        pallete : str
            Color pallete to use in levels of `attributes_to_plot`. Will be passed to
            `analysis.get_level_colors`.

        cmap : str
            Color map to use in numerical levels of `attributes_to_plot`.
            Will be passed to `analysis.get_level_colors`.

        rasterized : bool, optional
            Whether elements with many objects should be rasterized.
            Defaults to False.

        dpi : int, optional
            Definition of rasterized image in dots per inch (dpi).
            Defaults to 300.

        output_dir : str, optional
            Directory for generated files and plots.
            Defaults to "{results_dir}/unsupervised_analysis_{data_type}".

        **kwargs: optional
            kwargs are passed to ngs_toolkit.graphics.plot_projection
        """
        data_type = self._get_data_type(data_type)

        if data_type == "ATAC-seq":
            if plot_prefix is None:
                plot_prefix = "all_sites"
            if quant_matrix is None:
                quant_matrix = "accessibility"
        elif data_type == "ChIP-seq":
            if plot_prefix is None:
                plot_prefix = "all_sites"
            if quant_matrix is None:
                quant_matrix = "binding"
        elif data_type == "CNV":
            if plot_prefix is None:
                plot_prefix = "all_bins"
            if quant_matrix is None:
                quant_matrix = "cnv"
        elif data_type == "RNA-seq":
            if plot_prefix is None:
                plot_prefix = "all_genes"
            if quant_matrix is None:
                quant_matrix = "expression"
        else:
            raise ValueError("Data types can only be 'ATAC-seq', 'RNA-seq' or 'CNV'.")

        if ("{results_dir}" in output_dir) and ("{data_type}" in output_dir):
            output_dir = output_dir.format(results_dir=self.results_dir, data_type=data_type)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        matrix = getattr(self, quant_matrix)

        if not isinstance(matrix.columns, pd.core.indexes.multi.MultiIndex):
            msg = "Provided quantification matrix must have columns with MultiIndex."
            hint = " Use ngs_toolkit.general.annotate_with_sample_metadata to do that."
            _LOGGER.error(msg + hint)
            raise TypeError(msg)

        if samples is None:
            samples = [s for s in self.samples
                       if s.name in matrix.columns.get_level_values("sample_name")]
        else:
            samples = [s for s in samples
                       if s.name in matrix.columns.get_level_values("sample_name")]
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
        # remove attributes with all NaNs
        attributes_to_plot = [attr for attr in attributes_to_plot
                              if attr in matrix.columns.names]
        attributes_to_plot = [attr for attr in attributes_to_plot
                              if not pd.isnull(matrix.columns.get_level_values(attr)).all()]
        if len(attributes_to_plot) == 0:
            msg = ("None of the factors in `attributes_to_plot` could be found in the " +
                   "quantification matrix index or they are all NaN.")
            _LOGGER.error(msg)
            raise ValueError(msg)

        # This will always be a matrix for all samples
        color_dataframe = pd.DataFrame(
            self.get_level_colors(
                index=matrix.columns, levels=attributes_to_plot,
                pallete=pallete, cmap=cmap),
            index=attributes_to_plot, columns=matrix.columns)
        # will be filtered now by the requested samples if needed
        color_dataframe = color_dataframe[[s.name for s in samples]]

        # All regions, matching samples (provided samples in matrix)
        X = matrix.loc[:, matrix.columns.get_level_values("sample_name").isin([s.name for s in samples])]

        if standardize_matrix:
            std = StandardScaler()
            X = pd.DataFrame(std.fit_transform(X.T).T, index=X.index, columns=X.columns)

        if isinstance(X.columns, pd.MultiIndex):
            sample_display_names = X.columns.get_level_values("sample_name")
        else:
            sample_display_names = X.columns
        # TODO: Re-implement to accomodate multiindex
        # if prettier_sample_names:
        #     X.columns = (
        #         color_dataframe.columns
        #         .str.replace("ATAC-seq_", "")
        #         .str.replace("RNA-seq_", "")
        #         .str.replace("ChIP-seq_", ""))

        if "correlation" in steps:
            # Pairwise correlations
            for method in ["pearson", "spearman"]:
                _LOGGER.info("Plotting pairwise correlation with '{}' metric.".format(method))
                g = sns.clustermap(
                    X.astype(float).corr(method),
                    xticklabels=False, yticklabels=sample_display_names, annot=display_corr_values,
                    cmap="Spectral_r", figsize=(0.2 * X.shape[1], 0.2 * X.shape[1]),
                    cbar_kws={"label": "{} correlation".format(method.capitalize())},
                    row_colors=color_dataframe.T, col_colors=color_dataframe.T)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
                g.ax_heatmap.set_xlabel(None, visible=False)
                g.ax_heatmap.set_ylabel(None, visible=False)
                savefig(g.fig, os.path.join(
                    output_dir, "{}.{}.{}_correlation.clustermap.svg"
                    .format(self.name, plot_prefix, method)))

        if "manifold" in steps:
            # Manifolds
            # TODO: test usage of non default manifolds from sklearn
            params = defaultdict(dict)
            params.update({
                "MDS": {"n_jobs": -1},
                "Isomap": {"n_jobs": -1},
                "LocallyLinearEmbedding": {},
                "SpectralEmbedding": {"n_jobs": -1},
                "TSNE": {"init": "pca"},
            })
            for algo in manifold_algorithms:
                msg = "Learning manifold with '{}' algorithm".format(algo)
                _LOGGER.info(msg + ".")

                manif = getattr(manifold, algo)(**params[algo])
                try:
                    x_new = manif.fit_transform(X.T)
                except (TypeError, ValueError):
                    hint = " Number of samples might be too small to perform '{}'".format(algo)
                    _LOGGER.error(msg + " failed!" + hint)
                    continue

                x_new = pd.DataFrame(x_new, index=X.columns, columns=list(range(x_new.shape[1])))

                _LOGGER.info("Plotting projection of manifold with '{}' algorithm.".format(algo))
                plot_projection(
                    df=x_new, color_dataframe=color_dataframe, dims=1,
                    output_file=os.path.join(
                        output_dir, "{}.{}.{}.svg"
                        .format(self.name, plot_prefix, algo.lower())),
                    attributes_to_plot=attributes_to_plot, **kwargs)

        if "pca" in steps:
            # PCA
            pcs = min(*X.shape) - 1
            _LOGGER.info("Decomposing data with 'PCA' algorithm for {} dimensions.".format(pcs))
            pca = PCA(n_components=pcs, svd_solver="arpack")
            x_new = pca.fit_transform(X.T)

            pcs_order = range(pca.n_components_)
            x_new = pd.DataFrame(x_new, index=X.columns, columns=pcs_order)
            x_new.to_csv(
                os.path.join(output_dir, "{}.{}.pca.fit.csv".format(
                    self.name, plot_prefix)))
            comps = pd.DataFrame(pca.components_.T, index=X.index, columns=pcs_order)
            comps.to_csv(
                os.path.join(output_dir, "{}.{}.pca.loadings.csv".format(
                    self.name, plot_prefix)))

            # Write % variance expained to disk
            variance = pd.Series(
                pca.explained_variance_ratio_ * 100,
                name="percent_variance", index=pcs_order
                ).to_frame()
            variance["log_variance"] = np.log10(pca.explained_variance_)
            variance.index.name = "PC"
            variance.to_csv(
                os.path.join(output_dir, "{}.{}.pca.explained_variance.csv".format(
                    self.name, plot_prefix)))

            # plot % explained variance per PC
            _LOGGER.info("Plotting variance explained with PCA.")
            fig, axis = plt.subplots(1, 3, figsize=(4 * 3, 4))
            axis[0].plot(variance.index, variance["percent_variance"], "o-")
            axis[0].set_ylim((0, variance["percent_variance"].max() + variance["percent_variance"].max() * 0.1))
            axis[1].plot(variance.index, variance["log_variance"], "o-")
            axis[2].plot(variance.index, variance["percent_variance"].cumsum(), "o-")
            axis[2].set_ylim((0, 100))
            for ax in axis:
                ax.axvline(len(attributes_to_plot), linestyle="--")
                ax.set_xlabel("PC")
            axis[0].set_ylabel("% variance")
            axis[1].set_ylabel("log variance")
            axis[2].set_ylabel("Cumulative % variance")
            sns.despine(fig)
            savefig(fig, os.path.join(
                output_dir, "{}.{}.pca.explained_variance.svg"
                            .format(self.name, plot_prefix)))

            # plot pca
            pcs = min(x_new.shape[1] - 1, plot_max_pcs)
            _LOGGER.info("Plotting PCA up to {} dimensions.".format(pcs))
            plot_projection(
                df=x_new, color_dataframe=color_dataframe, dims=pcs,
                output_file=os.path.join(output_dir, "{}.{}.pca.svg".format(
                                         self.name, plot_prefix)),
                attributes_to_plot=attributes_to_plot, **kwargs)

        if ("pca" in steps) and ("pca_association" in steps):
            # Test association of PCs with attributes
            _LOGGER.info("Computing association of given attributes with principal components.")
            associations = list()
            for pc in pcs_order:
                for attr in attributes_to_plot:
                    _LOGGER.debug("PC {}; Attribute {}.".format(pc + 1, attr))

                    # Get all values of samples for this attr
                    groups = x_new.index.get_level_values(attr)

                    # Determine if attr is categorical or continuous
                    if all([isinstance(i, (str, bool)) for i in groups]) or len(groups) == 2:
                        variable_type = "categorical"
                    elif all([isinstance(i, (int, float, np.int64, np.float64)) for i in groups]):
                        variable_type = "numerical"
                    else:
                        _LOGGER.warning("attr %s cannot be tested." % attr)
                        associations.append([pc + 1, attr, variable_type, np.nan, np.nan, np.nan])
                        continue

                    if variable_type == "categorical":
                        # It categorical, test pairwise combinations of attributes
                        for group1, group2 in itertools.combinations(groups, 2):
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

            associations = pd.DataFrame(
                associations, columns=["pc", "attribute", "variable_type", "group_1", "group_2", "p_value"])

            if associations.empty:
                msg = "Couldn't test any associations between PCs and factors."
                hint = " Perhaps PCA produced only 1 PC?"
                _LOGGER.warning(msg + hint)
                return

            # correct p-values
            associations.loc[:, "adj_pvalue"] = multipletests(associations["p_value"], method="fdr_bh")[1]

            # write
            _LOGGER.info("Saving associations.")
            associations.to_csv(os.path.join(
                output_dir, "{}.{}.pca.variable_principle_components_association.csv"
                            .format(self.name, plot_prefix)), index=False)

            if len(attributes_to_plot) < 2:
                _LOGGER.info("Only one attribute given, can't plot associations.")
                return

            # Plot
            for var in ["p_value", "adj_pvalue"]:
                pivot = (
                    associations
                    .groupby(["pc", "attribute"])
                    .min()[var]
                    .reset_index()
                    .pivot(index="pc", columns="attribute", values=var)
                    .dropna(axis=1))

                # heatmap of -log p-values
                g = sns.clustermap(
                    -np.log10(pivot), row_cluster=False,
                    annot=True, cbar_kws={"label": "-log10(p_value) of association"},
                    square=True, rasterized=rasterized, vmin=0)
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
                savefig(g.fig, os.path.join(
                    output_dir, "{}.{}.pca.variable_principle_components_association.{}.svg"
                                .format(self.name, plot_prefix, var)))

                # heatmap of masked significant
                g = sns.clustermap(
                    (pivot < 0.05).astype(int),
                    row_cluster=False, cbar_kws={"label": "significant association"},
                    square=True, rasterized=rasterized, vmin=0, vmax=1, cmap="Paired")
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=45, ha="right")
                savefig(g.fig, os.path.join(
                    output_dir, "{}.{}.pca.variable_principle_components_association.{}.masked.svg"
                                .format(self.name, plot_prefix, var)))

    def differential_analysis(
            self,
            comparison_table=None,
            data_type=None,
            samples=None,
            covariates=None,
            output_dir="{results_dir}/differential_analysis_{data_type}",
            output_prefix="differential_analysis",
            alpha=0.05,
            overwrite=True,
            distributed=False,
            cpus=2,
            memory=16000):
        """
        Perform differential regions/genes across samples that are associated with a certain trait.
        Currently the only implementation is with DESeq2.
        This implies the rpy2 library and the respective R library are installed.

        Requires the R package "DESeq2" to be installed:
            >>> source("http://bioconductor.org/biocLite.R")
            >>> biocLite("DESeq2")

        For other implementations of differential analysis see `ngs_toolkit.general.least_squares_fit`
        and `ngs_toolkit.general.differential_from_bivariate_fit`.

        Parameters
        ----------
        comparison_table : pandas.DataFrame
            A dataframe with "comparison_name", "comparison_side" and "sample_name", "sample_group" columns.
            Defaults to the analysis' own "comparison_table" attribute.

        data_type : str, optional
            Type of data under analysis. One of "ATAC-seq" or "RNA-seq".
            Defaults to analysis' own data_type.

        samples : list, optional
            Samples to limit analysis to.
            Defaults to all samples in analysis object.

        covariates : list, optional
            Additional variables to take into account in the model fitting.
            Defaults to None.

        output_dir : str, optional
            Output directory for analysis.
            Defaults to "{results_dir}/differential_analysis_{data_type}".
            If containing "{data_type}", will format string with variable.

        output_prefix : str, optional
            Prefix for output files.
            Defaults to "differential_analysis".

        alpha : float, optional
            Significance level to use in differential analysis.
            Results for all features will be returned nonetheless. Defaults to 0.05.

        overwrite : bool, optional
            Whether results should be overwritten in case they already exist.
            Defaults to True.

        distributed : bool, optional
            Whether analysis should be distributed in a computing cluster for each comparison.
            Currently, only a SLURM implementation is available.
            If `True`, will not return results.
            Defaults to False.

        cpus : int, optional
            Number of CPUS to use when using distributed jobs.
            Default: 2.

        memory : int, optional
            Memory to use when using distributed jobs. Default: 16000 (16Gb).

        Returns
        -------
        pandas.DataFrame
            Results for all comparisons.
            Will be `None` if `distributed` is `True`.

        Attributes
        ----------
        differential_results : pandas.DataFrame
            Pandas dataframe with results.
        """
        if comparison_table is None:
            msg = "`comparison_table` was not given and is not set in analysis object."
            hint = "Add a `comparison_table` attribute to the analysis object."
            try:
                comparison_table = self.comparison_table
            except AttributeError as e:
                _LOGGER.error(msg)
                _LOGGER.info(hint)
                raise e

        data_type = self._get_data_type(data_type)

        # Check comparisons
        # check comparison table has required columns
        req_attrs = ["comparison_name", "comparison_side", "sample_name", "sample_group"]
        if not all([x in comparison_table.columns for x in req_attrs]):
            raise AssertionError("Given comparison table does not have all of '{}' columns."
                                 .format("".join(req_attrs)))
        # check all comparisons have samples in two sides
        if not all(comparison_table.groupby("comparison_name")["comparison_side"].nunique() == 2):
            msg = "All comparisons must have samples in each side of the comparison."
            raise AssertionError(msg)
        # check if any comparison and sample group has samples disagreeing in their side
        if not all(comparison_table.groupby(["comparison_name", "sample_group"])["comparison_side"].nunique() == 1):
            msg = "Samples in same comparison and group must agree on their side of the comparison."
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
        if data_type == "ATAC-seq":
            count_matrix = self.coverage
        elif data_type == "RNA-seq":
            count_matrix = self.expression_matrix_counts
        else:
            msg = "Differential analysis is only implemented for data types 'ATAC-seq' or 'RNA-seq'."
            raise ValueError(msg)

        if samples is None:
            samples = self.samples
        samples = [s for s in samples if
                   (s.name in comparison_table["sample_name"].tolist()) &
                   (s.name in count_matrix.columns)]
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
            comparison_table = comparison_table.set_index("sample_name").join(sample_table.set_index("sample_name")[covariates]).reset_index()

        # Make table for DESeq2
        experiment_matrix = comparison_table[
            ["sample_name", "sample_group"] + (covariates if covariates is not None else [])
        ].drop_duplicates()

        # Make formula for DESeq2
        formula = "~ {}sample_group".format(" + ".join(covariates) + " + " if covariates is not None else "")

        # Run DESeq2 analysis
        if not distributed:
            # TODO: for complex designs (one sample is in multiple groups/comparisons)
            #       implement running one comparison after the other
            results = deseq_analysis(
                count_matrix, experiment_matrix, comparison_table,
                formula, output_dir, output_prefix, alpha=alpha, overwrite=overwrite)
            try:
                results = results.set_index("index")
            except KeyError:
                pass
            _LOGGER.info("Setting results of differential analysis to a variable 'differential_results'.")
            self.differential_results = results
            return results

        else:
            tk = NGSTk()
            for comparison_name in comparison_table["comparison_name"].drop_duplicates():
                # make directory for comparison input/output
                out = os.path.join(os.path.abspath(output_dir), comparison_name)
                if not os.path.exists(out):
                    os.makedirs(out)

                comp = comparison_table[comparison_table["comparison_name"] == comparison_name]
                comp.to_csv(os.path.join(out, "comparison_table.csv"), index=False)

                exp = experiment_matrix[
                    experiment_matrix["sample_name"].isin(comp["sample_name"].tolist()) &
                    experiment_matrix["sample_group"].isin(comp["sample_group"].tolist())]
                exp.to_csv(os.path.join(out, "experiment_matrix.csv"), index=False)

                count = count_matrix[comp["sample_name"].drop_duplicates()]
                count.to_csv(os.path.join(out, "count_matrix.csv"), index=True)

                job_name = "deseq_job.{}".format(comparison_name)
                log_file = os.path.join(out, job_name + ".log")
                job_file = os.path.join(out, job_name + ".sh")

                cmd = tk.slurm_header(
                    job_name=job_name, output=log_file,
                    cpus_per_task=2, mem_per_cpu=16000)
                # TODO: add DESeq2 script to toolkit and make path configurable
                cmd += "python -u ~/deseq_parallel.py"
                cmd += " --output_prefix {}".format(output_prefix)
                cmd += " --formula '{}'".format(formula)
                cmd += " --alpha {}".format(alpha)
                if overwrite:
                    cmd += " --overwrite"
                cmd += " {}\n".format(out)
                cmd += tk.slurm_footer()

                with open(job_file, "w") as handle:
                    handle.write(textwrap.dedent(cmd).replace("\n ", "\n").replace("  ", ""))

                tk.slurm_submit_job(job_file)

    def collect_differential_analysis(
            self,
            comparison_table=None,
            data_type="ATAC-seq",
            output_dir="{results_dir}/differential_analysis_{data_type}",
            output_prefix="differential_analysis",
            permissive=True,
            assign=True, save=True,
            overwrite=False):
        """
        Collect results from DESeq2 differential analysis.
        Particularly useful when runing `differential_analysis` with in distributed mode.

        Parameters
        ----------
        comparison_table : pandas.DataFrame
            A dataframe with "comparison_name", "comparison_side" and "sample_name", "sample_group" columns.
            Defaults to the analysis's own "comparison_table" attribute.

        data_type : str, optional
            Type of data under analysis. One of "ATAC-seq" or "RNA-seq".
            Defaults to analysis' own data_type.

        output_dir : str, optional
            Output directory for analysis.
            Defaults to "{results_dir}/differential_analysis_{data_type}".
            If containing "{data_type}", will format string with variable.

        output_prefix : str, optional
            Prefix for output files.
            Defaults to "differential_analysis".

        permissive : bool, optional
            Whether non-existing files should be skipped or an error be thrown.
            Defaults to True.

        assign : bool, optional
            Whether to add results to a `differential_results` attribute.
            Defaults to True.

        save : bool, optional
            Whether to save results to disk.
            Defaults to True.

        overwrite : bool, optional
            Whether results should be overwritten in case they already exist.
            Defaults to False.

        Returns
        -------
        pandas.DataFrame
            Results for all comparisons.
            Will be `None` if `overwrite` is `False` and a results file already exists.

        Attributes
        ----------
        differential_results : pandas.DataFrame
            Pandas dataframe with results.
        """
        # TODO: Add "input_dir" and input_prefix"
        if comparison_table is None:
            msg = "`comparison_table` was not given and is not set in analysis object."
            hint = "Add a `comparison_table` attribute to the analysis object."
            try:
                comparison_table = self.comparison_table
            except AttributeError as e:
                _LOGGER.error(msg)
                _LOGGER.info(hint)
                raise e

        # Make output dir
        output_dir = self._format_string_with_attributes(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        results_file = os.path.join(output_dir, output_prefix + ".deseq_result.all_comparisons.csv")
        if not overwrite and os.path.exists(results_file):
            msg = "Differential analysis results '{}' already exist and argument `overwrite` is False."
            hint = " Will not do anything."
            _LOGGER.warning(msg.format(results_file) + hint)
            return

        comps = comparison_table["comparison_name"].drop_duplicates().sort_values()
        results = pd.DataFrame()
        for comp in tqdm(comps, total=len(comps), desc="Comparison"):
            out_file = os.path.join(output_dir, comp, output_prefix + ".deseq_result.{}.csv".format(comp))
            # print("Collecting comparison '{}'".format(comp))
            # read
            try:
                res2 = pd.read_csv(out_file, index_col=0)
            except IOError as e:
                if permissive:
                    _LOGGER.warning("Results file for comparison '{}' do not exist. Skipping.".format(comp))
                    continue
                else:
                    raise e
            # append
            results = results.append(res2.reset_index(), ignore_index=True)

        if save:
            # save all
            results.to_csv(results_file, index=False)

        # set index
        if "index" in results.columns:
            results = results.set_index("index")

        if assign:
            self.differential_results = results
        return results

    def plot_differential(
            self,
            results=None,
            comparison_table=None,
            samples=None,
            matrix=None,
            only_comparison_samples=False,
            data_type=None,
            alpha=0.05,
            corrected_p_value=True,
            fold_change=None,
            diff_based_on_rank=False,
            max_rank=1000,
            ranking_variable="padj",
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
            group_wise_colours=False,
            group_variables=None,
            pallete="tab20",
            cmap="RdBu_r"
            ):
        """
        Plot differential features (e.g. chromatin region, genes) discovered with supervised
        group comparisons by ``ngs_toolkit.general.differential_analysis``.
        This will plot number and direction of discovered features, scatter, MA and volcano
        plots for each comparison and joint heatmaps of log fold changes, normalized values
        or Z-scores of individual samples or groups in the differential features.

        Parameters
        ----------
        results : pandas.DataFrame
            Data frame with differential analysis results.
            See ``ngs_toolkit.general.differential_analysis`` for more information.

        comparison_table : pandas.DataFrame, optional
            Comparison table. If provided, group-wise plots will be produced.
            Defaults to the analysis' "comparison_table" attribute.

        samples : list, optional
            List of sample objects to restrict analysis to.
            Defaults to all samples in analysis.

        matrix : str, optional
            Matrix of quantification to use for plotting feature values across samples/groups.
           Defaults to either "accessibility" for ATAC-seq analysis or "expression" for RNA-seq.

        only_comparison_samples : bool, optional
            Whether to use only samples present in the `comparison_table`.
            Defaults to False.

        data_type : str, optional
            The data type being analyzed. Currently supported are "ATAC-seq" or "RNA-seq".
            Defaults to the analysis' own data_type.

        alpha : float, optional
            Significance level to consider a feature differential.
            Defaults to 0.05.

        corrected_p_value : bool, optional
            Whether to use a corrected p-valueto consider a feature differential.
            Defaults to True.

        fold_change : float, optional
            Effect size (log2 fold change) to consider a feature differential. Considers absolute values.
            Default is no log2 fold change threshold.

        diff_based_on_rank : bool, optional
            Whether a feature should be considered differential based on its rank.
            Defaults to False.

        max_rank : int, optional
            Rank to use when using `diff_based_on_rank`.
            Defaults to 1000.

        ranking_variable : str, optional
            Which variable to use for ranking when using `diff_based_on_rank`.
            Defaults to "padj".

        respect_stat_thresholds : bool, optional
            Whether the statistical thresholds from `alpha` and `fold_change` should still be
            respected when using `diff_based_on_rank`.
            Defaults to True

        output_dir : str, optional
            Directory to create output files.
            Defaults to "{results_dir}/differential_analysis_{data_type}"

        output_prefix : str, optional
            Prefix to use when creating output files.
            Defaults to "differential_analysis".

        plot_each_comparison : bool, optional
            Whether each comparison should be plotted in scatter, MA and volcano plots.
            Useful to turn off with many comparisons.
            Defaults to True.

        mean_column : str, optional
            Column  in `results` data frame containing values for mean values across samples.
            Defaults to "baseMean".

        log_fold_change_column : str, optional
            Column in `results` data frame containing values for log2FoldChange values across samples.
            Defaults to "log2FoldChange".

        p_value_column : str, optional
            Column  in `results` data frame containing values for p-values across samples.
            Defaults to "pvalue".

        adjusted_p_value_column : str, optional
            Column  in `results` data frame containing values for adjusted p-values across samples.
            Defaults to "padj".

        comparison_column : str, optional
            Column  in `results` data frame containing the name of the comparison.
            Defaults to "comparison_name".

        rasterized : bool, optional
            Whether plots with many objects should be rasterized.
            Defaults to True.

        robust : bool, optional
            Whether heatmap color scale ranges should be robust (using quantiles) rather than extreme values.
            Useful for noisy/extreme data.
            Defaults to False.

        feature_labels : bool, optional
            Whether features (regions/genes) should be labeled in heatmaps.
            Defaults to False.

        group_wise_colours : bool, optional
            Whether groups of samples should be coloured in heatmaps.
            Defaults to False.

        group_variables : list, optional
            Which variables to colour if `group_wise_colours` if True.
            Defaults to None (must be given).

        pallete : str
            Color pallete to use in levels of `group_variables`.
            Will be passed to `analysis.get_level_colors`.

        cmap : str
            Color map to use in numerical levels of `group_variables`.
            Will be passed to `analysis.get_level_colors`.
        """
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

        req_attrs = [mean_column, log_fold_change_column, p_value_column, adjusted_p_value_column, comparison_column]
        if not all([x in results.columns for x in req_attrs]):
            raise AssertionError("Results dataframe must have '{}' columns.".format(", ".join(req_attrs)))

        if data_type is None:
            msg = "Data type not defined and Analysis object does not have a `data_type` attribute."
            try:
                data_type = self.data_type
            except AttributeError as e:
                _LOGGER.error(msg)
                raise e
            if data_type is None:
                _LOGGER.error(msg)
                raise ValueError

        # Make output dir
        output_dir = self._format_string_with_attributes(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Get matrix and samples
        results = results.copy()
        if data_type == "ATAC-seq":
            if matrix is None:
                matrix = self.accessibility
            results.index.name = "region"
            var_name = "region"
            quantity = "Accessibility"
            unit = "RPM"
        elif data_type == "RNA-seq":
            if matrix is None:
                matrix = self.expression
            results.index.name = "gene_name"
            var_name = "gene"
            quantity = "Expression"
            unit = "TPM"
        else:
            msg = "Plot differential is only implemented for data types 'ATAC-seq' or 'RNA-seq'."
            raise AssertionError(msg)

        if samples is None:
            samples = self.samples
        samples = [s for s in samples if s.name in matrix.columns]

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
            samples = [s for s in samples if s.name in comparison_table["sample_name"].tolist()]
        matrix = matrix[[s.name for s in samples]]

        # Handle group colouring
        if group_wise_colours:
            if group_variables is None:
                msg = "If `group_wise_colours` is True, a list of `group_variables` must be passed."
                raise AssertionError(msg)

            # This will always be a matrix for all samples
            color_dataframe = pd.DataFrame(
                self.get_level_colors(
                    index=matrix.columns, levels=group_variables,
                    pallete=pallete, cmap=cmap),
                index=group_variables, columns=matrix.columns)
            # will be filtered now by the requested samples if needed
            color_dataframe = color_dataframe[[s.name for s in samples]]
            color_dataframe = color_dataframe.loc[:, matrix.columns]

        # Extract significant based on p-value and fold-change
        if fold_change is not None:
            fc = (results[log_fold_change_column].abs() > fold_change)
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
                    i = results.loc[results[comparison_column] == comparison, ranking_variable].abs().sort_values().tail(max_rank).index
                else:
                    i = results.loc[results[comparison_column] == comparison, ranking_variable].sort_values().head(max_rank).index
                results.loc[(results[comparison_column] == comparison) & results.index.isin(i), "diff_rank"] = True
            results.loc[:, "diff_rank"] = results.loc[:, "diff_rank"].fillna(False)
            if respect_stat_thresholds:
                results.loc[:, "diff"] = (results.loc[:, "diff"].isin([True])) & (results.loc[:, "diff_rank"].isin([True]))
            else:
                results.loc[:, "diff"] = results.loc[:, "diff_rank"]

        # Annotate direction of change
        results.loc[:, "direction"] = results.loc[
            :, log_fold_change_column].apply(lambda x: "up" if x >= 0 else "down")

        # PLOTS
        _LOGGER.info("Starting to generate plots for differential comparisons.")
        comparisons = sorted(results[comparison_column].drop_duplicates())
        n_side = int(np.ceil(np.sqrt(len(comparisons))))

        # P-value and Fold-change distributions
        for variable, label, axvline in [
                (p_value_column, "P-value", False),
                (adjusted_p_value_column, "Adjusted p-value", False),
                (log_fold_change_column, "log2(fold-change)", True)]:
            _LOGGER.info("Plotting distribution of {}.".format(label))

            # log fold-changes distributions
            fig, axis = plt.subplots(1, 1, figsize=(4, 4))
            sns.distplot(results[variable].dropna(), kde=False, ax=axis)
            if axvline:
                axis.axvline(0, color="black", alpha=0.5)
            axis.set_xlabel(label)
            axis.set_ylabel(var_name.capitalize() + "s (frequency)")
            sns.despine(fig)
            savefig(fig, os.path.join(output_dir, output_prefix + "." + variable + ".distribution.svg"))

            if plot_each_comparison:
                # per comparison
                g = sns.FacetGrid(data=results, col=comparison_column, col_wrap=n_side)
                g.map(sns.distplot, variable, kde=False)
                for ax in g.axes:
                    ax.set_yscale("log")
                    if axvline:
                        ax.axvline(0, color="black", alpha=0.5)
                    ax.set_xlabel(label)
                    ax.set_ylabel(var_name.capitalize() + "s (frequency)")
                sns.despine(g.fig)
                savefig(g.fig, os.path.join(output_dir, output_prefix + "." + variable + ".distribution.per_comparison.svg"))

        # Number of differential vars
        _LOGGER.info("Calculating number of differential {}s per comparison.".format(var_name))
        n_vars = float(matrix.shape[0])
        total_diff = results.groupby(
            [comparison_column])["diff"].sum().sort_values(ascending=False).reset_index()
        split_diff = results.groupby(
            [comparison_column, "direction"])["diff"].sum().sort_values(ascending=False).reset_index()
        split_diff.loc[split_diff["direction"] == "down", "diff"] *= -1
        split_diff["label"] = split_diff[comparison_column].astype(str) + ", " + split_diff["direction"]
        total_diff["diff_perc"] = (total_diff["diff"] / n_vars) * 100
        split_diff["diff_perc"] = (split_diff["diff"] / n_vars) * 100

        _LOGGER.info("Plotting number of differential {}s per comparison.".format(var_name))
        fig, axis = plt.subplots(2, 2, figsize=(4 * 2, 4 * 2))
        sns.barplot(data=total_diff, x="diff", y=comparison_column, orient="h", ax=axis[0, 0])
        sns.barplot(data=total_diff, x="diff_perc", y=comparison_column, orient="h", ax=axis[0, 1])
        sns.barplot(data=split_diff, x="diff", y=comparison_column, hue="direction", dodge=False, orient="h", ax=axis[1, 0])
        sns.barplot(data=split_diff, x="diff_perc", y=comparison_column, hue="direction", dodge=False, orient="h", ax=axis[1, 1])
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
        savefig(fig, os.path.join(output_dir, output_prefix + ".number_differential.directional.svg"))

        if plot_each_comparison:
            _LOGGER.info("Doing detailed plotting per comparison:")

            # Add same colour scale to all plots/comparisons
            smallest_p_value = -np.log10(np.nanpercentile(results[p_value_column], 1e-5))
            if smallest_p_value in [np.inf, np.nan]:
                smallest_p_value = 300
            _LOGGER.debug("Maximum -log10(p-value) across comparisons is {}".format(smallest_p_value))
            pval_cmap = "Reds"

            # Pairwise scatter plots
            if comparison_table is not None:
                _LOGGER.info(
                    "Plotting scatter of {} distribution for each group in each comparison.".format(var_name))
                fig, axes = plt.subplots(
                    n_side, n_side,
                    figsize=(n_side * 4, n_side * 4), sharex=True, sharey=True)
                if n_side > 1 or n_side > 1:
                    axes = iter(axes.flatten())
                else:
                    axes = iter([axes])
                for comparison in comparisons:
                    _LOGGER.debug("Comparison '{}'...".format(comparison))
                    c = comparison_table.loc[comparison_table[comparison_column] == comparison, :]
                    a = c.loc[c["comparison_side"] >= 1, "sample_name"]
                    b = c.loc[c["comparison_side"] <= 0, "sample_name"]

                    a = matrix.loc[:, [s.name for s in samples if s.name in a.tolist() and s.library == data_type]].mean(axis=1)
                    b = matrix.loc[:, [s.name for s in samples if s.name in b.tolist() and s.library == data_type]].mean(axis=1)

                    # Hexbin plot
                    ax = next(axes)
                    ax.hexbin(
                        b, a,
                        alpha=0.85, cmap="Greys", color="black", edgecolors="white",
                        linewidths=0, bins="log", mincnt=1, rasterized=True)

                    # Scatter for significant features
                    diff_vars = results.loc[
                        (results[comparison_column] == comparison) &
                        (results["diff"].isin([True])), :]
                    if diff_vars.shape[0] > 0:
                        # get color vector based on p-value
                        col = -np.log10(results.loc[
                            (results[comparison_column] == comparison) &
                            (results["diff"].isin([True])), p_value_column].squeeze())
                        _LOGGER.debug("Shapes: {} {} {}".format(a.shape, b.shape, diff_vars.shape))
                        # in case there's just one significant feature:
                        if isinstance(col, np.float_):
                            col = np.array([col])
                        collection = ax.scatter(
                            b.loc[diff_vars.index],
                            a.loc[diff_vars.index],
                            alpha=0.5, s=2, c=col, cmap=pval_cmap,
                            vmin=0, vmax=smallest_p_value)
                        add_colorbar_to_axis(collection, label="-log10(p-value)")
                    ax.set_title(comparison)
                    # Name groups
                    xl = c.loc[c["comparison_side"] <= 0, "sample_group"].drop_duplicates().squeeze()
                    yl = c.loc[c["comparison_side"] >= 1, "sample_group"].drop_duplicates().squeeze()
                    if not (isinstance(xl, str) and isinstance(yl, str)):
                        xl = "Down-regulated"
                        yl = "Up-regulated"
                    ax.set_xlabel(xl)
                    ax.set_ylabel(yl)

                    # x = y square
                    lims = [
                        np.min([ax.get_xlim(), ax.get_ylim()]),
                        np.max([ax.get_xlim(), ax.get_ylim()])]
                    ax.plot(
                        lims, lims,
                        linestyle="--", alpha=0.5, zorder=0, color="black")
                    ax.set_aspect("equal")
                    ax.set_xlim(lims)
                    ax.set_ylim(lims)
                for ax in axes:
                    ax.set_visible(False)
                sns.despine(fig)
                savefig(fig,  os.path.join(output_dir, output_prefix + ".scatter_plots.svg"))

            # Volcano plots
            _LOGGER.info("Plotting volcano plots for each comparison.")
            fig, axes = plt.subplots(n_side, n_side, figsize=(n_side * 4, n_side * 4), sharex=False, sharey=False)
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
                    t[log_fold_change_column], -np.log10(t[p_value_column]),
                    alpha=0.85, cmap="Greys", color="black", edgecolors="white",
                    linewidths=0, bins="log", mincnt=1, rasterized=True)

                # Scatter for significant
                diff_vars = t.loc[t["diff"].isin([True]), :]
                if diff_vars.shape[0] > 0:
                    collection = ax.scatter(
                        t.loc[diff_vars.index, log_fold_change_column],
                        -np.log10(t.loc[diff_vars.index, p_value_column]),
                        alpha=0.5, s=2, c=-np.log10(t.loc[diff_vars.index, p_value_column]), cmap=pval_cmap,
                        vmin=0, vmax=smallest_p_value)
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
                    linestyle="--", alpha=0.5, zorder=0, color="black")
                if fold_change is not None:
                    ax.axvline(-fold_change, linestyle="--", alpha=0.5, zorder=0, color="black")
                    ax.axvline(fold_change, linestyle="--", alpha=0.5, zorder=0, color="black")
            for ax in axes:
                ax.set_visible(False)
            sns.despine(fig)
            savefig(fig, os.path.join(output_dir, output_prefix + ".volcano_plots.svg"))

            # MA plots
            _LOGGER.info("Plotting MA plots for each comparison.")
            fig, axes = plt.subplots(n_side, n_side, figsize=(n_side * 4, n_side * 4), sharex=False, sharey=False)
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
                    np.log10(t[mean_column]), t[log_fold_change_column],
                    alpha=0.85, cmap="Greys", color="black", edgecolors="white",
                    linewidths=0, bins="log", mincnt=1, rasterized=True)

                # Scatter for significant
                diff_vars = t.loc[t["diff"].isin([True]), :]
                if diff_vars.shape[0] > 0:
                    collection = ax.scatter(
                        np.log10(t.loc[diff_vars.index, mean_column]),
                        t.loc[diff_vars.index, log_fold_change_column],
                        alpha=0.5, s=2, c=-np.log10(t.loc[diff_vars.index, p_value_column]), cmap=pval_cmap,
                        vmin=0, vmax=smallest_p_value)
                    add_colorbar_to_axis(collection, label="-log10(p-value)")
                ax.set_title(comparison)
                ax.set_xlabel("Mean {}".format(quantity.lower()))
                ax.set_ylabel("log2(fold-change)")
                ax.axhline(0, linestyle="--", alpha=0.5, zorder=0, color="black")
                ll = np.max([abs(ii) for ii in ax.get_ylim()])
                ax.set_ylim(-ll, ll)

                # Add lines of significance
                if fold_change is not None:
                    ax.axhline(-fold_change, linestyle="--", alpha=0.5, zorder=0, color="black")
                    ax.axhline(fold_change, linestyle="--", alpha=0.5, zorder=0, color="black")
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

        if comparison_table is not None:
            _LOGGER.info("A comparison table was given, will try to plot values per sample group.")
            if results[comparison_column].drop_duplicates().shape[0] > 1:
                _LOGGER.info("Getting per-group values for each comparison.")
                groups = pd.DataFrame()
                for sample_group in comparison_table["sample_group"].drop_duplicates():
                    c = comparison_table.loc[
                        comparison_table["sample_group"] == sample_group, "sample_name"].drop_duplicates()
                    if c.shape[0] > 0:
                        groups.loc[:, sample_group] = matrix[[d for d in c if d in sample_cols]].mean(axis=1)

                if groups.empty:
                    # It seems comparisons were not done in a all-versus-all fashion
                    for group in comparison_table["sample_group"].drop_duplicates():
                        c = comparison_table.loc[comparison_table["sample_group"] == group, "sample_name"].drop_duplicates()
                        if c.shape[0] > 0:
                            groups.loc[:, group] = matrix[c].mean(axis=1)

                # Select only differential regions from groups
                groups = groups.loc[all_diff, :].sort_index(axis=1)

                n = groups.isnull().sum().sum()
                if n > 0:
                    _LOGGER.warning(
                        "{} {}s (across all comparisons) were not found in quantification matrix!".format(n, var_name))
                    m = groups.columns[groups.isnull().sum() == groups.shape[0]]
                    if len(m) > 0:
                        _LOGGER.warning(
                            "{} comparison groups were not found in quantification matrix: '{}'!"
                            .format(len(m), ", ".join(m)) +
                            " Proceeding without those.")
                        groups = groups.loc[:, ~groups.columns.isin(m)]
                    f = groups.index[groups.isnull().sum(1) == groups.shape[1]]
                    if len(f) > 0:
                        _LOGGER.warning(
                            "{} {}s were not found in quantification matrix!"
                            .format(len(m), var_name) +
                            " Proceeding without those.")
                        groups = groups.dropna()
                    n = groups.isnull().sum().sum()
                    if n != 0:
                        _LOGGER.error(
                            "{} {}s (across all comparisons) still have NaNs. Cannot proceed!".format(n, var_name))
                    else:
                        _LOGGER.info("Plotting clustered heatmaps of sample groups in all differential {}s found.".format(var_name))
                        figsize = (max(5, 0.12 * groups.shape[1]), 5)
                        # Heatmaps
                        # Comparison level
                        g = sns.clustermap(
                            groups.corr(),
                            xticklabels=False, yticklabels=True, cbar_kws={"label": "Pearson correlation\non differential {}s".format(var_name)},
                            cmap="BuGn", metric="correlation", rasterized=True, figsize=(figsize[0], figsize[0]))
                        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
                        savefig(g.fig, os.path.join(output_dir, output_prefix + ".diff_{}.groups.clustermap.corr.svg".format(var_name)))

                        g = sns.clustermap(
                            groups,
                            xticklabels=True, yticklabels=feature_labels, cbar_kws={"label": "{} of\ndifferential {}s".format(quantity, var_name)},
                            cmap="BuGn", robust=robust, metric="correlation", rasterized=True, figsize=figsize)
                        g.ax_heatmap.set_ylabel("Differential {}s (n = {})".format(var_name, groups.shape[0]))
                        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
                        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
                        savefig(g.fig, os.path.join(output_dir, output_prefix + ".diff_{}.groups.clustermap.svg".format(var_name)))

                        g = sns.clustermap(
                            groups,
                            xticklabels=True, yticklabels=feature_labels, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}s".format(quantity, var_name)},
                            cmap="RdBu_r", robust=robust, metric="correlation", rasterized=True, figsize=figsize)
                        g.ax_heatmap.set_ylabel("Differential {}s (n = {})".format(var_name, groups.shape[0]))
                        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
                        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
                        savefig(g.fig, os.path.join(output_dir, output_prefix + ".diff_{}.groups.clustermap.z0.svg".format(var_name)))

                        # same without clustering
                        g = sns.clustermap(
                            groups,
                            col_cluster=False,
                            xticklabels=True, yticklabels=feature_labels, cbar_kws={"label": "{} of\ndifferential {}s".format(quantity, var_name)},
                            cmap="BuGn", robust=robust, metric="correlation", rasterized=True, figsize=figsize)
                        g.ax_heatmap.set_ylabel("Differential {}s (n = {})".format(var_name, groups.shape[0]))
                        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
                        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
                        savefig(g.fig, os.path.join(output_dir, output_prefix + ".diff_{}.groups.sorted.clustermap.svg".format(var_name)))

                        g = sns.clustermap(
                            groups,
                            col_cluster=False,
                            xticklabels=True, yticklabels=feature_labels, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}s".format(quantity, var_name)},
                            cmap="RdBu_r", center=0, robust=robust, metric="correlation", rasterized=True, figsize=figsize)
                        g.ax_heatmap.set_ylabel("Differential {}s (n = {})".format(var_name, groups.shape[0]))
                        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
                        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
                        savefig(g.fig, os.path.join(output_dir, output_prefix + ".diff_{}.groups.sorted.clustermap.z0.svg".format(var_name)))

        # Fold-changes and P-values
        # pivot table of genes vs comparisons
        _LOGGER.info("Getting fold-change and p-value values per comparison.")
        fold_changes = pd.pivot_table(
            results.loc[all_diff, :].reset_index(),
            index=results.index.name, columns=comparison_column,
            values=log_fold_change_column).fillna(0)
        p_values = -np.log10(pd.pivot_table(
            results.loc[all_diff, :].reset_index(),
            index=results.index.name, columns=comparison_column,
            values=adjusted_p_value_column))

        # get a signed p-value
        if fold_changes.shape == p_values.shape:
            p_values *= (fold_changes > 0).astype(int).replace(0, -1)

        for matrix_, label, desc in [
            (fold_changes, "log fold change", "log fold change"),
            (p_values, "p value", "-log10(signed p-value)"),
        ]:
            if matrix_.shape[1] > 1:
                _LOGGER.info("Plotting group-wise correlation of {}s per sample groups in all differential {}s found.".format(var_name, label))
                figsize = (max(5, 0.12 * matrix_.shape[1]), 5)
                grid = sns.clustermap(
                    matrix_.corr(),
                    xticklabels=False, yticklabels=True,
                    cbar_kws={"label": "Pearson correlation\non {}s".format(desc)},
                    cmap="BuGn", vmin=0, vmax=1, metric="correlation", rasterized=True, figsize=(figsize[0], figsize[0]))
                grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
                grid.ax_heatmap.set_xlabel("Comparison groups")
                grid.ax_heatmap.set_ylabel("Comparison groups")
                savefig(
                    grid.fig,
                    os.path.join(output_dir,
                                 output_prefix + ".diff_{}.groups.{}.clustermap.corr.svg".format(var_name, label.replace(" ", "_"))))
                _LOGGER.info("Plotting clustered heatmaps of {}s per sample groups in all differential {}s found.".format(var_name, label))
                try:
                    grid = sns.clustermap(
                        matrix_.loc[all_diff, :],
                        xticklabels=True, yticklabels=feature_labels,
                        cbar_kws={"label": "{} of\ndifferential {}s".format(desc, var_name)},
                        cmap="RdBu_r", center=0, robust=robust, metric="correlation", rasterized=True, figsize=figsize)
                    grid.ax_heatmap.set_ylabel("Differential {}s (n = {})".format(var_name, matrix_.loc[all_diff, :].shape[0]))
                    grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
                    grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
                    grid.ax_heatmap.set_xlabel("Comparison groups")
                    savefig(
                        grid.fig,
                        os.path.join(output_dir,
                                     output_prefix + ".diff_{}.groups.{}.clustermap.svg".format(var_name, label.replace(" ", "_"))))
                except FloatingPointError:
                    _LOGGER.error("{} contain null of infinite values. Cannot plot.".format(label))

        # Sample level
        _LOGGER.info("Getting per sample values of {} in all differential {}s found.".format(quantity, var_name))
        if isinstance(matrix.columns, pd.core.indexes.multi.MultiIndex):
            matrix.columns = matrix.columns.get_level_values("sample_name")

        matrix2 = matrix.loc[all_diff, :].sort_index(axis=1)

        n = matrix2.isnull().sum().sum()
        if n > 0:
            _LOGGER.warning(
                "WARNING! {} {} (across all comparisons) were not found in quantification matrix!".format(n, var_name) +
                " Proceeding without those.")
            matrix2 = matrix2.dropna()
        figsize = (max(5, 0.12 * matrix2.shape[1]), 5)
        if group_wise_colours:
            extra = {"col_colors": color_dataframe.T}
        else:
            extra = {}

        _LOGGER.info("Plotting sample-wise correlation heatmaps of {} values per sample in all differential {}s found.".format(quantity, var_name))
        grid = sns.clustermap(
            matrix2.corr(),
            yticklabels=True, xticklabels=False,
            cbar_kws={"label": "Pearson correlation\non differential {}s".format(var_name)},
            cmap="BuGn", metric="correlation", figsize=(figsize[0], figsize[0]), rasterized=rasterized, robust=robust, **extra)
        grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        savefig(grid.fig, os.path.join(output_dir, output_prefix + ".diff_{}.samples.clustermap.corr.svg".format(var_name)))

        _LOGGER.info("Plotting clustered heatmaps of {} values per sample in all differential {}s found.".format(quantity, var_name))
        grid = sns.clustermap(
            matrix2,
            yticklabels=feature_labels, cbar_kws={"label": "{} of\ndifferential {}s".format(quantity, var_name)},
            xticklabels=True, vmin=0, cmap="BuGn", metric="correlation", figsize=figsize, rasterized=rasterized, robust=robust, **extra)
        grid.ax_heatmap.set_ylabel("Differential {}s (n = {})".format(var_name, matrix2.shape[0]))
        grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
        grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        savefig(grid.fig, os.path.join(output_dir, output_prefix + ".diff_{}.samples.clustermap.svg".format(var_name)))

        grid = sns.clustermap(
            matrix2,
            yticklabels=feature_labels, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}s".format(quantity, var_name)},
            xticklabels=True, cmap="RdBu_r", center=0, metric="correlation", figsize=figsize, rasterized=rasterized, robust=robust, **extra)
        grid.ax_heatmap.set_ylabel("Differential {}s (n = {})".format(var_name, matrix2.shape[0]))
        grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
        grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        savefig(grid.fig, os.path.join(output_dir, output_prefix + ".diff_{}.samples.clustermap.z0.svg".format(var_name)))

        grid = sns.clustermap(
            matrix2,
            col_cluster=False,
            yticklabels=feature_labels, cbar_kws={"label": "{} of\ndifferential {}s".format(quantity, var_name)},
            xticklabels=True, vmin=0, cmap="BuGn", metric="correlation", figsize=figsize, rasterized=rasterized, robust=robust, **extra)
        grid.ax_heatmap.set_ylabel("Differential {}s (n = {})".format(var_name, matrix2.shape[0]))
        grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
        grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        savefig(grid.fig, os.path.join(output_dir, output_prefix + ".diff_{}.samples.sorted.clustermap.svg".format(var_name)))

        grid = sns.clustermap(
            matrix2,
            col_cluster=False,
            yticklabels=feature_labels, z_score=0, cbar_kws={"label": "Z-score of {}\non differential {}s".format(quantity, var_name)},
            xticklabels=True, cmap="RdBu_r", center=0, metric="correlation", figsize=figsize, rasterized=rasterized, robust=robust, **extra)
        grid.ax_heatmap.set_ylabel("Differential {}s (n = {})".format(var_name, matrix2.shape[0]))
        grid.ax_heatmap.set_xticklabels(grid.ax_heatmap.get_xticklabels(), rotation=90, fontsize="xx-small")
        grid.ax_heatmap.set_yticklabels(grid.ax_heatmap.get_yticklabels(), rotation=0, fontsize="xx-small")
        savefig(grid.fig, os.path.join(output_dir, output_prefix + ".diff_{}.samples.sorted.clustermap.z0.svg".format(var_name)))

    def differential_overlap(
            self,
            differential=None,
            data_type="ATAC-seq",
            output_dir="{results_dir}/differential_analysis_{data_type}",
            output_prefix="differential_analysis"):
        """
        Visualize intersection of sets of differential regions/genes.

        Parameters
        ----------
        differential : pandas.DataFrame
            DataFrame containing result of comparisons filtered for features considered as differential.

        data_type : str, optional
            Data type.
            Defaults to analysis' own data type.

        output_dir : str, optional
            Directory to create output files.
            Defaults to "{results_dir}/differential_analysis_{data_type}".

        output_prefix : str, optional
            Prefix to use when creating output files.
            Defaults to "differential_analysis".
        """
        data_type = self._get_data_type(data_type)
        # Make output dir
        output_dir = self._format_string_with_attributes(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        if data_type in ["ATAC-seq", "ChIP-seq"]:
            unit = "region"
            total = self.coverage.shape[0]
        if data_type == "CNV":
            unit = "bin"
        elif data_type == "RNA-seq":
            unit = "gene"
            total = self.expression.shape[0]
        else:
            _LOGGER.warning("Unknown data type. Will not use data-specific units.")
            unit = "feature"

        if "direction" not in differential.columns:
            differential.loc[:, "direction"] = differential["log2FoldChange"].apply(lambda x: "up" if x > 0 else "down")

        differential.index.name = "index"
        differential.loc[:, "intersect"] = 1
        piv = pd.pivot_table(
            differential.reset_index(),
            index="index", columns=["comparison_name", "direction"], values="intersect", fill_value=0)

        intersections = pd.DataFrame(columns=["group1", "group2", "dir1", "dir2", "size1", "size2", "intersection", "union"])
        perms = list(itertools.permutations(piv.T.groupby(level=["comparison_name", "direction"]).groups.items(), 2))
        for ((k1, dir1), i1), ((k2, dir2), i2) in tqdm(perms, total=len(perms), desc="Permutations"):
            i1 = set(piv[i1][piv[i1] == 1].dropna().index)
            i2 = set(piv[i2][piv[i2] == 1].dropna().index)
            intersections = intersections.append(
                pd.Series(
                    [k1, k2, dir1, dir2, len(i1), len(i2), len(i1.intersection(i2)), len(i1.union(i2))],
                    index=["group1", "group2", "dir1", "dir2", "size1", "size2", "intersection", "union"]
                ),
                ignore_index=True
            )
        # convert to %
        intersections.loc[:, "intersection"] = intersections["intersection"].astype(float)
        intersections.loc[:, "perc_1"] = intersections["intersection"] / intersections["size1"] * 100.
        intersections.loc[:, "perc_2"] = intersections["intersection"] / intersections["size2"] * 100.
        intersections.loc[:, "intersection_max_perc"] = intersections[["perc_1", "perc_2"]].max(axis=1)

        # calculate p-value from Fisher"s exact test
        intersections.loc[:, "a"] = total - intersections[["size1", "size2", "intersection"]].sum(axis=1)
        intersections.loc[:, "b"] = intersections["size1"] - intersections["intersection"]
        intersections.loc[:, "c"] = intersections["size2"] - intersections["intersection"]
        intersections.loc[:, "d"] = intersections["intersection"]

        for i, row in intersections[["d", "b", "c", "a"]].astype(int).iterrows():
            odds, p = fisher_exact(
                row
                .values
                .reshape((2, 2)),
                alternative="greater")
            intersections.loc[i, "odds_ratio"] = odds
            intersections.loc[i, "p_value"] = p
        # intersections["q_value"] = intersections["p_value"] * intersections.shape[0]
        intersections.loc[:, "q_value"] = multipletests(intersections["p_value"])[1]
        intersections.loc[:, "log_p_value"] = log_pvalues(intersections["p_value"])
        intersections.loc[:, "log_q_value"] = log_pvalues(intersections["q_value"])

        # save
        intersections.to_csv(os.path.join(output_dir, output_prefix + ".differential_overlap.csv"), index=False)
        intersections = pd.read_csv(os.path.join(output_dir, output_prefix + ".differential_overlap.csv"))

        for metric, label, description, fill_value in [
                ("intersection", "intersection", "total in intersection", 0),
                ("intersection_max_perc", "percentage_overlap", "max of intersection %", 0),
                ("log_p_value", "significance", "p-value", 0)]:
            _LOGGER.debug(metric)
            # make pivot tables
            piv_up = pd.pivot_table(
                intersections[(intersections["dir1"] == "up") & (intersections["dir2"] == "up")],
                index="group1", columns="group2", values=metric).fillna(fill_value)
            piv_down = pd.pivot_table(
                intersections[(intersections["dir1"] == "down") & (intersections["dir2"] == "down")],
                index="group1", columns="group2", values=metric).fillna(fill_value)
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
                square=True, cmap="Blues",
                cbar_kws={"label": "Concordant {}s ({})".format(unit, description)},
                ax=axis[0], **extra)
            sns.heatmap(
                piv_up,
                square=True, cmap="Reds",
                cbar_kws={"label": "Concordant {}s ({})".format(unit, description)},
                ax=axis[1], **extra)
            axis[0].set_title("Downregulated {}s".format(unit))
            axis[1].set_title("Upregulated {}s".format(unit))
            axis[0].set_xticklabels(axis[0].get_xticklabels(), rotation=90, ha="center")
            axis[1].set_xticklabels(axis[1].get_xticklabels(), rotation=90, ha="center")
            axis[0].set_yticklabels(axis[0].get_yticklabels(), rotation=0, ha="right")
            axis[1].set_yticklabels(axis[1].get_yticklabels(), rotation=0, ha="right")
            savefig(fig, os.path.join(
                    output_dir, output_prefix + ".differential_overlap.{}.up_down_split.svg"
                    .format(label)))

            # combined heatmap
            # with upregulated {}s in upper square matrix and downredulated in down square
            piv_combined = pd.DataFrame(np.triu(piv_up), index=piv_up.index, columns=piv_up.columns).replace(0, np.nan)
            piv_combined.update(pd.DataFrame(np.tril(-piv_down), index=piv_down.index, columns=piv_down.columns).replace(0, np.nan))
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
                square=True, cmap="RdBu_r", center=0,
                cbar_kws={"label": "Concordant {}s ({})".format(unit, description)},
                ax=axis, **extra)
            axis.set_xticklabels(axis.get_xticklabels(), rotation=90, ha="center")
            axis.set_yticklabels(axis.get_yticklabels(), rotation=0, ha="right")
            savefig(fig, os.path.join(
                    output_dir, output_prefix + ".differential_overlap.{}.up_down_together.svg"
                    .format(label)))

            # Rank plots
            if metric == "log_pvalue":
                r = pd.melt(
                    piv_combined.reset_index(),
                    id_vars=["group1"], var_name="group2", value_name="agreement")
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
                        row["rank"], row["agreement"],
                        s=row["group1"] + "\n" + row["group2"],
                        fontsize=5)
                for i, row in r.head(10).iterrows():
                    axis[2].text(
                        row["rank"], row["agreement"],
                        s=row["group1"] + "\n" + row["group2"],
                        fontsize=5)
                for ax in axis:
                    ax.set_ylabel("Agreement (-log(p-value))")
                    ax.set_xlabel("Rank")
                sns.despine(fig)
                savefig(fig, os.path.join(
                        output_dir, output_prefix + ".differential_overlap.{}.agreement.rank.svg"
                        .format(label)))

            # Observe disagreement
            # (overlap of down-regulated with up-regulated and vice-versa)
            piv_up = pd.pivot_table(
                intersections[(intersections["dir1"] == "up") & (intersections["dir2"] == "down")],
                index="group1", columns="group2", values=metric)
            piv_down = pd.pivot_table(
                intersections[(intersections["dir1"] == "down") & (intersections["dir2"] == "up")],
                index="group1", columns="group2", values=metric)

            piv_disagree = pd.concat([piv_up, piv_down]).groupby(level=0).max()
            if metric == "intersection":
                piv_disagree = np.log10(1 + piv_disagree)
            np.fill_diagonal(piv_disagree.values, np.nan)

            fig, axis = plt.subplots(1, 2, figsize=(16, 8), subplot_kw={"aspect": "equal"})
            sns.heatmap(
                piv_disagree, square=True, cmap="Greens",
                cbar_kws={"label": "Discordant {}s ({})".format(unit, description)}, ax=axis[0])

            norm = matplotlib.colors.Normalize(vmin=0, vmax=piv_disagree.max().max())
            cmap = plt.get_cmap("Greens")
            for j, g2 in enumerate(piv_disagree.index):
                for i, g1 in enumerate(piv_disagree.columns):
                    axis[1].scatter(
                        len(piv_disagree.index) - (j + 0.5), len(piv_disagree.index) - (i + 0.5),
                        s=(100 ** (norm(piv_disagree.loc[g1, g2]))) - 1,
                        color=cmap(norm(piv_disagree.loc[g1, g2])), marker="o")
                    axis[1].set_title("Rotate plot -90 degrees")
            axis[0].set_xticklabels(axis[0].get_xticklabels(), rotation=90, ha="center")
            axis[0].set_yticklabels(axis[0].get_yticklabels(), rotation=0, ha="right")
            axis[1].set_xlim((0, len(piv_disagree.index)))
            axis[1].set_ylim((0, len(piv_disagree.columns)))
            savefig(fig, os.path.join(
                    output_dir, output_prefix + ".differential_overlap.{}.disagreement.svg"
                    .format(label)))

            # Rank plots
            if metric == "log_pvalue":
                r = pd.melt(
                    piv_disagree.reset_index(),
                    id_vars=["group1"], var_name="group2", value_name="disagreement")
                r = r.dropna().sort_values("disagreement")
                r = r.iloc[range(0, r.shape[0], 2)]
                r["rank"] = r["disagreement"].rank(ascending=False)

                fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 4), subplot_kw={"aspect": "equal"})
                axis[0].scatter(r["rank"], r["disagreement"])
                axis[1].scatter(r["rank"].tail(10), r["disagreement"].tail(10))
                for i, row in r.tail(10).iterrows():
                    axis[1].text(
                        row["rank"], row["disagreement"],
                        s=row["group1"] + "\n" + row["group2"],
                        fontsize=5)
                for ax in axis:
                    ax.set_ylabel("Disagreement (-log(p-value))")
                    ax.set_xlabel("Rank")
                sns.despine(fig)
                savefig(fig, os.path.join(
                    output_dir, output_prefix + ".differential_overlap.{}.disagreement.rank.svg"
                    .format(label)))

    @check_organism_genome
    def differential_enrichment(
            self,
            differential=None,
            data_type=None,
            output_dir="{results_dir}/differential_analysis_{data_type}/enrichments",
            output_prefix="differential_analysis",
            genome=None,
            steps=["region", "lola", "meme", "homer", "enrichr"],
            directional=True,
            max_diff=1000,
            sort_var="pvalue",
            distributed=False,
            overwrite=False):
        """
        Perform various types of enrichment analysis given a dataframe of the results from differential analysis.
        Performs enrichment of gene sets (RNA-seq and ATAC-seq), genomic regions, chromatin states
        Location Overlap Analysis (LOLA) and TF motif enrichment (over-representation and de-novo search)
        (ATAC-seq only).

        Parameters
        ----------
        differential : pandas.DataFrame
            Data frame with differential results as produced by ``differential_analysis``.
            Must contain a "comparison_name" column.
            Defaults to `analysis.differential_results`.

        data_type : str, optional
            Data type. One of "ATAC-seq" and "RNA-seq".
            Defaults to the analysis' data_type attributes.

        output_dir : str, optional
            Directory to create output files.
            Defaults to "{results_dir}/differential_analysis_{data_type}".

        output_prefix : str, optional
            Prefix to use when creating output files.
            Defaults to "differential_analysis".

        genome : str, optional
            Genome assembly of the analysis.
            Defaults to Analysis's `genome` attribute.

        steps : list, optional
            Steps of the analysis to perform.
            Defaults to all possible: ["region", lola", "meme", "homer", "enrichr"].

        directional : bool, optional
            Whether enrichments should be performed in a direction-dependent way
            (up-regulated and down-regulated features separately).
            This requires a column named "log2FoldChange" to exist.
            Defaults to True.

        max_diff : int, optional
            Number of maximum features to perform enrichment for ranked by variable in `max_diff`.
            Defaults to 1000.

        sort_var : str, optional
            Variable to sort for when setting `max_diff`.
            Defaults to "pvalue".

        distributed : bool, optional
            Whether work should be submitted as jobs in a computing cluster.
            Defaults to False.

        overwrite : bool, optional
            Whether output files should be overwritten when `distributed` is True.
            Defaults to False.

        Attributes
        ----------
        enrichment_results : dict
            Dictionary with keys as in `steps` and values with pandas.DataFrame
            of enrichment results.
        """
        # TODO: separate and fix mouse TF ids
        # TODO: separate homer_consensus output processing
        # TODO: add overwrite function when distributed==False
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

        if data_type is None:
            msg = "Data type not given and Analysis object does not have a `data_type` attribute."
            try:
                data_type = self.data_type
            except AttributeError as e:
                _LOGGER.error(msg)
                raise e
            if data_type is None:
                _LOGGER.error(msg)
                raise ValueError

        if genome is None:
            genome = self.genome

        if data_type == "ATAC-seq":
            matrix = self.coverage_annotated
        elif data_type == "RNA-seq":
            matrix = self.expression
        else:
            msg = "Differential enrichment is only implemented for data types 'ATAC-seq' and 'RNA-seq'."
            raise ValueError(msg)

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
            ("region", region_enr, pd.read_csv, {}, "region_type_enrichment.csv", ".region_type_enrichment.csv"),
            ("meme", meme_enr, parse_ame, {}, "ame.txt", ".meme_ame.csv"),
            ("homer", homer_enr, parse_homer, {}, "homerResults", ".homer_motifs.csv"),
            ("lola", lola_enr, pd.read_csv, {"sep", "\t"}, "allEnrichments.tsv", ".lola.csv"),
            ("enrichr", pathway_enr, pd.read_csv, {"encoding": "utf-8"}, output_prefix + "_regions.enrichr.csv", ".enrichr.csv")]

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
                        (differential["comparison_name"] == comp) &
                        (f(differential["log2FoldChange"], arg)), :].index
                else:
                    diff = differential.loc[
                        (differential["comparison_name"] == comp), :].index

                # Handle extremes of regions
                if diff.shape[0] < 1:
                    continue
                if diff.shape[0] > max_diff:
                    if directional:
                        diff = (
                            getattr(
                                differential[
                                    (differential["comparison_name"] == comp) &
                                    (f(differential["log2FoldChange"], arg))]
                                [sort_var].sort_values(), top)
                            (max_diff).index)
                    else:
                        diff = (
                            getattr(
                                differential[
                                    (differential["comparison_name"] == comp)]
                                [sort_var].sort_values(), top)
                            (max_diff).index)

                # Add data_type specific info
                comparison_df = matrix.loc[diff, :]
                if comparison_df.shape != comparison_df.dropna().shape:
                    _LOGGER.warning(
                        "There are differential regions which are not in the set" +
                        " of annotated regions for comparison '{}'!".format(comp) +
                        " Continuing enrichment without those.")
                    comparison_df = comparison_df.dropna()

                # Prepare output dir
                comparison_dir = os.path.join(output_dir, "{}.{}".format(comp, direction))
                if not os.path.exists(comparison_dir):
                    os.makedirs(comparison_dir)

                # Prepare files and run (if not distributed)
                if data_type == "RNA-seq":
                    _LOGGER.info("Doing genes of comparison '{}', direction '{}'.".format(comp, direction))
                    comparison_df.index.name = "gene_name"
                    # write gene names to file
                    clean = comparison_df.reset_index()["gene_name"].drop_duplicates().sort_values()
                    clean.to_csv(
                        os.path.join(comparison_dir, output_prefix + ".gene_symbols.txt"),
                        header=None, index=False)

                    if "enrichr" in steps:
                        if serial:
                            if not os.path.exists(os.path.join(comparison_dir, output_prefix + ".enrichr.csv")):
                                enr = enrichr(comparison_df.reset_index())
                                enr.to_csv(
                                    os.path.join(comparison_dir, output_prefix + ".enrichr.csv"),
                                    index=False, encoding="utf-8")
                else:
                    _LOGGER.info("Doing regions of comparison '{}', direction '{}'.".format(comp, direction))
                    # do the suite of region enrichment analysis
                    self.characterize_regions_function(
                        comparison_df,
                        output_dir=comparison_dir, prefix=output_prefix, run=serial,
                        genome=genome, steps=steps)

                # collect enrichments
                if serial:
                    # read/parse, label and append
                    for name, df, function, kwargs, suffix, _ in possible_steps:
                        if name in steps:
                            enr = function(os.path.join(comparison_dir, suffix), **kwargs)
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
                        os.path.join(output_dir, output_prefix + output_suffix), index=False, encoding="utf-8")
        else:
            try:
                _LOGGER.info("Using background region set from analysis.sites")
                background = getattr(self, "sites").fn
            except AttributeError:
                _LOGGER.warning("Using no background region set because 'analysis.sites' is not set!")
                background = ""
            _LOGGER.info("Submitting enrichment jobs.")
            run_enrichment_jobs(
                results_dir=output_dir,
                genome=genome, background_bed=background,
                steps=steps, overwrite=overwrite,
                pickle_file=self.pickle_file)

    def collect_differential_enrichment(
            self,
            steps=["region", "lola", "motif", "homer", "homer_consensus", "enrichr"],
            directional=True,
            permissive=True,
            output_dir="{results_dir}/differential_analysis_{data_type}/enrichments",
            input_prefix="differential_analysis",
            output_prefix="differential_analysis",
            differential=None,
            data_type=None):
        """
        Collect the results of enrichment analysis ran after a differential analysis.

        Parameters
        ----------
        steps : list, optional
            Steps of the enrichment analysis to collect results for.
            Defaults to ["region", "lola", "meme", "homer", "enrichr"].

        directional : bool, optional
            Whether enrichments were made in a direction-dependent way
            (up-regulated and down-regulated features separately).
            This implies a column named "direction" exists".
            Defaults to True.

        differential : pandas.DataFrame, optional
            Data frame with differential results to select which comparisons to collect
            enrichments for. Usually produced by ``ngs_toolkit.general.differential_analysis``.
            Defaults to analysis `differential_results`.

        data_type : str, optional
            Data type. One of "ATAC-seq" and "RNA-seq".
            Defaults to Analysis' data_type.

        output_dir : str, optional
            Directory to create output files.
            Defaults to "{results_dir}/differential_analysis_{data_type}".

        output_prefix : str, optional
            Prefix to use when creating output files.
            Defaults to "differential_analysis".

        permissive : bool, optional
            Whether to skip non-existing files, giving a warning.
            Defaults to True.

        Attributes
        ----------
        enrichment_results : dict
            Dictionary with keys as in `steps` and values with pandas.DataFrame
            of enrichment results.
        """
        # TODO: separate and fix mouse TF ids
        # TODO: separate homer_consensus output processing
        data_type = self._get_data_type(data_type)
        if data_type not in ["ATAC-seq", "RNA-seq"]:
            raise ValueError("`data_type` must match one of 'ATAC-seq' or 'RNA-seq'.")

        output_dir = self._format_string_with_attributes(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        data_type_steps = {
            "ATAC-seq": ["region", "lola", "motif", "homer", "homer_consensus", "enrichr"],
            "ChIP-seq": ["region", "lola", "motif", "homer", "homer_consensus", "enrichr"],
            "RNA-seq": ["enrichr"]}
        if steps is None:
            steps = [s for s in steps if s in data_type_steps[data_type]]
        for step in steps:
            if step not in data_type_steps[data_type]:
                steps.pop(step)
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
            ("region", region_enr, pd.read_csv, {}, "region_type_enrichment.csv", ".region_type_enrichment.csv"),
            ("meme", meme_enr, parse_ame, {}, "ame.txt", ".meme_ame.csv"),
            ("homer", homer_enr, parse_homer, {}, "homerResults", ".homer_motifs.csv"),
            ("homer_consensus", homer_consensus, pd.read_csv, {"sep": "\t"}, "knownResults.txt", ".homer_consensus.csv"),
            ("lola", lola_enr, pd.read_csv, {"sep": "\t"}, "allEnrichments.tsv", ".lola.csv"),
            ("enrichr", pathway_enr, pd.read_csv, {"encoding": "utf-8"}, output_prefix + "_genes.enrichr.csv", ".enrichr.csv")]

        # Examine each region cluster
        comps = differential["comparison_name"].drop_duplicates()
        for comp in tqdm(comps, total=len(comps), desc="Comparison"):
            if directional:
                # Separate in up/down-regulated genes
                params = list()
                if differential[
                        (differential["comparison_name"] == comp) &
                        (differential["log2FoldChange"] > 0)].shape[0] > 0:
                    params.append("up")
                if differential[
                        (differential["comparison_name"] == comp) &
                        (differential["log2FoldChange"] < 0)].shape[0] > 0:
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
                        except IOError as e:
                            if permissive:
                                _LOGGER.warn(error_msg.format(name, comp, direction))
                            else:
                                raise e
                        else:
                            if not enr.empty:
                                enr.loc[:, "comparison_name"] = comp
                                enr.loc[:, "direction"] = direction
                                enr.loc[:, "label"] = "{}.{}".format(comp, direction)
                                df.append(enr)
                            else:
                                _LOGGER.warning("Comparison '{}' {} results are empty!".format(comp, name))

        # write combined enrichments
        _LOGGER.info("Saving combined enrichments for all comparisons.")
        self.enrichment_results = dict()

        for name, df, function, kwargs, suffix, output_suffix in possible_steps:
            if name in steps:
                if len(df) == 0:
                    msg = "No comparison has {} results. Saving empty dataframe!".format(name)
                    _LOGGER.warning(msg)
                    self.enrichment_results[name] = pd.DataFrame()
                else:
                    self.enrichment_results[name] = pd.concat(df, axis=0)
                self.enrichment_results[name].to_csv(
                    os.path.join(output_dir, output_prefix + output_suffix), index=False, encoding="utf-8")

    def plot_differential_enrichment(
            self,
            steps=None,
            enrichment_type=None,
            enrichment_table=None,
            direction_dependent=True,
            output_dir="{results_dir}/differential_analysis_{data_type}/enrichments",
            comp_variable="comparison_name",
            output_prefix="differential_analysis",
            rasterized=True,
            barplots=True,
            correlation_plots=True,
            clustermap_metric="correlation",
            top_n=5,
            z_score=0,
            cmap=None):
        """
        Make plots illustrating enrichment of features for various comparisons.

        Input can be the dictionary under `analysis.enrichment_results` or
        a single dataframe of enrichment terms across several comparisons for a given type of enrichment.
        In the later case both `enrichment_table` and `enrichment_type` must be given.

        Parameters
        ----------
        steps : list, optional
            Types of the enrichment analysis to plot.
            One of {"region", "lola", "motif", "great", "enrichr"}.
            Defaults (None) to all keys present in analysis.enrichment_results.

        enrichment_type : str, optional
            Type of enrichment if run for a single type of enrichment.
            In this case `enrichment_table` must be given.
            One of {"region", "lola", "motif", "great", "enrichr"}.
            Default (None) is to run all keys present in analysis.enrichment_results.

        enrichment_table : pandas.DataFrame, optional
            Data frame with enrichment results as produced by
            ``differential_enrichment`` or ``collect_differential_enrichment``.
            If given, `enrichment_type` must be given too.
            Default (None) is the dataframes in all values present in analysis.enrichment_results.

        direction_dependent : bool, optional
            Whether enrichments were made in a direction-dependent way (up-regulated and down-regulated features separately).
            This implies a column named "direction" exists".
            Defaults to True.

        output_dir : str, optional
            Directory to create output files.
            Defaults to "{results_dir}/differential_analysis_{data_type}/enrichments".

        comp_variable : str, optional
            Column defining which comparison enrichment terms belong to.
            Defaults to "comparison_name".

        output_prefix : str, optional
            Prefix to use when creating output files.
            Defaults to "differential_analysis".

        rasterized : bool, optional
            Whether or not to rasterize heatmaps for efficient plotting.
            Defaults to True.

        barplots : bool, optional
            Whether barplots with top enriched terms per comparison should be produced.
            Defaults to True.

        correlation_plots : bool, optional
            Whether correlation plots of comparisons across enriched terms should be produced.
            Defaults to True.

        clustermap_metric : str, optional
            Distance metric to use for clustermap clustering,
            Default to "correlation" (Pearson's).

        top_n : int, optional
            Top terms to be used to make barplots.
            Defaults to 5

        z_score : bool, optional
            Which dimention/axis to perform Z-score transformation for.
            Numpy/Pandas conventions are used:
            `0` is row-wise (in this case across comparisons) and `1` is column-wise (across terms).
            Defaults to 0.

        cmap : str, optional
            Colormap to use in heatmaps.
            Default None.
        """
        # TODO: split function in its smaller parts and call them appropriately.
        def enrichment_barplot(
                input_df,
                x,
                y,
                group_variable,
                top_n,
                output_file):
            n = len(input_df[group_variable].drop_duplicates())
            n_side = int(np.ceil(np.sqrt(n)))

            top_data = enrichment_table.set_index(x).groupby(group_variable)[y].nlargest(top_n).reset_index()

            fig, axis = plt.subplots(n_side, n_side, figsize=(
                4 * n_side, n_side * max(5, 0.12 * top_n)), sharex=False, sharey=False)
            if isinstance(axis, np.ndarray):
                axis = iter(axis.flatten())
            else:
                axis = iter(np.array([axis]))
            for comp in top_data[group_variable].drop_duplicates().sort_values():
                df2 = top_data.loc[top_data[group_variable] == comp, :]
                ax = next(axis)
                sns.barplot(
                    df2[y], df2[x], estimator=max,
                    orient="horizontal", ax=ax, color=sns.color_palette("colorblind")[0])
                ax.set_title(comp)
            for ax in axis:
                ax.set_visible(False)
            sns.despine(fig)
            savefig(fig, output_file)

        def enrichment_correlation_plot(
                input_df,
                output_file,
                label="Pearson correlation of enrichment"):
            try:
                g = sns.clustermap(
                    input_df.T.corr(),
                    rasterized=rasterized,
                    xticklabels=True,
                    yticklabels=True,
                    cbar_kws={"label": label})
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
                savefig(g.fig, output_file)
            except FloatingPointError:
                msg = "Plotting of correlation matrix failed: {}".format(output_file)
                _LOGGER.warn(msg)

        def enrichment_clustermap(
                input_df,
                output_file,
                label="Enrichment\nof differential regions",
                z_score=None, params={}):
            # plot clustered heatmap
            shape = input_df.shape
            if z_score is not None:
                params.update({"cmap": "RdBu_r", "center": 0, "z_score": z_score})
            try:
                g = sns.clustermap(
                    input_df, figsize=(
                        max(6, 0.12 * shape[1]), max(6, 0.12 * shape[0])),
                    metric=clustermap_metric,
                    xticklabels=True, yticklabels=True, rasterized=rasterized,
                    cbar_kws={"label": label}, **params)
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(),
                                             rotation=90, ha="right", fontsize="xx-small")
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(),
                                             rotation=0, fontsize="xx-small")
                savefig(g.fig, output_file)
            except FloatingPointError:
                msg = "Plotting of correlation matrix failed: {}".format(output_file)
                _LOGGER.warn(msg)

        if steps is None:
            steps = ["region", "lola", "motif", "great", "enrichr"]

        if (enrichment_table is None) and (enrichment_type is None):
            if not hasattr(self, "enrichment_results"):
                msg = "'enrichment_table' and 'enrichment_type' were not given"
                msg += "but analysis also does not have a 'enrichment_results' attribute."
                _LOGGER.error(msg)
                raise ValueError(msg)
            else:
                for enrichment_name, enrichment_table in self.enrichment_results.items():
                    if enrichment_name in steps:
                        self.plot_differential_enrichment(
                            steps=[enrichment_name],
                            enrichment_table=enrichment_table,
                            enrichment_type=enrichment_name,
                            direction_dependent=direction_dependent,
                            output_dir=output_dir,
                            comp_variable=comp_variable,
                            output_prefix=output_prefix,
                            rasterized=rasterized,
                            barplots=barplots,
                            correlation_plots=correlation_plots,
                            clustermap_metric=clustermap_metric,
                            top_n=top_n if enrichment_name != "motif" else 300,
                            z_score=z_score,
                            cmap=cmap)
                return

        if enrichment_type not in ["region", "lola", "enrichr", "motif", "homer_consensus", "great"]:
            raise AssertionError("`enrichment_type` must be one of 'lola', 'enrichr', 'motif', 'homer_consensus', 'great'.")

        output_dir = self._format_string_with_attributes(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        if z_score == 0:
            z_score_label = "Row"
        elif z_score == 1:
            z_score_label = "Column"
        elif z_score is None:
            pass
        else:
            raise ValueError("Argument 'z_score' must be on of 0, 1 or None.")

        enrichment_table = enrichment_table.copy()
        if "direction" in enrichment_table.columns and direction_dependent:
            enrichment_table.loc[:, comp_variable] = enrichment_table[comp_variable].astype(
                str) + " " + enrichment_table["direction"].astype(str)

        if enrichment_type == "region":
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
                    pvalue=0.05, top_n=top_n)

            # pivot table
            region_pivot = pd.pivot_table(
                enrichment_table, values="log2_odds_ratio", columns=comp_variable, index="region").fillna(0)

            # plot correlation
            if correlation_plots:
                enrichment_correlation_plot(
                    input_df=region_pivot, label="Correlation of enrichment\nof differential regions",
                    output_file=os.path.join(output_dir, output_prefix + ".region_type_enrichment.correlation.svg"))

            # plot clustered heatmaps
            enrichment_clustermap(
                region_pivot,
                output_file=os.path.join(output_dir, output_prefix + ".region_type_enrichment.cluster_specific.svg"),
                label="log2(odd ratio) of enrichment\nof differential regions", params={"cmap": "RdBu_r", "center": 0})

        if enrichment_type == "lola":
            # get a unique label for each lola region set
            enrichment_table.loc[:, "label"] = (
                enrichment_table["description"].astype(str) + ", " +
                enrichment_table["cellType"].astype(str) + ", " +
                enrichment_table["tissue"].astype(str) + ", " +
                enrichment_table["antibody"].astype(str) + ", " +
                enrichment_table["treatment"].astype(str))
            enrichment_table.loc[:, "label"] = (
                enrichment_table["label"]
                .str.replace("nan", "").str.replace("None", "")
                .str.replace(", , ", "").str.replace(", $", "")
                .str.decode("unicode_escape").str.encode("ascii", "ignore"))

            # Replace inf values with maximum non-inf p-value observed
            r = enrichment_table.loc[enrichment_table["pValueLog"] != np.inf, "pValueLog"]
            r += r * 0.1
            enrichment_table.loc[:, "pValueLog"] = enrichment_table["pValueLog"].replace(np.inf, r)

            # Plot top_n terms of each comparison in barplots
            if barplots:
                enrichment_barplot(
                    enrichment_table, x="label", y="pValueLog",
                    group_variable=comp_variable, top_n=top_n,
                    output_file=os.path.join(output_dir, output_prefix + ".lola.barplot.top_{}.svg".format(top_n)))

            # Plot heatmaps of terms for each comparison
            if len(enrichment_table[comp_variable].drop_duplicates()) < 2:
                return

            # pivot table
            lola_pivot = pd.pivot_table(
                enrichment_table, values="pValueLog", columns=comp_variable, index="label").fillna(0)
            lola_pivot = lola_pivot.replace(np.inf, lola_pivot[lola_pivot != np.inf].max().max())

            # plot correlation
            if correlation_plots:
                enrichment_correlation_plot(
                    input_df=lola_pivot, label="Correlation of enrichment\nof differential regions",
                    output_file=os.path.join(output_dir, output_prefix + ".lola.correlation.svg"))

            # plot clustered heatmaps of top terms
            top = enrichment_table.set_index("label").groupby(comp_variable)["pValueLog"].nlargest(top_n)
            top_terms = top.index.get_level_values("label").unique()
            enrichment_clustermap(
                lola_pivot.loc[top_terms, :],
                output_file=os.path.join(output_dir, output_prefix + ".lola.cluster_specific.svg"),
                label="-log10(p-value) of enrichment\nof differential regions")
            if z_score is not None:
                enrichment_clustermap(
                    lola_pivot.loc[top_terms, :],
                    output_file=os.path.join(output_dir, output_prefix + ".lola.cluster_specific.{}_z_score.svg".format(z_score_label)),
                    label="{} Z-score of enrichment\nof differential regions".format(z_score_label), z_score=z_score)

        if enrichment_type == "motif":
            enrichment_table.loc[:, "log_p_value"] = log_pvalues(enrichment_table["p_value"])

            # Plot top_n terms of each comparison in barplots
            if barplots:
                enrichment_barplot(
                    enrichment_table, x="TF", y="log_p_value",
                    group_variable=comp_variable, top_n=top_n,
                    output_file=os.path.join(output_dir, output_prefix + ".lola.barplot.top_{}.svg".format(top_n)))
            if len(enrichment_table[comp_variable].drop_duplicates()) < 2:
                return
            # Plot heatmaps of terms for each comparison
            motifs_pivot = pd.pivot_table(
                enrichment_table, values="log_p_value", columns="TF", index=comp_variable).fillna(0)
            # plot correlation
            if correlation_plots:
                enrichment_correlation_plot(
                    input_df=motifs_pivot, label="Correlation of enrichment\nof differential regions",
                    output_file=os.path.join(output_dir, output_prefix + ".motifs.correlation.svg"))

            # plot clustered heatmaps of top terms
            top = enrichment_table.set_index("TF").groupby(comp_variable)["log_p_value"].nlargest(top_n)
            top_terms = top.index.get_level_values("TF").unique()
            enrichment_clustermap(
                motifs_pivot.loc[:, top_terms],
                output_file=os.path.join(output_dir, output_prefix + ".motifs.cluster_specific.svg"),
                label="-log10(p-value) of enrichment\nof differential regions")
            if z_score is not None:
                enrichment_clustermap(
                    motifs_pivot.loc[:, top_terms],
                    output_file=os.path.join(output_dir, output_prefix + ".motifs.cluster_specific.{}_z_score.svg".format(z_score_label)),
                    label="{} Z-score of enrichment\nof differential regions".format(z_score_label), z_score=z_score)

        if enrichment_type == "homer_consensus":
            enrichment_table.loc[:, "enrichment_over_background"] = (
                enrichment_table["% of Target Sequences with Motif"] /
                enrichment_table["% of Background Sequences with Motif"])
            enrichment_table.loc[:, "log_p_value"] = log_pvalues(enrichment_table["P-value"])

            # Plot top_n terms of each comparison in barplots
            top_n = min(top_n, enrichment_table.set_index("Motif Name").groupby(comp_variable)["log_p_value"].count().min() - 1)
            if barplots:
                enrichment_barplot(
                    enrichment_table, x="Motif Name", y="log_p_value",
                    group_variable=comp_variable, top_n=top_n,
                    output_file=os.path.join(output_dir, output_prefix + ".lola.barplot.top_{}.svg".format(top_n)))

            # Significance vs fold enrichment over background
            n = len(enrichment_table[comp_variable].drop_duplicates())
            n_side = int(np.ceil(np.sqrt(n)))
            fig, axis = plt.subplots(n_side, n_side, figsize=(3 * n_side, 3 * n_side), sharex=False, sharey=False)
            axis = axis.flatten()
            for i, comp in enumerate(enrichment_table[comp_variable].drop_duplicates().sort_values()):
                enr = enrichment_table[enrichment_table[comp_variable] == comp]
                enr.loc[:, "Motif Name"] = enr["Motif Name"].str.replace(".*BestGuess:", "").str.replace(r"-ChIP-Seq.*", "")

                enr.loc[:, "combined"] = enr[["enrichment_over_background", "log_p_value"]].apply(zscore).mean(axis=1)
                axis[i].scatter(
                    enr["enrichment_over_background"],
                    enr["log_p_value"],
                    c=enr["combined"],
                    s=8, alpha=0.75)

                # label top points
                for j in enr["combined"].sort_values().tail(5).index:
                    axis[i].text(
                        enr.loc[j, "enrichment_over_background"],
                        enr.loc[j, "log_p_value"],
                        s=enr.loc[j, "Motif Name"], ha="right", fontsize=5)
                axis[i].set_title(comp)

            for ax in axis.reshape((n_side, n_side))[:, 0]:
                ax.set_ylabel("-log10(p-value)")
            for ax in axis.reshape((n_side, n_side))[-1, :]:
                ax.set_xlabel("Enrichment over background")
            sns.despine(fig)
            savefig(fig, os.path.join(output_dir, output_prefix + ".homer_consensus.scatterplot.svg"))

            # Plot heatmaps of terms for each comparison
            if len(enrichment_table[comp_variable].drop_duplicates()) < 2:
                return

            for label, metric in [
                    ("-log10(p-value) of enrichment\nin", "log_p_value"),
                    ("Enrichment over background of\n", "enrichment_over_background")]:
                # pivot table
                motifs_pivot = pd.pivot_table(
                    enrichment_table, values=metric, columns="Motif Name", index=comp_variable).fillna(0)

                # plot correlation
                if correlation_plots:
                    enrichment_correlation_plot(
                        input_df=motifs_pivot, label="Correlation of enrichment\nof differential regions",
                        output_file=os.path.join(output_dir, output_prefix + ".homer_consensus.correlation.svg"))

                top = enrichment_table.set_index("Motif Name").groupby(comp_variable)[metric].nlargest(top_n)
                top_terms = top.index.get_level_values("Motif Name").unique()

                # plot clustered heatmaps of top terms
                enrichment_clustermap(
                    motifs_pivot.loc[:, top_terms].T,
                    output_file=os.path.join(output_dir, output_prefix + ".homer_consensus.cluster_specific.svg"),
                    label=label + " differential regions")
                if z_score is not None:
                    enrichment_clustermap(
                        motifs_pivot.loc[:, top_terms].T,
                        output_file=os.path.join(
                            output_dir, output_prefix + ".homer_consensus.cluster_specific.{}_z_score.svg"
                            .format(z_score_label)),
                        label="{} Z-score of {} differential regions".format(z_score_label, label), z_score=z_score)

        if enrichment_type == "enrichr":
            enrichment_table["log_p_value"] = log_pvalues(enrichment_table["p_value"])

            for gene_set_library in enrichment_table["gene_set_library"].unique():
                _LOGGER.info(gene_set_library)

                # Plot top_n terms of each comparison in barplots
                n = len(enrichment_table[comp_variable].drop_duplicates())
                n_side = int(np.ceil(np.sqrt(n)))

                if barplots:
                    top_data = (
                        enrichment_table[enrichment_table["gene_set_library"] == gene_set_library]
                        .set_index("description")
                        .groupby(comp_variable)
                        ["log_p_value"]
                        .nlargest(top_n)
                        .reset_index())

                    enrichment_barplot(
                        top_data, x="description", y="log_p_value",
                        group_variable=comp_variable, top_n=top_n,
                        output_file=os.path.join(
                            output_dir, output_prefix + ".enrichr.{}.barplot.top_{}.svg"
                            .format(gene_set_library, top_n)))

                    # # ^^ possible replacement
                    # grid = sns.catplot(
                    #     data=top_data, x="log_p_value", y="description",
                    #     order=top_data.groupby("description")["log_p_value"].mean().sort_values(ascending=False).index,
                    #     kind="bar", orient="horiz", col=comp_variable, col_wrap=n_side, palette="magma_r")
                    # grid.savefig(os.path.join(output_dir, output_prefix + ".enrichr.{}.barplot.top_{}.joint_comparisons.svg".format(
                    #         gene_set_library, top_n)), bbox_inches="tight", dpi=300)

                # Scatter plots of Z-score vs p-value vs combined score
                fig, axis = plt.subplots(
                    n_side, n_side,
                    figsize=(4 * n_side, 4 * n_side),
                    sharex=True, sharey=True)
                axis = axis.flatten()
                # normalize color across comparisons
                d = enrichment_table.loc[
                    (enrichment_table["gene_set_library"] == gene_set_library),
                    "combined_score"].describe()
                norm = matplotlib.colors.Normalize(vmin=d["min"], vmax=d["max"])
                for i, comparison in enumerate(enrichment_table[comp_variable].unique()):
                    enr = enrichment_table[
                            (enrichment_table["gene_set_library"] == gene_set_library) &
                            (enrichment_table[comp_variable] == comparison)]
                    sns.scatterplot(
                        data=enr,
                        x="z_score", y="log_p_value", size="combined_score", hue="combined_score",
                        hue_norm=norm,
                        ax=axis[i], rasterized=rasterized, palette="magma")
                    axis[i].set_title(comparison)

                    done = list()
                    for metric in ["log_p_value", "z_score", "combined_score"]:
                        f = pd.DataFrame.head if metric == "z_score" else pd.DataFrame.tail
                        for s in f(enr.sort_values(metric), top_n).index:
                            if enr.loc[s, "description"] not in done:
                                axis[i].text(
                                    enr.loc[s, "z_score"], enr.loc[s, "log_p_value"],
                                    s=enr.loc[s, "description"])
                                done.append(enr.loc[s, "description"])
                sns.despine(fig)
                savefig(fig, os.path.join(output_dir, output_prefix + ".enrichr.{}.zscore_vs_pvalue.scatterplot.svg".format(gene_set_library)))

                # Plot heatmaps of terms for each comparison
                if len(enrichment_table[comp_variable].drop_duplicates()) < 2:
                    continue

                # pivot table
                enrichr_pivot = pd.pivot_table(
                    enrichment_table[enrichment_table["gene_set_library"] == gene_set_library],
                    values="log_p_value", columns="description", index=comp_variable).fillna(0)

                # plot correlation
                if correlation_plots:
                    enrichment_correlation_plot(
                        input_df=enrichr_pivot, label="Correlation of enrichment\nof differential gene sets",
                        output_file=os.path.join(output_dir, output_prefix + ".enrichr.{}.correlation.svg".format(gene_set_library)))

                top = enrichment_table[enrichment_table["gene_set_library"] == gene_set_library].set_index(
                    "description").groupby(comp_variable)["p_value"].nsmallest(top_n)
                top_terms = top.index.get_level_values("description").unique()
                # top_terms = top_terms[top_terms.isin(lola_pivot.columns[lola_pivot.sum() > 5])]

                # plot clustered heatmap
                enrichment_clustermap(
                    enrichr_pivot[list(set(top_terms))].T,
                    output_file=os.path.join(
                        output_dir, output_prefix + ".enrichr.{}.cluster_specific.svg"
                        .format(gene_set_library)),
                    label="-log10(p-value) of enrichment\nof differential genes")
                if z_score is not None:
                    enrichment_clustermap(
                        enrichr_pivot[list(set(top_terms))].T,
                        output_file=os.path.join(
                            output_dir, output_prefix + ".enrichr.{}.cluster_specific.{}_z_score.svg"
                            .format(gene_set_library, z_score_label)),
                        label="{} Z-score of enrichment\nof differential regions".format(z_score_label), z_score=z_score)

        if enrichment_type == "great":
            enrichment_table["log_q_value"] = log_pvalues(enrichment_table["HyperFdrQ"])

            for gene_set_library in enrichment_table["Ontology"].unique():
                _LOGGER.info(gene_set_library)

                # Plot top_n terms of each comparison in barplots
                top_data = (
                    enrichment_table[enrichment_table["Ontology"] == gene_set_library]
                    .set_index("Desc")
                    .groupby(comp_variable)
                    ["log_q_value"]
                    .nlargest(top_n)
                    .reset_index())

                n = len(enrichment_table[comp_variable].drop_duplicates())
                n_side = int(np.ceil(np.sqrt(n)))

                if barplots:
                    fig, axis = plt.subplots(
                        n_side, n_side,
                        figsize=(4 * n_side, n_side * max(5, 0.12 * top_n)),
                        sharex=False, sharey=False)
                    axis = iter(axis.flatten())
                    for i, comp in enumerate(top_data[comp_variable].drop_duplicates().sort_values()):
                        df2 = top_data.loc[top_data[comp_variable] == comp, :]
                        ax = next(axis)
                        sns.barplot(
                            df2["log_q_value"], df2["Desc"], estimator=max,
                            orient="horizontal", ax=ax, color=sns.color_palette("colorblind")[0])
                        ax.set_title(comp)
                    for ax in axis:
                        ax.set_visible(False)
                    sns.despine(fig)
                    savefig(fig, os.path.join(output_dir, output_prefix + ".great.{}.barplot.top_{}.svg".format(gene_set_library, top_n)))

                # Plot heatmaps of terms for each comparison
                if len(enrichment_table[comp_variable].drop_duplicates()) < 2:
                    return

                # pivot table
                great_pivot = pd.pivot_table(
                    enrichment_table[enrichment_table["Ontology"] == gene_set_library],
                    values="log_q_value", columns="Desc", index=comp_variable).fillna(0)

                # plot correlation
                if correlation_plots:
                    enrichment_correlation_plot(
                        input_df=great_pivot, label="Correlation of enrichment\nof differential gene sets",
                        output_file=os.path.join(output_dir, output_prefix + ".great.{}.correlation.svg".format(gene_set_library)))

                top = enrichment_table[enrichment_table["Ontology"] == gene_set_library].set_index(
                    "Desc").groupby(comp_variable)["HyperP"].nsmallest(top_n)
                top_terms = top.index.get_level_values("Desc").unique()

                # plot clustered heatmaps
                enrichment_clustermap(
                    great_pivot[list(set(top_terms))].T,
                    output_file=os.path.join(
                        output_dir, output_prefix + ".great.{}.cluster_specific.svg"
                        .format(gene_set_library)),
                    label="-log10(p-value) of enrichment\nof differential genes")
                if z_score is not None:
                    enrichment_clustermap(
                        great_pivot[list(set(top_terms))].T,
                        output_file=os.path.join(
                            output_dir, output_prefix + ".great.{}.cluster_specific.{}_z_score.svg"
                            .format(gene_set_library, z_score_label)),
                        label="{} Z-score of enrichment\nof differential regions".format(z_score_label), z_score=z_score)
