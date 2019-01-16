#!/usr/bin/env python

import os

import numpy as np
import pandas as pd

from ngs_toolkit import _LOGGER
from ngs_toolkit import _CONFIG
from ngs_toolkit.decorators import check_organism_genome


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

    :param name: Name of the analysis. Defaults to ``analysis``.
    :type name: str, optional
    :param samples: List of ``peppy.Sample`` objects that this analysis is tied to.
                    Defaults to ``None``.
    :type samples: list, optional
    :param prj: A ``peppy.Project`` object that this analysis is tied to.
                Defaults to ``None``.
    :type prj: peppy.Project, optional
    :param data_dir: Directory containing processed data (e.g. by looper) that will
                     be input to the analysis. This is in principle not required.
                     Defaults to ``data``.
    :type data_dir: str, optional
    :param results_dir: Directory to contain outputs produced by the analysis.
                        Defaults to ``results``.
    :type results_dir: str, optional
    :param pickle_file: A pickle file to serialize the object.
                        Defaults to "`name`.pickle".
    :type pickle_file: str, optional
    :param from_pickle: Whether the analysis should be loaded from an existing
                        serialized analysis object in ``pickle_file``.
                        Defaults to False.
    :type from_pickle: bool, optional

    Additional keyword arguments will simply be stored as object attributes.
    """
    def __init__(
            self,
            name="analysis",
            samples=None,
            prj=None,
            root_dir=None,
            data_dir="data",
            results_dir="results",
            pickle_file=None,
            from_pickle=False,
            **kwargs):
        # parse kwargs with default
        self.name = name
        if root_dir is None:
            self.root_dir = os.curdir
        self.root_dir = os.path.abspath(self.root_dir)

        # # if given absolute paths, keep them, otherwise append to root directory
        for dir_, attr in [(data_dir, 'data_dir'), (results_dir, 'results_dir')]:
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

        _LOGGER.debug("Setting data type-specific attributes to None.")
        attrs = [
            "data_type", "var_names",
            "quantity", "norm_units", "raw_matrix_name",
            "norm_matrix_name", "annot_matrix_name"]
        for attr in attrs:
            if not hasattr(self, attr):
                setattr(self, attr, None)

    def __repr__(self):
        t = "'{}' analysis".format(self.data_type) if self.data_type is not None else "Analysis"
        samples = " with {} samples".format(len(self.samples)) if self.samples is not None else ""
        organism = " of organism '{}'".format(self.organism) if self.organism is not None else ""
        genome = " ({})".format(self.genome) if self.genome is not None else ""
        suffix = "."
        return t + " object '{}'".format(self.name) + samples + organism + genome + suffix

    @staticmethod
    def _overwride_sample_representation():
        from peppy import Sample

        def r(self): return self.name
        Sample.__repr__ = r

    @staticmethod
    def _format_string_with_attributes(obj, string):
        """
        Detect whether a string should be formatted with attributes from obj.
        """
        # TODO: test
        to_format = pd.Series(string).str.extractall(r"{(.*?)}")[0].values
        attrs = obj.__dict__.keys()
        if not all([x in attrs for x in to_format]):
            msg = "Not all required patterns were found as attributes of object '{}'.".format(obj)
            _LOGGER.error(msg)
            raise ValueError(msg)
        return string.format(**obj.__dict__)

    @staticmethod
    def _check_data_type_is_supported(data_type):
        # TODO: test
        supported = _CONFIG['supported_data_types']
        return data_type in supported

    def _get_data_type(self, data_type=None):
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
                raise ValueError
        else:
            msg = "Data type is not supported."
            hint = " Check which data types are supported in the 'supported_data_types'"
            hint += " section of the configuration file."
            if not self._check_data_type_is_supported(data_type):
                raise ValueError(msg + hint)
        return data_type

    def _check_samples_have_file(self, attr, f=all):
        return f([os.path.exists(str(getattr(sample, attr))) for sample in self.samples])

    def _get_samples_have_file(self, attr):
        return [sample for sample in self.samples if os.path.exists(str(getattr(sample, attr)))]

    def _get_samples_missing_file(self, attr):
        return [sample for sample in self.samples if not os.path.exists(str(getattr(sample, attr)))]

    def update(self, pickle_file=None):
        """
        Update all of the object's attributes with the attributes from a serialized
        object (i.e. object stored in a file) object.

        :param pickle_file: Pickle file to load. By default this is the object's attribute
                            `pickle_file`.
        :type pickle_file: str, optional
        """
        self.__dict__.update(self.from_pickle(pickle_file=pickle_file).__dict__)

    def set_organism_genome(self):
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

        :param overwrite: Whether to overwrite attribute values if existing. Defaults to True
        :type overwrite: bool, optional
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
                if isinstance(self.comparison_table, str):
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
        under 'sample_input_files:<data type>:<attribute>'.

        :param overwrite: Whether to overwrite attribute values if existing. Defaults to True
        :type overwrite: bool, optional
        """
        if self.samples is None:
            _LOGGER.error("Analysis object does not have attached Samples. " +
                          "Will not add special attributes to samples such as " +
                          "input file locations.")
            return

        msg = "Setting '{}' in sample {} as '{}'."
        for data_type in _CONFIG['sample_input_files']:
            for attr, value in _CONFIG['sample_input_files'][data_type].items():
                for s in [s for s in self.samples if s.protocol == data_type]:
                    if value is None:
                        pass
                    elif ("{data_dir}" in value) and ("{sample_name}" in value):
                        value = value.format(data_dir=self.data_dir, sample_name=s.name)
                    elif "{data_dir}" in value:
                        value = value.format(data_dir=self.data_dir)
                    elif "{sample_name}" in value:
                        value = value.format(sample_name=s.name)
                    if overwrite:
                        _LOGGER.info(msg.format(attr, s.name, value))
                        setattr(s, attr, value)
                    else:
                        if not hasattr(s, attr):
                            _LOGGER.info(msg.format(attr, s.name, value))
                            setattr(s, attr, value)
                        else:
                            if getattr(s, attr) is None:
                                _LOGGER.info(msg.format(attr, s.name, value))
                                setattr(s, attr, value)
                            else:
                                _LOGGER.debug("{} already exists in sample, not overwriting."
                                              .format(attr.replace("_", " ").capitalize()))

    def to_pickle(self, timestamp=False):
        """
        Serialize object (i.e. save to disk) to hickle format.
        """
        import pickle
        if timestamp:
            import time
            import datetime
            ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d-%H%M%S')
            p = self.pickle_file.replace(".pickle", ".{}.pickle".format(ts))
        else:
            p = self.pickle_file
        pickle.dump(self, open(p, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

    def from_pickle(self, pickle_file=None):
        """
        Load object from pickle file.

        :param pickle_file: Pickle file to load. By default this is the object's
                            attribute `pickle_file`.
        :type pickle_file: str, optional
        """
        import pickle
        if pickle_file is None:
            pickle_file = self.pickle_file
        return pickle.load(open(pickle_file, 'rb'))

    @check_organism_genome
    def get_annotations(
            self,
            organism=None, genome_assembly=None, output_dir=None,
            steps=['blacklist', 'tss', 'genomic_context']):
        """
        Get genome annotations and other resources for several ngs_toolkit analysis.

        :param organism: Organism to get for. Currently supported are 'human' and 'mouse'.
                         Defaults to analysis' own organism.
        :type organism: str, optional
        :param genome_assembly: Genome assembly to get for.
                                Currently supported are 'hg19', 'hg38' and 'mm10'.
                                Defaults to analysis' own genome assembly.
        :type genome_assembly: str, optional
        :param output_dir: Directory to save results to.
                           Defaults to 'reference' in analysis root directory.
        :type output_dir: str, optional
        :param steps: Which steps to get.
                      Options are:
                        - 'fasta': Genome sequence in FASTA format
                        - 'blacklist': Locations of blacklisted regions for genome
                        - 'tss': Locations of gene's TSSs
                        - 'genomic_context': Genomic context of genome
                      Defaults to ['blacklist', 'tss', 'genomic_context']
        :type steps: list, optional
        :returns: Dictionary with keys same as the options as steps, containing
                  paths to the requested files.
        :rtype: dict
        """
        from ngs_toolkit.general import (
            get_genome_reference,
            get_blacklist_annotations,
            get_tss_annotations,
            get_genomic_context)

        if organism is None:
            organism = self.organism
        if genome_assembly is None:
            genome_assembly = self.genome
        if output_dir is None:
            output_dir = os.path.join(self.root_dir, "reference")

        args = {'organism': organism,
                'genome_assembly': genome_assembly,
                'output_dir': output_dir}
        mapping = {"hg19": "grch37", "hg38": "grch38", "mm10": "grcm38"}
        output = dict()

        if 'fasta' in steps:
            output['fasta_file'] = get_genome_reference(**args)
        if 'blacklist' in steps:
            output['blacklist_file'] = get_blacklist_annotations(**args)
        if 'tss' in steps:
            get_tss_annotations(**args)
            output['tss_file'] = os.path.join(
                output_dir, "{}.{}.gene_annotation.protein_coding.tss.bed"
                .format(self.organism, mapping[self.genome]))
        if 'genomic_context' in steps:
            get_genomic_context(**args)
            output['genomic_context_file'] = os.path.join(
                output_dir, "{}.{}.genomic_context.bed"
                .format(self.organism, mapping[self.genome]))

        return output

    def get_matrix(self, matrix=None, matrix_name=None, samples=None):
        """
        Return a matrix that is an attribute of self subsetted for the requested samples.

        :param pandas.DataFrame matrix: Pandas DataFrame.
        :param list samples: Iterable of peppy.Sample objects to restrict matrix to.
                             If not provided (`None` is passed) the matrix will not be subsetted.
        :param str matrix_name: Name of the matrix that is an attribute of the object
                                with values for samples in `samples`.
        :returns pandas.DataFrame: Requested DataFrame.
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

        :param str quant_matrix: Attribute name of matrix to annotate. By default this will
                                 be infered from the analysis data_type in the following way:
                                 ATAC-seq or ChIP-seq: ``coverage_annotated``;
                                 RNA-seq: ``expression_annotated``.
        :param list attributes: Desired attributes to be annotated. This defaults
                                to all attributes in the original sample annotation sheet
                                of the analysis Project.
        :param list numerical_attributes: Attributes which are numeric even though they
                                          might be so in the samples' attributes. Will attempt
                                          to convert values to numeric.
        :param bool save: Whether to write normalized DataFrame to disk.
        :param bool assign: Whether to assign the normalized DataFrame to an attribute
                            ``accessibility`` if ``data_type`` is "ATAC-seq,
                            ``binding`` if ``data_type`` is "ChIP-seq, or
                            ``expression`` if ``data_type`` is "RNA-seq.
        :var pd.DataFrame {accessibility,binding,expression}: A pandas DataFrame with
                                                              MultiIndex column index
                                                              containing the sample's
                                                              attributes specified.
        :returns pd.DataFrame: Annotated dataframe with requested sample attributes.
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
            # TODO: implement dataframe return
            as_dataframe=False):
        """
        Get tuples of floats representing a colour for a sample in a given variable in a
        dataframe's index (particularly useful with MultiIndex dataframes).

        If given, will use the provieded ``index`` argument, otherwise, the the columns
        and its levels of an attribute of self named ``matrix``.
        ``levels`` can be passed to subset the levels of the index.

        Will try to guess if each variable is categorical or numerical and return either colours
        from a colour ``pallete`` or a ``cmap``, respectively with null values set to ``nan_color``
        (a 4-value tuple of floats).

        :param index: Pandas Index to use. If not provided (default == None), this will be
        the column Index of the provided ``matrix``.
        :type index: pandas.Index, optional
        :param matrix: Name of analysis attribute containing a dataframe with pandas.MultiIndex
                       columns to use.
        :type matrix: str, optional
        :param levels: Levels of multiindex to restrict to. Defaults to all in index.
        :type levels: list, optional
        :param pallete: Name of matplotlib color palete to use with categorical levels.
                        See matplotlib.org/examples/color/colormaps_reference.html.
                        Defaults to ``Paired``.
        :type pallete: str, optional
        :param cmap: Name of matplotlib color palete to use with numerical levels.
                     See matplotlib.org/examples/color/colormaps_reference.html.
                     Defaults to ``RdBu_r``.
        :type cmap: str, optional
        :param nan_color: Color for missing (i.e. NA) values.
                          Defaults to ``(0.662745, 0.662745, 0.662745, 0.5)`` == ``grey``.
        :type nan_color: tuple, optional
        :param as_dataframe: Whether a dataframe should be return. Defaults to False.
                             Not implemented yet.
        :type as_dataframe: bool, optional
        :returns: List of list tuples (matrix) of shape (level, sample) with rgb values of
                  each of the variable.
        :rtype: {list}
        """
        import matplotlib
        import matplotlib.pyplot as plt
        from collections import Counter

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
        if isinstance(index, pd.core.indexes.base.Index):
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

        return colors
