from ._version import __version__


def setup_logger(level="INFO", logfile=None):
    """
    Set up a logger for the library.

    :param level: Level of logging to display. Defaults to "INFO". See possible levels here: https://docs.python.org/2/library/logging.html#levels
    :type level: str, optional
    :param logfile: File to write log to. Defaults to ``~/ngs_toolkit.log.txt``
    :type logfile: str, optional
    :returns: The logger
    :rtype: logging.Logger
    """
    import logging
    import os

    global _LOGGER

    _LOGGER = logging.getLogger('ngs_toolkit')
    _LOGGER.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    if logfile is None:
        logfile = os.path.join(os.path.expanduser("~"), ".ngs_toolkit.log.txt")
    fh = logging.FileHandler(logfile)
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.getLevelName(level))
    # create formatter and add it to the handlers
    fmt = ("ngs_toolkit.v{}:%(module)s:L%(lineno)d (%(funcName)s) [%(levelname)s] %(asctime)s > %(message)s"
           .format(__version__))
    formatter = logging.Formatter(fmt, datefmt='%Y-%m-%d %H:%M:%S')
    fh.setFormatter(formatter)
    fmt = ("ngs_toolkit.v{}:%(module)s:L%(lineno)d (%(funcName)s) [%(levelname)s] > %(message)s"
           .format(__version__))
    formatter = logging.Formatter(fmt)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    _LOGGER.addHandler(fh)
    _LOGGER.addHandler(ch)

    _LOGGER.info("This is ngs_toolkit (http://ngs-toolkit.rtfd.io), version: {}"
                 .format(__version__))
    return _LOGGER


def setup_config(custom_yaml_config=None):
    """
    Set up global library configuration.

    It reads ngs_toolkit's package data to load a default configuration,
    tries to update it by reading a file in ``~/.ngs_toolkit.config.yaml`` if present,
    and lastly, updates it by reading a possible passed yaml file ``custom_yaml_config``.

    :returns: Dictionary with configurations
    :rtype: dict
    """
    import pkg_resources
    import os
    import yaml

    global _CONFIG

    default_config_path = 'config/default.yaml'
    default_config_path = pkg_resources.resource_filename(__name__, default_config_path)
    _LOGGER.info("Reading default configuration file distributed with package from '{}'.".format(default_config_path))
    try:
        _CONFIG = yaml.load(open(default_config_path, 'r'))
        _LOGGER.debug("Default config: {}".format(_CONFIG))
    except IOError:
        _LOGGER.error("Couldn't read configuration file from '{}'.".format(default_config_path))
        _CONFIG = dict()

    user_config_path = os.path.join(os.path.expanduser("~"), ".ngs_toolkit.config.yaml")
    if os.path.exists(user_config_path):
        # Read up
        _LOGGER.info("Found custom user config: {}".format(user_config_path))
        try:
            custom_config = yaml.load(open(user_config_path, "r"))
            _LOGGER.debug("Custom user config: {}".format(custom_config))
            # Update config
            _LOGGER.info("Updating configuration with custom file from '{}'.".format(user_config_path))
            _CONFIG.update(custom_config)
            _LOGGER.debug("Current config: {}".format(custom_config))
        except IOError:
            _LOGGER.error("Configuration file from '{}' exists but is not readable. Ignoring.".format(user_config_path))
    else:
        _LOGGER.info("To use custom configurations including paths to static files, create a '{}' file.".format(user_config_path))

    if custom_yaml_config is not None:
        # Read up
        try:
            custom_config = yaml.load(open(custom_yaml_config, "r"))
            _LOGGER.debug("Custom passed config: {}".format(custom_config))
            # Update config
            _LOGGER.info("Updating configuration with custom file from '{}'.".format(custom_yaml_config))
            _CONFIG.update(custom_config)
            _LOGGER.debug("Current config: {}".format(custom_config))
        except IOError as e:
            _LOGGER.error("Pased configuration file from '{}' exists but is not readable.".format(custom_yaml_config))
            raise e

    return _CONFIG


setup_logger()
setup_config()
