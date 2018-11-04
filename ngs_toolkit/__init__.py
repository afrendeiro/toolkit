from ._version import __version__
import logging
import sys

def setup_logger(level=1):
	global _LOGGER
	_LOGGER = logging.getLogger("ngs_toolkit")
	handler = logging.StreamHandler(sys.stdout)
	handler.setLevel(level)
	if level != 1:
		fmt = "%(module)s:%(lineno)d (%(funcName)s) [%(levelname)s] > %(message)s"
	else:
		fmt = "[%(levelname)s] > %(message)s"
	formatter = logging.Formatter(fmt=fmt)
	handler.setFormatter(formatter)
	_LOGGER.addHandler(handler)
	return _LOGGER


__version__

setup_logger()
