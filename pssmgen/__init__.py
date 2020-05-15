from .pssm import PSSM
from .__version__ import __version__
import logging

LOG_FORMAT = "%(message)s"
logging.basicConfig(level=logging.DEBUG, format=LOG_FORMAT)