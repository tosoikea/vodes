"""Top-level package for vodes."""

__author__ = """Torben Soennecken"""
__email__ = 'soennecken@rootitup.de'
__version__ = '0.1.1'

import os

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


# Initialize logging handlers
import logging
import logging.config
import yaml

with open(os.path.join(__location__, 'logging.yml'), 'r') as f:
    config = yaml.safe_load(f.read())
    logging.config.dictConfig(config)

logger = logging.getLogger(__name__)

logging.debug('Package successfully loaded.')