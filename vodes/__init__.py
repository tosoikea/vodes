"""Top-level package for vodes."""

__author__ = """Torben Soennecken"""
__email__ = 'soennecken@rootitup.de'
__version__ = '0.2.0'

import os

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


# Initialize logging handlers
import logging
import logging.config
import yaml

path = os.path.join(__location__, 'logging.yml')

if os.path.exists(path):
    with open(path, 'r') as f:
        config = yaml.safe_load(f.read())
        logging.config.dictConfig(config)

logger = logging.getLogger(__name__)
logger.debug('Package successfully loaded.')