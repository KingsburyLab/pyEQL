'''
pyEQL salt matching library

This file contains functions that allow a pyEQL Solution object composed of
individual species (usually ions) to be mapped to a solution of one or more
salts. This mapping is necessary because some parameters (such as activity
coefficient data) can only be determined for salts (e.g. NaCl) and not individual
species (e.g. Na+)

'''
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


