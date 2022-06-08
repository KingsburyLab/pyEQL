import sys

if sys.version_info[:2] >= (3, 8):
    # TODO: Import directly (no need for conditional) when `python_requires = >= 3.8`
    from importlib.metadata import PackageNotFoundError, version  # pragma: no cover
else:
    from importlib_metadata import PackageNotFoundError, version  # pragma: no cover

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = version(dist_name)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError

"""
pyEQL
=====

pyEQL is a python package for calculating the properties of aqueous solutions
and performing chemical thermodynamics computations.

:copyright: 2013-2018 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

"""
from pyEQL.database import Paramsdb

# initialize the parameters database
paramsDB = Paramsdb()

from pyEQL.parameter import unit
from pyEQL.functions import *
from pyEQL.solution import Solution
