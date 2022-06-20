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
## Units handling
# per the pint documentation, it's important that pint and its associated Unit
# Registry are only imported once.
from pint import UnitRegistry

# here we assign the identifier 'unit' to the UnitRegistry
unit = UnitRegistry()

# use this to enable legacy handling of offset units
# TODO fix this to handle offsets the way pint wants us to since 0.7
unit.autoconvert_offset_to_baseunit = True

# append custom unit definitions and contexts
from pkg_resources import resource_filename

fname = resource_filename("pyEQL", "pint_custom_units.txt")
unit.load_definitions(fname)
# activate the "chemistry" context globally
unit.enable_contexts("chem")
# set the default string formatting for pint quantities
unit.default_format = "P~"

from pyEQL.database import Paramsdb

# initialize the parameters database
paramsDB = Paramsdb()

from pyEQL.functions import *
from pyEQL.solution import Solution
