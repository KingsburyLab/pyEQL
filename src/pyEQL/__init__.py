"""
pyEQL
=====

pyEQL is a python package for calculating the properties of aqueous solutions
and performing chemical thermodynamics computations.

:copyright: 2013-2023 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.
"""
from importlib.metadata import PackageNotFoundError, version  # pragma: no cover

from pint import UnitRegistry
from pkg_resources import resource_filename

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = version(dist_name)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError

# Units handling
# per the pint documentation, it's important that pint and its associated Unit
# Registry are only instantiated once.
# here we assign the identifier 'unit' to the UnitRegistry
# the cache_folder arg is added to speed up registry instantiation
ureg = UnitRegistry(cache_folder=":auto:")
# convert "offset units" so that, e.g. Quantity('25 degC') works without error
# see https://pint.readthedocs.io/en/0.22/user/nonmult.html?highlight=offset#temperature-conversion
ureg.autoconvert_offset_to_baseunit = True
# append custom unit definitions and contexts
fname = resource_filename("pyEQL", "pint_custom_units.txt")
ureg.load_definitions(fname)
# activate the "chemistry" context globally
ureg.enable_contexts("chem")
# set the default string formatting for pint quantities
ureg.default_format = "P~"


from pyEQL.functions import *  # noqa: E402, F403
from pyEQL.solution import Solution  # noqa: E402
