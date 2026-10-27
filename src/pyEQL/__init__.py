"""
pyEQL is a python package for calculating the properties of aqueous solutions
and performing chemical thermodynamics computations.

:copyright: 2013-2024 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.
"""

# The _version.py file is managed by setuptools-scm
#   and is not in version control.
import logging
from importlib.resources import files

from maggma.stores import JSONStore
from pint import UnitRegistry

try:
    from pyEQL._version import version as __version__
except ImportError:
    __version__ = "src"

# logging
logger = logging.getLogger("pyEQL")
logger.setLevel(logging.WARNING)
logger.addHandler(logging.NullHandler())


# Units handling
# per the pint documentation, it's important that pint and its associated Unit
# Registry are only instantiated once.
# here we assign the identifier 'unit' to the UnitRegistry
# the cache_folder arg is added to speed up registry instantiation
ureg = UnitRegistry(cache_folder=":auto:")
# convert "offset units" so that, e.g. Quantity('25 degC') works without error
# see https://pint.readthedocs.io/en/0.22/user/nonmult.html?highlight=offset#temperature-conversion
ureg.autoconvert_offset_to_baseunit = True
# NOTE: the "chemistry" context (which enables mole <-> mass conversions using a
# molecular weight) is NOT enabled globally. Starting with pint 0.25, Quantity.check()
# and is_compatible_with() honor the enabled contexts, so a globally active "chemistry"
# context makes [substance] and [mass] appear mutually compatible. That breaks the unit
# classification in Solution.get_amount(), where e.g. Quantity("kg").check("[substance]")
# would return True and send a mass request down the substance branch. Every conversion
# that actually needs the context passes it explicitly, e.g. qty.to("kg", "chem", mw=mw).
# set the default string formatting for pint quantities
ureg.formatter.default_format = "P~"

# create a Store for the default database
json_db_file = files("pyEQL") / "database" / "pyeql_db.json"
IonDB = JSONStore(str(json_db_file), key="formula", encoding="utf8")
# By calling connect on init, we get the expensive JSON reading operation out
# of the way. Subsequent calls to connect will bypass this and access the already-
# instantiated Store in memory, which should speed up instantiation of Solution objects.
IonDB.connect()

from pyEQL.solution import Solution  # noqa: E402
