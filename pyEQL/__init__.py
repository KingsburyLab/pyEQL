'''
pyEQL
=====

pyEQL is a python package for calculating the properties of aqueous solutions
and performing chemical thermodynamics computations.

:copyright: 2013-2015 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

'''
# initialize the parameters database
from pyEQL.database import Paramsdb
paramsDB = database.Paramsdb()

from pyEQL.parameter import unit
from pyEQL.functions import *
from pyEQL.solution import Solution
