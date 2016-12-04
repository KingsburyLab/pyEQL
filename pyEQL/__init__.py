'''
pyEQL
=====

pyEQL is a python package for calculating the properties of aqueous solutions
and performing chemical thermodynamics computations.

:copyright: 2013-2016 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

'''
# initialize the parameters database
from pyEQL.database import Paramsdb
paramsDB = database.Paramsdb()

from pyEQL.parameter import unit
from pyEQL.functions import *
from pyEQL.solution import Solution

# define custom assertion functions to compare model output with experimental
# data when running units tests.
# See https://stackoverflow.com/questions/6655724/how-to-write-a-custom-assertfoo-method-in-python
# for the method I'm using here. 
class CustomAssertions:
    def assertWithinExperimentalError(self,result,expected,tol=0.05):
        '''
        Test whether 'result' is within 'tol' relative error of
        'expected'
        '''
        rel_error = abs(result-expected)/expected
        if not rel_error < tol:
            raise AssertionError('Result {:} differs from expected value by {:.2f}%'.format(result,rel_error*100))

# enable easy testing
def test():
    """Run all tests.
    :return: a :class:`unittest.TestResult` object
    """
    from .tests import run
    return run()
