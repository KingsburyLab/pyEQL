import unittest
import os
import logging

def testsuite():
    """A testsuite that includes all the pyEQL tests.
    """
    # turn off log warning messages during tests
    logging.disable(logging.WARNING)
    
    # include all files in the directory of this script
    suite = unittest.TestLoader().discover(os.path.dirname(__file__))

    return suite


def main():
    """Runs the testsuite as command line application.
    """
    try:
        unittest.main()
    except Exception as e:
        print('Error: %s' % e)


def run():
    """Run all tests.
    :return: a :class:`unittest.TestResult` object
    """
    test_runner = unittest.TextTestRunner()
    return test_runner.run(testsuite())
