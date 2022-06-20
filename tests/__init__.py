import logging
import os
import unittest


def testsuite():
    """A testsuite that includes all the pyEQL tests."""
    # turn off log warning messages during tests
    logging.disable(logging.WARNING)

    # include all files in the directory of this script
    suite = unittest.TestLoader().discover(os.path.dirname(__file__))

    return suite


def main():
    """Runs the testsuite as command line application."""
    try:
        unittest.main()
    except Exception as e:
        print("Error: %s" % e)


def run():
    """Run all tests.
    :return: a :class:`unittest.TestResult` object
    """
    test_runner = unittest.TextTestRunner()
    return test_runner.run(testsuite())


# define custom assertion functions to compare model output with experimental
# data when running units tests.
# See https://stackoverflow.com/questions/6655724/how-to-write-a-custom-assertfoo-method-in-python
# for the method I'm using here.
class CustomAssertions:
    def assertWithinExperimentalError(self, result, expected, tol=0.05):
        """
        Test whether 'result' is within 'tol' relative error of
        'expected'
        """
        rel_error = abs(result - expected) / expected
        if not rel_error < tol:
            raise AssertionError(
                "Result {:} differs from expected value by {:.2f}%".format(
                    result, rel_error * 100
                )
            )
