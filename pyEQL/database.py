"""
This module contains classes, functions, and methods for reading input files
and assembling database entries for use by pyEQL.

By default, pyEQL searches all files in the /database subdirectory for parameters.

:copyright: 2013-2018 by Ryan S. Kingsbury
:license: LGPL, see LICENSE for more details.

"""

# logging system
import logging

logger = logging.getLogger(__name__)

# add a filter to emit only unique log messages to the handler
from pyEQL.logging_system import Unique

unique = Unique()
logger.addFilter(unique)

# add a handler for console output, since pyEQL is meant to be used interactively
ch = logging.StreamHandler()

# create formatter for the log
formatter = logging.Formatter("(%(name)s) - %(levelname)s - %(message)s")

# add formatter to the handler
ch.setFormatter(formatter)
logger.addHandler(ch)

# for parameter creation functions
import pyEQL.parameter as pm

# for file input/output functions
import os


## Database Management Functions


class Paramsdb:
    """
    create a global dictionary to contain a dynamically-generated list of Parameters
    for solute species. The dictionary keys are the individual chemical species
    formulas. The dictionary's values are a python set object containing all parameters
    that apply to the species.
    """

    def __init__(self):

        self.parameters_database = {}

        # set the directory containing database files
        self.database_dir = [os.path.dirname(__file__) + "/database"]

    def add_path(self, path):
        """
        Add a user-defined directory to the database search path
        """
        self.database_dir.append(path)

    def list_path(self):
        """
        List all search paths for database files
        """
        for path in self.database_dir:
            print(path)

    def search_parameters(self, formula):
        """Each time a new solute species is created in a solution, this function:
        
        1) searches to see whether a list of parameters for the species has already been
        compiled from the database
        2) searches all files in the specified database directory(ies) for the species
        3) creates a Parameter object for each value found 
        4) compiles these objects into a set
        5) adds the set to a dictionary indexed by species name (formula)
        6) points the new solute object to the dictionary
        
        formula : str
                String representing the chemical formula of the species.
        """
        # get the hill_order() and is_valid_formula() methods from the chemistry module
        import pyEQL.chemical_formula as chem

        # if the formula is already in the database, then we've already searched
        # and compiled parameters, so there is no need to do it again.
        if formula in self.parameters_database:
            pass
        else:
            # add an entry to the parameters database
            self.parameters_database[formula] = set()

            # search all the files in each database directory
            for directory in self.database_dir:
                for file in os.listdir(directory):
                    # reset the line number count
                    line_num = 0

                    # ignore the template file
                    if file == "template.tsv":
                        continue

                    # look at only .csv files
                    if ".tsv" in file:

                        # open each file
                        current_file = open(directory + "/" + file, "r")

                        # read each line of the file, looking for the formula
                        try:
                            for line in current_file:

                                line_num += 1

                                try:
                                    # look for keywords in the first column of each file. If found,
                                    # store the entry from the 2nd column
                                    if "Name" in line:
                                        param_name = _parse_line(line)[1]

                                    elif "Description" in line:
                                        param_desc = _parse_line(line)[1]

                                    elif "Unit" in line:
                                        param_unit = _parse_line(line)[1]

                                    elif "Reference" in line:
                                        param_ref = _parse_line(line)[1]

                                    elif "Temperature" in line:
                                        param_temp = _parse_line(line)[1]

                                    elif "Pressure" in line:
                                        param_press = _parse_line(line)[1]

                                    elif "Ionic Strength" in line:
                                        param_ionic = _parse_line(line)[1]

                                    elif "Comment" in line:
                                        param_comment = _parse_line(line)[1]

                                    # use the hill_order() function to standardize the
                                    # supplied formula. Then standardize teh formulas in the
                                    # database and see if they match.
                                    # this allows a database entry for 'MgCl2' to be matched
                                    # even if the user enters 'Mg(Cl)2', for example
                                    elif chem.is_valid_formula(_parse_line(line)[0]):
                                        if chem.hill_order(formula) == chem.hill_order(
                                            _parse_line(line)[0]
                                        ):
                                            # if there are multiple columns, pass the values as a list.
                                            # If a single column, then just pass the value
                                            if len(_parse_line(line)) > 2:
                                                param_value = _parse_line(line)[1:]
                                            else:
                                                param_value = _parse_line(line)[1]

                                            # Create a new parameter object
                                            parameter = pm.Parameter(
                                                param_name,
                                                param_value,
                                                param_unit,
                                                reference=param_ref,
                                                pressure=param_press,
                                                temperature=param_temp,
                                                ionic_strength=param_ionic,
                                                description=param_desc,
                                                comment=param_comment,
                                            )

                                            # Add the parameter to the set for this species
                                            self.parameters_database[formula].add(
                                                parameter
                                            )

                                except ValueError:
                                    logger.warning(
                                        "Error encountered when reading line %s in %s"
                                        % (line_num, file)
                                    )
                                    continue

                        # log a warning if an invalid character prevents reading a line
                        except UnicodeDecodeError:
                            logger.warning(
                                "Invalid character found when reading %s. File skipped."
                                % file
                            )

                        current_file.close()

    def has_parameter(self, formula, name):
        """
        Boolean test to determine whether a parameter exists in the database for a given species
        """
        # search the available database files if the species is not currently
        # in the database
        if self.has_species(formula) is False:
            self.search_parameters(formula)

        found = False

        try:
            for item in self.parameters_database[formula]:
                if item.get_name() == name:
                    found = True
                else:
                    continue

            return found

        except KeyError:
            logger.error("Species %s not found in database" % formula)
            return None

    def get_parameter(self, formula, name):
        """
        Retrieve a parameter from the database
        """
        found = False

        try:
            for item in self.parameters_database[formula]:
                if item.get_name() == name:
                    found = True
                    return item

            if found == False:
                logger.error(
                    "Parameter %s for species %s not found in database"
                    % (name, formula)
                )

        except KeyError:
            logger.error("Species %s not found in database" % formula)
            return None

    def add_parameter(self, formula, parameter):
        """
        Add a parameter to the database
        """
        self.parameters_database[formula].add(parameter)

    def has_species(self, formula):
        """
        Boolean test to determine whether a species is present in the database
        """
        if formula in self.parameters_database:
            return True
        else:
            return False

    def print_database(self, solute=None):
        """ Function to generate a human-friendly summary of all the database parameters
        that are actually used in the simulation
        
        Parameters
        ----------
        solute : str, optional
                The chemical formula for a species. If this argument of supplied, the output
                will contain only the database entries for this species. Otherwise,
                all database entries will be printed.
      
        """
        if solute is not None:
            try:
                key = solute
                print("Parameters for species %s:" % key)
                print("--------------------------\n")
                for item in self.parameters_database[key]:
                    print(item)
            except KeyError:
                print("Species %s not found in database." % solute)
        else:
            for key in self.parameters_database.keys():
                print("Parameters for species %s:" % key)
                print("--------------------------\n")
                for item in self.parameters_database[key]:
                    print(item)


def _parse_line(line):
    """
    Function to parse lines in a tab-seprated value file format.
    
    This function accepts a string (a line read from a tab-separated
    input file). It removes the newline character and splits the string
    at each tab stop, returning a list of the remaining substrings in which each
    list entry corresponds to the contents of one cell in the file.
    
    """
    # remove the newline character
    line = line.replace("\n", "")

    # separate the string at every tab stop
    str_list = line.split("\t")

    # return the list of string entries
    return str_list
