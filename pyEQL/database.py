# -*- coding: utf-8 -*-
"""
This module contains classes, functions, and methods for reading input files
and assembling database entries for use by pyEQL.

By default, pyEQL searches all files in the /database subdirectory for parameters.

File Format:
-----------

Database files are tab-separated text files. The first column must contain a key (usually a chemical formula), 
while the actual parameter values are in the 2nd and following columns. Each database
file is intended to contain a single parameter type from a single source for multiple chemical species.
I.e., we have one database file that contains diffusion coefficients for a bunch ofions as found
in the CRC Handbook of Chemistry and Physics.

The top of each database file must, at a minimum, contain rows for 'Name' and 'Units'. Preferably,
other information such as conditions, notes and a reference are also supplied. See 'template.csv' in the
/database directory for an example.

Multi-member parameters (e.g. coefficients for an equation) can be defined by using mutiple
columns. Multiple columns are only searched for the actual parameter, not for the header
rows for Name, Units, etc.

"""

## Logging System
''' Create a logging system using Python's built-in module. 
Add the null handler to avoid errors in case the calling application doesn't configure any handlers.

NOTE: make sure to set the disable_existing_loggers option in the log configuration
options of the calling application in order to avoid disabling the pyEQL module's log
 
The default logging levels are mapped to pyEQL events as follows:
 
DEBUG       -   detailed messages about function execution including methods used, data sources,
                temperature adjustments, etc.
INFO        -   Messages indicating calculation steps, function calls, etc.
WARNING     -   assumptions or limitations of module output
ERROR       -   Module could not complete a task due to invalid input or other problem
CRITICAL    -   not used

'''
import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# for parameter creation functions
import parameter as pm
# for file input/output functions
import os


## Database Management Functions
''''
create a global dictionary to contain a dynamically-generated list of Parameters
for solute species. The dictionary keys are the individual chemical species
formulas. The dictionary's values are a python set object containing all parameters
that apply to the species.
'''
parameters_database={}

# set the directory containing database files
database_dir = os.getcwd()+'/database'

        
def search_parameters(formula):
    '''Each time a new solute species is created in a solution, this function:
    
    1) searches to see whether a list of parameters for the species has already been
    compiled from the database
    2) searches all files in the specified database directory(ies) for the species
    3) creates a Parameter object for each value found 
    4) compiles these objects into a set
    5) adds the set to a dictionary indexed by species name (formula)
    6) points the new solute object to the dictionary
    
    formula : str
            String representing the chemical formula of the species.
    '''
    # get the hill_order() and is_valid_formula() methods from the chemistry module
    import chemical_formula as chem
    
    # if the formula is already in the database, then we've already searched
    # and compiled parameters, so there is no need to do it again.    
    if formula in parameters_database:
        pass
    else:
        # add an entry to the parameters database
        parameters_database[formula] = set()
        
        # search all the files in the database directory
        for file in os.listdir(database_dir):
            # reset the line number count
            line_num = 0
            
            # ignore the template file
            if file == 'template.csv':
                continue
            
            # open each file
            current_file = open(database_dir+'/'+file,'r')
            
            # read each line of the file, looking for the formula
            try:            
                for line in current_file:
                    
                    line_num += 1
                    
                    try:
                        # look for keywords in the first column of each file. If found,
                        # store the entry from the 2nd column
                        if 'Name' in line:
                            param_name = _parse_line(line)[1]
                            
                        elif 'Description' in line:
                            param_desc = _parse_line(line)[1]
                            
                        elif 'Unit' in line:
                            param_unit = _parse_line(line)[1]
                            
                        elif 'Reference' in line:
                            param_ref = _parse_line(line)[1]
        
                        elif 'Temperature' in line:
                            param_temp = _parse_line(line)[1]
                            
                        elif 'Pressure' in line:
                            param_press = _parse_line(line)[1]
                            
                        elif 'Ionic Strength' in line:
                            param_ionic = _parse_line(line)[1]
                            
                        elif 'Comment' in line:
                            param_comment = _parse_line(line)[1]
                            
                        
                        # use the hill_order() function to standardize the 
                        # supplied formula. Then standardize teh formulas in the
                        # database and see if they match.
                        # this allows a database entry for 'MgCl2' to be matched
                        # even if the user enters 'Mg(Cl)2', for example
                        elif chem.is_valid_formula(_parse_line(line)[0]):
                            if chem.hill_order(formula) == chem.hill_order(_parse_line(line)[0]):                                                
                                # if there are multiple columns, pass the values as a list. 
                                # If a single column, then just pass the value
                                if len(_parse_line(line)) >2:
                                    param_value = _parse_line(line)[1:]
                                else:
                                    param_value = _parse_line(line)[1]
                                                                            
                                # Create a new parameter object
                                parameter = pm.Parameter(param_name,param_value,param_unit, \
                                reference=param_ref,pressure=param_press,temperature=param_temp,ionic_strength=param_ionic,description=param_desc,comment=param_comment)
                                
                                # Add the parameter to the set for this species
                                parameters_database[formula].add(parameter)
                        
                    except ValueError:                   
                        logger.warning('Error encountered when reading line %s in %s' % (line_num,file))
                        continue
            
            # log a warning if an invalid character prevents reading a line
            except UnicodeDecodeError:
                    logger.warning('Invalid character found when reading %s. File skipped.' % file)
                
            current_file.close()

def _parse_line(line):
    '''
    Function to parse lines in a tab-seprated value file format.
    
    This function accepts a string (a line read from a tab-separated
    input file). It removes the newline character and splits the string
    at each tab stop, returning a list of the remaining substrings in which each
    list entry corresponds to the contents of one cell in the file.
    
    '''
    # remove the newline character
    line = line.replace('\n','')                        
    
    # separate the string at every tab stop
    str_list = line.split('\t')
    
    # return the list of string entries
    return str_list
    

def print_database():
    ''' Function to generate a human-friendly summary of all the database parameters
    that are actually used in the simulation
  
    '''
    for key in parameters_database.keys():
        print('Parameters for species %s:' % key)
        print('--------------------------')
        for item in parameters_database[key]:
            print(item)
    
    return None
    

    