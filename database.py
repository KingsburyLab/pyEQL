# -*- coding: utf-8 -*-
"""
This module contains classes, functions, and methods for reading input files
and assembling database entries for use by pyEQL.

@author: ryan
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
# create a global dictionary to contain a dynamically-generated list of Parameters
# for solute species. The dictionary keys are the individual chemical species
# formulas. The dictionary's values are a python set object containing all parameters
# that apply to the species.
parameters_database={}

# set the directory containing database files
database_dir = os.getcwd()+'/database'

        
def search_parameters(formula):
    '''Each time a new solute species is created in a solution, this function:
    
    formula : str
            String representing the chemical formula of the species.
    
    1) searches to see whether a list of parameters for the species has already been
    compiled from the database
    2) searches all files in the specified database directory(ies) for the species
    3) creates a Parameter object for each value found 
    4) compiles these objects into a list
    5) adds the list to a dictionary indexed by species name (formula)
    6) points the new solute object to the dictionary
    
    
    pseudo code:
    
    if formula = 'H2O':
            # call some special methods to calculate parameters
    
    if parameters_database[self.name]:
        # point self to the database for parameters
    else:
        # Add  formula:[] entry to parameters_database dictionary
        # Search all the files in the database
        parameter_list = []
        for file in database_dir:
            current_file = open(file,r)
            if formula is in the file:
            #       create a parameter object
                new_parameter = Parameter(FILL IN)
            #      append object to list
                parameter_list.append(new_parameter)
                parameters_database.update({formula:parameter_list})

    '''
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
            for line in current_file:
                
                line_num += 1
                
                try:
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
                        
                    # make sure to only read values when the formula is followed by a tab
                    # stop. Otherwise a search for 'NaCl' will also read parameters for
                    # 'NaClO4', for example.
                    elif formula+'\t' in line:                                                
                        # convert values into a tuple (Excluding the formula) and pass
                        # to the Parameter class
                        param_value = tuple(_parse_line(line)[1:])
                                                                    
                        # Create a new parameter object
                        parameter = pm.Parameter(param_name,param_value,param_unit, \
                        reference=param_ref,pressure=param_press,temperature=param_temp,ionic_strength=param_ionic,description=param_desc,comments=param_comment)
                        
                        # Add the parameter to the set for this species
                        parameters_database[formula].add(parameter)
                        
                except ValueError:                   
                    logger.warning('Error encountered when reading line %s in %s' %(line_num,file))
            
            current_file.close()

def _parse_line(line):
    '''
    Function to parse lines in a tab-seprated value file format
    
    '''
    # remove the newline character
    line = line.replace('\n','')                        
    
    # separate the string at every tab stop
    str_list = line.split('\t')
    
    # return the list of string entries
    return str_list
    

def print_database(sort_type):
    ''' Function to generate a human-friendly summary of all the database parameters
    that are actually used in the simulation
    
    sort_type : str
                String that specifies whether to group the output by parameter name
                e.g. 'diffusion coefficient' or by species. Acceptable values are 
                'name' and 'species'.
    
    1. Loop through the parameters_database dictionary
    2. group parameters either by type or species
    3. generate a list, one line each, with value, units, reference, conditions
    
    
    '''
    # TODO
    if sort_type == 'name':
        pass
    elif sort_type == 'species':
        pass
    pass

    