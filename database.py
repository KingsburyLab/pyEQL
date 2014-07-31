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
            # ignore the template file
            if file == 'template.csv':
                continue
            
            # TODO - ignore the SMILES file b/c it is not numeric
            # Figure out how to handle non-numeric parameter values
            if file == 'SMILES.csv':
                continue
            
            # open each file
            current_file = open(database_dir+'/'+file,'r')
            
            # read each line of the file, looking for the formula
            for line in current_file:
                if 'Name' in line:
                    # identify the index of the tab separator
                    start_slice = line.index('\t')+1
                    # get the string for the name
                    param_name = line[start_slice:]
                elif 'Description' in line:
                     # identify the index of the tab separator
                    start_slice = line.index('\t')+1
                    # get the string for the Description
                    param_desc = line[start_slice:]
                    
                elif 'Unit' in line:
                     # identify the index of the tab separator
                    start_slice = line.index('\t')+1
                    # get the string for the unit
                    param_unit = line[start_slice:]
                    
                elif 'Reference' in line:
                     # identify the index of the tab separator
                    start_slice = line.index('\t')+1
                    # get the string for the reference, 
                    param_ref = line[start_slice:]

                elif 'Temperature' in line:
                     # identify the index of the tab separator
                    start_slice = line.index('\t')+1
                    # get the string for the temperature,
                    param_temp = line[start_slice:]
                    
                elif 'Pressure' in line:
                     # identify the index of the tab separator
                    start_slice = line.index('\t')+1
                    # get the string for the pressure, 
                    param_press = line[start_slice:]
                    
                elif 'Ionic Strength' in line:
                     # identify the index of the tab separator
                    start_slice = line.index('\t')+1
                    # get the string for the reference, 
                    param_ionic = line[start_slice:]
                    
                elif 'Comment' in line:
                     # identify the index of the tab separator
                    start_slice = line.index('\t')+1
                    # get the string for the reference, 
                    param_comment = line[start_slice:]
                    
                elif formula in line:
                     # identify the index of the tab separator
                    start_slice = line.index('\t')+1
                    # get the string for the value, 
                    try:
                        param_value = float(line[start_slice:])
                    
                        # Create a new parameter object
                        parameter = pm.Parameter(param_name,param_value,param_unit,
                                                 reference=param_ref,pressure=param_press,temperature=param_temp,ionic_strength=param_ionic,description=param_desc,comments=param_comment)
                        
                        # Add the parameter to the set for this species
                        parameters_database[formula].add(parameter)
                    except ValueError:                   
                        pass
                
            current_file.close()

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

    