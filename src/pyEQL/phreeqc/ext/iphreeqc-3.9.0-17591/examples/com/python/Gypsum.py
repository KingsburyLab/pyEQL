"""Compares gypsum solubility for WATEQ4F and Pitzer databases.
"""
# Import standard library modules first.
import os
# Then get third party modules.
from win32com.client import Dispatch
import matplotlib.pyplot as plt

def selected_array(db_path, input_string):
    """Load database via COM and run input string.
    """
    dbase = Dispatch('IPhreeqcCOM.Object')
    dbase.LoadDatabase(db_path)
    dbase.RunString(input_string)
    return dbase.GetSelectedOutputArray()

def show_results(input_string):
    """Get results for different databases
    """
    wateq4f_result = selected_array('wateq4f.dat', input_string)
    pitzer_result  = selected_array('pitzer.dat', input_string)
    # Get data from the arrays.
    nacl_conc      = [entry[0] for entry in wateq4f_result][1:]
    wateq4f_values = [entry[1] for entry in wateq4f_result][1:]
    pitzer_values  = [entry[1] for entry in pitzer_result][1:]
    # Plot
    plt.plot(nacl_conc, pitzer_values, 'k', nacl_conc, wateq4f_values,'k--')
    plt.axis([0, 6, 0, .06])
    plt.legend(('PITZER','WATEQ4F'), loc = (0.4, 0.4))
    plt.ylabel('GYPSUM SOLUBILITY, MOLES PER KILOGRAM WATER')
    plt.xlabel('NaCl, MOLES PER KILOGRAM WATER')
    plt.savefig("Figure2.png")
    plt.show()
    
if __name__ == '__main__':
    # This will only run when called as script from the command line
    # and not when imported from another script.
    INPUT_STRING = """
    SOLUTION 1
    END
    INCREMENTAL_REACTIONS
    REACTION
    	NaCl 1.0
    	0 60*0.1 moles
    EQUILIBRIUM_PHASES
    	Gypsum
    USE solution 1
    SELECTED_OUTPUT
    	-reset false
    	-total Na S(6)
    END"""
    show_results(INPUT_STRING)