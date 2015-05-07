.. _database:


Database System
***************

pyEQL creates a database to collect various parameters needed to perform
it's calculations. pyEQL's default database includes a collection of the
following parameters for some common electrolytes:

 - Diffusion coefficients for 104 ions
 - Pitzer model activity correction coefficients for 157 salts
 - Pitzer model partial molar volume coefficients for 120 salts
 - Partial molar volumes for 10 ions (see note)
 - Viscosity model coefficients for 6 salts (see note)

.. note:: Due to copyright restrictions, pyEQL's built-in databases contain only a small selection of partial molar volumes and viscosity coefficients for some common ions like H+, Na+, Cl-, and OH-. We are working on securing permission to distribute a more complete dataset. In the mean time, see the references in the example databases for good data sources. Alternatively, you can provide your own parameters in a custom database (see below). pyEQL does already contain a fairly large collection of Pitzer parameters for both activity correction and partial molar volume; and this will be expanded in the future.

Basics
======

The Paramsdb class creates a container for parameters. Each paramter
is an object which contains not only the value, but also information about
the units, the reference, and the conditions of measurement. paramsdb() also
defines several methods that are helpful for retrieving parameters.

pyEQL automatically initializes an instance of Paramsdb under the name 'db'.
You can access database methods like this:

.. doctest::

   >>> import pyEQL
   >>> pyEQL.db
   <pyEQL.database.Paramsdb at 0x7fead183f240>
   >>> pyEQL.db.has_species('H+')
   True

Anytime a new solute is added to a solution, the search_parameters() method
is called. This method searches every database file within the search path
(by default, only pyEQL's built-in databases) for any parameters associated
with that solute, and adds them to the database.

Adding your own Database Files
==============================

Custom Search Paths
-------------------

The database system is meant to be easily extensible. To include your own
parameters, first you need to add a directory of your choosing to the
search path.

.. doctest::

   >>> pyEQL.db.add_path('/home/user')

You can always check to see which paths pyEQL is searching by using list_path():

.. doctest::

   >>> pyEQL.db.list_path()
   <default installation directory>/database
   /home/user

Then, place your custom database file inside that directory. **NOTE: custom
database files are searched IN ADDITION TO the default databases.** You don't
need to re-create the information from the built-in files. Custom databases
only need to contain extra parameters that are not included already.

File Format
-----------

Databases are formatted as TAB-SEPARATED text files and carry the .tsv extension.
The intent of this format is to make database files easy to edit with common
spreadsheet software. 

.. warning:: If you open an existing or template database file for editing, some spreadsheet software will try to replace the tabs with commas when you save it again. pyEQL does NOT read comma-separated files.

Since pyEQL compiles the database from multiple files, the intent is for each
file to contain values for one type of parameter (such as a diffusion coefficient) 
from one source. The file can then list values of that parameter for a number of
different solutes.

The upper section of each file contains information about the source of the
data, the units, the name of the parameter, and the conditions of measurement.
The top of each database file must, at a minimum, contain rows for 'Name' and 'Units'. 
Preferably, other information such as conditions, notes and a reference are also supplied.
See `template.tsv` in the \database subdirectory for an example.

The remainder of the file contains solute formulas in the first column (see
:ref:`chemistry`) and corresponding values of the parameter in the following columns.
Sets of parameters (such as activity correction coefficients) can be specified
by using more than one column.

.. warning:: Currently there is no way to handle duplicated parameters. So if you supply a parameter with the same name as a built-in one, unexpected behavior may result.

Special Names
-------------
The name of a parameter is used as a kind of index within pyEQL. Certain methods
expect certain parameter names. The following are the currently-used internal
names:

 - 'diffusion_coefficient' - diffusion coefficient
 - 'pitzer_parameters_activity' - coefficients for the Pitzer model for activity correction
 - 'pitzer_parameters_volume'- coefficients for the Pitzer model for partial molar volume
 - 'erying_viscosity_coefficients' - coefficients for an Erying-type viscosity correction model
 - 'partial_molar_volume'- the partial molar volume (used if Pitzer parameters are not available)

If you wish to supply these paramaters for a custom solute not included in the built-in
database, make sure to format the name exactly the same way.

You can also specify a custom parameter name, and retrieve it using the get_parameter()
method. If the solute is 'Na+'

.. doctest::

   >>> pyEQL.db.get_parameter('Na+','my_parameter_name')

Viewing the Database
====================

You can view the entire contents of the database using the print_database() method.
Since pyEQL searches for parameters as they are added, the database will only
contain parameters for solutes that have actually been used during the execution
of your script. The output is organized by solute.

.. doctest::
   
   >>> pyEQL.db.print_database()
   
   >>> s1 = pyEQL.Solution([['Na+','0.5 mol/kg'],['Cl-','0.5 mol/kg']])
   >>> pyEQL.db.print_database()
   Parameters for species Cl-:
   --------------------------
   Parameter diffusion_coefficient
   Diffusion Coefficient
   -------------------------------------------
   Value: 2.032e-05 cm²/s
   Conditions (T,P,Ionic Strength): 25 celsius, 1 atm, 0
   Notes: For most ions, increases 2-3% per degree above 25C
   Reference: CRC Handbook of Chemistry and Physics, 92nd Ed., pp. 5-77 to 5-79
   
   Parameter partial_molar_volume
   Partial molar volume
   -------------------------------------------
   Value: 21.6 cm³/mol
   Conditions (T,P,Ionic Strength): 25 celsius, 1 atm, 0
   Notes: correction factor 5e-4 cm3/g-K
   Reference: Durchschlag, H., Zipper, P., 1994. "Calculation of the Partial Molal Volume of Organic Compounds and Polymers." Progress in Colloid & Polymer Science (94), 20-39.
   ...

API Documentation (database.py)
===============================

.. automodule:: pyEQL.database
   :members:

API Documentation (parameter.py)
================================

.. automodule:: pyEQL.parameter
   :members:
