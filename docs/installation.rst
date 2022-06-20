.. _installation:


Installation
************

Dependencies
============

pyEQL requires Python 3.0 or greater. We highly recommend the Anaconda distribution of Python, which bundles many other
scientific computing packages and is easier to install (especially on Windows). You can download it at
https://www.continuum.io/downloads.

pyEQL also requires the following packages:

     * `pint <https://github.com/hgrecco/pint>`_
     * `scipy <http://scipy.org/>`_

If you use pip to install pyEQL (recommended), they should be installed automatically.

Automatically install via pip and PyPI
======================================

Once Python is installed, The `Python Package Index <https://pypi.python.org/pypi>`_ repository will allow installation
to be done easily from the command line as follows::

    pip install pyEQL

This should automatically pull in the required dependencies as well.

.. note:: You may have to run 'pip3' rather than 'pip' if you intend to use your system's default Python installation
    rather than Anaconda. For example, on many Linux and Mac systems Python 2.x and Python 3.x are installed side-by-side.
    You can tell if this is the case on your system by going to a command line and typing 'python' like so::

      $ python
      Python 2.7.12 (default, Jul  1 2016, 15:12:24)
      [GCC 5.4.0 20160609] on linux2
      Type "help", "copyright", "credits" or "license" for more information.
      >>>

    This means Python 2.x is installed. If you run 'pip install' it will point to the Python 2.7 installation, but pyEQL
    only works on Python 3. So, try this::

      $ python3
      Python 3.5.1+ (default, Mar 30 2016, 22:46:26)
      [GCC 5.3.1 20160330] on linux
      Type "help", "copyright", "credits" or "license" for more information.
      >>>

    To get to Python 3.x, you have to type 'python3'. In this case, you would run 'pip3 install'

Manually install via Git
========================
Simply navigate to a directory of your choice on your computer and clone the repository by executing the following terminal command::

    git clone https://github.com/rkingsbury/pyEQL

Then install by executing::

    pip install -e pyEQL

.. note:: You may have to run 'pip3' rather than 'pip'. See the note in the Automatic installation section.
