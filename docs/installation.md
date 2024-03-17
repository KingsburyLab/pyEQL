(installation)=

# Installing

## Use a conda environment

We highly recommend installing python in an isolated environment using [`conda`](https://docs.conda.io/en/latest/) (or its speedier, backward-compatible successor, [mamba](https://mamba.readthedocs.io/en/latest/)). In particular, we recommend the [miniforge](https://github.com/conda-forge/miniforge#miniforge3) or [mambaforge](https://github.com/conda-forge/miniforge#mambaforge) distributions of Python, which are lightweight distributions of conda that automatically activate the `conda-forge` channel for up-to-date scientific packages.

```{note}
If you are on a Windows machine, we recommend you install the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install) and set up your conda environments inside the WSL environment.
```

After installing `conda` / `mamba`, follow their instructions to create an environment. The steps should be similar to the following:

1. Open your terminal (or "Anaconda prompt" or "Miniforge prompt" on Windows)
2. Pick a name for your environment (note: you can create many environments if you want)
3. type `conda create -n <name-you-picked> python=3.10` (if you install miniforge) or `mamba create -n <name-you-picked> python=3.10` (if you installed mambaforge) and press enter
4. After the environment is installed, type `conda activate <name-you-picked>` / `mamba activate <name-you-picked>` and press enter

## pip install

Once Python is installed and your environment is activated you can install `pyEQL` from [PyPi](https://pypi.python.org/pypi) by typing the following command:

```
pip install pyEQL
```

This should automatically pull in the required [dependencies](#other-dependencies) as well.

```{important}
If you are NOT using a `conda` environment, may have to run 'pip3' rather than 'pip'. This will be the case if Python 2.x and Python 3.x are installed side-by-side on your system.
You can tell if this is the case by typing the following command:

```

```
$ python --version
Python 2.7.12
```

This means Python 2.x is installed. If you run 'pip install' it will point to the Python 2.7 installation, but pyEQL
only works on Python 3. So, try this:

```
$ python3 --version
Python 3.9.7
```

To get to Python 3.x, you have to type 'python3'. In this case, you would run 'pip3 install'

```{warning}
If you are using a Mac with an Apple M1, M2, etc. chip (i.e., Arm64 architecture), some features of `pyEQL` will be
unavailable. Specifically, anything which depends on PHREEQC (e.g., the
`equilibrate` method in the native engine and the entire
 [`phreeqc` engine](engines.md)) will not work. This is because `phreeqpython` is currently
 not available for this platform. All other functions of `pyEQL` should work as expected.

Feel free to post your experiences or proposed solutions at https://github.com/KingsburyLab/pyEQL/issues/109

NOTE: Some users have reported being able to use `phreeqpython` functions by installing
an x86 version of `conda`/`miniconda`. See the issue report for more details.
```

## Other dependencies

pyEQL also requires the following packages:

- [pint](https://github.com/hgrecco/pint) - for automated unit conversion
- [pymatgen](https://github.com/materialsproject/pymatgen/) - used to interpret chemical formulas
- [iapws](https://github.com/jjgomera/iapws/) - used to calculate the properties of water
- [monty](https://github.com/materialsvirtuallab/monty) - used for saving and loading `Solution` objects to files
- [maggma](https://materialsproject.github.io/maggma/) - used by the internal property database
- [scipy](http://scipy.org/)
- [numpy](http://numpy.org/)

If you use pip to install pyEQL (recommended), they should be installed automatically.

## Installing the development branch

If you want to use the bleeding edge version before it is released to PyPi instead of the latest stable release, you can substitute the following for the above 'pip install' command:

```
pip install git+https://github.com/KingsburyLab/pyEQL.git@main
```

## Manually install via Git

Simply navigate to a directory of your choice on your computer and clone the repository by executing the following terminal command:

```
git clone https://github.com/KingsburyLab/pyEQL
```

Then install by executing:

```
pip install -e pyEQL
```

```{note}
You may have to run 'pip3' rather than 'pip'. See the note in the [pip install](#pip-install) section.
```
