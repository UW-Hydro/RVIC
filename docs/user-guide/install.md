# Installing RVIC

## Dependencies
- [Python 2.7 or later including Python3](http://www.python.org)
- [NumPy](http://www.numpy.org)
- [SciPy](http://www.scipy.org)
- [netcdf4-python](https://code.google.com/p/netcdf4-python)
- [pandas](http://pandas.pydata.org)

If using `REMAP=True`:

- [CDO](https://code.zmaw.de/projects/cdo)
- [cdo.py](https://github.com/Try2Code/cdo-bindings)

## Installing using a package manager

RVIC is available via [PyPi](https://pypi.python.org/pypi/rvic):

```
pip install rvic
```

or Anaconda via the UW-Hydro channel:

```
conda install --channel https://conda.anaconda.org/UW-Hydro rvic
```

## Building RVIC

### Option 1:  Using Anaconda

The easiest way to install RVIC and its dependencies is to use the [Anaconda Python distribution](https://store.continuum.io/cshop/anaconda/).

To install Anaconda, follow these two simple steps (check to make sure the installer version is the most current)

1.  download and run the Anaconda installer:  [http://continuum.io/downloads](http://continuum.io/downloads)

2.  setup a virtual environment for RVIC

```shell
conda create -n rvic anaconda
source activate rvic
```

*Note:  you'll need to do the `source activate rvic` to activate the RVIC virtual environment in any new shells.*

Now, download the RVIC source code:

```shell
git clone git@github.com:UW-Hydro/RVIC.git
```

From the RVIC source code repository, RVIC can be installed using Python's `distutils`:

```shell
python setup.py install
```

This installs a top level script, `rvic`, into your bin/ directory and the `rvic` package into your Python path.

If you don't want to use the Anaconda installation I've shown above, you can build the package in your local python installation using:
```python
python setup.py develop
```

### Option 2a:  Using a local Python Install (With Write Permissions)

If you have write permissions to the location of your Python distribution, you can just run

```shell
python setup.py install
```

from the top level RVIC directory.  This will install RVIC into your `$PYTHONPATH`.

### Option 2b:  Using a local Python Install (Without Write Permissions)

If you do not have write permissions, you can install RVIC in your local `$PYTHONPATH` by following these steps:

Create a `lib/python` directory in your `$HOME` directory:

```shell
mkdir -p $HOME/lib/python/
```

Add this library path to your `$PYTHONPATH` in your `.bashrc`:

```shell
export PYTHONPATH=$HOME/lib/python:$PYTHONPATH
```

Run `setup.py`:

```shell
python setup.py install --home=$HOME
```

## Testing your install

From the following command line:

```shell
rvic -h
python -c 'import rvic'
```

If you don't get any errors, you should be ready to use RVIC.
