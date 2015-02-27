# Installing RVIC

## Dependencies
- [Python 2.7](http://www.python.org/)
- [Numpy](http://www.numpy.org)
- [Scipy](http://www.scipy.org/)
- [netcdf4-python](https://code.google.com/p/netcdf4-python/)

If using `REMAP=True`:

- [cdo](https://code.zmaw.de/projects/cdo)
- [cdo.py](https://github.com/Try2Code/cdo-bindings)

## Building RVIC
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
git clone https://github.com/jhamman/RVIC.git
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
