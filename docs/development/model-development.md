## Model Development

Most of the development of RVIC has been done in the UW Hydro | Computational Hydrology Group at the University of Washington.  There are others interested in expanding the RVIC model to include alternate forms of unit hydrograph development, parrallization of the convolution routines, and coupling to VIC and other land surface models.

For a list of planned features, or to suggest a new feature, go to the [RVIC github issues page](https://github.com/UW-Hydro/RVIC/issues).


## Contribuiting

RVIC is a [Open Source](development/open-source) model.  As such, we want encourage other researchers to contribute to the source code.  Please fork the github repository and add new features.

## Testing
A three tiered test platform has been developed for the RVIC model.

These tests are executed and controlled by `RVIC/tests/run_tests.py`.  This script controls which tests are executed and where they read/write to.

1.  **Unit Testing** - Using the [`pytest`](http://pytest.org/latest/) module, small pieces of code are tested using a standard [unit testing approach](http://en.wikipedia.org/wiki/Unit_testing).  The existing unit testing coverage is lower than the ideal level so further development in this area would be useful.

2.  **System Tests** - System tests are basically short examples that test the function of a large piece of the RVIC model.  A set of tests has been developed, based on 3 real life applications of the RVIC model.  The data and configuration files for these tests are located [here](ftp://ftp.hydro.washington.edu/pub/UW-Hydro/RVIC/rvic_tests.tar.gz).  To use these tests, follow these simple steps:
    - `cd RVIC/tests`
    - `curl -O ftp://ftp.hydro.washington.edu/pub/jhamman/RVIC/rvic_tests.tar.gz`
    - `tar zxvf rvic_tests.tar.gz`
    - `./run_tests.py {all, examples}`

    To run a different set of system tests, simply make a copy of the test configuration file `RVIC/tests/examples/examples.cfg` and follow the instructions in the header of file on how to setup a new test.  Then run the test using `run_tests.py {all, examples} --examples=path_to_new_test_config_file`

3.  **Continuous Integration** - Continuous integration is done using the [`travis-ci.org`](https://travis-ci.org/UW-Hydro/RVIC) platform.  After each commit or merge to the RVIC repository, the model is rebuilt using a fresh install of the Anaconda Python package and the *unit tests* are run.  The status of these tests are shown on the RVIC [README.md](https://github.com/UW-Hydro/RVIC/blob/master/README.md) page.  The current status of the master branch build is: [![Build Status](https://travis-ci.org/UW-Hydro/RVIC.svg?branch=master)](https://travis-ci.org/UW-Hydro/RVIC)
