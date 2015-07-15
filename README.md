# RVIC Streamflow Routing Model

[![Documentation Status](https://readthedocs.org/projects/rvic/badge/?version=latest)](https://readthedocs.org/projects/rvic/?badge=latest) [![Build Status](https://travis-ci.org/UW-Hydro/RVIC.svg?branch=master)](https://travis-ci.org/UW-Hydro/RVIC) [![Code Health](https://landscape.io/github/UW-Hydro/RVIC/master/landscape.svg?style=flat)](https://landscape.io/github/UW-Hydro/RVIC/master) [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.17740.svg)](http://dx.doi.org/10.5281/zenodo.17740)



The RVIC streamflow routing model is an adapted version of the model the model typically used as a post-processor with the Variable Infiltration Capacity (VIC) hydrology model. The routing model is a source-to-sink model that solves a linearized version of the Saint-Venant equations. This model, developed by Lohmann et al. (1996, 1998a, 1998b), has been used in many offline studies at a variety of spatial scales. Furthermore, the development of the impulse response functions (IRFs) is done as a preprocessing step, which considerably reduces the computation time in subsequent routing steps.

### Usage
See the [RVIC Documentation Site](http://rvic.readthedocs.org/en/latest/).

### License
RVIC is available under the GNU GPL v3.0 License.  See LICENSE.txt for more information.
