# RVIC Streamflow Routing Model

**This is just a quick example.  A second change.**

[![Build Status](https://travis-ci.org/UW-Hydro/RVIC.svg?branch=master)](https://travis-ci.org/UW-Hydro/RVIC) [![Documentation Status](https://readthedocs.org/projects/rvic/badge/?version=latest)](https://readthedocs.org/projects/rvic/?badge=latest)

The RVIC streamflow routing model is an adapted version of the model the model typically used as a post-processor with the Variable Infiltration Capacity (VIC) hydrology model. The routing model is a source-to-sink model that solves a linearized version of the Saint-Venant equations. This model, developed by Lohmann et al. (1996, 1998a, 1998b), has been used in many offline studies at a variety of spatial scales. Furthermore, the development of the impulse response functions (IRFs) is done as a preprocessing step, which considerably reduces the computation time in subsequent routing steps.

### Usage
See the [RVIC Documentation Site](http://rvic.readthedocs.org/en/latest/).

### License
RVIC is available under the GNU GPL v3.0 License.  See LICENSE.txt for more information.
