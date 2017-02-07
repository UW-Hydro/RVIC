# RVIC Streamflow Routing Model

| VIC Links & Badges              |                                                                             |
|------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| VIC Documentation      | [![Documentation Status](https://readthedocs.org/projects/rvic/badge/?version=latest)](https://readthedocs.org/projects/rvic/?badge=latest)                                                              |
| Travis Build           | [![Build Status](https://travis-ci.org/UW-Hydro/RVIC.svg?branch=master)](https://travis-ci.org/UW-Hydro/RVIC)                                                                                            |
| Code Health            | [![Code Health](https://landscape.io/github/UW-Hydro/RVIC/master/landscape.svg?style=flat)](https://landscape.io/github/UW-Hydro/RVIC/master)                                                            |
| License                | [![GitHub license](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://raw.githubusercontent.com/UW-Hydro/RVIC/master/LICENSE.txt)                                                              |
| Current Release DOI    | [![DOI](https://zenodo.org/badge/11590212.svg)](https://zenodo.org/badge/latestdoi/11590212) |

The RVIC streamflow routing model is an adapted version of the model the model typically used as a post-processor with the Variable Infiltration Capacity (VIC) hydrology model. The routing model is a source-to-sink model that solves a linearized version of the Saint-Venant equations. This model, developed by Lohmann et al. (1996, 1998a, 1998b), has been used in many offline studies at a variety of spatial scales. Furthermore, the development of the impulse response functions (IRFs) is done as a preprocessing step, which considerably reduces the computation time in subsequent routing steps.

### Usage
See the [RVIC Documentation Site](http://rvic.readthedocs.org/en/latest/).

### License
RVIC is available under the GNU GPL v3.0 License.  See [LICENSE.txt](https://raw.githubusercontent.com/UW-Hydro/RVIC/master/LICENSE.txt) for more information.
