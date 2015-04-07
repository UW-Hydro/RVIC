# Welcome to the RVIC User Guide

The RVIC Streamflow Routing model is a simple source to sink routing model. The model represents each grid cell by a node in the channel network. Unit hydrographs are developed that described the time distribution of flow from each source grid cell to a corresponding sink grid cell. The development of the unit hydrographs is done as a pre-processing step (i.e. `parameters` or `convert`). The final step is the `convolution` of the unit hydrographs with fluxes from a land surface model, typically VIC.

## Using RVIC from the command line
From the command line, the `rvic` executable can be used giving one of the {[`parameters`](user-guide/parameters), [`convert`](user-guide/conversion), [`convolution`](user-guide/convolution)} subcommands and a configuration file:

```shell
rvic {parameters, convert, convolution} configuration_file.cfg
```

Other command line options can be found using `rvic -h`

## Using RVIC inside a Python interpreter

```python
from rvic.parameters import parameters
from rvic.convolution import convolution

# Run the parameter generation routine
parameters.parameters(config_file, np=1)

# Or run the convolution routine
convolution.convolution(config_file)
```

For more information on the `RVIC` api, see [this](user-guide/api) page.

## The RVIC Model Workflow

There are three main utilities in the RVIC Model:

1.  [Parameter Development](user-guide/parameters): Development of impulse response functions using input datasets such as a flow direction grid, outlet locations, etc.
2.  [Flux Convolution](user-guide/convolution): This utility generates streamflows by convolving the impulse response function developed in the previous setp with runoff fluxes from a land model (e.g. VIC).
4.  [Conversion from Older VIC Routing Model Versions](user-guide/conversion): A simple conversion utility to provide users with the ability to convert old routing model setups into RVIC parameters.
