# RVIC API

Almost all of RVIC is available via a public API.

```Python
import rvic
```

`rvic.parameters.parameters` and `rvic.convolution.convolution` both support either a path (string) to an INI style configuration file or a dictionary including configuration options. This feature is quite useful when generating ensembles of RVIC simulations.

Most of the internals of RVIC are stored int the `rvic.core` module:

```Python
from rvic.core import latlon2yx
```

More coming soon!  For now, take a look at the [RVIC github repository](https://github.com/UW-Hydro/RVIC).
