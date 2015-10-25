# What's New

v1.1.0 (25 October 2015)

This release contains a number of bug and compatibility fixes.

### Enhancements

- Compatibility with Python 3.4 and 3.5 ([GH38](https://github.com/UW-Hydro/RVIC/pull/38)).
- Simplified multiprocessing in `parameters` which should improve stability ([GH60](https://github.com/UW-Hydro/RVIC/pull/60)).

### Bug Fixes

- Fixed bug that caused `REMAP=False` to fail ([GH60](https://github.com/UW-Hydro/RVIC/pull/60)).
- Improvements were made to the `SEARCH_FOR_CHANNEL` option in `parameters`.
- Fixed bug where last timestep of history output was not populated ([GH71](https://github.com/UW-Hydro/RVIC/pull/71)).
