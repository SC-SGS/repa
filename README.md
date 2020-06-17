# Librepa ("Repa")

Librepa is a library for parallel, distributed-memory, load-balanced, regular,
three-dimensional grids.

It provides these under a common interface, that can be used to query the
structure of the grid, map positions to cells and ranks and repartition it.
Librepa does not handle numeric payload.

**For the documentation, see [DOC.md](DOC.md).**

## Dependencies

For CI purposes, etc. a working Docker container with all dependencies as well as its Dockerfile can be found on [Docker Hub](https://hub.docker.com/r/hirschsn/repa).

### Required

- MPI distribution (e.g. [OpenMPI](https://www.open-mpi.org/))
- [Boost](https://www.boost.org/) (mpi, serialization) v1.67, 1.68 or 1.72 or newer
- [CMake](https://cmake.org/) and Make/Ninja/... for building
- C++14 compatible compiler and standard library

### Optional

- [Boost](https://www.boost.org/) (unit_test) for tests
- [KDPart](https://github.com/hirschsn/kdpart) for kd-tree-based load-balancing method
- [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview) for graph-partitioning-based load-balancing method
- [lahnerml's p4est](https://github.com/lahnerml/p4est/tree/p4est-ESPResSo-integration) for SFC-based grid

## Install

[Spack](https://github.com/spack/spack) is the recommended way to install repa.
First, add the [spack repo spack-hirschsn](https://github.com/hirschsn/spack-hirschsn/) with the receipe for repa and its dependencies.
Then simply:

```sh
spack install librepa
```

## Development

For development purposes, you can also install repa's dependencies with spack (you, again, need the additional spack repo, see above):

```sh
# Install dependencies in a spack environment called 'repa-dev'
spack env create repa-dev
spack env activate repa-dev
spack install --only dependencies librepa
```

### Manual Build

If you need to build repa by hand:

```sh
git clone https://github.com/hirschsn/repa && cd repa
spack env activate repa-dev
mkdir build && cd build
cmake ..
make
make test # Optional
```

## License

Copyright 2017-2020 The Repa Authors

Repa is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Repa is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Repa.  If not, see <https://www.gnu.org/licenses/>.

## Contributors

Core:

- Steffen Hirschmann

Contributors:

- Adriaan Nieß
- Maximilian Wildbrett
- Malte Brunn
- Michael Lahnert
- Simon Hauser
- Benjamin Vier
