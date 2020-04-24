# Librepa ("Repa")

Repa is a library for *load-balanced*, regular grids.

## Dependencies

For CI purposes, etc. a working Docker container with all dependencies as well as its Dockerfile can be found on [Docker Hub](https://hub.docker.com/r/hirschsn/repa).

### Required

- MPI distribution (e.g. [OpenMPI](https://www.open-mpi.org/))
- [Boost](https://www.boost.org/) (mpi, serialization) [1] [2]
- [CMake](https://cmake.org/) and Make/Ninja/... for building
- C++14 compatible compiler and standard library [3]

[1] Note, that Boost::MPI must be compiled with the chosen MPI distribution. <br/>
[2] Currently, only Boost 1.67.0 and Boost 1.68.0 are supported. Previous and later versions (<1.67.0 and 1.69.0–1.72.0) contain Boost::MPI bugs that we cannot work around.<br/>
[3] If you have to use an older standard library, have a look at the "libstdc++-4.8.5" branch.

### Optional

- [Boost](https://www.boost.org/) (unit_test) for tests
- [KDPart](https://github.com/hirschsn/kdpart) for kd-tree-based load-balancing method
- [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview) for graph-partitioning-based load-balancing method
- [lahnerml's p4est](https://github.com/lahnerml/p4est/tree/p4est-ESPResSo-integration) for SFC-based grid

Same note as above goes for KDPart and ParMETIS.

### Spack

Using [spack](https://github.com/spack/spack), you can install repa easily.
First, add the [spack repo spack-hirschsn](https://github.com/hirschsn/spack-hirschsn/) with the receipe for repa and its dependencies.
Then simply:

```sh
spack install librepa
```

For development purposes, you can also install repa's dependencies (you, again, need the additional spack repo) and load them as a spack environment:

```sh
git clone https://github.com/hirschsn/repa
cd repa
spack install
# Wait...
spack env activate .
```

### Linux Distributions

MPI, Boost, ParMETIS and CGAL can be installed e.g. on Debian/Ubuntu using:

```sh
apt-get install openmpi-bin libboost-all-dev libparmetis-dev cmake
```

Note for manual ParMETIS installations: ParMETIS's `make install` does *not* install metis.h from metis/include/. Copy this file manually to the appropriate prefix/include directory.

Install dependencies:

```sh
mkdir -p dep/src && cd dep/src
sh -c 'git clone https://github.com/hirschsn/kdpart \
          && cd kdpart \
          && make \
          && make install PREFIX="$(pwd)/../.."'
sh -c 'git clone --recursive https://github.com/lahnerml/p4est --branch p4est-ESPResSo-integration \
          && cd p4est \
          && ./bootstrap \
          && ./configure --prefix="$(pwd)/../.." --enable-mpi \
          && make \
          && make install'
DEP_DIR="$(pwd)/.."
```

## Build

```sh
git clone https://github.com/hirschsn/repa && cd repa
mkdir build && cd build
cmake .. -DKDPART_DIR="${DEP_DIR}" -DP4EST_DIR="${DEP_DIR}"
make
# If necessary, don't forget: export LD_LIBRARY_PATH="${DEP_DIR}/lib:$LD_LIBRARY_PATH"
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
