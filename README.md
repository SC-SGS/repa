# Librepa ("Repa")

Repa is a library for *load-balanced*, regular grids.

## Dependencies

### Required

- MPI distribution (e.g. [OpenMPI](https://www.open-mpi.org/))
- [Boost](https://www.boost.org/) (mpi, serialization) [1] [2]
- [CMake](https://cmake.org/) and Make/Ninja/... for building
- C++14 compatible compiler

[1] Note, that Boost::MPI must be compiled with the chosen MPI distribution. <br/>
[2] Currently, only Boost 1.67.0 and Boost 1.68.0 are supported. Previous and later versions (<1.67.0 and 1.69.0–1.72.0) contain Boost::MPI bugs that we cannot work around.

### Optional

- [Boost](https://www.boost.org/) (unit_test) for tests
- [KDPart](https://github.com/hirschsn/kdpart) for kd-tree-based load-balancing method
- [CGAL](https://www.cgal.org/) for grid-based load-balancing method
- [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview) for graph-partitioning-based load-balancing method
- [lahnerml's p4est](https://github.com/lahnerml/p4est/tree/p4est-ESPResSo-integration) for SFC-based grid

Same note as above goes for KDPart and ParMETIS.

### Linux Distributions

MPI, Boost, ParMETIS and CGAL can be installed e.g. on Debian/Ubuntu using:

```sh
apt-get install openmpi-bin libboost-all-dev libparmetis-dev libcgal-dev cmake
```

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

Copyright 2017-2019 The Repa Authors

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
