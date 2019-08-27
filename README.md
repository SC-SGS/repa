# Librepa ("Repa")

Repa is a library for *load-balanced*, regular grids.

## Dependencies

### Required

- MPI distribution (e.g. [OpenMPI](https://www.open-mpi.org/))
- [Boost](https://www.boost.org/) (mpi, serialization)
- [CMake](https://cmake.org/) and Make/Ninja/... for building
- C++14 compatible compiler

Note, that Boost::MPI must be compiled with the chosen MPI distribution.

### Optional

- [Boost](https://www.boost.org/) (unit_test) for tests
- [KDPart](https://github.com/hirschsn/kdpart) for kd-tree-based load-balancing method
- [CGAL](https://www.cgal.org/) for grid-based load-balancing method
- [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview) for graph-partitioning-based load-balancing method

Same note as above goes for KDPart and ParMETIS.

### Linux Distributions

MPI, Boost, ParMETIS and CGAL can be installed e.g. on Debian/Ubuntu using:
```sh
$ apt-get install openmpi-bin libboost-all-dev libparmetis-dev libcgal-dev cmake
```

Install dependencies:
```sh
$ mkdir dep && cd dep
$ sh -c 'git clone https://github.com/hirschsn/kdpart && cd kdpart && make'
$ DEP_DIR="$(pwd)"
```

## Build

```sh
$ git clone https://github.com/hirschsn/repa && cd repa
$ mkdir build && cd build
$ cmake .. -GNinja -DKDPART_DIR=$(DEP_DIR)/kdpart
$ ninja
$ ninja test
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
