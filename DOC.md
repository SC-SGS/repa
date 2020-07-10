# Introduction

Librepa is a library for parallel, distributed-memory, load-balanced, regular,
three-dimensional grids.

It provides these under a common interface, that can be used to query the
structure of the grid, map positions to cells and ranks and repartition it.
Librepa does not handle numeric payload.

# Examples

The following examples are stripped of unnecessary details. Almost all examples
are available fully compilable and runnable in the `examples/` folder, though
not built by default. Run `make examples` to build them.

## Generation of a distributed grid

Creating a fully-periodic grid over the unit cube
![\Omega = \[0, 1)^3](https://render.githubusercontent.com/render/math?math=%5COmega%20%3D%20%5B0%2C%201%29%5E3)
with a mesh width of *at least*
![h = 10^{-1}](https://render.githubusercontent.com/render/math?math=h%20%3D%2010%5E%7B-1%7D):

```c++
// example01.cpp

#include <repa/repa.hpp>
// boost::mpi::communicator comm;
// Create a 3d Cartesian grid distributed across "comm"
auto grid = repa::make_pargrid(repa::GridType::CART, comm, {1., 1., 1.}, 1e-1);
// "grid" is of type std::unique_ptr<repa::ParallelLCGrid>
```

- The lower bound of the grid is *always* implicitly the origin.
- The mesh width *can* be chosen larger by the implementation.
- The call to `repa::make_pargrid` is collective, i.e. has to be done by all
  processes in the communicator passed to it.
- Internally, the state is *automatically* distributed over all MPI processes
- Non-periodic grids are currently not supported.
- To use a different partitioning strategy than Cartesian, just change the
  `repa::GridTypes::CART` to some other supported value, e.g.
  `repa::GridTypes::GRAPH`.

## Querying the grid size and the mesh width

```c++
// example02.cpp

#include <repa/repa.hpp>
auto grid = repa::make_pargrid(repa::GridType::CART, comm, {1., 1., 1.}, 1e-1);
repa::Vec3i grid_size = grid->grid_size();
repa::Vec3d cell_size = grid->cell_size();
std::cout << "Grid size: " << grid_size << std::endl;
std::cout << "Cell size: " << cell_size << std::endl;
```

Output:

```sh
Grid size: {10, 10, 10}
Cell size: {0.1, 0.1, 0.1}
```

- The "grid size" is the global number of cells in each dimension. The ouput is
  the same on all processes.
- For every dimension, the actual mesh width ("cell size")
  ![s](https://render.githubusercontent.com/render/math?math=s)
  is guaranteed to be:
  ![h \leq s < 2\,h](https://render.githubusercontent.com/render/math?math=h%20%5Cleq%20s%20%3C%202%5C%2Ch)
- The grid size is *always* integer. Thus the variable cell size.
- Let ![b](https://render.githubusercontent.com/render/math?math=b) be the box
  size a specific dimension. If
  ![b = k\,h](https://render.githubusercontent.com/render/math?math=b%20%3D%20k%5C%2Ch)
  for some integer
  ![k](https://render.githubusercontent.com/render/math?math=k), then
  ![s = h](https://render.githubusercontent.com/render/math?math=s%20%3D%20h)
  (as in the example).

## Getting the number of local and ghost cells

```c++
// example03.cpp

#include <repa/repa.hpp>
auto grid = repa::make_pargrid(repa::GridType::CART, comm, {1., 1., 1.}, 1e-1);
std::cout << "Number of local cells: "
          << grid->local_cells().size()
          << std::endl;
std::cout << "Number of ghost cells: "
          << grid->ghost_cells().size()
          << std::endl;
```

The output will depend on the number of MPI processes the program runs on.
It is, however, guaranteed, that:

- The sum of the number of local cells across all processes is always:
  `grid_size[0] * grid_size[1] * grid_size[2]`
- There are only necessary ghost cells (we call this "minimal ghost layer"
  instead of a "full halo")

```sh
$ mpiexec -n 1 ./example03
Number of local cells: 1000
Number of ghost cells: 0

$ mpiexec -n 2 ./example03 # Ordering of output unspecified
Number of local cells: 500
Number of ghost cells: 200
Number of local cells: 500
Number of ghost cells: 200
```

## Querying the rank and cell of a position

```c++
// example04.cpp

#include <repa/repa.hpp>
auto grid = repa::make_pargrid(repa::GridType::CART, comm, {1., 1., 1.}, 1e-1);
const auto pos = repa::Vec3d{.75, .75, .75};

const int rank = grid->position_to_rank(pos);
if (rank == comm.rank()) {
    // On ranks that do not own "pos" the following call will throw
    const auto cellidx = grid->position_to_cell_index(pos);
    std::cout << pos
              << " belongs to rank " << rank
              // "cellidx" is of type repa::local_cell_index_type and can be
              // (implicitly) converted to an integer representing the cell
              // index.
              << " into cell: " << static_cast<int>(cellidx)
              << std::endl;
}
```

Output:

```sh
$ mpiexec -n 1 ./example04
{0.75, 0.75, 0.75} belongs to rank 0 into cell: 777

$ mpiexec -n 2 ./example04
{0.75, 0.75, 0.75} belongs to rank 1 into cell: 277
```

Note:

- `position_to_X` functions may only be called with valid positions in
  ![\Omega](https://render.githubusercontent.com/render/math?math=%5COmega).
- `position_to_cell_index(pos)` will throw if called on processes with rank
  `r` where `r != position_to_rank(pos)`
- Before the first call to `repartition`, `position_to_rank` is guaranteed to
  resolve the whole domain on *all* processes.
- After the first call to `repartition`, `position_to_rank` might only resolve
  a subset of the domain. This subset is guaranteed to include at least the own
  subdomain and the ghost layers.
- In any case, `position_to_rank` throws a `std::domain_error` if the position
  cannot be resolved.

## Cell ordering

Cell indices are guaranteed to start from 0 and be continuously numbered:

- Local cell indices
  ![l](https://render.githubusercontent.com/render/math?math=l):
  ![l \in \{0, \ldots, L - 1\}](https://render.githubusercontent.com/render/math?math=l%20%5Cin%20%5C%7B0%2C%20%5Cldots%2C%20L%20-%201%5C%7D),
  where `L = grid->local_cells().size()`.
- Ghost cell indices
  ![g](https://render.githubusercontent.com/render/math?math=g):
  ![g \in \{0, \ldots, G - 1\}](https://render.githubusercontent.com/render/math?math=g%20%5Cin%20%5C%7B0%2C%20%5Cldots%2C%20G%20-%201%5C%7D),
  where `G = grid->ghost_cells().size()`.

This applies to each process individually.
To conveniently iterate over these ranges of cell, `grid->local_cells()` and
`grid->ghost_cells()` return range objects. Elements of these ranges can be
converted to integers to get a raw index.

```c++
// example05.cpp

#include <repa/repa.hpp>
auto grid = repa::make_pargrid(repa::GridType::CART, comm, {1., 1., 1.}, 1e-1);

// Use the range objects in conjunction with range-based for loops
for (const auto local_idx : grid->local_cells()) {
    // Index types cast implicitly to int
    std::cout << "Rank " << comm.rank() << " local cell index: " << local_idx;
    // On debug builds the function "global_hash" will return the associated
    // global cell index to a local index:
    std::cout << " global index: " << grid->global_hash(local_idx) << std::endl;
}
```

Output:

```sh
$ mpiexec -n 1 ./example05
Rank 0 local cell index: 0 global index: 0
Rank 0 local cell index: 1 global index: 1
Rank 0 local cell index: 2 global index: 2
Rank 0 local cell index: 3 global index: 3
[...]

$ mpiexec -n 2 ./example05 # Output of the processes might interleave
Rank 0 local cell index: 0 global index: 0
Rank 0 local cell index: 1 global index: 1
Rank 0 local cell index: 2 global index: 2
Rank 0 local cell index: 3 global index: 3
[...]
Rank 1 local cell index: 0 global index: 500
Rank 1 local cell index: 1 global index: 501
Rank 1 local cell index: 2 global index: 502
Rank 1 local cell index: 3 global index: 503
[...]
```

Note:

- Local cell indices are unique only on the same rank. So are ghost indices.
  Ghost indices are, however, associated to a different process via the ghost
  exchange volume (see below).
- Global cell indices are unique among all processes, they should, however, not
  be relied on. The function `global_hash()` might always return `0` on
  non-debug builds. (This is mostly the case if keeping track of the global
  indices causes overhead.)

## Finding neighbor cells of a local cell

```c++
// example06.cpp

#include <repa/repa.hpp>
auto grid = repa::make_pargrid(repa::GridType::CART, comm, {1., 1., 1.}, 1e-1);
const auto c = repa::local_cell_index_type{0};
const int which_neighbor = 0; // 0 - 26

const auto neigh = grid->cell_neighbor_index(c, which_neighbor);
std::cout << "Neighbor " << which_neighbor << " of cell " << c << " is ";
if (neigh.is<repa::local_cell_index_type>())
    std::cout << " local cell " << neigh.as<repa::local_cell_index_type>();
else
    std::cout << " ghost cell " << neigh.as<repa::ghost_cell_index_type>();
std::cout << std::endl;
```

Output:

```sh
$ mpiexec -n 1 ./example06
Neighbor 0 of cell 0 is local cell 0
$ mpiexec -n 2 ./example06
Neighbor 0 of cell 0 is local cell 0
Neighbor 0 of cell 0 is local cell 0
```

- The neighbor index (`which_neighbor`) is encoded as follows:
  - `0` is the cell itself, i.e. `c == grid->cell_neighbor_index(c, 0)` is
    always true for valid `c`.
  - `1-13` are neighbors such that traversing all local cells and their
    neighbors `1-13` guarantees that each pair of neighboring cells has been
    traversed exactly once (see the following example).
  - `14-26` are the remaining neighbors.
- The ouput of `grid->cell_neighbor_index` is of type
  `local_or_ghost_cell_index_type`, which is a
  `repa::util::simple_variant<local_cell_index_type, ghost_cell_index_type>`.
  See its documentation for the provided methods. The most relevant are:
  - `neigh.is<T>() -> bool` where `T` is either `local_cell_index_type` or
    `ghost_cell_index_type`.
  - `neigh.as<T>() -> T`. If `neigh.is<T>() != true`, triggers an assertion on
    debug and an abort on non-debug builds.
  - `neigh.fmap(F)` and `neigh.visit(F1, F2)` for advanced uses.

## Iterating over all pairs of neighboring cells exactly once

```c++
// example07.cpp

#include <repa/repa.hpp>
auto grid = repa::make_pargrid(repa::GridType::CART, comm, {1., 1., 1.}, 1e-1);
// Combination of example05 and example06
for (const auto &c : grid->local_cells()) {
    for (int which_neigh = 0; which_neigh <= 13; ++which_neigh) {
        const auto d = grid->cell_neighbor_index(c, which_neigh);

        // Gets called exactly once for all unordered pairs {c, d} (which
        // includes the pair {c, c}). We use global cells here for
        // demonstration purposes.
        f(grid->global_hash(c), grid->global_hash(d));
    }
}

// Let's see how often we come by global cell index "c_want"
const int c_want = 120;
void f(repa::global_cell_index_type c, repa::global_cell_index_type d)
{
    static int ncalls = 0;
    if (c == c_want || d == c_want) {
        ncalls++;
        std::cout << ncalls << " ";
        std::cout << c << ", " << d << std::endl;
    }
}
```

Output:

```sh
1 10, 120
2 19, 120
3 20, 120
4 29, 120
5 39, 120
6 110, 120
7 119, 120
8 120, 120
9 120, 220
10 120, 30
11 120, 130
12 120, 230
13 120, 11
14 120, 111
15 120, 211
16 120, 21
17 120, 121
18 120, 221
19 120, 31
20 120, 131
21 120, 231
22 129, 120
23 139, 120
24 210, 120
25 219, 120
26 229, 120
27 239, 120
```

Note:

- Exactly the expected `27 = 3*3*3` pairs.
- `c_want = 120` is the first index in the pairs in exactly `14 = 13+1` cases.
- Use this pattern to implement "half-shell" neighborhood traversal in MD codes.

## Ghost Exchange Volume

Since librepa doesn't handle numeric payload for you, you have to perform
the ghost exchange yourself. Librepa, however, tells you which data needs
to be exchanged via a format called `GhostExchangeDesc` defined in
`pargrid.hpp`:

```c++
struct GhostExchangeDesc {
    rank_type dest; // Destination rank
    std::vector<ghost_cell_index_type>
        recv; // Ghost cell indices which are to be received
    std::vector<local_cell_index_type>
        send; // Local cell indices which are to be sent

    //...
};
```

The function `ParallelLCGrid::get_boundary_info()` returns a range of
`GhostExchangeDesc`.

```c++
// example08.cpp

#include <repa/repa.hpp>
auto grid = repa::make_pargrid(repa::GridType::CART, comm, {1., 1., 1.}, 1e-1);

const auto gexds = grid->get_boundary_info();

    for (const auto &gexd : gexds) {
        // Print the number of cells to be send and received for each ghost
        // exchange
        std::cout << "Communication " << comm.rank() << " -> " << gexd.dest
                  << " send: " << gexd.send.size()
                  << " recv: " << gexd.recv.size() << std::endl;
    }
```

Output:

```sh
$ mpiexec -n 1 ./example08
[No output]

$ mpiexec -n 2 ./example08
Communication 1 -> 0 send: 200 recv: 200
Communication 0 -> 1 send: 200 recv: 200
```

- Ghost layer is restricted to *one* layer of cells only.
- Each neighboring process is included at most *once* in the ghost exchange
  volume and there is *no* self-communication.
- The exchange descriptors on either side have the same ordering, e.g.:
  Let `gexd_p` on process `p` correspond to `gexd_q` on process `q` (that means
  `gexd_p.dest == q` and `gexd_q.dest == p`), then `gexd_p.send[i]` corresponds
  to `gexd_q.recv[i]`, meaning these local/ghost indices will represent the
  same global cell.
- On each process holds: the number of ghost cells is the sum of `recv.size()`
  across all `GhostExchangeDesc`.

## Combined Example 1: Sequential Short-Range Molecular Dynamics Example

This example contains a couple of stubs and partially implemented functions
to show how librepa can be used in an MD application.

```c++
// This example is not available under the examples/ directory.

#include <repa/repa.hpp>
// Linked-cell grid with maximum short-range interaction range of 2.5
// over a simulation box of [0, 80)^3
auto grid = repa::make_pargrid(repa::GridType::CART, {80., 80., 80.}, 2.5);

struct Particle {
    //...
};
// For the sake of this example, assume this function
repa::Vec3d position(Particle &);

using Particles = std::vector<Particle>;
// Storage of particles - a "cell" is just a vector of particles
std::vector<Particles> linked_cells;

// Sort Particles into their corresponding cell
void insert_particles(Particles &&ps)
{
    linked_cells.clear();
    // Allocate the particle storage (ordering according to "grid")
    linked_cells.resize(grid->local_cells().size());

    for (auto &p : ps) {
        linked_cells.at(
            grid->position_to_cell_index(position(p))
        ).push_back(std::move(p));
    }
}

// Half-shell force calculation of all particles
void linked_cell_algorithm()
{
    set_all_forces_zero();
    fold_particles_back_into_primary_simulation_box();
    insert_particles(extract_particles(linked_cells()));

    // Traverse local cells
    for (const auto &c : grid->local_cells()) {
        // Calculates interactions of particles in "c"
        calculate_forces(linked_cells[c]);

        // Iterate over "half-shell" neighborhood of "c"
        for (int which_neigh = 1; which_neigh <= 13; ++which_neigh) {
            const auto d = grid->cell_neighbor_index(c, which_neigh);
            // Calculate interactions of particles in "c" with particles in "d"
            calculate_forces(linked_cells[c], linked_cells[d]);
        }
    }
}

// MD stuff
void calculate_forces(Particles &a, Particles &b)
{
    for (Particle &p : a)
        for (Particle &q : b)
            add_forces(p, q);
}

void calculate_forces(Particles &a)
{
    for (auto it_p = a.begin(); it_p != a.end(); ++it_p)
        for (auto it_q = std::next(it_p); it_q != a.end(); ++it_q)
            add_forces(*it_p, *it_q);
}

// Adds force of interaction between "p" and "q" to both, "p" and "q".
void add_forces(Particle &p, Particle &q);
```

Note: This example neither contains MD details like `add_forces` nor
driver code, that actually calls `linked_cell_algorithm` and integrates the
particles, etc.

## Combined Example 2: Parallel Short-Range Molecular Dynamics Example

This example contains a couple of stubs and partially implemented functions
to show how librepa can be used in a *parallel* MD application.
The following code extends the one from above.

```c++
// This example is not available under the examples/ directory.

#include <boost/mpi.hpp>

std::vector<Particles> ghost_cells;

// grid->get_boundary_info() as well as grid->neighbors() is guaranteed to list
// a only a single communication per neighbor rank, so the communication
// tag is irrelevant.
const int tag = 0;

// Communicator (and boost::mpi::environment) defined somewhere
boost::mpi::communicator comm;

// The following assumes "Particle" is serializable

void ghost_exchange()
{
    // Ghost communication sends around vectors of cells to neighbors
    using CommVol = decltype(linked_cells);
    std::map<int, CommVol> send_data; // Using map for ease of demonstration
    std::map<int, CommVol> recv_data;

    // Iterate over the "boundary info"
    for (const auto &ghost_exchange : grid->get_boundary_info()) {
        const auto dest = ghost_exchange.dest;
        // Post receives
        comm.irecv(dest, tag, recv_data[dest]);

        // Collect boundary cells and send them
        CommVol &send_buf = send_data[dest];
        for (const auto &c : ghost_exchange.send)
            send_buf.push_back(linked_cells[c]);
        comm.isend(dest, tag, send_buf);
    }

    // To not clutter this example, the boost::mpi::request objects
    // returned by comm.irecv and comm.isend are not collected.
    // If we did, they would need to be waited for here.
    boost::mpi::wait_all(...);

    // Put the received data into ghost cells.
    ghost_cells.clear();
    ghost_cells.resize(grid->ghost_cells().size());
    for (const auto &ghost_exchange : grid->get_boundary_info()) {
        CommVol &recv_buf = send_data[ghost_exchange.dest];
        for (const auto &c : ghost_exchange.recv)
            ghost_cells[c] = std::move(recv_buf[c]);
    }
}

void exchange_particles()
{
    std::map<int, Particles> send_vol;
    std::map<int, Particles> recv_vol;

    // Sort out outliers (particles that don't belong on this process anymore)
    for (auto &c : linked_cells) {
        for (auto it = c.begin(); it != c.end(); ++it) {
            if (const auto rank = grid->position_to_rank(position(p));
                rank != comm.rank()) {
                    // Swapping particle to the end, so it can be moved from.
                    std::swap(p, c.back());
                    send_vol[rank].push_back(std::move(c.back()));
                    c.pop_back();
                    --it; // Not advancing iterator, because we swapped a new
                          // particle here.
            }
        }
    }

    for (const auto dest : grid->neighbors()) {
        comm.irecv(dest, 0, recv_vol[dest]);
        comm.isend(dest, 0, send_vol[dest]);
    }

    // See above
    boost::mpi::wait_all(...);

    // Insert particles into "linked_cells" using the
    // basic linked-cell insertion function from above
    for (const auto dest : grid->neighbors()) {
        insert_particles(std::move(recv_vol[dest]));
    }
}

// Minimal change required in:
void linked_cell_algorithm()
{
    set_all_forces_zero();

    for (const auto &c : grid->local_cells()) {
        calculate_forces(linked_cells[c]);

        for (int which_neigh = 1; which_neigh <= 13; ++which_neigh) {
            const auto d = grid->cell_neighbor_index(d, which_neigh);

            // Note the change here:
            using local_type = repa::local_cell_index_type; // Make the following
            using ghost_type = repa::ghost_cell_index_type; // lines shorter
            Particles &neigh_cell = d.is<local_type>()
                                    ? linked_cells[d.as<local_type>()]
                                    : ghost_cells[d.as<ghost_type>()]
            calculate_forces(linked_cells[c], neigh_cell);
        }
    }
}
```

## Initial Partitioning

```c++
#include <repa/repa.hpp>

// C++20 designated initializers
auto grid = repa::make_pargrid(repa::GridType::CART, comm, {1., 1., 1.}, 1e-3, repa::ExtraParams{ .init_part = "cart3d" });

// pre-C++20
repa::ExtraParams ep;
ep.init_part = "cart3d"; // or "linear", etc.
auto grid = repa::make_pargrid(repa::GridType::CART, comm, {1., 1., 1.}, 1e-3, ep);
```

TODO...

## Repartition

TODO...

## Commands & Variants

TODO...

# Notes

See the individual notes under each example.

# Rationale

This section list miscellaneous reasons or choices made in the library.
Though they are ordered for easy referencing, the order is arbitrary.

1. `local_cells()` and `ghost_cells()` return range object since they are
   more convenient to use with range-based for loops. In previous versions,
   lines like `for (int i = 0; i < grid->n_local_cells();  ++i)` were
   omnipresent. `for (const auto cell : grid->local_cells())` is not only
   easier to comprehend in my opinion but also leaves the ordering/domain of
   the local cells to the grid implementation.
1. Historically, there has been a function called `position_to_neighidx`. It
   returned an index `i` to `grid->neighbors()` to easily manage data sent to
   neighbors (allowing e.g. the examples to use a `std::vector` instead a
   `std::map`). This function was removed. If you need it use:

    ```c++
    std::ptrdiff_t position_to_neighidx(ParallelLCGrid *grid, repa::Vec3i position)
    {
        const auto rank = grid->position_to_rank(position);
        const auto it = std::find(grid->neighbors().begin(), grid->neighbors().end(), rank);

        if (it != grid->neighbors.end())
            return std::distance(grid->neighbors().begin(), it);
        else
            throw no_such_neighbor();
    }
    ```

1. Ghost-Layers are restricted to one cell width, because this code was and is
   intended to be used for short-range MD.
1. Regular grids only, because this code was and is intended to be used for
   short-range MD, i.e. linked-cell grids, which are mostly chosen regular.
   If you want software that supports adaptive grids, librepa is (currently)
   of no use for you.

# Reference

The exposed functions are:

```c++
namespace repa {

// Returns the GridType associated with a descriptive string for the grid type.
GridType parse_grid_type(const std::string &desc);

// Inverse of 'parse_grid_type'.
std::string grid_type_to_string(GridType gt);

// Returns true if support for a certain grid type is compiled in. 
bool has_grid_type(GridType gt);

// Returns a set of all supported grid types.
std::set< GridType > supported_grid_types();

// Returns a reference to a variant setter object if the grid supports it.
// TOOD:
VariantSetter & variants (grids::ParallelLCGrid *grid);

// Grid factory method.
std::unique_ptr< grids::ParallelLCGrid >
make_pargrid(GridType gt,
             const boost::mpi::communicator &comm,
             Vec3d box_size,
             double min_cell_size,
             ExtraParams ep = ExtraParams{});

// VariantSetter interface, see below
struct VariantSetter {
    virtual std::set<std::string> get_supported_variants() const;
    virtual void set_variant(const std::string &);
};

namespace grids {
// See below
struct ParallelLCGrid;
} // namespace grids
} // namespace repa
```

## ParallelLCGrid

The `ParallelLCGrid` interface is defined as follows:

```c++
struct ParallelLCGrid {
    /** Returns the range of local cells.
     */
    cell_range<local_cell_index_type> local_cells() const;

    /** Returns the range of ghost cells.
     */
    cell_range<ghost_cell_index_type> ghost_cells() const;

    /** Returns a span/range of ranks of all neighbor processes.
     */
    virtual util::const_span<rank_type> neighbor_ranks() const;

    /** Returns the cell sizes of Linked Cell grid.
     */
    virtual Vec3d cell_size() const;

    /** Returns the number of grid cells in total in each direction.
     */
    virtual Vec3i grid_size() const;

    /** Returns the index of a cell neighboring a given cell (by index).
     *
     * The neighbor can either be a local cell or a ghost cell.
     *
     * @throws std::domain_error if cellidx is not a valid local cell.
     *
     * Neighbor 0 is the cells itself, i.e. "cell_neighbor_index(c, 0) == c"
     * Neighbors 1-13: Half shell neighborhood
     * Neighbors 14-26: Rest of full shell neighborhood
     *
     * @param cellidx Base cell
     * @param neigh Neighbor
     */
    virtual local_or_ghost_cell_index_type
    cell_neighbor_index(local_cell_index_type cellidx, fs_neighidx neigh);

    /** Returns the ghost exchange info.
     * @see GhostExchangeDesc
     */
    virtual util::const_span<GhostExchangeDesc> get_boundary_info();

    /** Returns the index of a local cell at position "pos".
     * @throws std::domain_error if position is not in the local subdomain.
     */
    virtual local_cell_index_type position_to_cell_index(Vec3d pos);

    /** Returns the rank of the process which is responsible for the cell at
     * position "pos". Before the first call to repartition() is guaranteed to
     * work for the whole domain! After the first repartition() might only work
     * for the process itself and its neighbors or its ghost layer.
     *
     * @throws std::runtime_error if position cannot be resolved because the
     * specific class supports resolving only its subdomain and ghost layer (see
     * above).
     */
    virtual rank_type position_to_rank(Vec3d pos);

    /** *Maybe* repartitions the grid. Returns true if grid has been changed
     * (repartitioned). This means all data of this class is invalidated.
     * If false is returned, *no* data returned since the last call to
     * repartition() is invalidated.
     *
     * The data invalidation includes cell indices. These silently get a new
     * meaning (underlying global cell index).
     *
     * @param exchange_start_callback is a function with no arguments which
     * starts the data migration. This function is only called if the return
     * value is "true". Also, it is called as soon as "position_to_rank" can
     * safely be called.
     */
    virtual bool
    repartition(CellMetric m,
                CellCellMetric ccm,
                Thunk exchange_start_callback);

    /** Deliver implementation-defined commands to the partitioner.
     *
     * @throws UnknownCommandError if command cannot be interpreted.
     */
    virtual void command(std::string s);

    /** Returns a globally unique id for a local cell.
     * This id is uniquely assigned to the global cell corresponding to a local
     * one, i.e. two different processes will return the same global_hash if the
     * (most likely different) local cellidxs correspond to the same global
     * cell. If NDEBUG is set, additionally to the above stated semantics, this
     * function is allowed to return constant 0.
     *
     * This function is useful for testing purposes only.
     * Use *only* if NDEBUG is *not* set.
     *
     */
    virtual global_cell_index_type
    global_hash(local_or_ghost_cell_index_type cellidx);
};
```

## Variants

Some grids have different variants (possibilities) of used algorithms.
Query them like so:

```c++
auto grid = ...;
const auto supported_variants = variants(grid.get()).get_supported_variants();
```

Set them:

```c++
variants(grid.get()).set_variant("...");
```

A call to `get_supported_variants` is always well-defined and might return an
empty set if the grid does not support any variants.
A call to `set_variants(X)` is only well-defined if `X` is in
`get_supported_variants()`. Otherwise, this function will throw an error.

## Doxygen Documentation

Doxygen-generated documentation of all the interfaces, classes and functions
can be found at: [https://librepa.github.io](https://librepa.github.io)
