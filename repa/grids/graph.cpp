/**
 * Copyright 2017-2019 Steffen Hirschmann
 * Copyright 2017-2018 Maximilian Wildbrett
 *
 * This file is part of Repa.
 *
 * Repa is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Repa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Repa.  If not, see <https://www.gnu.org/licenses/>.
 */

//#ifdef HAVE_METIS

#include "graph.hpp"
#include <algorithm>
#include <boost/mpi.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>
#include <mpi.h>

#include "util/ensure.hpp"
#include "util/mpi_type.hpp"
#include "util/push_back_unique.hpp"

#define MPI_IDX_T util::mpi_type<idx_t>::type

#ifndef NDEBUG
#define GRAPH_DEBUG
#endif

namespace repa {
namespace grids {

lidx Graph::n_local_cells()
{
    return localCells;
}

gidx Graph::n_ghost_cells()
{
    return ghostCells;
}

nidx Graph::n_neighbors()
{
    return neighbors.size();
}

rank Graph::neighbor_rank(nidx i)
{
    return neighbors[i];
}

Vec3d Graph::cell_size()
{
    return gbox.cell_size();
}

Vec3i Graph::grid_size()
{
    return gbox.grid_size();
}

lgidx Graph::cell_neighbor_index(lidx cellidx, int neigh)
{
    return global_to_local[gbox.neighbor(cells[cellidx], neigh)];
}

std::vector<GhostExchangeDesc> Graph::get_boundary_info()
{
    return exchangeVector;
}

lidx Graph::position_to_cell_index(double pos[3])
{
    if (position_to_rank(pos) != comm_cart.rank())
        throw std::domain_error("Particle not in local box");

    return global_to_local[gbox.cell_at_pos(pos)];
}

rank Graph::position_to_rank(double pos[3])
{
    return static_cast<rank>(partition[gbox.cell_at_pos(pos)]);
}

nidx Graph::position_to_neighidx(double pos[3])
{
    rank rank = position_to_rank(pos);
    auto ni = std::find(std::begin(neighbors), std::end(neighbors), rank);

    if (ni != std::end(neighbors))
        return std::distance(std::begin(neighbors), ni);
    else
        throw std::domain_error("Position not within a neighbor process.");
}

/*
 * Repartition.
 * Every node is responsible for a certain range of cells along the
 * linearization defined by gbox. Every process evaluates the weights for its
 * cells and sends them to the process that is responsible for the graph node
 * that corresponds to the cell.
 * Partitioning is performed in parallel via ParMETIS.
 */
bool Graph::repartition(CellMetric m,
                        CellCellMetric ccm,
                        Thunk exchange_start_callback)
{
    static constexpr idx_t w_fac = 100;
    auto vertex_weights = m();
    if (vertex_weights.size() != n_local_cells()) {
        throw std::runtime_error(
            "Metric only supplied " + std::to_string(vertex_weights.size())
            + "weights. Necessary: " + std::to_string(n_local_cells()));
    }

    idx_t nglocells = static_cast<idx_t>(gbox.ncells());
    idx_t ncells_per_proc = static_cast<idx_t>(
        std::ceil(static_cast<double>(nglocells) / comm_cart.size()));

#ifdef GRAPH_DEBUG
    // if (comm_cart.rank() == 0)
    //  std::cout << "Graph::repartition with comm_cart.size() = " <<
    //  comm_cart.size() << std::endl;
#endif

    // Vertex ranges per process
    std::vector<idx_t> vtxdist(comm_cart.size() + 1);
    for (int i = 0; i < comm_cart.size(); ++i)
        vtxdist[i] = i * ncells_per_proc;
    vtxdist[comm_cart.size()] = nglocells;

#ifdef GRAPH_DEBUG
    // std::cout << "Vtxdist: ";
    // std::copy(std::begin(vtxdist), std::end(vtxdist),
    // std::ostream_iterator<idx_t>(std::cout, " ")); std::cout << std::endl;
#endif

#ifdef GRAPH_DEBUG
    ENSURE(vtxdist.size() == comm_cart.size() + 1);
    for (int i = 0; i < comm_cart.size(); ++i) {
        ENSURE(0 <= vtxdist[i] && vtxdist[i] < nglocells);
    }
    ENSURE(vtxdist[comm_cart.size()] == nglocells);
#endif

    // Receive vertex and edge weights
    std::vector<int> recvranks;
    // Determine ranks from which to receive
    // (Via old "partition" field and the range of graph nodes this
    // process is responsible for.)
    for (idx_t i = vtxdist[comm_cart.rank()]; i < vtxdist[comm_cart.rank() + 1];
         ++i)
        util::push_back_unique(recvranks, static_cast<int>(partition[i]));

#ifdef GRAPH_DEBUG
    for (int r : recvranks) {
        ENSURE(0 <= r && r < comm_cart.size());
        ENSURE(std::count(std::begin(recvranks), std::end(recvranks), r) == 1);
    }
    ENSURE(recvranks.size() <= comm_cart.size());
#endif

    // [0]: vertex weight
    // [1-26]: edge weights
    using Weights = std::array<idx_t, 27>;

    std::vector<boost::mpi::request> rreq(comm_cart.size());
    std::vector<std::vector<Weights>> their_weights(comm_cart.size());
    idx_t wsum = 0; // Check for possible overflow
    for (int n : recvranks)
        rreq[n] = comm_cart.irecv(n, 0, their_weights[n]);

    // Sending vertex weights
    std::vector<boost::mpi::request> sreq(comm_cart.size());
    std::vector<std::vector<Weights>> my_weights(comm_cart.size());
    for (int i = 0; i < localCells; ++i) {
        // "Rank" is responsible for cell "gidx" / "i" (local)
        // during graph parititioning
        int gidx = cells[i];
        rank rank = gidx / ncells_per_proc;

        Weights w;
        w[0] = static_cast<idx_t>(vertex_weights[i]);

        // Edge weights must not be strictly positive. Add 1 and scale the true
        // weight. Also leave out the loop (i, i), i.e. start at n = 1.
        for (int n = 1; n < 27; ++n) {
            auto neigh = cell_neighbor_index(i, n);
            w[n] = static_cast<idx_t>(ccm(i, neigh)) * w_fac + 1;
#ifdef GRAPH_DEBUG
            ENSURE(w[n] > 0);
            if (neigh < n_local_cells()) {
                // Local symmetry -- only ensures local symmetry, however,
                // symmetry is also required for cross-boundary edges.
                ENSURE(w[n] == ccm(neigh, i) * w_fac + 1);
            }
#endif
        }

        my_weights[rank].push_back(std::move(w));

        // Catch total weight sum overlow
        wsum += std::accumulate(std::begin(w), std::end(w),
                                static_cast<idx_t>(0));
        if (wsum > std::numeric_limits<idx_t>::max() / comm_cart.size()) {
            std::fprintf(stderr,
                         "Warning: Graph weights are too large for chosen "
                         "index type width.\n");
            if (w_fac > 1)
                std::fprintf(stderr,
                             "- Try to reduce the weight factor in graph.cpp.");
            if (IDXTYPEWIDTH == 32)
                std::fprintf(stderr, "- Recompile a 64bit ParMETIS as this "
                                     "version only uses 32bit.");
            errexit();
        }
    }

    for (int i = 0; i < comm_cart.size(); ++i) {
        if (my_weights[i].size() > 0)
            sreq[i] = comm_cart.isend(i, 0, my_weights[i]);
    }

    // Prepare graph

    idx_t nvtx = vtxdist[comm_cart.rank() + 1] - vtxdist[comm_cart.rank()];

    // Regular grid as graph
    std::vector<idx_t> xadj(nvtx + 1), adjncy(26 * nvtx);
    for (idx_t i = 0; i < nvtx; ++i) {
        xadj[i] = 26 * i;
        for (int n = 0; n < 26; ++n) {
            adjncy[26 * i + n]
                = gbox.neighbor(vtxdist[comm_cart.rank()] + i, n + 1);
        }
    }
    xadj[nvtx] = 26 * nvtx;

#ifdef GRAPH_DEBUG
    // std::cout << "xadj: ";
    // std::copy(std::begin(xadj), std::end(xadj),
    // std::ostream_iterator<idx_t>(std::cout, " ")); std::cout << std::endl;
#endif

#ifdef GRAPH_DEBUG
    ENSURE(nvtx == vtxdist[comm_cart.rank() + 1] - vtxdist[comm_cart.rank()]);
    ENSURE(nvtx <= nglocells);
    ENSURE(xadj.size() == nvtx + 1);
    for (int i = 0; i < nvtx; ++i) {
        ENSURE(xadj[i] < xadj[i + 1]);
        ENSURE(xadj[i + 1] - xadj[i] == 26);
    }

    for (int i = 0; i < adjncy.size(); ++i) {
        ENSURE(adjncy[i] >= 0 && adjncy[i] < nglocells);
    }
#endif

    std::vector<idx_t> vwgt(nvtx), adjwgt(26 * nvtx);
#ifdef GRAPH_DEBUG
    std::fill(std::begin(vwgt), std::end(vwgt), -1);
    std::fill(std::begin(adjwgt), std::end(adjwgt), -1);
#endif

    boost::mpi::wait_all(std::begin(rreq), std::end(rreq));

    std::vector<size_t> ii(comm_cart.size(), 0); // Iteration indices
    size_t li = 0;
    for (idx_t i = vtxdist[comm_cart.rank()]; i < vtxdist[comm_cart.rank() + 1];
         ++i) {
        auto rank = partition[i];
        auto &w = their_weights[rank][ii[rank]++];
        vwgt[li] = w[0];
        for (int n = 0; n < 26; ++n) {
            adjwgt[26 * li + n] = w[n + 1];
        }
        li++;
    }

#ifdef GRAPH_DEBUG
    for (int i = 0; i < comm_cart.size(); ++i) {
        ENSURE(ii[i] == their_weights[i].size());
    }
    ENSURE(li == nvtx);

    for (idx_t w : vwgt) {
        ENSURE(w != -1);
        ENSURE(w >= 0 && w < 10000000);
    }
    for (idx_t w : adjwgt) {
        ENSURE(w != -1);
        ENSURE(w >= 0 && w < 10000000);
        //if (!m.has_cell_cell_metric())
        //    ENSURE(w == 1.0 * w_fac + 1);
    }
#endif

    boost::mpi::wait_all(std::begin(sreq), std::end(sreq));

    // Weights -- 0: no, 1: edges only, 2: vertices only, 3: both
    idx_t wgtflag = 3;
    // Numbering -- 0: C-style, 1: Fortran-style
    idx_t numflag = 0;
    // Number of weights per vertex
    idx_t ncon = 1;
    // Number of wanted subdomains
    idx_t nparts = static_cast<idx_t>(comm_cart.size());
    // Target weights -- must add up to 1.0
    std::vector<real_t> tpwgts(nparts, 1.0 / nparts);
    // Imbalance tolerance
    real_t ubvec = 1.05; // Recommended by ParMETIS doc

    // METIS options -- ParMETIS does not use METIS's options format
    // idx_t options[METIS_NOPTIONS];
    // METIS_SetDefaultOptions(options);
    // options[METIS_OPTION_MINCONN] = 1;
    // options[METIS_OPTION_CONTIG] = 1;
    // options[METIS_OPTION_NOOUTPUT] = 1;
    idx_t options[3] = {0, 0, comm_cart.rank()};

    // Comm
    MPI_Comm communicator = comm_cart;

    // Result parameters
    idx_t edgecut;
    std::vector<idx_t> part(nvtx, -1);

    // if (comm_cart.rank() == 0)
    //  std::cout << "Calling ParMETIS_V3_PartKway with " << nglocells << "
    //  total vertices for " << nparts << " subdomains." << std::endl;

    if (ParMETIS_V3_PartKway(vtxdist.data(), xadj.data(), adjncy.data(),
                             vwgt.data(), adjwgt.data(), &wgtflag, &numflag,
                             &ncon, &nparts, tpwgts.data(), &ubvec, options,
                             &edgecut, part.data(), &communicator)
        != METIS_OK) {
        if (comm_cart.rank() == 0) {
            std::fprintf(
                stderr,
                "Error when graph partitioning: ParMETIS returned error.\n");
            errexit();
        }
    }

#ifdef GRAPH_DEBUG
    ENSURE(part.size() == nvtx);
    for (int r : part) {
        ENSURE(r != -1);
        ENSURE(0 <= r && r < comm_cart.size());
    }
#endif

#ifdef GRAPH_DEBUG
    std::fill(std::begin(partition), std::end(partition),
              static_cast<idx_t>(-1));
#endif
    // int li = 0;
    // for (idx_t i = vtxdist[comm_cart.rank()]; i < vtxdist[comm_cart.rank() +
    // 1]; ++i) {
    //  partition[i] = part[li++];
    //}

    // MPI expects integer recvcounts and displacements for MPI_Allgatherv.
    std::vector<int> recvcount(comm_cart.size()), displ(comm_cart.size());
    for (size_t i = 0; i < comm_cart.size(); ++i) {
        recvcount[i] = static_cast<int>(vtxdist[i + 1] - vtxdist[i]);
        displ[i] = static_cast<int>(vtxdist[i]);
    }
    MPI_Allgatherv(part.data(), nvtx, MPI_IDX_T, partition.data(),
                   recvcount.data(), displ.data(), MPI_IDX_T, comm_cart);
#ifdef GRAPH_DEBUG
    ENSURE(partition.size() == nglocells);
    for (int r : partition) {
        ENSURE(r != -1);
        ENSURE(0 <= r && r < comm_cart.size());
    }
#endif

#ifdef GRAPH_DEBUG
    int nlc = static_cast<int>(std::count(
        std::begin(partition), std::end(partition), comm_cart.rank()));

    std::vector<int> nlcs(comm_cart.size());
    MPI_Allgather(&nlc, 1, MPI_INT, nlcs.data(), 1, MPI_INT, comm_cart);

    if (comm_cart.rank() == 0) {
        std::cout << "Partitioning result:";
        std::copy(std::begin(nlcs), std::end(nlcs),
                  std::ostream_iterator<int>(std::cout, " "));
        std::cout << std::endl;
    }

    ENSURE(std::accumulate(std::begin(nlcs), std::end(nlcs), 0) == nglocells);
#endif

    // Position to rank is answered solely via the "partition" vector.
    // So the particle migration is ready to go.
    exchange_start_callback();

    init();

    return true;
}

Graph::Graph(const boost::mpi::communicator &comm,
             Vec3d box_size,
             double min_cell_size)
    : ParallelLCGrid(comm, box_size, min_cell_size),
      gbox(box_size, min_cell_size)
{
    int nglocells = gbox.ncells();
    int ncells_per_proc = static_cast<int>(
        std::ceil(static_cast<double>(nglocells) / comm_cart.size()));

    // Initial partitioning
    partition.resize(nglocells);

    // Line-wise init
    for (int i = 0; i < nglocells; ++i) {
        partition[i] = i / ncells_per_proc;
    }

    //// Init to equally sized boxes on Cartesian grid
    // int dims[3] = {0, 0, 0};
    // MPI_Dims_create(comm_cart.size(), 3, dims);

    // auto cellgrid = gbox.grid_size();
    // Vec3i cells_per_proc = {{
    //    static_cast<int>(std::ceil(static_cast<double>(cellgrid[0]) /
    //    dims[0])), static_cast<int>(std::ceil(static_cast<double>(cellgrid[1])
    //    / dims[1])),
    //    static_cast<int>(std::ceil(static_cast<double>(cellgrid[2]) /
    //    dims[2])),
    //}};

    // for (int i = 0; i < nglocells; ++i) {
    //  auto cellidx = gbox.unlinearize(i);
    //  // Transform cellidx to 3d proc coord
    //  for (int i = 0; i < 3; ++i)
    //    cellidx[i] /= cells_per_proc[i];
    //  int rank;
    //  MPI_Cart_rank(comm_cart, cellidx.data(), &rank);
    //  partition[i] = rank;
    //}

    init();
}

Graph::~Graph()
{
}

/*
 * Rebuild the data structures describing subdomain and communication.
 */
void Graph::init()
{
    const int nglocells = partition.size();

    localCells = 0;
    ghostCells = 0;
    cells.clear();
    global_to_local.clear();
    neighbors.clear();

    // Extract the local cells from "partition".
    for (int i = 0; i < nglocells; i++) {
        if (partition[i] == comm_cart.rank()) {
            // Vector of own cells
            cells.push_back(i);
            // Index mapping from global to local
            global_to_local[i] = localCells;
            // Number of own cells
            localCells++;
        }
    }

    // Temporary storage for exchange descriptors.
    // Will be filled only for neighbors
    // and moved from later.
    std::vector<GhostExchangeDesc> tmp_ex_descs(comm_cart.size());

    // Determine ghost cells and communication volume
    for (int i = 0; i < localCells; i++) {
        for (int neighborIndex :
             gbox.full_shell_neigh_without_center(cells[i])) {
            rank owner = static_cast<rank>(partition[neighborIndex]);
            if (owner == comm_cart.rank())
                continue;

            // Find ghost cells. Add only once to "cells" vector.
            if (global_to_local.find(neighborIndex)
                == std::end(global_to_local)) {
                // Add ghost cell to cells vector
                cells.push_back(neighborIndex);
                // Index mapping from global to ghost
                global_to_local[neighborIndex] = localCells + ghostCells;
                // Number of ghost cells
                ghostCells++;
            }

            // Initialize exdesc and add "rank" as neighbor if unknown.
            if (tmp_ex_descs[owner].dest == -1) {
                neighbors.push_back(owner);
                tmp_ex_descs[owner].dest = owner;
            }

            util::push_back_unique(tmp_ex_descs[owner].recv, neighborIndex);
            util::push_back_unique(tmp_ex_descs[owner].send, cells[i]);
        }
    }

    // Move all existent exchange descriptors from "tmp_ex_descs" to
    // "exchangeVector".
    exchangeVector.clear();
    for (int i = 0; i < comm_cart.size(); ++i) {
        if (tmp_ex_descs[i].dest != -1) {
            auto ed = std::move(tmp_ex_descs[i]);

            // Make sure, index ordering is the same on every process
            // and global to local index conversion
            std::sort(std::begin(ed.recv), std::end(ed.recv));
            std::transform(std::begin(ed.recv), std::end(ed.recv),
                           std::begin(ed.recv),
                           [this](int i) { return global_to_local[i]; });
            std::sort(std::begin(ed.send), std::end(ed.send));
            std::transform(std::begin(ed.send), std::end(ed.send),
                           std::begin(ed.send),
                           [this](int i) { return global_to_local[i]; });

            exchangeVector.push_back(std::move(ed));
        }
    }

#ifdef GRAPH_DEBUG
    for (int i = 0; i < comm_cart.size(); ++i) {
        if (tmp_ex_descs[i].dest != -1)
            ENSURE(tmp_ex_descs[i].recv.size() == 0
                   && tmp_ex_descs[i].send.size() == 0);
    }
#endif
}

int Graph::global_hash(lgidx cellidx)
{
    return cells[cellidx];
}

} // namespace grids
} // namespace repa

//#endif
