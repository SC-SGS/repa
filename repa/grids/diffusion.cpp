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

#include "diffusion.hpp"
#include <algorithm>
#include <boost/mpi/nonblocking.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>
#include <numeric>

#include "util/ensure.hpp"
#include "util/fill.hpp"
#include "util/mpi_graph.hpp"
#include "util/push_back_unique.hpp"

#ifndef NDEBUG
#define DIFFUSION_DEBUG
#endif

namespace boost {
namespace serialization {
template <typename Archive>
void load(Archive &ar,
          repa::grids::Diffusion::CellNeighborhood &n,
          const unsigned int /* file_version */)
{
    ar >> n.basecell;
    ar >> n.neighranks;
}

template <typename Archive>
void save(Archive &ar,
          const repa::grids::Diffusion::CellNeighborhood &n,
          const unsigned int /* file_version */)
{
    ar << n.basecell;
    ar << n.neighranks;
}

template <class Archive>
void serialize(Archive &ar,
               repa::grids::Diffusion::CellNeighborhood &n,
               const unsigned int file_version)
{
    split_free(ar, n, file_version);
}
} // namespace serialization
} // namespace boost

/*
 * Determines the status of each process (underloaded, overloaded)
 * in the neighborhood given the local load and returns the volume of load to
 * send to each neighbor. On underloaded processes, returns a vector of zeros.
 *
 * This call is collective on neighcomm.
 *
 * See Willebeek Le Mair and Reeves, IEEE Tr. Par. Distr. Sys. 4(9), Sep 1993
 *
 * @param neighcomm Graph communicator which reflects the neighbor relationship
 *                  amongst processes (undirected edges), without edges to the
 *                  process itself.
 * @param load The load of the calling process.
 * @returns Vector of load values ordered according to the neighborhood
 *          ordering in neighcomm.
 */
static std::vector<double> compute_send_volume(MPI_Comm neighcomm, double load)
{
    int nneigh = repa::util::mpi_undirected_neighbor_count(neighcomm);
    // Exchange load in local neighborhood
    std::vector<double> neighloads(nneigh);
    MPI_Neighbor_allgather(&load, 1, MPI_DOUBLE, neighloads.data(), 1,
                           MPI_DOUBLE, neighcomm);

    double avgload
        = std::accumulate(std::begin(neighloads), std::end(neighloads), load)
          / (nneigh + 1);

    // Return empty send volume if this process is underloaded
    if (load < avgload)
        return std::vector<double>(neighloads.size(), 0.0);

    std::vector<double> deficiency(neighloads.size());

    // Calculate deficiency
    for (size_t i = 0; i < neighloads.size(); ++i) {
        deficiency[i] = std::max(avgload - neighloads[i], 0.0);
    }

    auto total_deficiency
        = std::accumulate(std::begin(deficiency), std::end(deficiency), 0.0);
    double overload = load - avgload;

    // Make "deficiency" relative and then scale it to be an
    // absolute part of this process's overload
    for (size_t i = 0; i < neighloads.size(); ++i) {
        deficiency[i] = overload * deficiency[i] / total_deficiency;
    }

    return deficiency;
}

namespace repa {
namespace grids {

void Diffusion::clear_unknown_cell_ownership()
{
    auto is_my_cell = [this](local_or_ghost_cell_index_type neighcell) {
        return partition[neighcell] == comm_cart.rank();
    };

    fill_if_index(std::begin(partition), std::end(partition), -1,
                  [this, is_my_cell](size_t glocellidx) {
                      auto neighborhood = gbox.full_shell_neigh(glocellidx);
                      return std::none_of(std::begin(neighborhood),
                                          std::end(neighborhood), is_my_cell);
                  });
}

bool Diffusion::sub_repartition(CellMetric m, CellCellMetric ccm)
{
    auto cellweights = m();
    if (cellweights.size() != n_local_cells()) {
        throw std::runtime_error(
            "Metric only supplied " + std::to_string(cellweights.size())
            + "weights. Necessary: " + std::to_string(n_local_cells()));
    }

    clear_unknown_cell_ownership();

    // compute local, estimated load
    double local_load
        = std::accumulate(std::begin(cellweights), std::end(cellweights), 0.0);

    std::vector<double> send_volume
        = compute_send_volume(neighcomm, local_load);
#ifdef DIFFUSION_DEBUG
    ENSURE(send_volume.size() == neighbors.size());
#endif

    PerNeighbor<GlobalCellIndices> toSend(neighbors.size());

    if (std::any_of(std::begin(send_volume), std::end(send_volume),
                    [](double d) { return d > 0.0; })) {

        // Create list of border cells which can be send to neighbors.
        // First element of each vector is the rank to which the sells should
        // be send. The other elements in the vectors are the global cell IDs
        // of the cells.
        toSend = compute_send_list(std::move(send_volume), cellweights);

        // Update partition array
        for (size_t i = 0; i < toSend.size(); ++i) {
            fill_index_range(partition, std::begin(toSend[i]),
                             std::end(toSend[i]), neighbors[i]);
        }
    }

    //
    // First communication step
    // Send *all* vectors in "toSend" to *all* neighbors.
    // (Not only their respective receive volumes.)
    // This is used to avoid inconsistencies, especially at newly created
    // neighborhood relationships
    //
    std::vector<boost::mpi::request> sreq_cells(neighbors.size());
    std::vector<boost::mpi::request> rreq_cells(neighbors.size());

    for (rank_index_type i = 0; i < neighbors.size(); ++i) {
        // Push back the rank, that is save an extra communication of
        // "neighbors" and interleave it into "toSend"
        toSend[i].push_back(static_cast<global_cell_index_type>(neighbors[i]));
    }

    // Extra loop as all ranks need to be added before sending
    for (rank_index_type i = 0; i < neighbors.size(); ++i) {
        sreq_cells[i] = comm_cart.isend(neighbors[i], 2, toSend);
    }

    // All send volumes from all processes
    PerNeighbor<PerNeighbor<GlobalCellIndices>>
    //          ^^^^^^^^^^^ this "PerNeighbor" is actually "PerNeighborsNeighbor"
        received_cells(neighbors.size());
    for (rank_index_type i = 0; i < neighbors.size(); ++i) {
        rreq_cells[i] = comm_cart.irecv(neighbors[i], 2, received_cells[i]);
    }

    boost::mpi::wait_all(std::begin(rreq_cells), std::end(rreq_cells));

    // Update the partition entry for all received cells.
    for (size_t from = 0; from < received_cells.size(); ++from) {
        for (size_t to = 0; to < received_cells[from].size(); ++to) {
            // Extract target rank, again.
            rank_type target_rank = static_cast<rank_type>(received_cells[from][to].back());
            received_cells[from][to].pop_back();

            fill_index_range(partition, std::begin(received_cells[from][to]),
                             std::end(received_cells[from][to]), target_rank);
        }
    }

    boost::mpi::wait_all(std::begin(sreq_cells), std::end(sreq_cells));

    //
    // END of first communication step
    //
#ifdef DIFFUSION_DEBUG
    std::vector<int> p2 = partition;
    for (auto &el : p2)
        if (el != comm_cart.rank())
            el = -1;

    MPI_Allreduce(MPI_IN_PLACE, p2.data(), p2.size(), MPI_INT, MPI_MAX,
                  comm_cart);
    for (auto el : p2)
        ENSURE(el > -1);
#endif

    // Remove ranks from "toSend", again.
    for (rank_index_type i = 0; i < neighbors.size(); ++i)
        toSend[i].pop_back();

    //
    // Second communication Step
    // Send neighbourhood of sent cells.
    //
    std::vector<boost::mpi::request> rreq_neigh(neighbors.size());
    std::vector<boost::mpi::request> sreq_neigh(neighbors.size());

    auto sendVectors = sendNeighbourhood(toSend);

    for (rank_index_type i = 0; i < neighbors.size(); ++i) {
        sreq_neigh[i] = comm_cart.isend(neighbors[i], 2, sendVectors[i]);
    }

    // All send volumes from all processes
    PerNeighbor<CellNeighborhoodPerCell> received_neighborhood(neighbors.size());
    for (rank_index_type i = 0; i < neighbors.size(); ++i) {
        rreq_neigh[i]
            = comm_cart.irecv(neighbors[i], 2, received_neighborhood[i]);
    }

    boost::mpi::wait_all(std::begin(rreq_neigh), std::end(rreq_neigh));
    updateReceivedNeighbourhood(received_neighborhood);
    boost::mpi::wait_all(std::begin(sreq_neigh), std::end(sreq_neigh));

#ifdef DIFFUSION_DEBUG
    for (global_cell_index_type i = 0; i < partition.size(); ++i) {
        if (partition[i] != comm_cart.rank())
            continue;

        for (int j = 0; j < 27; ++j) {
            global_cell_index_type n = gbox.neighbor(i, j);
            ENSURE(partition[n] > -1);
        }
    }
#endif

    return true;
}
/*
 * Initialization
 */
Diffusion::Diffusion(const boost::mpi::communicator &comm,
                     Vec3d box_size,
                     double min_cell_size)
    : GloMethod(comm, box_size, min_cell_size), neighcomm(MPI_COMM_NULL)
{
}

Diffusion::~Diffusion()
{
    int finalized = 0;
    MPI_Finalized(&finalized);
    if (neighcomm != MPI_COMM_NULL && !finalized)
        MPI_Comm_free(&neighcomm);
}

/*
 * Computes a vector of vectors. The inner vectors contain a rank of the
 * process where the cells shall send and the cellids of this cells.
 */
Diffusion::PerNeighbor<Diffusion::GlobalCellIndices>
Diffusion::compute_send_list(std::vector<double> &&send_loads,
                             const std::vector<double> &weights)
{
    std::vector<std::tuple<int, double, local_cell_index_type>> plist;
    for (size_t i = 0; i < borderCells.size(); i++) {
        // Profit when sending this cell away
        double profit = weights[borderCells[i]];

        // Additional cell communication induced if this cell is sent away
        int nadditional_comm = 0;
        for (global_cell_index_type neighCell :
             gbox.full_shell_neigh_without_center(cells[borderCells[i]])) {
            if (partition[neighCell] == comm_cart.rank()
                && std::find(std::begin(borderCells), std::end(borderCells),
                             global_to_local[neighCell])
                       != std::end(borderCells)) {
                nadditional_comm++;
            }
        }
#ifdef DIFFUSION_DEBUG
        ENSURE(nadditional_comm < 27);
#endif

        if (profit > 0)
            plist.emplace_back(27 - nadditional_comm, profit, borderCells[i]);
    }

    PerNeighbor<GlobalCellIndices> to_send(send_loads.size());

    // Use a maxheap: Always draw the maximum element
    // (1. least new border cells, 2. most profit)
    // and find a process that can take this cell.
    std::make_heap(std::begin(plist), std::end(plist));
    while (!plist.empty()) {
        std::pop_heap(std::begin(plist), std::end(plist));
        local_cell_index_type cidx = std::get<2>(plist.back());
        plist.pop_back();

        for (auto neighrank : borderCellsNeighbors[cidx]) {
            auto neighidx
                = std::distance(std::begin(neighbors),
                                std::find(std::begin(neighbors),
                                          std::end(neighbors), neighrank));

            if (weights[cidx] <= send_loads[neighidx]) {
                to_send[neighidx].push_back(cells[cidx]);
                send_loads[neighidx] -= weights[cidx];
                // This cell is done. Continue with the next.
                break;
            }
        }
    }

    return to_send;
}

Diffusion::PerNeighbor<Diffusion::CellNeighborhoodPerCell>
Diffusion::sendNeighbourhood(
    const PerNeighbor<GlobalCellIndices> &toSend)
{
    PerNeighbor<CellNeighborhoodPerCell> sendVectors(toSend.size());
    for (size_t i = 0; i < toSend.size(); ++i) {
        sendVectors[i].resize(toSend[i].size());
        for (size_t j = 0; j < toSend[i].size(); ++j) {
            sendVectors[i][j].basecell = toSend[i][j];
            int k = 0;
            for (global_cell_index_type n :
                 gbox.full_shell_neigh_without_center(
                     sendVectors[i][j].basecell)) {
                sendVectors[i][j].neighranks[k] = partition[n];
                k++;
            }
        }
    }

    return sendVectors;
}

/*
 * Based on neighbourhood, received in function "receiveNeighbourhood",
 * partition array is updated. (Only neighbourhood is changed)
 */
void Diffusion::updateReceivedNeighbourhood(
    const PerNeighbor<CellNeighborhoodPerCell> &neighs)
{
    for (size_t i = 0; i < neighs.size(); ++i) {
        for (size_t j = 0; j < neighs[i].size(); ++j) {
            global_cell_index_type basecell = neighs[i][j].basecell;
            int k = 0;
            for (global_cell_index_type n :
                 gbox.full_shell_neigh_without_center(basecell)) {
                partition[n] = neighs[i][j].neighranks[k++];
            }
        }
    }
}

void Diffusion::pre_init(bool firstcall)
{
    borderCells.clear();
    borderCellsNeighbors.clear();

    if (!firstcall)
        clear_unknown_cell_ownership();
}

void Diffusion::post_init(bool firstcall)
{
    // Create graph comm with current process structure
    if (neighcomm != MPI_COMM_NULL)
        MPI_Comm_free(&neighcomm);

    // Edges to all processes in "neighbors"
    MPI_Dist_graph_create_adjacent(
        comm_cart, neighbors.size(), neighbors.data(),
        static_cast<const int *>(MPI_UNWEIGHTED), neighbors.size(),
        neighbors.data(), static_cast<const int *>(MPI_UNWEIGHTED),
        MPI_INFO_NULL, 0, &neighcomm);

#ifdef DIFFUSION_DEBUG
    int indegree = 0, outdegree = 0, weighted = 0;
    MPI_Dist_graph_neighbors_count(neighcomm, &indegree, &outdegree, &weighted);
    ENSURE(!weighted);
    ENSURE(static_cast<size_t>(indegree) == neighbors.size());
    ENSURE(static_cast<size_t>(outdegree) == neighbors.size());
    std::vector<int> __ineighs(indegree, -1), __iw(indegree, -1);
    std::vector<int> __oneighs(outdegree, -1), __ow(outdegree, -1);
    MPI_Dist_graph_neighbors(neighcomm, indegree, __ineighs.data(), __iw.data(),
                             outdegree, __oneighs.data(), __ow.data());
    for (size_t i = 0; i < neighbors.size(); ++i) {
        ENSURE(__ineighs[i] == neighbors[i]);
        ENSURE(__oneighs[i] == neighbors[i]);
    }
#endif
}

void Diffusion::init_new_foreign_cell(local_cell_index_type localcell,
                                      global_cell_index_type foreigncell,
                                      rank_type owner)
{
    // First cell identifying "localcell" as border cell?
    if (borderCells.empty() || borderCells.back() != localcell)
        borderCells.push_back(localcell);

    util::push_back_unique(borderCellsNeighbors[localcell], owner);
}

} // namespace grids
} // namespace repa
