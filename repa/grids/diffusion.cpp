/**
 * Copyright 2017-2019 The repa authors
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
#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/utility.hpp> // std::pair
#include <boost/serialization/vector.hpp>
#include <numeric>
#include <regex>

#include "util/get_keys.hpp"
#include "util/mpi_graph.hpp"
#include "util/mpi_neighbor_allgather.hpp"
#include "util/mpi_neighbor_alltoall.hpp"
#include "util/push_back_unique.hpp"
#include "util/range.hpp"

#ifndef NDEBUG
#define DIFFUSION_DEBUG
#endif

namespace _impl {

using GlobalIndices = std::vector<repa::global_cell_index_type>;
using SendVolIndices = std::vector<GlobalIndices>;
using RankVector = std::vector<repa::rank_type>;

template <typename F>
void for_each_reassignment(
    const std::pair<SendVolIndices, RankVector> &sendvolume, F &&f)
{
    const auto &indicess = std::get<0>(sendvolume);
    const auto &new_values = std::get<1>(sendvolume);
    assert(indicess.size() == new_values.size());
    // Do *not* loop over new_values, "indicess" can be empty.
    for (size_t set_num = 0; set_num < indicess.size(); ++set_num) {
        const auto &indices = indicess[set_num];
        const auto &n_value = new_values[set_num];

        for (const auto &i : indices) {
            f(i, n_value);
        }
    }
}

#ifndef NDEBUG
template <typename PartitionEntryType>
static bool is_correct_distributed_partitioning(
    const std::vector<PartitionEntryType> &partition,
    const boost::mpi::communicator &comm)
{
    // Check that every cell has exactly one owner
    std::vector<int> nowners(partition.size(), 0);
    for (size_t i = 0; i < partition.size(); ++i)
        if (partition[i] == comm.rank())
            nowners[i]++;

    MPI_Allreduce(MPI_IN_PLACE, nowners.data(), nowners.size(), MPI_INT,
                  MPI_SUM, comm);
    auto x = std::all_of(std::begin(nowners), std::end(nowners),
                         [](int el) { return el == 1; });
    assert(x);
    return x;
}

template <typename GBox, typename PartitionEntryType>
bool is_ghost_layer_fully_known(
    const std::vector<PartitionEntryType> &partition,
    const boost::mpi::communicator &comm,
    const GBox &gbox)
{
    // Check that the neighborhood of every owned cell is known.
    for (const auto i :
         repa::util::range(repa::global_cell_index_type{partition.size()})) {
        if (partition[i] != comm.rank())
            continue;

        for (const auto ni : gbox.full_shell_neigh(i)) {
            if (!partition[ni].has_value())
                return false;
        }
    }
    return true;
}
#endif

} // namespace _impl

namespace boost {
namespace serialization {
template <typename Archive>
void load(Archive &ar,
          repa::grids::__diff_impl::CellNeighborhood &n,
          const unsigned int /* file_version */)
{
    ar >> n.basecell;
    ar >> n.neighranks;
}

template <typename Archive>
void save(Archive &ar,
          const repa::grids::__diff_impl::CellNeighborhood &n,
          const unsigned int /* file_version */)
{
    ar << n.basecell;
    ar << n.neighranks;
}

template <class Archive>
void serialize(Archive &ar,
               repa::grids::__diff_impl::CellNeighborhood &n,
               const unsigned int file_version)
{
    split_free(ar, n, file_version);
}
} // namespace serialization
} // namespace boost

namespace repa {
namespace grids {

static const std::unordered_map<std::string, diff_variants::FlowCalcKind>
    supported_default_diffusion_variants
    = {{"willebeek", diff_variants::FlowCalcKind::WILLEBEEK},
       {"schornbaum", diff_variants::FlowCalcKind::SCHORN}};

template <typename Rng, typename Pred>
bool none_of(const Rng &rng, Pred &&p)
{
    return std::none_of(rng.begin(), rng.end(), std::forward<Pred>(p));
}

void Diffusion::invalidate_if_unknown(global_cell_index_type cellidx)
{
    auto is_my_cell = [this](global_cell_index_type neighcell) {
        return partition[neighcell] == comm_cart.rank();
    };

    if (const auto neighborhood = gbox.full_shell_neigh(cellidx);
        none_of(neighborhood, is_my_cell)) {
        partition[cellidx] = {}; // Declare the owner to be unknown.
    }
}

void Diffusion::clear_unknown_cell_ownership()
{
    for (const global_cell_index_type i : gbox.global_cells()) {
        if (!partition[i].has_value())
            continue;
        invalidate_if_unknown(i);
    }
}

bool Diffusion::sub_repartition(CellMetric m, CellCellMetric ccm)
{
    const auto cellweights = m();
    assert(cellweights.size() == local_cells().size());

    if (stores_full_partitioning) {
        clear_unknown_cell_ownership();
        stores_full_partitioning = false;
    }

    // compute local, estimated load
    double local_load
        = std::accumulate(std::begin(cellweights), std::end(cellweights), 0.0);

    std::vector<double> send_volume
        = flow_calc->compute_flow(neighcomm, neighbors, local_load);
    assert(send_volume.size() == neighbors.size());

    const PerNeighbor<GlobalCellIndices> cells_to_send
        = compute_send_list(std::move(send_volume), cellweights);
    const auto send_information
        = std::make_pair(std::cref(cells_to_send), std::cref(neighbors));

    // Update partition array and invalidate unknown neighbors
    _impl::for_each_reassignment(
        send_information, [this](global_cell_index_type i, rank_type r) {
            partition[i] = r;
            // DO NOT DO THIS. WE NEED TO SEND THIS INFO TO NEIGHBORS LATER.
            //for (const auto neigh : gbox.full_shell_neigh_without_center(i)) {
            //    invalidate_if_unknown(neigh);
            //}
        });


    std::set<global_cell_index_type> _relevant_cells;
    // Add all ghost layer cells - not taking into account the reassignment "partition[i] = r" in line 200.
    for (const auto lidx: local_cells()) {
        for (const auto neighidx: gbox.full_shell_neigh_without_center(cell_store.as_global_index(lidx))) {
            if (cell_store.as_ghost_index(neighidx)) {
                _relevant_cells.emplace(neighidx);
            }
        }
    }


    //
    // First communication step
    // Send *all* vectors in "cells_to_send" to *all* neighbors.
    // (Not only their respective receive volumes.)
    // This is used to avoid inconsistencies, especially at newly created
    // neighborhood relationships
    //
    {
        // All send volumes from all processes
        const auto neighbor_sendinformation
            = util::mpi_neighbor_allgather(neighcomm, send_information);

        // Which cells do we need?
        // - Our own ghost layer
        // - Neighbors of "cells_to_send"

        // Which cells are these?
        // - Ghost layer: clear.
        // - Neighbors of "cells_to_send": local cells or ghost layer cells
        //   (ghost layer as per definition before l.200 reassignment of partition[i]).
        //   should still be present in "local_cells" somewhere?
        // => Only update ghost layer cells.
        //    local cells CANNOT be updated here because other processes do not own them and hence cannot reassign them.

        auto is_relevant_cell = [&_relevant_cells](global_cell_index_type i){
            return _relevant_cells.find(i) != _relevant_cells.end();
        };

        // Ganz am Ende dann "cells_to_send" rausnehmen. //
        // Moment. cells_to_send sind schon als "jemand anderem gehörend" markiert.
        std::for_each(std::begin(neighbor_sendinformation),
                      std::end(neighbor_sendinformation),
                      [this, &is_relevant_cell](const auto &neighbor_info) {
                          _impl::for_each_reassignment(
                              neighbor_info,
                              [this, &is_relevant_cell](global_cell_index_type i, rank_type r) {
                                  // Note: We must also update the neighbors of cells in "cells_to_send", here.
                                  // Because these neighborhoods are sent to their new owners in the next step.
                                  if (is_relevant_cell(i))
                                    partition[i] = r;
                              });
                      });
    }

    assert(_impl::is_correct_distributed_partitioning(partition, comm_cart));

    //
    // Second communication Step
    // Send neighbourhood of sent cells.
    //
    {
        const auto received_neighborhood_info = util::mpi_neighbor_alltoall(
            neighcomm, get_neighborhood_information(cells_to_send));

        // TODO: hier sicherstellen, dass keine falschen eingefügt werden.
        update_partitioning_from_received_neighbourhood(
            received_neighborhood_info);
    }

    // Remove unnecessary entries from "partition".
    for (const auto &i: _relevant_cells) {
        bool is_needed = false;
        for (const auto neigh: gbox.full_shell_neigh(i)) {
            if (partition[neigh] == comm.rank())
                is_needed = true;
        }
        if (!is_needed)
            partition[i] = {};
    }

    assert(_impl::is_correct_distributed_partitioning(partition, comm_cart));
    assert(_impl::is_ghost_layer_fully_known(partition, comm_cart, gbox));

    return true;
}

Diffusion::Diffusion(const boost::mpi::communicator &comm,
                     Vec3d box_size,
                     double min_cell_size,
                     ExtraParams ep)
    : GloMethod(comm, box_size, min_cell_size, ep),
      stores_full_partitioning(true),
      flow_calc(diff_variants::create_flow_calc(
          diff_variants::FlowCalcKind::WILLEBEEK))
{
    // Initial partitioning
    partition.reserve(gbox.global_cells().size());
    boost::push_back(partition, util::StaticRankAssigner{initial_partitioning,
                                                         gbox, comm_cart}
                                    .partitioning());
    assert(partition.size() == gbox.global_cells().size());
}

Diffusion::~Diffusion()
{
}

std::set<std::string> Diffusion::get_supported_variants() const
{
    return util::get_keys(supported_default_diffusion_variants);
}

void Diffusion::set_variant(const std::string &var)
{
    // Hacky
    command(std::string{"set flow "} + var);
}

Diffusion::PerNeighbor<Diffusion::GlobalCellIndices>
Diffusion::compute_send_list(std::vector<double> &&send_loads,
                             const std::vector<double> &weights) const
{
    std::vector<std::tuple<int, double, local_cell_index_type>> plist;

    // Return empty vector, if nothing to send
    if (std::none_of(send_loads.begin(), send_loads.end(),
                     [](double v) { return v > 0.; }))
        return PerNeighbor<GlobalCellIndices>(send_loads.size());

    for (size_t i = 0; i < borderCells.size(); i++) {
        // Profit when sending this cell away
        double profit = weights[borderCells[i]];

        // Additional cell communication induced if this cell is sent away
        int nadditional_comm = 0;
        for (global_cell_index_type neighCell :
             gbox.full_shell_neigh_without_center(
                 cell_store.as_global_index(borderCells[i]))) {
            // Count all local, non-border cells.
            if (partition[neighCell] != comm_cart.rank())
                continue;
            if (const auto local_index = cell_store.as_local_index(neighCell);
                local_index // Is a local and not a ghost cell
                && std::find(std::begin(borderCells), std::end(borderCells),
                             *local_index)
                       != std::end(borderCells)) {
                nadditional_comm++;
            }
        }
#ifdef DIFFUSION_DEBUG
        assert(nadditional_comm < 27);
#endif

        plist.emplace_back(27 - nadditional_comm, profit, borderCells[i]);
    }

    // A node is required to keep at least one cell.
    // Normally, the flow count should not require a process to hand away
    // all its cells.^* But we enfore this here by removing one cell
    // if all local cells are candidates.
    // * It can happen on very heterogeneous scenarios and flow count > 1.
    //   For flow_count == 1 the following code is irrelevant because a
    //   subdomain will never be required to send off *all* load.
    if (borderCells.size() == static_cast<size_t>(n_local_cells()))
        plist.pop_back();

    PerNeighbor<GlobalCellIndices> to_send(send_loads.size());
    double load = std::accumulate(weights.begin(), weights.end(), 0.0);
    std::vector<double> send_loads_copy = send_loads;

    // Use a maxheap: Always draw the maximum element
    // (1. least new border cells, 2. most profit)
    // and find a process that can take this cell.
    std::make_heap(std::begin(plist), std::end(plist));
    while (!plist.empty()) {
        std::pop_heap(std::begin(plist), std::end(plist));
        local_cell_index_type cidx = std::get<2>(plist.back());
        plist.pop_back();

        for (auto neighrank : borderCellsNeighbors.at(cidx)) {
            if (!accept_transfer(cidx, neighrank))
                continue;
            auto neighidx
                = std::distance(std::begin(neighbors),
                                std::find(std::begin(neighbors),
                                          std::end(neighbors), neighrank));

            if (weights[cidx] <= 0
                && send_loads_copy[neighidx]
                       < profit_percentage_pass_through * load)
                continue;

            if (weights[cidx] <= send_loads[neighidx]) {
                to_send[neighidx].push_back(cell_store.as_global_index(cidx));
                send_loads[neighidx] -= weights[cidx];
                // This cell is done. Continue with the next.
                break;
            }
        }
    }

    return to_send;
}

Diffusion::PerNeighbor<__diff_impl::CellNeighborhoodPerCell>
Diffusion::get_neighborhood_information(
    const PerNeighbor<GlobalCellIndices> &cells_to_send) const
{
    PerNeighbor<__diff_impl::CellNeighborhoodPerCell> sendVectors(
        cells_to_send.size());
    for (size_t i = 0; i < cells_to_send.size(); ++i) {
        sendVectors[i].resize(cells_to_send[i].size());
        for (size_t j = 0; j < cells_to_send[i].size(); ++j) {
            sendVectors[i][j].basecell = cells_to_send[i][j];
            int k = 0;
            for (global_cell_index_type n :
                 gbox.full_shell_neigh_without_center(
                     sendVectors[i][j].basecell)) {
                sendVectors[i][j].neighranks[k] = partition[n].value();
                k++;
            }
        }
    }

    return sendVectors;
}

void Diffusion::update_partitioning_from_received_neighbourhood(
    const PerNeighbor<__diff_impl::CellNeighborhoodPerCell> &neighs)
{
    for (size_t i = 0; i < neighs.size(); ++i) {
        for (size_t j = 0; j < neighs[i].size(); ++j) {
            global_cell_index_type basecell = neighs[i][j].basecell;
            assert(partition[basecell] == comm.rank());
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
    neighcomm = util::undirected_graph_communicator(comm_cart, neighbors);
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

void Diffusion::command(std::string s)
{
    std::smatch m;

    static const std::regex profit_re(
        "(set) (profit) (pass) (through) (([[:digit:]]*[.])?[[:digit:]]+)");
    if (std::regex_match(s, m, profit_re)) {
        double ppt = std::stod(m[5].str().c_str(), NULL);
#ifndef NDEBUG
        if (comm_cart.rank() == 0)
            std::cout << "Setting profit pass through = " << ppt << std::endl;
#endif
        profit_percentage_pass_through = ppt;
        return;
    }

    static const std::regex iter_re("(set) (flow_count) ([[:digit:]]+)");
    if (std::regex_match(s, m, iter_re)) {
        uint32_t flow_count = std::stoul(m[3].str().c_str(), NULL);

        const bool r = diff_variants::diffusion_maybe_set_nflow_iter(
            flow_calc.get(), flow_count);
        if (comm_cart.rank() == 0) {
            if (r) {
#ifndef NDEBUG
                std::cout << "Setting flow_count = " << flow_count << std::endl;
#endif
            }
            else {
                std::cerr << "Cannot set nflow iter."
                          << "Not supported by flow implementation."
                          << std::endl;
            }
        }
        return;
    }

    static const std::regex flow_re("(set) (flow) (.*)");
    if (std::regex_match(s, m, flow_re)) {
        const std::string &impl = m[3].str();
        try {
            flow_calc = diff_variants::create_flow_calc(
                supported_default_diffusion_variants.at(impl));
#ifndef NDEBUG
            if (comm_cart.rank() == 0)
                std::cout << "Setting implementation to: " << impl << std::endl;
#endif
        }
        catch (const std::out_of_range &) {
            if (comm_cart.rank() == 0) {
                std::cerr << "Cannot set implementation! Implementation \""
                          << impl << "\" is not found or not allowed.\n"
                          << "Note: \"so\" or \"sof\" are not supported by "
                          << "grid \"diff\". Use \"ps_diff\" instead.\n\n"
                          << "Setting default implementation: \"willebeek\"."
                          << std::endl;
            }
            assert(false);
        }
    }

    if (s == "synchronize") {
        // Make partition array globally known
        std::vector<rank_type> buf(partition.size());
        for (size_t i = 0; i < partition.size(); ++i)
            buf[i] = partition[i].value_or(-1);
        MPI_Allreduce(MPI_IN_PLACE, buf.data(), buf.size(),
                      boost::mpi::get_mpi_datatype<rank_type>(), MPI_MAX,
                      comm_cart);
        for (size_t i = 0; i < partition.size(); ++i)
            partition[i] = buf[i];
        stores_full_partitioning = true;
    }
}
} // namespace grids
} // namespace repa
