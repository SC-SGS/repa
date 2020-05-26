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
#include <boost/serialization/array.hpp>
#include <boost/serialization/utility.hpp> // std::pair
#include <boost/serialization/vector.hpp>
#include <numeric>
#include <regex>

#include "util/fill.hpp"
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

/** Does a number of "set partition[i] = v for all i in i-vector".
 *
 * The different operations are given as vectors and the values and
 * i-vector-vector are passed as std::pair.
 */
template <typename PartitionEntryType>
void mark_new_owners_from_sendvolume(
    std::vector<PartitionEntryType> &partition,
    const std::pair<SendVolIndices, RankVector> &sendvolume)
{
    const auto &indicess = std::get<0>(sendvolume);
    const auto &new_values = std::get<1>(sendvolume);
    assert(indicess.size() == new_values.size());
    // Do *not* loop over new_values, "indicess" can be empty.
    for (size_t set_num = 0; set_num < indicess.size(); ++set_num) {
        const auto &indices = indicess[set_num];
        const auto &n_value = new_values[set_num];

        for (const auto &i : indices)
            partition[i] = n_value;
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

void Diffusion::clear_unknown_cell_ownership()
{
    auto is_my_cell = [this](global_cell_index_type neighcell) {
        return partition[neighcell] == comm_cart.rank();
    };

    fill_if_index(
        std::begin(partition), std::end(partition),
        util::ioptional<rank_type>{}, [this, is_my_cell](size_t glocellidx) {
            auto neighborhood
                = gbox.full_shell_neigh(global_cell_index_type{glocellidx});
            return std::none_of(std::begin(neighborhood),
                                std::end(neighborhood), is_my_cell);
        });
}

bool Diffusion::sub_repartition(CellMetric m, CellCellMetric ccm)
{
    const auto cellweights = m();
    assert(cellweights.size() == n_local_cells());

    clear_unknown_cell_ownership();

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

    // Update partition array
    _impl::mark_new_owners_from_sendvolume(partition, send_information);

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

        // Update the partition entry for all received cells.
        using namespace std::placeholders;
        std::for_each(std::begin(neighbor_sendinformation),
                      std::end(neighbor_sendinformation),
                      [this](const auto &neighbor_info) {
                          _impl::mark_new_owners_from_sendvolume(partition,
                                                                 neighbor_info);
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

        update_partitioning_from_received_neighbourhood(
            received_neighborhood_info);
    }
    assert(_impl::is_ghost_layer_fully_known(partition, comm_cart, gbox));

    return true;
}

Diffusion::Diffusion(const boost::mpi::communicator &comm,
                     Vec3d box_size,
                     double min_cell_size,
                     ExtraParams ep)
    : GloMethod(comm, box_size, min_cell_size, ep),
      flow_calc(diff_variants::create_flow_calc(
          diff_variants::FlowCalcKind::WILLEBEEK))
{
    // Initial partitioning
    partition.resize(gbox.ncells());
    util::make_initial_partitioner(gbox, comm_cart)
        .apply(initial_partitioning,
               [this](global_cell_index_type idx, rank_type r) {
                   assert(r >= 0 && r < this->comm.size());
                   this->partition[idx] = r;
               });
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
             gbox.full_shell_neigh_without_center(cells[borderCells[i]])) {
            if (partition[neighCell] == comm_cart.rank()
                && std::find(std::begin(borderCells), std::end(borderCells),
                             global_to_local.at(neighCell))
                       != std::end(borderCells)) {
                nadditional_comm++;
            }
        }
#ifdef DIFFUSION_DEBUG
        assert(nadditional_comm < 27);
#endif

        plist.emplace_back(27 - nadditional_comm, profit, borderCells[i]);
    }

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
                to_send[neighidx].push_back(cells[cidx]);
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
        if (comm_cart.rank() == 0)
            std::cout << "Setting profit pass through = " << ppt << std::endl;
        profit_percentage_pass_through = ppt;
        return;
    }

    static const std::regex iter_re("(set) (flow_count) ([[:digit:]]+)");
    if (std::regex_match(s, m, iter_re)) {
        uint32_t flow_count = std::stoul(m[3].str().c_str(), NULL);
        if (diff_variants::diffusion_maybe_set_nflow_iter(flow_calc.get(),
                                                          flow_count)
            && comm_cart.rank() == 0)
            std::cout << "Setting flow_count = " << flow_count << std::endl;
        else if (comm_cart.rank() == 0)
            std::cerr << "Cannot set nflow iter. Not supported by your "
                         "selected flow calculation."
                      << std::endl;
        return;
    }

    static const std::regex flow_re("(set) (flow) (.*)");
    if (std::regex_match(s, m, flow_re)) {
        const std::string &impl = m[3].str();
        try {
            flow_calc = diff_variants::create_flow_calc(
                supported_default_diffusion_variants.at(impl));
            if (comm_cart.rank() == 0)
                std::cout << "Setting implementation to: " << impl << std::endl;
        }
        catch (const std::out_of_range &) {
            if (comm_cart.rank() == 0) {
                std::cerr
                    << "Cannot set implementation! Implementation \"" << impl
                    << "\" is not found or not allowed to be used with the "
                       "default Diffusion. If you wanna use \"so\" or \"sof\" "
                       "use \"ps_diff\" instead. Now using default "
                       "implementation: \"willebeek\"."
                    << std::endl;
            }
            assert(false);
        }
    }
}
} // namespace grids
} // namespace repa
