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

#include "psdiffusion.hpp"

#include <algorithm>
#include <boost/mpi/collectives.hpp>
#include <regex>
#include <set>

#include "grids/util/mpi_cart.hpp"
#include "grids/util/mpi_graph.hpp"
#include "util/get_keys.hpp"
#include "util/mpi_cart.hpp"
#include "util/set_union.hpp"
#include "util/vec_arith.hpp"

#ifndef NDEBUG
#define PSDIFFUSION_DEBUG
#endif

namespace repa {
namespace grids {

static const std::unordered_map<std::string, diff_variants::FlowCalcKind>
    supported_ps_diffusion_variants
    = {{"so", diff_variants::FlowCalcKind::SO}};

/** Returns the new Coords of a Rank after it beeing maped to opposite site.
 * This is a necessary step for the metric because of the periodic edge.
 */
static Vec3i map_coords_to_opposite_side(const Vec3i &c0,
                                         const Vec3i &c2,
                                         const Vec3i &gsize)
{
    using namespace util::vector_arithmetic;
    Vec3i te, ts;
    for (int i = 0; i < te.size(); i++) {
        te[i] = c2[i] - c0[i] > 1;
        ts[i] = c0[i] - c2[i] > 1;
    }
    return c2 + ts * gsize - te * gsize;
}

/*
 * Initialization
 */
PSDiffusion::PSDiffusion(const boost::mpi::communicator &comm,
                         Vec3d box_size,
                         double min_cell_size,
                         ExtraParams ep)
    : Diffusion(comm, box_size, min_cell_size, ep),
      init_topology_comm(
          util::make_init_part_communicator(comm, initial_partitioning))
{
    ensure(initial_partitioning != util::InitialPartitionType::LINEAR,
           "PSDiffusion does not support initial linear partitioning."
           " Please use Cartesian1D or Cartesian3D");
}

PSDiffusion::~PSDiffusion()
{
}

std::set<std::string> PSDiffusion::get_supported_variants() const
{
    return util::set_union(Diffusion::get_supported_variants(),
                           util::get_keys(supported_ps_diffusion_variants));
}

void PSDiffusion::post_init(bool firstcall)
{
    Diffusion::post_init(firstcall);

#ifdef PSDIFFUSION_DEBUG
    // Copy initial neighborhood, so that the neighborhood can be checked for
    // consistency in later iterations.
    if (firstcall) {
        const int nneigh = neighbors.size();
        const int value
            = boost::mpi::all_reduce(comm_cart, nneigh, std::bit_or<void>{});
        ensure(nneigh == value,
               "Not all processes have the same amount of neighbors.");

        initial_neighborhood
            = std::vector<rank_type>(neighbors.begin(), neighbors.end());
    }

    // Neighborhood consistency check
    std::set<rank_type> a(initial_neighborhood.begin(),
                          initial_neighborhood.end());
    std::set<rank_type> b(neighbors.begin(), neighbors.end());

    assert(std::includes(a.begin(), a.end(), b.begin(), b.end()));
#endif
}

bool PSDiffusion::accept_transfer(local_cell_index_type cidx,
                                  rank_type neighrank) const
{
    return coords_based_allow_sending(cidx, neighrank);
}

bool PSDiffusion::coords_based_allow_sending(local_cell_index_type c,
                                             rank_type neighrank) const
{
    const Vec3i comm_dims = util::mpi_cart_get_dims(init_topology_comm);
    const Vec3i c0 = util::mpi_cart_get_coords(init_topology_comm);
    const Vec3i cn = map_coords_to_opposite_side(
        c0, util::mpi_cart_get_coords(init_topology_comm, neighrank),
        comm_dims);
    // Determine any of the processes owning neighboring cells of "c" would
    // be new neighbors to "neighrank".
    for (const global_cell_index_type &d :
         gbox.full_shell_neigh_without_center(cells[c])) {
        rank_type rank_d = rank_of_cell(d);
        if (rank_d == rank_of_cell(cells[c]) || rank_d == neighrank)
            continue;
        Vec3i c2 = map_coords_to_opposite_side(
            c0, util::mpi_cart_get_coords(init_topology_comm, rank_d),
            comm_dims);

        if (std::abs(cn[0] - c2[0]) >= 2 || std::abs(cn[1] - c2[1]) >= 2
            || std::abs(cn[2] - c2[2]) >= 2)
            return false;
    }
    return true;
}

void PSDiffusion::command(std::string s)
{
    std::smatch m;
    static const std::regex beta_re(
        "(set) (beta) (([[:digit:]]*[.])?[[:digit:]]+)");
    if (std::regex_match(s, m, beta_re)) {
        double beta_value = std::stod(m[3].str().c_str(), NULL);
        if (diff_variants::diffusion_maybe_set_beta(flow_calc.get(), beta_value)
            && comm_cart.rank() == 0)
            std::cout << "Setting beta = " << beta_value << std::endl;
        else if (comm_cart.rank() == 0)
            std::cerr << "Cannot set beta value. Not supported by your "
                         "selected flow calculation."
                      << std::endl;
        return;
    }

    static const std::regex flow_re("(set) (flow) (.*)");
    if (std::regex_match(s, m, flow_re)) {
        const std::string &impl = m[3].str();
        try {
            flow_calc = diff_variants::create_flow_calc(
                supported_ps_diffusion_variants.at(impl));
            if (comm_cart.rank() == 0)
                std::cout << "Setting implementation to: " << impl << std::endl;
            return;
        }
        catch (const std::out_of_range &) {
        }
    }

    Diffusion::command(s);
}

} // namespace grids
} // namespace repa
