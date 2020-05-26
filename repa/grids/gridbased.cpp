
/**
 * Copyright 2017-2020 Steffen Hirschmann
 * Copyright 2020 Benjamin Vier
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

#include "gridbased.hpp"

#include <algorithm>
#include <boost/mpi/collectives.hpp>
#include <map>
#include <regex>

#include "util/mpi_cart.hpp"
#include "util/mpi_cart_coloring.hpp"
#include "util/mpi_graph.hpp"
#include "util/push_back_unique.hpp"
#include "util/range.hpp"
#include "util/vdist.hpp"
#include "util/vec_arith.hpp"

#ifndef NDEBUG
#define GRID_DEBUG
#endif

namespace repa {
namespace grids {

std::array<Vec3d, 8> GridBasedGrid::bounding_box(rank_type r) const
{
    const Vec3i coord = util::mpi_cart_get_coords(comm_cart, r);
    const Vec3i dims = util::mpi_cart_get_dims(comm_cart);

    std::array<Vec3d, 8> result;
    size_t i = 0;
    // Ranks holding the bounding box grid points of "r" = (c0, c1, c2) are:
    // (c0,     c1,     c2) upper right back corner,
    // (c0 - 1, c1,     c2) upper left back corner,
    // (c0,     c1 - 1, c2) lower right back corner,
    // (c0,     c1,     c2 - 1) upper right front corner,
    // (c0 - 1, c1 - 1, c2) lower left back corner
    // ... 2 more ...
    // (c0 - 1, c1 - 1, c2 - 1) lower left front corner
    // In total the set: {c0, c0 - 1} x {c1, c1 - 1} x {c2, c2 - 1}
    Vec3i off;
    for (off[2] = 0; off[2] <= 1; ++off[2]) {
        for (off[1] = 0; off[1] <= 1; ++off[1]) {
            for (off[0] = 0; off[0] <= 1; ++off[0]) {
                using namespace util::vector_arithmetic;
                Vec3i nc = (coord - off) % dims;
                rank_type proc = util::mpi_cart_rank(comm_cart, nc);

                // Mirror the gridpoint back to where this subdomain is
                // expecting it.
                const Vec3i mirror
                    = -static_cast_vec<Vec3i>((coord == 0) && (off == 1));
                result[i] = gridpoints[proc] + mirror * box_l;
                i++;
            }
        }
    }
    return result;
}

rank_index_type GridBasedGrid::n_neighbors() const
{
    return rank_index_type{const_neighborhood.size()};
}

rank_type GridBasedGrid::neighbor_rank(rank_index_type i) const
{
    return const_neighborhood[i];
}

void GridBasedGrid::pre_init(bool firstcall)
{
    if (firstcall) {
        auto dims = util::mpi_cart_get_dims(comm_cart);
        if (dims[0] % 2 == 1 || dims[1] % 2 == 1 || dims[2] % 2 == 1) {
            if (comm_cart.rank() == 0)
                std::cerr
                    << "Warning: There is an odd number of processes in at "
                       "least one "
                       "dimension. "
                    << " The nodes on the domain boundary in these "
                       "dimensions are *not* "
                       "shifted. "
                    << " This *only* affects load-balancing quality, not "
                       "functionality ."
                    << std::endl;
        }
        init_regular_partitioning();
    }
}

void GridBasedGrid::post_init(bool firstcall)
{
}

void GridBasedGrid::init_regular_partitioning()
{
    is_regular_grid = true;

    using namespace util::vector_arithmetic;
    gridpoint = (node_pos + 1) * (box_l / node_grid);

    create_cartesian_neighborhood(); // Necessary only once. Communication
                                     // structure does not change
    init_octagons();
}

void GridBasedGrid::create_cartesian_neighborhood()
{
    const Vec3i coord = util::mpi_cart_get_coords(comm_cart);
    const Vec3i dims = util::mpi_cart_get_dims(comm_cart);

    std::vector<rank_type> source_neigh,
        dest_neigh; // Send and receive neighborhood for repart
    Vec3i off;
    for (off[0] = -1; off[0] <= 1; ++off[0]) {
        for (off[1] = -1; off[1] <= 1; ++off[1]) {
            for (off[2] = -1; off[2] <= 1; ++off[2]) {
                using namespace util::vector_arithmetic;
                const Vec3i nc = (coord + off) % dims;
                rank_type r = util::mpi_cart_rank(comm_cart, nc);

                // Insert "r" as a new neighbor if yet unseen.
                if (r == comm_cart.rank())
                    continue;

                util::push_back_unique(const_neighborhood, r);

                if (all(off >= 0))
                    util::push_back_unique(source_neigh, r);
                if (all(off <= 0))
                    util::push_back_unique(dest_neigh, r);
            }
        }
    }

    source_neigh.push_back(comm_cart.rank());
    dest_neigh.push_back(comm_cart.rank());
    neighcomm = util::directed_graph_communicator(comm_cart, source_neigh,
                                                  dest_neigh);
}

void GridBasedGrid::init_octagons()
{
    boost::mpi::all_gather(comm_cart, gridpoint, gridpoints);

    my_dom = util::tetra::Octagon(bounding_box(comm_cart.rank()));

    neighbor_doms.clear();
    neighbor_doms.reserve(n_neighbors());

    for (const auto nidx : util::range(n_neighbors())) {
        neighbor_doms.push_back(
            util::tetra::Octagon(bounding_box(neighbor_rank(nidx))));
    }
}

GridBasedGrid::GridBasedGrid(const boost::mpi::communicator &comm,
                             Vec3d box_size,
                             double min_cell_size,
                             ExtraParams ep)
    : GloMethod(comm, box_size, min_cell_size, ep),
      mu(1.0),
      subdomain_midpoint(ep.subdomain_midpoint
                             ? ep.subdomain_midpoint
                             : decltype(subdomain_midpoint){std::bind(
                                 &GridBasedGrid::get_subdomain_center, this)})
{
    util::tetra::precision = static_cast<int16_t>(10. / min_cell_size);
}

GridBasedGrid::~GridBasedGrid()
{
}

rank_type GridBasedGrid::cart_topology_position_to_rank(Vec3d pos) const
{
    // Cache the octagons for each process.
    // They become invalid after the first repartitioning.
    // DO NOT CALL this function then.
    assert(is_regular_grid);

    static std::map<rank_type, util::tetra::Octagon> all_octs;
    for (rank_type i = 0; i < comm_cart.size(); ++i) {
        if (all_octs.find(i) == std::end(all_octs))
            all_octs[i] = util::tetra::Octagon(bounding_box(i));

        if (all_octs[i].contains(pos))
            return i;
    }

    throw std::runtime_error(
        "Position globally unknown. This is a bug, please report it.");
}

rank_type GridBasedGrid::rank_of_cell(global_cell_index_type cellidx) const
{
    // Cell ownership is based on the cell midpoint.
    const auto mp = gbox.midpoint(cellidx);

    if (is_regular_grid)
        return cart_topology_position_to_rank(mp);

    // .contains() is mutually exclusive. The expectation is that most
    // queried positions belong to this node, so check it first. The order
    // of the neighbors is irrelevant.
    if (my_dom.contains(mp))
        return comm.rank();

    for (const auto i : util::range(rank_index_type{neighbor_doms.size()})) {
        if (neighbor_doms[i].contains(mp))
            return neighbor_rank(i);
    }

    return UNKNOWN_RANK;
}

Vec3d GridBasedGrid::get_subdomain_center()
{
    using namespace util::vector_arithmetic;
    Vec3d c{0., 0., 0.};

    for (const auto i : util::range(n_local_cells()))
        c += gbox.midpoint(cells[i]);
    return c / static_cast<double>(n_local_cells());
}

static Vec3d calc_shift(double local_load,
                        Vec3d subdomain_midpoint,
                        Vec3d cur_gridpoint,
                        MPI_Comm neighcomm)
{
    using namespace util::vector_arithmetic;

    // The node displacement is calculated according to
    // C. Begau, G. Sutmann, Comp. Phys. Comm. 190 (2015), p. 51 - 61
    const auto nneigh
        = rank_index_type{util::mpi_undirected_neighbor_count(neighcomm)};
    const double lambda_p = local_load;
    const auto r_p = subdomain_midpoint;

    std::vector<double> lambda(nneigh);
    MPI_Neighbor_allgather(&lambda_p, 1, MPI_DOUBLE, lambda.data(), 1,
                           MPI_DOUBLE, neighcomm);

    const double lnormalizer
        = std::accumulate(lambda.begin(), lambda.end(), 0.0) / nneigh;

    std::vector<double> lambda_hat(nneigh);
    for (int i = 0; i < nneigh; ++i)
        lambda_hat[i] = lambda[i] / lnormalizer;

    std::vector<Vec3d> r(nneigh);
    MPI_Neighbor_allgather(r_p.data(), sizeof(Vec3d), MPI_BYTE, r.data(),
                           sizeof(Vec3d), MPI_BYTE, neighcomm);

    for (const auto i : util::range(nneigh)) {
        // Form "u"
        r[i] -= cur_gridpoint;
        const double len = util::norm2(r[i]);
        // Form "f"
        r[i] *= (lambda_hat[i] - 1) / len;
    }

    // Note: We do not need to consider neighbors multiple times even
    // if two processes neighbor themselves along multiple boundaries.
    // We have a Cartesian grid. That means that if a process
    // appears twice in the neighborhood, all do.
    // So we can safely neglect multiple neighbors.

    Vec3d shift_vector = {0., 0., 0.};
    for (int i = 0; i < nneigh; i++) {
        shift_vector += r[i];
    }

    return shift_vector;
}

static Vec3d shift_gridpoint(Vec3d gp,
                             Vec3d shift_vector,
                             double factor,
                             const boost::mpi::communicator &comm_cart)
{
    const Vec3i coords = util::mpi_cart_get_coords(comm_cart);
    const Vec3i dims = util::mpi_cart_get_dims(comm_cart);

    // Note: Since we do not shift gridpoints over periodic boundaries,
    // f values from periodic neighbors are not considered.
    // Therefore, they do not need periodic mirroring.
    for (int d = 0; d < 3; ++d) {
        // Shift only non-boundary coordinates
        if (coords[d] == dims[d] - 1)
            continue;
        gp[d] += shift_vector[d];
    }

    return gp;
}

bool GridBasedGrid::sub_repartition(CellMetric m, CellCellMetric ccm)
{
    using namespace util::vector_arithmetic;

    const auto weights = m();
    assert(weights.size() == n_local_cells());

    const double lambda_p
        = std::accumulate(std::begin(weights), std::end(weights), 0.0);
    const auto r_p = this->subdomain_midpoint();

    auto shift_vector = calc_shift(lambda_p, r_p, gridpoint, neighcomm);

    // const Vec3i coords = util::mpi_cart_get_coords(comm_cart);
    const std::vector<rank_type> domains_to_check
        = util::mpi_directed_neighbors(neighcomm).first;

    // Colored shifting scheme to avoid multiple node conflicts at once,
    // according to:
    // C. Begau, G. Sutmann, Comp. Phys. Comm. 190 (2015), p. 51 - 61
    util::independent_process_sets(comm_cart)
        .for_each([&]() {
            double neighborhood_valid = false;
            for (double factor = 1.0; !neighborhood_valid && factor > .2;
                 factor /= 2.) {
                gridpoints[comm_cart.rank()] = shift_gridpoint(
                    gridpoint, shift_vector, mu * factor, comm_cart);
                neighborhood_valid
                    = check_validity_of_subdomains(domains_to_check);
            }
            // Restore old info in "gridpoints" vector or accept new gridpoint
            if (neighborhood_valid)
                gridpoint = gridpoints[comm_cart.rank()];
            else
                gridpoints[comm_cart.rank()] = gridpoint;
        })
        .for_all_after_each_round([&]() {
            // Update gridpoint and gridpoints
            // Currently allgather. Structly, only the changed gridpoints need
            // to be communicated
            gridpoints.clear();
            boost::mpi::all_gather(comm_cart, gridpoint, gridpoints);
            assert(gridpoints.size() == comm_cart.size());
        })();

    is_regular_grid = false;
    init_octagons(); // Necessary for position_to_rank queries
    return true;
}

bool GridBasedGrid::check_validity_of_subdomains(
    const std::vector<rank_type> &ranks) const
{
    const auto cs = cell_size();
    const auto max_cs = std::max(std::max(cs[0], cs[1]), cs[2]);

    return std::all_of(
        std::begin(ranks), std::end(ranks), [max_cs, this](rank_type r) {
            return util::tetra::Octagon(bounding_box(r), max_cs).is_valid();
        });
}

void GridBasedGrid::command(std::string s)
{
    static const std::regex mure("\\s*mu\\s*=\\s*(\\d+\\.|\\.\\d+|\\d+.\\d+)");
    std::smatch m;

    if (std::regex_match(s, m, mure)) {
        mu = std::strtod(m[1].str().c_str(), NULL);
        if (comm_cart.rank() == 0)
            std::cout << "Setting mu = " << mu << std::endl;
    }
}

} // namespace grids
} // namespace repa