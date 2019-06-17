#pragma once

//#ifdef HAVE_METIS

#include "diffusion.hpp"
#include "globox.hpp"
#include "graph.hpp"
#include "pargrid.hpp"
#include <array>

namespace repa {
namespace grids {

struct HybridGPDiff : public ParallelLCGrid {
    HybridGPDiff(const boost::mpi::communicator &comm,
                 Vec3d box_size,
                 double min_cell_size);
    lidx n_local_cells() override;
    gidx n_ghost_cells() override;
    nidx n_neighbors() override;
    rank neighbor_rank(nidx i) override;
    Vec3d cell_size() override;
    Vec3i grid_size() override;
    lgidx cell_neighbor_index(lidx cellidx, int neigh) override;
    std::vector<GhostExchangeDesc> get_boundary_info() override;
    lidx position_to_cell_index(double pos[3]) override;
    rank position_to_rank(double pos[3]) override;
    nidx position_to_neighidx(double pos[3]) override;
    bool repartition(const repart::Metric &m,
                     std::function<void()> exchange_start_callback) override;

    void command(std::string s) override;

private:
    /** Underlying implementations
     * Note that currently the "partition" vector is copied between the two
     * sub-partitionert is the partitioning method is changed.
     * Both partitioners could be made to operate on a joint "partition" vector,
     * however currently they both can(!) use different value types.
     */
    Diffusion diff_impl;
    Graph graph_impl;

    enum class State { DIFF, GRAPH };

    /** Stores the state of the partitioner for switching purpose */
    State state;
    /** Stores if the state should be switched before the next repartition call
     */
    State switch_to_state;
    /** Reference to the implementation that is currently in use. */
    ParallelLCGrid *active_implementation;

    /** Switches between graph partitioning and diffusion. Activates
     * the partitioner that is currently not active.
     */
    void switch_implementation();
};
} // namespace grids
} // namespace repa

//#endif // HAVE_METIS
