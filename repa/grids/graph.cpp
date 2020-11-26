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

#include "graph.hpp"
#include <algorithm>
#include <boost/mpi.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/datatype.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>
#include <mpi.h>
#include <parmetis.h>

#include "util/all_gatherv.hpp"
#include "util/initial_partitioning.hpp"
#include "util/push_back_unique.hpp"
#include "util/vector_coerce.hpp"

#define MPI_IDX_T boost::mpi::get_mpi_datatype(static_cast<idx_t>(0))

#ifndef NDEBUG
#define GRAPH_DEBUG
#endif

namespace repa {
namespace grids {

Graph::Graph(const boost::mpi::communicator &comm,
             Vec3d box_size,
             double min_cell_size,
             ExtraParams ep)
    : GloMethod(comm, box_size, min_cell_size, ep)
{
    // Initial partitioning
    partition.reserve(gbox.global_cells().size());
    boost::push_back(partition, util::StaticRankAssigner{initial_partitioning,
                                                         gbox, comm_cart}
                                    .partitioning());
    assert(partition.size() == gbox.global_cells().size());
}

Graph::~Graph()
{
}

constexpr inline double clamp(double min, double val, double max)
{
    if (val < min)
        return min;
    else if (val > max)
        return max;
    else
        return val;
}

namespace {
/** Searches in range [begin, end) for an element greater than the specified.
 * Behaves like a call so std::upper_bound. However, if called successively,
 * this class assumes that the elements searched for monotonously increase.
 */
template <typename It>
struct FwdSearchUpperBound {
    using T = typename It::value_type;
    FwdSearchUpperBound(It begin, It end) : _begin(begin), _end(end)
    {
    }

    // Returns the iterator to the first element greater than "element".
    // If no such exists, return "end".
    It upper_bound(const T &element)
    {
#ifndef NDEBUG
        // Check pre-condition
        if (!_firstcall) {
            assert(element >= _last_el);
            _last_el = element;
        }
#endif
        for (; _begin != _end; _begin++) {
            if (*_begin > element)
                return _begin;
        }
        return _end;
    }

private:
    It _begin, _end;

#ifndef NDEBUG
    T _last_el = T{0};
    bool _firstcall = true;
#endif
};

template <typename Cont>
FwdSearchUpperBound<typename Cont::const_iterator>
new_fwd_searcher(const Cont &c)
{
    return FwdSearchUpperBound<typename Cont::const_iterator>{std::cbegin(c),
                                                              std::cend(c)};
}
} // namespace

/*
 * Repartition.
 * Every node is responsible for a certain range of cells along the
 * linearization defined by gbox. Every process evaluates the weights for its
 * cells and sends them to the process that is responsible for the graph node
 * that corresponds to the cell.
 * Partitioning is performed in parallel via ParMETIS.
 */
bool Graph::sub_repartition(CellMetric m, CellCellMetric ccm)
{
    static constexpr idx_t w_fac = 100;
    const auto vertex_weights = m();
    assert(vertex_weights.size() == local_cells().size());

    std::vector<idx_t> vtxdist(comm_cart.size() + 1);

    // Calculate "vtxdist". vtxdist[i] = first cell index belonging to rank i
    {
        const idx_t nglocells = static_cast<idx_t>(gbox.global_cells().size());
        const double ncells_per_proc
            = static_cast<double>(nglocells) / comm_cart.size();

        // Vertex ranges per process
        for (rank_type i = 0; i < comm_cart.size(); ++i) {
            const idx_t first_index
                = static_cast<idx_t>(std::round(i * ncells_per_proc));
            vtxdist[i] = clamp(0, first_index, nglocells - 1);
        }
        vtxdist[comm_cart.size()] = nglocells;

#ifdef GRAPH_DEBUG
        assert(vtxdist.size() == static_cast<size_t>(comm_cart.size()) + 1);
        for (int i = 0; i < comm_cart.size(); ++i) {
            assert(0 <= vtxdist[i] && vtxdist[i] < nglocells);
        }
        assert(vtxdist[comm_cart.size()] == nglocells);
#endif
    }

    // Receive vertex and edge weights
    std::vector<rank_type> recvranks;
    // Determine ranks from which to receive
    // (Via old "partition" field and the range of graph nodes this
    // process is responsible for.)
    for (idx_t i = vtxdist[comm_cart.rank()]; i < vtxdist[comm_cart.rank() + 1];
         ++i)
        util::push_back_unique(recvranks, partition[i]);

#ifdef GRAPH_DEBUG
    for (int r : recvranks) {
        assert(0 <= r && r < comm_cart.size());
        assert(std::count(std::begin(recvranks), std::end(recvranks), r) == 1);
    }
    assert(recvranks.size() <= static_cast<size_t>(comm_cart.size()));
#endif

    // [0]: vertex weight
    // [1-26]: edge weights
    using Weights = std::array<idx_t, 27>;

    std::vector<boost::mpi::request> rreq;
    std::vector<std::vector<Weights>> their_weights(comm_cart.size());
    idx_t wsum = 0; // Check for possible overflow
    rreq.reserve(comm_cart.size());
    for (rank_type n : recvranks)
        rreq.push_back(comm_cart.irecv(n, 0, their_weights[n]));

    // Sending vertex weights
    std::vector<boost::mpi::request> sreq;
    std::vector<std::vector<Weights>> my_weights(comm_cart.size());

    auto gp_part = new_fwd_searcher(vtxdist);
    for (const auto i : cell_store.local_cells()) {
        // "Rank" is responsible for cell "gidx" / "i" (local)
        // during graph parititioning
        const global_cell_index_type gidx = cell_store.as_global_index(i);

        // The inverse of cell-rank-assignment.
        // Note that due to rounding and clamping this is not straightforward
        // a simple "gidx / n_cells_per_proc". Therefore, we search for it
        // in vtxdist. We don't use a binary search because of its logarithmic
        // cost but rather use a linear forward search because the inputs
        // ("gidx") are ordered/monotonously increasing.
        // This makes the whole process (all upper_bound-calls) O(N+P)
        // because "gp_part" traverses the vtxdist array concurrently to the
        // traversal of cell_store.local_cells().
        const auto it = gp_part.upper_bound(gidx);
        assert(it != vtxdist.end());
        assert(it != vtxdist.begin());
        const rank_type rank = std::distance(vtxdist.cbegin(), it) - 1;

        Weights w;
        w[0] = static_cast<idx_t>(vertex_weights[i]);

        // Edge weights must not be strictly positive. Add 1 and scale the true
        // weight. Also leave out the loop (i, i), i.e. start at n = 1.
        for (int n = 1; n < 27; ++n) {
            auto neigh = cell_neighbor_index(i, n);
            w[n] = static_cast<idx_t>(ccm(i, neigh)) * w_fac + 1;
#ifdef GRAPH_DEBUG
            assert(w[n] > 0);
            if (neigh.is<local_cell_index_type>()) {
                // Local symmetry -- only ensures local symmetry, however,
                // symmetry is also required for cross-boundary edges.
                assert(w[n]
                       == ccm(neigh.as<local_cell_index_type>(), i) * w_fac
                              + 1);
            }
#endif
        }

        my_weights[rank].push_back(std::move(w));

        // Catch total weight sum overlow
        wsum += std::accumulate(std::begin(w), std::end(w),
                                static_cast<idx_t>(0));
        ensure(wsum < std::numeric_limits<idx_t>::max() / comm_cart.size(),
               "Graph weights too large for chosen Metis IDX_TYPE_WIDTH");
    }

    sreq.reserve(comm_cart.size());
    for (rank_type i = 0; i < comm_cart.size(); ++i) {
        if (my_weights[i].size() > 0)
            sreq.push_back(comm_cart.isend(i, 0, my_weights[i]));
    }

    // Prepare graph

    idx_t nvtx = vtxdist[comm_cart.rank() + 1] - vtxdist[comm_cart.rank()];

    // Regular grid as graph
    std::vector<idx_t> xadj(nvtx + 1), adjncy(26 * nvtx);
    for (idx_t i = 0; i < nvtx; ++i) {
        xadj[i] = 26 * i;
        for (int n = 0; n < 26; ++n) {
            adjncy[26 * i + n] = gbox.neighbor(
                global_cell_index_type{vtxdist[comm_cart.rank()] + i}, n + 1);
        }
    }
    xadj[nvtx] = 26 * nvtx;

#ifdef GRAPH_DEBUG
    assert(nvtx == vtxdist[comm_cart.rank() + 1] - vtxdist[comm_cart.rank()]);
    assert(nvtx <= vtxdist.back()); // Last element is equal to the number of
                                    // global cells
    assert(xadj.size() == static_cast<size_t>(nvtx) + 1);
    for (int i = 0; i < nvtx; ++i) {
        assert(xadj[i] < xadj[i + 1]);
        assert(xadj[i + 1] - xadj[i] == 26);
    }

    for (size_t i = 0; i < adjncy.size(); ++i) {
        assert(adjncy[i] >= 0 && adjncy[i] < vtxdist.back());
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
        assert(ii[i] == their_weights[i].size());
    }
    assert(li == static_cast<size_t>(nvtx));

    for (idx_t w : vwgt) {
        assert(w != -1);
        assert(w >= 0 && w < 10000000);
    }
    for (idx_t w : adjwgt) {
        assert(w != -1);
        assert(w >= 0 && w < 10000000);
        // if (!m.has_cell_cell_metric())
        //    assert(w == 1.0 * w_fac + 1);
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
    std::vector<real_t> tpwgts(nparts, real_t{1.0} / nparts);
    // Imbalance tolerance
    real_t ubvec = real_t{1.05}; // Recommended by ParMETIS doc

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
    std::vector<idx_t> part(nvtx, static_cast<idx_t>(-1));

    auto metis_ret = ParMETIS_V3_PartKway(
        vtxdist.data(), xadj.data(), adjncy.data(), vwgt.data(), adjwgt.data(),
        &wgtflag, &numflag, &ncon, &nparts, tpwgts.data(), &ubvec, options,
        &edgecut, part.data(), &communicator);
    ensure(metis_ret == METIS_OK, "ParMETIS_V3_PartKway returned error.");

    // Copy idx_t to rank. Avoid copying if idx_t and rank are the same types.
    auto parti = util::coerce_vector_to<rank_type>(part);

#ifdef GRAPH_DEBUG
    assert(parti.size() == static_cast<size_t>(nvtx));
    for (auto r : parti) {
        assert(r != static_cast<idx_t>(-1));
        assert(0 <= r && r < comm_cart.size());
    }
#endif

#ifdef GRAPH_DEBUG
    std::fill(std::begin(partition), std::end(partition),
              static_cast<rank_type>(-1));
#endif

    util::all_gatherv_displ(comm_cart, parti.cref(), vtxdist, partition);

#ifdef GRAPH_DEBUG
    assert(partition.size() == gbox.global_cells().size());
    for (int r : partition) {
        assert(0 <= r && r < comm_cart.size());
    }
#endif

#ifdef GRAPH_DEBUG
    size_t nlc = std::count(std::begin(partition), std::end(partition),
                            comm_cart.rank());
    std::vector<decltype(nlc)> nlcs;
    boost::mpi::all_gather(comm, nlc, nlcs);

    if (comm_cart.rank() == 0) {
        std::cout << "Partitioning result:";
        std::copy(std::begin(nlcs), std::end(nlcs),
                  std::ostream_iterator<int>(std::cout, " "));
        std::cout << std::endl;
    }

    assert(std::accumulate(std::begin(nlcs), std::end(nlcs), size_t{0})
           == gbox.global_cells().size());
#endif

    return true;
}

} // namespace grids
} // namespace repa
