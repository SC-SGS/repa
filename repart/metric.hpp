
#pragma once

#include <boost/mpi/communicator.hpp>
#include <functional>
#include <string>
#include <vector>

namespace repa {
namespace repart {

/** Represents a linear combination of single metric functions.
 */
struct Metric {
    using metric_func = std::function<void(std::vector<double> &)>;
    using cc_metric_func = std::function<double(int, int)>;

    Metric(const boost::mpi::communicator &comm_cart) : comm_cart(comm_cart)
    {
    }
    Metric(const boost::mpi::communicator &comm_cart, const std::string &desc)
        : Metric(comm_cart)
    {
        set_metric(desc);
    }

    /** Set "desc" as metric. Might throw a std::invalid_argument exception if
     * desc is not understood. Metric description strings are linear
     * combinations of single metrics. E.g. "2.0*ncells +1.7*nghostpart" The
     * space after the metric name ("ncell") is mandatory. Factor,
     * multiplication and addition sign are mandatory. Negative constants are
     * only allowed for the first factor. For further use subtraction instead of
     * addition, e.g. "-1.0*ncells -1.7*nghostpart". Single metric names are
     * also acceptable and interpreted as "1.0<name>". Valid metrics are:
     * ncells, npart, ndistpairs, nforcepairs, nbondedia, nghostcells,
     * nghostpart, runtime and rand. \param desc string to describe the metric
     */
    void set_metric(const std::string &desc);

    /** Returns cell weights.
     * \return vector of weights. Length: local_cells.n
     */
    std::vector<double> operator()() const;

    /** Returns true is this metric is set to calculate cell-cell weights
     * (e.g. transfer counts).
     */
    bool has_cell_cell_metric() const
    {
        return ccmdesc.size() > 0;
    }

    /** Returns cell-cell metric weights if has_cell_cell_metric() == true. Else
     * returns 1.0 always.
     */
    double cell_cell_weight(int i, int j) const;

    /** Returns the load accumulated over all cells for this process.
     */
    double curload() const;

    /** Returns the average of curload() over all processes. Needs to be called
     * by all processes.
     */
    double paverage() const;

    /** Returns the maximum of curload() over all processes. Needs to be called
     * by all processes.
     */
    double pmax() const;

    /** Returns the imbalance of curload() over all processes. Needs to be
     * called by all processes.
     */
    double pimbalance() const;

private:
    void parse_metric_desc(const std::string &desc);
    void parse_cell_metric_desc(const std::string &desc);
    void parse_cell_cell_metric_desc(const std::string &desc);

    const boost::mpi::communicator &comm_cart;
    std::vector<std::pair<double, metric_func>> mdesc;
    std::vector<std::pair<double, cc_metric_func>> ccmdesc;
};

} // namespace repart
} // namespace repa
