
#include "metric.hpp"
#include "cells.hpp"
#include "domain_decomposition.hpp"
#include <algorithm>
#include <boost/iterator/iterator_facade.hpp>
#include <cctype>
#include <chrono>
#include <stdexcept>
#include <random>
#include <regex>
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "short_range_loop.hpp"

// Fills weights with a constant.
static void metric_ncells(std::vector<double> &weights) {
  std::fill(weights.begin(), weights.end(), 1.0);
}

// Fills weights with the number of particles per cell.
static void metric_npart(std::vector<double> &weights) {
  std::transform(local_cells.cell, local_cells.cell + local_cells.n,
                 weights.begin(), [](const Cell *c) { return c->n; });
}

int cell_ndistpairs(Cell *c) {
  int nnp =
      std::accumulate(c->m_neighbors.red().begin(), c->m_neighbors.red().begin(), 0,
                      [](int acc, const Cell *neigh) { return acc + neigh->n; });
  return c->n * nnp;
}

// Fills weights with the number of distance pairs per cell.
static void metric_ndistpairs(std::vector<double> &weights) {
  std::transform(local_cells.begin(), local_cells.end(), weights.begin(),
                 cell_ndistpairs);
}

static void metric_nforcepairs(std::vector<double> &weights) {
#ifdef LENNARD_JONES
  // Ugly hack: Use pairkernel to advance the current cell number
  int cellno = 0;
  for(;cellno < (local_cells.n - 1) && local_cells.cell[cellno]->n == 0; cellno++);

  auto particlekernel = [&cellno](Particle &p) {
    const Cell *c = local_cells.cell[cellno];
    // At the end of a cell? -> Increment cellno consistently to
    // short_range_loop. I.e. next non-empty cell
    if (p.p.identity == c->part[c->n - 1].p.identity) {
      if (cellno < (local_cells.n - 1))
        cellno++;
      // Cell "cellno" might be empty, skip ahead if so
      while (cellno < (local_cells.n - 1) && local_cells.cell[cellno]->n == 0)
        cellno++;
    }
  };

  auto pairkernel = [&weights, &cellno](Particle &p1, Particle &p2,
                                       Distance &d) {
#ifdef EXCLUSIONS
    if (do_nonbonded(&p1, &p2))
#endif
    {
      IA_parameters *ia_params = get_ia_param(p1.p.type, p2.p.type);
      double dist = std::sqrt(d.dist2);
      if (dist < ia_params->LJ_cut + ia_params->LJ_offset &&
          dist > ia_params->LJ_min + ia_params->LJ_offset)
        weights[cellno]++;
    }
  };

  short_range_loop(particlekernel, pairkernel);
#else
  fprintf(stderr, "Error: metric nforcepairs only available with LENNARD_JONES feature enabled.\n");
  errexit();
#endif
}

int cell_nbondedia(Cell *cell) {
  int nbondedia = 0;
  for (Particle *p = cell->part; p < cell->part + cell->n; p++) {
    for (int i = 0; i < p->bl.n;) {
      int type_num = p->bl.e[i++];
      Bonded_ia_parameters *iaparams = &bonded_ia_params[type_num];
      // int type = iaparams->type;

      // This could be incremented conditionally if "type" has a specific value
      // to only count bonded_ia of a certain type.
      nbondedia++;
      i += iaparams->num; // Skip the all partner particle numbers
    }
  }
  return nbondedia;
}

static void metric_nbondedia(std::vector<double> &weights) {
  std::transform(local_cells.begin(), local_cells.end(), weights.begin(),
                 cell_nbondedia);
}

static void metric_nghostcells(std::vector<double> &weights) {
  // Reminder: Weights is a vector over local cells.
  // We simply count the number of ghost cells and
  // add this cost in equal parts to all local cells
  double nghostfrac = static_cast<double>(ghost_cells.n) / local_cells.n;
  std::fill(weights.begin(), weights.end(), nghostfrac);
}

static void metric_nghostpart(std::vector<double> &weights) {
  // Reminder: Weights is a vector over local cells.
  // We simply count the number of ghost particles and
  // add this cost in equal parts to all local cells
  int nghostpart =
      std::accumulate(ghost_cells.begin(), ghost_cells.end(), 0,
                      [](int acc, const Cell *c) { return acc + c->n; });
  double nghostfrac = static_cast<double>(nghostpart) / local_cells.n;
  std::fill(weights.begin(), weights.end(), nghostfrac);
}

static void metric_runtime(std::vector<double> &weights) {
  throw std::runtime_error("Currently not implemented.");
}

// Generator for random integers
struct Randintgen {
  Randintgen()
      : mt(std::chrono::high_resolution_clock::now()
               .time_since_epoch()
               .count()),
        d(1, 1000) {}
  Randintgen(Randintgen &&other)
      : mt(std::move(other.mt)), d(std::move(other.d)) {}
  int operator()() { return d(mt); }

private:
  std::mt19937 mt;
  std::uniform_int_distribution<int> d;
};

// Fills weights with random integers
static void metric_rand(std::vector<double> &weights) {
  std::generate(weights.begin(), weights.end(), Randintgen());
}

static double cc_metric_uniform(int i, int j) {
  (void) i;
  (void) j;
  return 1.0;
}

static double cc_metric_multiply_npart(int i, int j) {
  return cells[i].n * cells[j].n;
}

static double cc_metric_add_npart(int i, int j) {
  return cells[i].n + cells[j].n;
}

static double cc_metric_forcepairs(int i, int j) {
#ifdef LENNARD_JONES
  auto& c1 = cells[i];
  auto& c2 = cells[j];

  int nfp = 0;
  double r[3];

  for (int p1 = 0; p1 < c1.n; ++p1) {
    for (int p2 = 0; p2 < c2.n; ++p2) {
      IA_parameters *ia_params = get_ia_param(c1.part[p1].p.type, c2.part[p2].p.type);
      get_mi_vector(r, c1.part[p1].r.p, c2.part[p2].r.p);
      double dist = std::sqrt(sqrlen(r));
      if (dist < ia_params->LJ_cut + ia_params->LJ_offset &&
          dist > ia_params->LJ_min + ia_params->LJ_offset)
        nfp++;
    }
  }
  return static_cast<double>(nfp);
#else
  fprintf(stderr, "Error: metric nforcepairs only available with LENNARD_JONES feature enabled.\n");
  errexit();
  // not reached.
  return 0;
#endif
}

namespace generic_dd {
namespace repart {

// Get the appropriate metric function described in string "desc".
// Throws a std::invalid_argument exception if no such metric is available
static Metric::metric_func get_single_metric_func(const std::string &desc) {
  using desc_func_pair = std::pair<std::string, Metric::metric_func>;
  static const std::vector<desc_func_pair> mets = {
      {"ncells", metric_ncells},
      {"npart", metric_npart},
      {"ndistpairs", metric_ndistpairs},
      {"nforcepairs", metric_nforcepairs},
      {"nbondedia", metric_nbondedia},
      {"nghostcells", metric_nghostcells},
      {"nghostpart", metric_nghostpart},
      {"runtime", metric_runtime},
      {"rand", metric_rand}};

  for (const auto &t : mets) {
    if (desc == std::get<0>(t)) {
      return std::get<1>(t);
    }
  }

  throw std::invalid_argument(std::string("No such metric available: ") + desc);
}

static Metric::cc_metric_func get_single_cc_metric_func(const std::string &desc) {
  using desc_func_pair = std::pair<std::string, Metric::cc_metric_func>;
  static const std::vector<desc_func_pair> mets = {
      {"uniform", cc_metric_uniform},
      {"add_npart", cc_metric_add_npart},
      {"multiply_npart", cc_metric_multiply_npart},
      {"forcepairs", cc_metric_forcepairs}};
  for (const auto &t : mets) {
    if (desc == std::get<0>(t)) {
      return std::get<1>(t);
    }
  }

  throw std::invalid_argument(std::string("No such metric available: ") + desc);
}

void Metric::set_metric(const std::string &desc) {
  mdesc.clear();
  ccmdesc.clear();
  parse_metric_desc(desc);
}

static std::string remove_whitespace(const std::string &s) {
  static const std::regex ws("\\s");
  return std::regex_replace(s, ws, "");
}

void Metric::parse_metric_desc(const std::string &desc) {
  std::string::size_type pos;
  // Cell cell metric denoted in square brackets at the end.
  if ((pos = desc.find('[')) != std::string::npos) {
    auto cdesc = desc.substr(0, pos - 1);
    auto ccdesc = desc.substr(pos + 1, desc.find(']') - pos - 1);
    //std::cout << "Parsing '" << cdesc << "' as metric and '" << ccdesc << "' as cell cell metric" << std::endl;
    parse_cell_metric_desc(cdesc);
    parse_cell_cell_metric_desc(ccdesc);
  } else {
    parse_cell_metric_desc(desc);
  }
}

// Parses a string representing a linear combination
template <typename MetricVec, typename StringParseFunc>
static void parse_linear_combination(const std::string &desc,
                                     MetricVec& mvec,
                                     StringParseFunc parse_metric) {
  static const std::regex term_re(
      "\\s*([\\+-]?\\s*\\d*(\\.\\d*)?)?\\s*\\*?\\s*(\\w+)");
  static const auto match_to_term = [parse_metric](const std::smatch &match) {
    // There is whitespace in the numbers (first group in term_re), so we have
    // to remove it before calling atof.
    std::string num = remove_whitespace(match.str(1));
    // Number might be omitted, i.e. "foo+bar". In this case 'num' only holds a
    // sign or nothing at all (for the very first term).
    if (num == "+" || num == "-" || num == "")
      num += "1.0";
    try {
      return std::make_pair(stod(num), parse_metric(match.str(3)));
    } catch (...) {
      std::cout << "Error in metric parsing at: " << match.str(0) << std::endl;
      throw;
    }
  };

  auto t_begin =
      std::sregex_iterator(std::begin(desc), std::end(desc), term_re);
  auto t_end = decltype(t_begin){};

  std::transform(t_begin, t_end, std::back_inserter(mvec), match_to_term);
}

void Metric::parse_cell_cell_metric_desc(const std::string &desc) {
  parse_linear_combination(desc, ccmdesc, get_single_cc_metric_func);
}

void Metric::parse_cell_metric_desc(const std::string &desc) {
  parse_linear_combination(desc, mdesc, get_single_metric_func);
}

std::vector<double> Metric::operator()() const {
  std::vector<double> w(local_cells.n, 0.0), tmp(local_cells.n);

  for (const auto &t : mdesc) {
    double factor = std::get<0>(t);
    metric_func func = std::get<1>(t);

    func(tmp);
    // Multiply add to w.
    std::transform(
        w.begin(), w.end(), tmp.begin(), w.begin(),
        [factor](double acc, double val) { return acc + factor * val; });
  }
  return w;
}

double Metric::curload() const {
  auto ws = (*this)();
  return std::accumulate(std::begin(ws), std::end(ws), 0.0);
}

double Metric::paverage() const {
  double l = curload();
  double tot;
  MPI_Allreduce(&l, &tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return tot / n_nodes;
}

double Metric::pmax() const {
  double l = curload();
  double max;
  MPI_Allreduce(&l, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return max;
}

double Metric::pimbalance() const {
  double l = curload();

  MPI_Request req[2];
  double tot, max;

  MPI_Iallreduce(&l, &tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &req[0]);
  MPI_Iallreduce(&l, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD, &req[1]);
  MPI_Waitall(2, req, MPI_STATUS_IGNORE);

  double avg = tot / n_nodes;
  return max / avg;
}

double Metric::cell_cell_weight(int i, int j) const {
  if (!has_cell_cell_metric()) {
    return 1.0;
  }

  double res = 0.0;
  for (const auto &t : ccmdesc) {
    double factor = std::get<0>(t);
    cc_metric_func func = std::get<1>(t);
    res += factor * func(i, j);
  }

  return res;
}

}
}
