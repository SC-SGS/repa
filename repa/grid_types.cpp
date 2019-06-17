
#include "grid_types.hpp"

namespace repa {

namespace {

static const std::unordered_map<std::string, GridType> grid_type_registry
    = {{"p4est", GridType::P4EST},
       {"cart", GridType::CART},
       {"graph", GridType::GRAPH},
       {"diff", GridType::DIFF},
       {"hybrid_gp_diff", GridType::HYB_GP_DIFF},
       {"kd_tree", GridType::KD_TREE},
       {"gridbased", GridType::GRIDBASED}};

#ifdef HAVE_KDPART
#define KDPART_AVAIL true
#else
#define KDPART_AVAIL false
#endif

#ifdef HAVE_P4EST
#define P4EST_AVAIL true
#else
#define P4EST_AVAIL false
#endif

#ifdef HAVE_PARMETIS
#define PARMETIS_AVAIL true
#else
#define PARMETIS_AVAIL false
#endif

#ifdef HAVE_TETRA
#define TETRA_AVAIL true
#else
#define TETRA_AVAIL false
#endif

// Awaiting C++20 which will finally have designated initializers -.-
static const std::unordered_map<GridType, bool> grid_type_availability
    = {{GridType::P4EST, P4EST_AVAIL},
       {GridType::CART, true},
       {GridType::GRAPH, PARMETIS_AVAIL},
       {GridType::DIFF, true},
       {GridType::HYB_GP_DIFF, PARMETIS_AVAIL},
       {GridType::KD_TREE, KDPART_AVAIL},
       {GridType::GRIDBASED, true}};

} // namespace

inline GridType parse_grid_type(const std::string &desc)
{
    try {
        return grid_type_registry.at(desc);
    }
    catch (const std::out_of_range &) {
        throw UnknownGridTypeError(desc);
    }
}

std::string grid_type_to_string(GridType gt)
{
    for (const auto &p : grid_type_registry) {
        if (gt == p.second)
            return p.first;
    }

    throw UnknownGridTypeError();
}

bool has_grid_type(GridType gt)
{
    try {
        return grid_type_availability.at(gt);
    }
    catch (const std::out_of_range &) {
        throw UnknownGridTypeError();
    }
}

std::set<GridType> supported_grid_types()
{
    std::set<GridType> gts;
    for (const auto &p : grid_type_availability) {
        if (p.second)
            gts.insert(p.first);
    }
    return gts;
}

} // namespace repa