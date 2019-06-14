#pragma once

#include <string>
#include <unordered_map>

namespace repa {

/** Enum of supported grid types.
 * Caution: Might not all be compiled in.
 */
enum class GridType {
    NONE,
    P4EST,
    CART,
    GRAPH,
    DIFF,
    KD_TREE,
    HYB_GP_DIFF,
    GRIDBASED
};

struct UnknownGridTypeError {
    UnknownGridTypeError(std::string s)
        : w(std::string("Unknown grid type: `") + s + std::string("'"))
    {
    }
    virtual const char *what() const
    {
        return w.c_str();
    }

private:
    std::string w;
};

/** Returns the GridType associated with a descriptive string for the grid
 * type.
 */
inline GridType parse_grid_type(const std::string &desc)
{
    static const std::unordered_map<std::string, GridType> grid_type_registry
        = {{"p4est", GridType::P4EST},
           {"cart", GridType::CART},
           {"graph", GridType::GRAPH},
           {"diff", GridType::DIFF},
           {"hybrid_gp_diff", GridType::HYB_GP_DIFF},
           {"kd_tree", GridType::KD_TREE},
           {"gridbased", GridType::GRIDBASED}};

    try {
        return grid_type_registry.at(desc);
    }
    catch (const std::out_of_range &) {
        throw UnknownGridTypeError(desc);
    }
}

} // namespace repa
