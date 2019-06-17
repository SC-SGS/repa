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
    UnknownGridTypeError() : w(std::string("Unknown grid type."))
    {
    }
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
GridType parse_grid_type(const std::string &desc);

/** Returns true if support for a certain grid type is compiled in.
 */
bool has_grid_type(GridType gt);

} // namespace repa
