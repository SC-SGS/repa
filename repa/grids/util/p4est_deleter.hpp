/**
 * Copyright 2017-2019 Steffen Hirschmann
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

#pragma once

#include <memory>
#include <p8est.h>
#include <p8est_ghost.h>
#include <p8est_mesh.h>

namespace std {
template <>
struct default_delete<p8est_t> {
    void operator()(p8est_t *p) const
    {
        if (p != nullptr)
            p8est_destroy(p);
    }
};
template <>
struct default_delete<p8est_ghost_t> {
    void operator()(p8est_ghost_t *p) const
    {
        if (p != nullptr)
            p8est_ghost_destroy(p);
    }
};
template <>
struct default_delete<p8est_mesh_t> {
    void operator()(p8est_mesh_t *p) const
    {
        if (p != nullptr)
            p8est_mesh_destroy(p);
    }
};
template <>
struct default_delete<p8est_connectivity_t> {
    void operator()(p8est_connectivity_t *p) const
    {
        if (p != nullptr)
            p8est_connectivity_destroy(p);
    }
};
template <>
struct default_delete<sc_array_t> {
    void operator()(sc_array_t *p) const
    {
        if (p != nullptr)
            sc_array_destroy(p);
    }
};
} // namespace std