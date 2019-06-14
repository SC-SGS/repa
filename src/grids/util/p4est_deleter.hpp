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