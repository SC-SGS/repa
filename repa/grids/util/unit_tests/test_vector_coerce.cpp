/**
 * Copyright 2017-2020 Steffen Hirschmann
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

#include <algorithm>
#include <doctest/doctest.h>

#include "../vector_coerce.hpp"

TEST_CASE("vector_coerce different types")
{
    std::vector<int> v{1, 2, 3, 4};

    const auto w = repa::util::coerce_vector_to<unsigned>(v);

    CHECK(v.size() == w.size());
    CHECK(std::equal(v.begin(), v.end(), w.begin()));
    CHECK(static_cast<const void *>(w.cref().data())
          != static_cast<void *>(v.data()));
}

TEST_CASE("vector_coerce same types")
{
    std::vector<int> v{1, 2, 3, 4};

    const auto w = repa::util::coerce_vector_to<int>(v);

    CHECK(v.size() == w.size());
    CHECK(std::equal(v.begin(), v.end(), w.begin()));
    CHECK(w.cref().data() == v.data());
}
