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

/** Assigns a new given value to the elements in the index range
 * [first, last) if index satisfies a predicate.
 * Here, "first" corresponds to index 0 and "last-1" to index (last-first-1).
 */
template <typename FwdIt, typename T, typename Pred>
static void fill_if_index(FwdIt first, FwdIt last, const T &val, Pred p)
{
    for (size_t i = 0; first != last; ++first) {
        if (*first != val && p(i))
            *first = val;
        ++i;
    }
}

/** Sets "data[i]" to "val" for all i in the half-open interval
 * ["first_index", "last_index").
 */
template <typename T, typename FwdIt>
static void fill_index_range(std::vector<T> &data,
                             FwdIt first_index,
                             FwdIt last_index,
                             const T &val)
{
    using idx_type = typename FwdIt::value_type;
    std::for_each(first_index, last_index, [&data, &val](idx_type i) {
#ifdef NDEBUG
        data[i] = val;
#else
        data.at(i) = val;
#endif
    });
}