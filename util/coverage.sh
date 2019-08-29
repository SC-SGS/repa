#!/bin/sh

# Copyright 2017-2019 Steffen Hirschmann
#
# This file is part of Repa.
#
# Repa is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Repa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Repa.  If not, see <https://www.gnu.org/licenses/>.
#

#
# This script runs lcov and genhtml.
# ONLY call this script from a build directory.
#

eprint() {
    printf "Error in coverage.sh: %s\\n" "$@" >&2
    exit 1
}

dirflags() {
    printf " -d %s" $(dirname $(find . -iname "*.gcda") | sort | uniq)
}

[ -z "`which lcov`" ] && eprint "Cound not find lcov."
[ -z "`which genhtml`" ] && eprint "Cound not find genhtml."

[ -f "rules.ninja" ] || eprint "Call this script only from build directories."

[ -e "cov" ] && eprint "Please remove the old 'cov' directory first."

# Find all directories with gcda files and prepend "-d"

lcov -z $(dirflags)
ninja test
# Note: The following dirflags might be different, if it is a new build.
lcov -t "repa_coverage" -o all_cov.info -c $(dirflags)

# Select all entries with the current working directory (strip "/build")
lcov -e all_cov.info "$(dirname $(pwd))/*" -o all_cov.info

genhtml all_cov.info -o cov
