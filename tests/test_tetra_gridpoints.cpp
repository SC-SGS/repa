
#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tetra
#include <boost/test/unit_test.hpp>

#include <array>
#include <iostream>
#include <vector>

#include <repa/grids/util/tetra.hpp>

using repa::Vec3d;
using namespace repa::util;

struct PointArray {
    static const int size = 3;
    std::array<std::array<std::array<Vec3d, size>, size>, size> point = {{{}}};

    std::array<Vec3d, 8> getVerticesAtPosition(int id) const;
    PointArray();
};

static inline int extract_bit(int value, int bitno)
{
    return !!(value & (1 << bitno));
}

PointArray::PointArray()
{
    std::array<Vec3d, 8> gridpoints{{{35.567, 38.4418, 45.0975},
                                     {35.7365, 38.3409, 75.8135},
                                     {35.9562, 82.1775, 44.5612},
                                     {35.9879, 82.2762, 76.3616},
                                     {83.8035, 38.3232, 45.2142},
                                     {83.641, 38.3249, 75.6815},
                                     {83.4199, 82.2701, 44.6563},
                                     {83.3893, 82.2873, 76.2623}}};

    repa::Vec3<bool> minus{true, true, true};
    for (int x = 0; x < size; x++) {
        minus[0] = x == 0;
        int use_x = std::abs(x - 1) * 4;
        for (int y = 0; y < size; y++) {
            minus[1] = y == 0;
            int use_y = std::abs(y - 1) * 2;
            for (int z = 0; z < size; z++) {
                minus[2] = z == 0;
                int use_z = std::abs(z - 1);

                point[x][y][z] = gridpoints[use_x + use_y + use_z];
                for (int d = 0; d < 3; d++) {
                    if (minus[d]) {
                        point[x][y][z][d] -= 80.;
                    }
                }
            }
        }
    }
}

std::array<Vec3d, 8> PointArray::getVerticesAtPosition(int id) const
{
    const int x = extract_bit(id, 0);
    const int y = extract_bit(id, 1);
    const int z = extract_bit(id, 2);
    return {
        point[1 + x][1 + y][1 + z], point[0 + x][1 + y][1 + z],
        point[1 + x][0 + y][1 + z], point[0 + x][0 + y][1 + z],
        point[1 + x][1 + y][0 + z], point[0 + x][1 + y][0 + z],
        point[1 + x][0 + y][0 + z], point[0 + x][0 + y][0 + z],
    };
}

BOOST_AUTO_TEST_CASE(gridpoints_1)
{
    tetra::init_tetra(1, {80., 80., 80.});

    auto points = PointArray{};
    std::array<tetra::Octagon, 8> octas;
    for (int i = 0; i < 8; i++) {
        octas[i] = tetra::Octagon(points.getVerticesAtPosition(i), 0.1);
        BOOST_CHECK(octas[i].is_valid());
    }

    for (int x = 0; x < 80; x++) {
        for (int y = 0; y < 80; y++) {
            for (int z = 0; z < 80; z++) {
                const Vec3d p = {double(x), double(y), double(z)};
                auto n_accept
                    = std::count_if(std::begin(octas), std::end(octas),
                                    [p](const tetra::Octagon &octa) {
                                        return octa.contains(p);
                                    });
                BOOST_CHECK(n_accept == 1);
            }
        }
    }
}

int main(int argc, char **argv)
{
    return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
