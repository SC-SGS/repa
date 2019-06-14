
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <cstdio>

#include "pargrid.hpp"
#include "pargrid_factory.hpp"

int main()
{
    boost::mpi::environment env;
    boost::mpi::communicator comm;

    repa::Vec3d box = {{10., 10., 10.}};
    double min_cell_size = 1.;

    auto pg = repa::grids::make_pargrid(repa::GridType::CART, comm, box,
                                        min_cell_size);
    auto cs = pg->cell_size();

    std::printf("%lf %lf %lf\n", cs[0], cs[1], cs[2]);
}