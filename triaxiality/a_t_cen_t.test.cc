#include "catch2/catch.hpp"
#include "a_cen_t.hh"
#include "t_cen_t.hh"

using y3_cluster::A_CEN_t;
using y3_cluster::T_CEN_t;

TEST_CASE("a_cen_t and t_cen_t works")
{
  //Unit test with the files: DS_M_2e+14_4e+14_z_0_0.34.dat and DS_M_2e+14_4e+14_z_0_0.34_cosi_0.8_1.dat
  double rp[12] = {0.125893, 0.199526, 0.316228, 0.501187, 0.794328, 1.25893, 1.99526, 3.16228, 5.01187, 7.94328, 12.5893, 19.9526}; //units of [Mpc/h]
  double ds[12] = {167.484, 141.745, 109.689, 79.222, 53.9142, 34.3401, 19.7749, 10.0391, 4.9866, 2.79118, 1.71323, 1.09489};        //units of [h Msun/pc^2]  
  double dsmu[12] = {215.694, 179.938, 136.966, 97.6729, 65.1235, 39.7234, 22.1475, 11.6005, 6.06163, 3.45574, 2.09941, 1.32538};    //units of [h Msun/pc^2]
  cosmosis::DataBlock mu=0.9; //for the bin 0.8<mu<1

  A_CEN_t acen(mu);

  for (std::size_t i = 0, sz = 12; i != sz; ++i)
  {
    double const fz = (acen(mu) * T_CEN_t(rp[i]) );
    double constexpr epsrel = 1.0e-1;
    CHECK(fz == Approx(dsmu[i]/ds[i]).epsilon(epsrel));
  }
}

