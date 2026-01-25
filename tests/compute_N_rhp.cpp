#include "../src/compute_N_rhp.hpp"

#include "convert_reader_data.hpp"

#include <gtest/gtest.h>
#include <tdsreader/reader.hpp>
#include <tdsreader/system.hpp>

GTEST_TEST(test_compute_n_rhp, non_commensurate)
{
    tdsreader::tds data = tdsreader::read("TDSControl-examples/data/VerheydenEtAl2008");
    tds::tds sys = create_tds_from_data(std::move(data));
    const std::size_t N = tds::compute_N_rhp(sys, -2.5, 20);
    sys.A(1) *= std::exp(2.5 * sys.hA(1));
    sys.A(0) += 2.5 * Eigen::MatrixXd::Identity(sys.n(), sys.n());
    EXPECT_EQ(tds::compute_N_rhp(sys, 0., 20), N);
}