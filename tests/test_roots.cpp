#include <Eigen/Core>
#include <gtest/gtest.h>
#include <tdscontrol/roots.hpp>
#include <tdscontrol/tds.hpp>

GTEST_TEST(roots, scalar)
{
    auto A0 = Eigen::MatrixXd(1, 1);
    A0 << 1.;
    auto A1 = Eigen::MatrixXd(1, 1);
    A1 << -1.;
    const tds::tds sys({A0, A1}, {0., 1.});
    const auto roots = tds::roots(sys, 10);
    EXPECT_EQ(roots.size(), 11);
}
GTEST_TEST(roots, two_dimensional)
{
    auto A0 = Eigen::MatrixXd(2, 2);
    A0 << 1., -1., 2., 3.;
    auto A1 = Eigen::MatrixXd(2, 2);
    A1 << -1., 1., 6., -7.;
    const tds::tds sys({A0, A1}, {0., 3.});
    const auto roots = tds::roots(sys, 10);
    EXPECT_EQ(roots.size(), 22);
}
