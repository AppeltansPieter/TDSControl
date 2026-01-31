#include <Eigen/Core>
#include <Eigen/SVD>
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
    const auto roots = tds::roots(sys, -5.);
    for (const auto root : roots) {
        const std::complex<double> M = root - 1. + std::exp(-1. * root);
        EXPECT_LE(std::abs(M), 1e-7);
    }
}
GTEST_TEST(roots, two_dimensional)
{
    auto A0 = Eigen::MatrixXd(2, 2);
    A0 << 1., -1., 2., 3.;
    auto A1 = Eigen::MatrixXd(2, 2);
    A1 << -1., 1., 6., -7.;
    const tds::tds sys({A0, A1}, {0., 3.});
    const auto roots = tds::roots(sys, -0.5);
    Eigen::MatrixXcd M(2, 2);
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(2, 2);
    for (const auto root : roots) {
        M = root * Eigen::MatrixXcd::Identity(2, 2) - A0 - A1 * std::exp(-3. * root);
        svd.compute(M);
        EXPECT_LE(svd.singularValues()[1], 1e-7);
    }
}
