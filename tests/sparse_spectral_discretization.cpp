#include "../src/discretization/spectral.hpp"

#include <Eigen/Core>
#include <gtest/gtest.h>
#include <tdscontrol/tds.hpp>

GTEST_TEST(sparse_spectral_discretization, scalar)
{
    Eigen::MatrixXd A0(1, 1);
    A0 << 1.;
    Eigen::MatrixXd A1(1, 1);
    A1 << -1.;
    Eigen::MatrixXd A2(1, 1);
    A2 << 2;
    const tds::tds sys({A0, A1, A2}, {0., 1., 2.});
    auto [Sigma, Pi] = sparse_spectral_discretization(sys, 5);
    Eigen::MatrixXd Sigma_expected(6, 6);
    Sigma_expected << 2., -1., 4., -1., 2., -1., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 0.,
        0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1.;
    EXPECT_LE((Sigma - Sigma_expected).norm(), 1e-12);
    Eigen::MatrixXd Pi_expected(6, 6);
    Pi_expected << 1., 1., 1., 1., 1., 1., 1., 0., -0.5, 0., 0., 0., 0., 0.25, 0., -0.25, 0., 0., 0., 0., 1. / 6., 0.,
        -1. / 6., 0., 0., 0., 0., 0.125, 0., -0.125, 0., 0., 0., 0., 0.1, 0.;
    EXPECT_LE((Pi - Pi_expected).norm(), 1e-12);
}

GTEST_TEST(sparse_spectral_discretization, two_dimensional)
{
    Eigen::MatrixXd A0(2, 2);
    A0 << 1., 2., 3., 4.;
    Eigen::MatrixXd A1(2, 2);
    A1 << -1., -2., -3., -4.;
    Eigen::MatrixXd A2(2, 2);
    A2 << 6., 7., 8., 9.;
    const tds::tds sys({A0, A1, A2}, {0., 1., 3.});
    auto [Sigma, Pi] = sparse_spectral_discretization(sys, 3);
    Eigen::MatrixXd Sigma_expected(8, 8);
    Sigma_expected << 6., 7., -16. / 3., -17. / 3., 70. / 9., 95. / 9., -112. / 27., -89. / 27., 8., 9., -6., -19. / 3.,
        120. / 9., 145. / 9., -66. / 27., -43. / 27., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,
        0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,
        0., 0., 0., 1.;
    EXPECT_LE((Sigma - Sigma_expected).norm(), 1e-12);
    Eigen::MatrixXd Pi_expected(8, 8);
    Pi_expected << 1., 0., 1., 0., 1., 0., 1., 0., 0., 1., 0., 1., 0., 1., 0., 1., 3. / 2., 0., 0., 0., -3. / 4., 0.,
        0., 0., 0., 3. / 2., 0., 0., 0., -3. / 4., 0., 0., 0., 0., 3. / 8., 0., 0., 0., -3. / 8., 0., 0., 0., 0.,
        3. / 8., 0., 0., 0., -3. / 8., 0., 0., 0., 0., 0.25, 0., 0., 0., 0., 0., 0., 0., 0., 0.25, 0., 0.;
    EXPECT_LE((Pi - Pi_expected).norm(), 1e-12);
}
