#include "../src/newtons_method.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <tdsreader/reader.hpp>
#include <tdsreader/system.hpp>
#include <vector>

using namespace tds::internal;

void check_quadratic_convergence(const std::vector<double> &residuals)
{
    const std::size_t nb_iterations = residuals.size();
    const double e_next = residuals[nb_iterations - 1];
    const double e_curr = residuals[nb_iterations - 2];
    const double e_prev = residuals[nb_iterations - 3];
    EXPECT_GE(std::log(e_next / e_curr) / std::log(e_curr / e_prev), 1.8);
}

GTEST_TEST(Newtons_method, scalar)
{
    NewtonsMethod solver(1);
    Eigen::VectorXcd res(1);
    res << 1.2;
    solver(res, [](const Eigen::VectorXcd &x, Eigen::VectorXcd &fun, Eigen::MatrixXcd &Jac) {
        fun(0) = std::cos(x(0));
        Jac(0, 0) = -std::sin(x(0));
    });
    EXPECT_NEAR(res(0).real(), std::numbers::pi / 2.0, 1e-6);
    EXPECT_LE(std::abs(std::cos(res(0))), 1e-10);
    check_quadratic_convergence(solver.get_residuals());
}

GTEST_TEST(Newtons_method, linear)
{
    NewtonsMethod solver(3);
    Eigen::VectorXcd res(3);
    res << 0.5, -0.5, 0;
    Eigen::MatrixXcd A(3, 3);
    A << 1., 2., 3., -2., 3., 4., -6., 7., -8;
    Eigen::VectorXcd b(3);
    b << 4., 3, 2.;
    solver(
        res,
        [&A, &b](const Eigen::VectorXcd &x, Eigen::VectorXcd &fun, Eigen::MatrixXcd &Jac) {
            fun = A * x - b;
            Jac = A;
        },
        Tolerance(1e-10));

    // Check backward error
    EXPECT_LE((A * res - b).norm(), 1e-10);

    // Check forward error
    EXPECT_LE((res - A.colPivHouseholderQr().solve(b)).norm(), 1e-10);
    EXPECT_EQ(solver.get_residuals().size(), 2);
}

void example_two_dimensional(const Eigen::VectorXcd &x, Eigen::VectorXcd &fun, Eigen::MatrixXcd &Jac)
{
    fun << std::sin(x(0)), std::cos(x(1));
    Jac << std::cos(x(0)), 0., 0., -std::sin(x(1));
}

GTEST_TEST(Newtons_method, two_dimensional)
{
    constexpr Tolerance tol(1e-12);
    NewtonsMethod solver(2);
    Eigen::VectorXcd res(2);
    res << 1.1, -0.5;
    solver(res, &example_two_dimensional, tol);

    // Check backward error
    Eigen::VectorXcd residual(2);
    Eigen::MatrixXcd Jac(2, 2);
    example_two_dimensional(res, residual, Jac);
    EXPECT_LE(residual.norm(), tol.value);

    // Check forward error
    Eigen::VectorXcd expected(2);
    expected << 0., -std::numbers::pi / 2.;
    EXPECT_LE((res - expected).norm(), 1e-8);
    check_quadratic_convergence(solver.get_residuals());
}

GTEST_TEST(Newtons_method, diverge)
{
    constexpr MaxIter nb_iterations{20};
    NewtonsMethod solver(1);
    Eigen::VectorXcd res(1);
    res << 0.5;
    solver(
        res,
        [](const Eigen::VectorXcd &x, Eigen::VectorXcd &fun, Eigen::MatrixXcd &Jac) {
            fun(0) = std::pow(x(0), 1. / 3.);
            Jac(0, 0) = 1. / 3. * std::pow(x(0), -2. / 3.);
        },
        Tolerance(1e-10), nb_iterations);
    EXPECT_EQ(solver.get_residuals().size(), nb_iterations.value + 1);
}

GTEST_TEST(Newtons_method, no_step)
{
    NewtonsMethod solver(2);
    Eigen::VectorXcd res(2);
    Eigen::VectorXcd orig(2);
    res << 0., -std::numbers::pi / 2.;
    orig = res;
    solver(res, &example_two_dimensional);

    EXPECT_EQ(solver.get_residuals().size(), 1);
    EXPECT_EQ((orig - res).norm(), 0.);
}

GTEST_TEST(Newtons_method, multiple_root)
{
    constexpr MaxIter nb_iterations{100};
    constexpr Tolerance tol{1e-12};
    NewtonsMethod solver(3);
    Eigen::VectorXcd res(3);
    res << 0.5, -0.5, 0.5;
    auto fun = [](const Eigen::VectorXcd &x, Eigen::VectorXcd &fun, Eigen::MatrixXcd &Jac) {
        fun << std::sin(x(0)), std::pow(std::sin(x(1)), 2) * (x(0) - 1.), x(2) * (x(0) - 1.);
        Jac << std::cos(x(0)), 0., 0., std::pow(std::sin(x(1)), 2), 2. * (x(0) - 1.) * std::sin(x(1)) * std::cos(x(1)),
            0., x(2), 0., x(0) - 1.;
    };
    solver(res, fun, tol, nb_iterations);

    // Check backward error
    Eigen::VectorXcd residual(3);
    Eigen::MatrixXcd Jac(3, 3);
    fun(res, residual, Jac);
    EXPECT_LE(residual.norm(), tol.value);

    // Check forward error
    EXPECT_LE(res.norm(), 1e-6);
}
using namespace std::complex_literals;

GTEST_TEST(Newtons_method, evp)
{
    // TODO: How best run relative to root directory => check AURA?
    auto data = tdsreader::read("TDSControl-examples/data/VerheydenEtAl2008");
    Eigen::Index n = data.A.at(0).rows();
    auto f = [&data, n](const Eigen::VectorXcd &x, Eigen::VectorXcd &fun, Eigen::MatrixXcd &Jac) {
        const auto lambda = x(0);
        const auto v = x.tail(n);
        fun.head(n) = (lambda * Eigen::MatrixXcd::Identity(n, n) - data.A.at(0) - data.A.at(1) * std::exp(-lambda)) * v;
        fun(n) = v.squaredNorm() - 1;
        Jac.block(0, 0, n, 1) = (Eigen::MatrixXcd::Identity(n, n) + data.A.at(1) * std::exp(-x(0))) * v;
        Jac.block(0, 1, n, n) =
            lambda * Eigen::MatrixXcd::Identity(n, n) - data.A.at(0) - data.A.at(1) * std::exp(-x(0));
        Jac(n, 0) = 0.;
        Jac.block(n, 1, 1, n) = 2 * v.adjoint();
    };

    NewtonsMethod solver(n + 1);
    Eigen::VectorXcd res(n + 1);
    res(0) = 0.6;
    res.tail(n) = Eigen::VectorXcd::Random(n);
    res.tail(n) /= res.tail(n).norm();
    solver(res, f);
    const std::complex<double> lambda = res(0);
    const auto v = res.tail(n);
    // Backward error
    EXPECT_LE(
        ((lambda * Eigen::MatrixXcd::Identity(n, n) - data.A.at(0) - data.A.at(1) * std::exp(-lambda)) * v).norm(),
        1e-10);
    EXPECT_NEAR(v.norm(), 1., 1e-10);
}
