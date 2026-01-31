#include "compute_N_rhp.hpp"
#include "discretization/spectral.hpp"
#include "newtons_method.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <tdscontrol/roots.hpp>
#include <tdscontrol/tds.hpp>

namespace tds {

/**
 * @brief Helper function whose roots correspond to the sought for eigenvalue-eigenvector pair and its Jacobian.
 *
 * @details F(l, v) = [(l*I - sum_{k=1}^m A_k exp(-l * tau_k)) * v; v'*v - 1]
 *          J(l, v) = [I + sum_{k=1}^m tau_k A_k exp(-l *tau_k) & (l*I - sum_{k=1}^m A_k exp(-l * tau_k))]
 *                    [                0                        &                 2*adj(v)               ]
 * @param sys The time-delay system
 * @param x Vector containing the eigenvalue followed by the eigenvector
 * @param fun The function to find the root of
 * @param Jac Its Jacobian
 */
void fun_evp(const tds &sys, const Eigen::VectorXcd &x, Eigen::VectorXcd &fun, Eigen::MatrixXcd &Jac)
{
    const Eigen::Index n = sys.n();
    const std::size_t m = sys.mA();
    const auto lambda = x(0);
    const auto v = x.tail(n);
    auto M = Jac.block(0, 1, n, n);
    Jac.block(0, 0, n, 1) = v;
    M = lambda * Eigen::MatrixXcd::Identity(n, n);
    for (std::size_t i = 0; i < m; ++i) {
        M -= sys.A(i) * std::exp(-lambda * sys.hA(i));
        Jac.block(0, 0, n, 1) += sys.A(i) * v * sys.hA(i) * std::exp(-lambda * sys.hA(i));
    }
    fun.head(n) = M * v;
    fun(n) = v.squaredNorm() - 1;
    Jac(n, 0) = 0.;
    Jac.block(n, 1, 1, n) = 2 * v.adjoint();
};

std::vector<std::complex<double>> roots(const tds &system, const double r)
{
    std::size_t N = compute_N_rhp(system, r);
    // TODO: make max size eigenvalue problem configurable
    N = std::min(N, 100UL);
    Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;

    auto [Sigma, Pi] = sparse_spectral_discretization(system, N);

    ges.compute(Sigma, Pi);
    std::vector<std::complex<double>> res;
    res.reserve(ges.eigenvalues().size());
    const Eigen::Index n = system.n();
    const size_t m = system.mA();
    internal::NewtonsMethod newtons_method(n + 1);
    Eigen::VectorXcd x0(n + 1);
    Eigen::MatrixXcd M(n, n);
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(M);
    for (const auto &ev : ges.eigenvalues()) {
        x0(0) = ev;
        M = ev * Eigen::MatrixXcd::Identity(n, n);
        for (std::size_t i = 0; i < m; ++i) {
            M -= system.A(i) * std::exp(-ev * system.hA(i));
        }
        svd.compute(M, Eigen::ComputeThinV);
        x0.tail(n) = svd.matrixV().col(n - 1);
        const double residual = newtons_method(x0, [&system](const Eigen::VectorXcd &x, Eigen::VectorXcd &fun,
                                                             Eigen::MatrixXcd &Jac) { fun_evp(system, x, fun, Jac); });
        if (residual < 1e-5 && std::real(x0(0)) >= r) {
            res.push_back(x0(0));
        }
    }
    return res;
}

} // namespace tds
