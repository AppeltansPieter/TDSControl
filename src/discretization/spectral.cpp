#include "spectral.hpp"

#include "cheb.hpp"

#include <algorithm>

namespace tds {
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> sparse_spectral_discretization(const tds &system, const std::size_t N)
{
    // TODO: Examine effect of scaling largest delay to one and shifting origin with r?
    /* Based on Michiels, Wim, and Silviu-Iulian Niculescu,Stability, control, and computation for time-delay systems:
     * an eigenvalue-based approach. Society for Industrial and Applied Mathematics, 2014. Chapter 2.2, Computing all
     * characteristic roots in a right half plane.
     */
    const Eigen::Index n = system.n();
    const auto N_ = static_cast<Eigen::Index>(N);
    const double tau_m = *std::ranges::max_element(system.hA());
    // Set-up Sigma (See equation 2.13)
    Eigen::MatrixXd Sigma = Eigen::MatrixXd::Identity(n * (N_ + 1), n * (N_ + 1));
    Sigma.block(0, 0, n, n) = Eigen::MatrixXd::Zero(n, n);
    for (Eigen::Index i = 0; i <= N_; i++) {
        for (std::size_t k = 0; k < system.mA(); k++) {
            Sigma.block(0, i * n, n, n) += system.A()[k] * cheb(-2.0 * system.hA()[k] / tau_m + 1.0, i);
        }
    }

    // Set-up Pi (See equation 2.12)
    Eigen::MatrixXd Pi = Eigen::MatrixXd::Zero(n * (N_ + 1), n * (N_ + 1));
    for (Eigen::Index i = 1; i <= N_; i++) {
        Pi.block(0, (i - 1) * n, n, n).noalias() = 4.0 / tau_m * Eigen::MatrixXd::Identity(n, n);
        Pi.block(i * n, (i - 1) * n, n, n) = 1.0 / static_cast<double>(i) * Eigen::MatrixXd::Identity(n, n);
        if (i < N_) {
            Pi.block(i * n, (i + 1) * n, n, n) = -1.0 / static_cast<double>(i) * Eigen::MatrixXd::Identity(n, n);
        }
    }
    Pi.block(0, N_ * n, n, n).noalias() = 4.0 / tau_m * Eigen::MatrixXd::Identity(n, n);
    Pi.block(n, 0, n, n) = 2.0 * Eigen::MatrixXd::Identity(n, n);
    Pi = tau_m / 4.0 * Pi;

    return {Sigma, Pi};
}
} // namespace tds