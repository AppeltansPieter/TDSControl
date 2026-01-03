#include "discretization/spectral.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <tdscontrol/roots.hpp>
#include <tdscontrol/tds.hpp>

namespace tds {

std::vector<std::complex<double>> roots(const tds &system, const unsigned int N)
{
    Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges;

    auto [Sigma, Pi] = sparse_spectral_discretization(system, N);

    ges.compute(Sigma, Pi);
    std::vector<std::complex<double>> res;
    Eigen::VectorXcd ev = ges.eigenvalues();
    for (Eigen::Index i = 0; i < ges.eigenvalues().size(); i++) {
        res.push_back(ges.eigenvalues()[i]);
    }
    return res;
}

} // namespace tds
