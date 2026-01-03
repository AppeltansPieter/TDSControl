#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <tdscontrol/roots.hpp>
#include <tdscontrol/tds.hpp>

int main(int /*argc*/, char const ** /*argv*/)
{
    auto A0 = Eigen::MatrixXd(2, 2);
    A0 << 1., 3., -7., 1.;
    auto A1 = Eigen::MatrixXd(2, 2);
    A1 << -1., 2., 8., -1.;
    constexpr double tau = 3.0;
    const tds::tds sys({A0, A1}, {0., tau});
    const auto roots = tds::roots(sys, 15);
    std::cout << std::fixed << std::setprecision(16);
    for (const auto &root : roots) {
        Eigen::MatrixXcd M = (root * Eigen::MatrixXd::Identity(2, 2) - A0 - A1 * std::exp(-tau * root));
        std::cout << "Eigenvalue: ";
        std::cout << root;
        std::cout << ", Residue: ";
        std::cout << M.determinant();
        std::cout << "\n";
    }
    std::cout << "\n";

    return 0;
}
