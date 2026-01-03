#pragma once

#include <Eigen/Core>
#include <tdscontrol/tds.hpp>
#include <utility>

namespace tds {

/**
 * @brief Compute a spectral discretization of a given degree for a time-delay system.
 *
 * @param system The time-delay system to discretize
 * @param N The desired degree of the spectral discretization
 * @return The matrix pencil corresponding to the spectral discretization of the tds.
 */
std::pair<Eigen::MatrixXd, Eigen::MatrixXd> sparse_spectral_discretization(const tds &system, std::size_t N);

} // namespace tds