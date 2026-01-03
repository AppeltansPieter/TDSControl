#pragma once
#include <complex>
#include <tdscontrol/tds.hpp>
#include <vector>

namespace tds {
/**
 * @brief Compute the characteristic roots of a time-delay system based on a spectral discretization of a given degree
 *
 * @param system The system of which to compute the characteristic roots
 * @param N The desired degree of the spectral discretization
 * @return The computed characteristic roots
 */
std::vector<std::complex<double>> roots(const tds &system, unsigned int N);
} // namespace tds
