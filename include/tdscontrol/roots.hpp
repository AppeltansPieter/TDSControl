#pragma once
#include <complex>
#include <tdscontrol/tds.hpp>
#include <vector>

namespace tds {
/**
 * @brief Compute the characteristic roots of a time-delay system in a given right half-plane.
 *
 * @param system The system of which to compute the characteristic roots
 * @param r Defines the right half-plane
 * @return The computed characteristic roots
 */
std::vector<std::complex<double>> roots(const tds &system, double r);
} // namespace tds
