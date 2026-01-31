#pragma once

#include <tdscontrol/tds.hpp>

namespace tds {

/**
 * @brief Compute the necessary number of discretization points such that the characteristic function is sufficiently
 *        well approximated in all characteristic roots in a given right half-plane.
 *
 * @param sys The time delay-system
 * @param r Defines the right half-plane { \lambda \in \C: Re(\lambda) >= r}
 * @param p Number of discretization points for \Psi
 * @return std::size_t required The number of discretization points
 */
std::size_t compute_N_rhp(const tds &sys, double r, std::size_t p = 20);
// TODO: add possibility to specify tolerance
} // namespace tds
