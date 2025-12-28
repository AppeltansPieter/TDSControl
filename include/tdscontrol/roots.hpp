#pragma once
#include <complex>
#include <tdscontrol/tds.hpp>
#include <vector>

namespace tds {
std::vector<std::complex<double>> roots(const tds &system, unsigned int N);
} // namespace tds
