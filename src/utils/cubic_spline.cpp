#include "cubic_spline.hpp"

#include "thomas_algorithm.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace tds::utils {

double CubicSpline::operator()(const double x) const
{
    assert(x >= x0_);
    assert(x <= x0_ + (dx_ * a_.size()));

    size_t ind = std::floor((x - x0_) / dx_);
    if (ind >= a_.size()) {
        ind = a_.size() - 1;
    }
    const double dx = x - x0_ - (static_cast<double>(ind) * dx_);
    return a_.at(ind) + (b_.at(ind) * dx) + (c_.at(ind) * std::pow(dx, 2)) + (d_.at(ind) * std::pow(dx, 3));
}

void CubicSpline::compute_coefficients(const std::span<const double> &y)
{
    // Create system
    // [4*dx  dx   0 0 ... 0]
    // [dx   4*dx dx 0 ... 0]
    // [        ...         ]
    // [   0  ...  0 dx 4*dx]
    std::vector<double> sub(y.size() - 3);
    std::ranges::fill(sub, dx_);
    std::vector<double> diag(y.size() - 2);
    std::ranges::fill(diag, 4 * dx_);
    std::vector<double> super(y.size() - 3);
    std::ranges::fill(super, dx_);
    std::vector<double> rhs(y.size() - 2);

    for (size_t i = 0; i < rhs.size(); ++i) {
        rhs[i] = 3.0 * (y[i] - (2 * y[i + 1]) + y[i + 2]) / dx_;
    }
    thomas_algorithm(sub, diag, super, rhs);
    // Extract coefficients of the different polynomials.
    c_[0] = 0.0;
    std::ranges::copy(rhs, std::next(c_.begin()));
    std::copy(y.begin(), std::prev(y.end()), a_.begin());
    for (size_t i = 0; i < c_.size() - 1; ++i) {
        d_[i] = (c_[i + 1] - c_[i]) / (3. * dx_);
        b_[i] = (y[i + 1] - y[i]) / dx_ - (c_[i + 1] + 2. * c_[i]) / 3. * dx_;
    }
    d_.back() = -c_.back() / (3. * dx_);
    b_.back() = ((y.back() - *std::prev(y.end(), 2UL)) / dx_) - (2. * c_.back() / 3. * dx_);
}

} // namespace tds::utils