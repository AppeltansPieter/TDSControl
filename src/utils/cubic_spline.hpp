#pragma once

#include <span>
#include <vector>

namespace tds::utils {
/**
 * @brief Class representing cubic spline with equidistant nodes and natural boundary conditions.
 */
class CubicSpline {
  public:
    /**
     * @brief Cubic spline with equidistant nodes and natural boundary conditions.
     *
     * @param x0 First node
     * @param dx Distant between subsequent nodes
     * @param y Function values at the nodes
     */
    CubicSpline(const double x0, const double dx, const std::span<const double> &y)
        : x0_(x0), dx_(dx), a_(y.size() - 1), b_(y.size() - 1), c_(y.size() - 1), d_(y.size() - 1)
    {
        compute_coefficients(y);
    }

    /**
     * @brief Interpolate the spline at a given position.
     *
     * @param x The position at which to evaluate the spline.
     * @return double The value of the spline at the given position.
     */
    double operator()(double x) const;

  private:
    // TODO: Make constexpr? Create polynomial at compile time!
    /**
     * @brief Compute the coefficients of the polynomials.
     *
     * @param y Function values at the nodes.
     */
    void compute_coefficients(const std::span<const double> &y);

    double x0_;
    double dx_;
    std::vector<double> a_;
    std::vector<double> b_;
    std::vector<double> c_;
    std::vector<double> d_;
};

} // namespace tds::utils
