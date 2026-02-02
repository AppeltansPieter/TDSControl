#include "compute_N_rhp.hpp"

#include "utils/cubic_spline.hpp"
#include "utils/grid_iterator.hpp"
#include "utils/ipow.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <complex>
#include <numbers>
#include <ranges>

using namespace std::complex_literals;

namespace tds {

static constexpr std::array<double, 33> b{
    -0.2462, -0.2266, -0.2053, -0.2302, -0.2326, -0.2335, -0.2362, -0.2421, -0.2463, -0.2604, -0.2656,
    -0.2749, -0.2919, -0.3030, -0.3140, -0.3265, -0.3491, -0.3664, -0.3892, -0.4204, -0.4548, -0.4855,
    -0.5339, -0.5872, -0.6491, -0.7354, -0.8478, -0.9930, -1.1800, -1.4448, -1.9414, -2.7149, -3.0292};

static constexpr std::array<double, 33> a{0.9124, 0.9123, 0.9136, 0.9165, 0.9195, 0.9234, 0.9285, 0.9345, 0.9416,
                                          0.9501, 0.9592, 0.9698, 0.9818, 0.9947, 1.0090, 1.0249, 1.0427, 1.0620,
                                          1.0833, 1.1069, 1.1331, 1.1614, 1.1936, 1.2289, 1.2685, 1.3132, 1.3642,
                                          1.4231, 1.4913, 1.5731, 1.6783, 1.7867, 1.8183};

/**
 * @brief Constructs the matrix sum_{k=0}^{m} A_k * exp(-(r + Xi) * n_k * tau) * exp(1j * theta) - r * I
 *        for a system with commensurate delays in which tau the base delay.
 *
 * @param sys [in] The time-delay system
 * @param A [out] Matrix storing the result
 * @param r [in] Defines right half-plane in which to compute eigenvalues.
 * @param Xi [in] Shift
 * @param i [in]
 * @param p [in]
 */
void create_matrix(const tds &sys, Eigen::MatrixXcd &A, const double r, const double Xi, const std::size_t i,
                   std::size_t p)
{
    // Assume commensurate delays with base delay 1/20 * max delay
    const double tau = *std::ranges::max_element(sys.hA()) / 20.;
    A = sys.A(0) * std::exp(-(r + Xi) * std::round(sys.hA(0) / tau)) - r * Eigen::MatrixXcd::Identity(sys.n(), sys.n());
    const double theta = 2. * std::numbers::pi * static_cast<double>(i) / static_cast<double>(p);
    for (std::size_t j = 0; j < sys.mA() - 1; ++j) {
        const double n = std::round(sys.hA(j + 1) / tau);
        // use sys.hA(j+1) instead of n * tau?
        A += sys.A(j + 1) * std::exp(-(r + Xi) * n * tau) * std::exp(theta * 1.0i * n);
    }
}

/**
 * @brief Constructs the matrix sum_{k=0}^{m} A_k * exp(-(r + Xi) * h_k) * exp(1j * theta_k) - r * I
 *        for a system with commensurate delays in which tau the base delay.
 *
 * @param sys [in] The time-delay system
 * @param A [out] Matrix storing the result
 * @param r [in] Defines right half-plane in which to compute eigenvalues.
 * @param Xi [in] Shift
 * @param i_grid [in] I-grid vector
 * @param p [in]
 */
void create_matrix(const tds &sys, Eigen::MatrixXcd &A, const double r, const double Xi,
                   const std::vector<std::size_t> &i_grid, const std::size_t p)
{
    A = sys.A(0) * std::exp(-(r + Xi) * sys.hA(0)) - r * Eigen::MatrixXcd::Identity(sys.n(), sys.n());
    for (std::size_t j = 0; j < sys.mA() - 1; ++j) {
        const double theta = static_cast<double>(i_grid[j]) * 2. * std::numbers::pi / static_cast<double>(p + 1);
        A += sys.A(j + 1) * std::exp(-(r + Xi) * sys.hA(j + 1)) * std::exp(theta * 1.0i);
    }
}

/**
 * @brief Get the correct iterator
 *
 * @tparam commensurate Whether the delays are commensurate
 * @param m The number of state-terms
 * @param p The number of discretization points
 * @return The correct iterator
 */
template <bool commensurate> auto get_grid_iterator(const std::size_t m, const std::size_t p)
{
    if constexpr (commensurate) {
        return std::views::iota(0UL, p);
    }
    else {
        return utils::ThetaRange(m - 1, p);
    }
}

/**
 * @brief Helper function to distinguish between commensurate and non-commensurate delays.
 *
 * @tparam commensurate Whether or not the system should be interpreted as having commensurate delays
 * @param sys The time delay-system
 * @param r The considered right half-plane
 * @param p The number of discretization points
 * @return The required degree of spectral discretization
 */
template <bool commensurate> std::size_t compute_N_rhp_impl(const tds &sys, const double r, const std::size_t p)
{
    const std::size_t m = sys.mA();
    Eigen::MatrixXcd A(sys.n(), sys.n());
    // TODO: advantage in place or does it only allocate once?
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> eslvr;
    std::vector<std::complex<double>> evs;
    if constexpr (commensurate) {
        evs.reserve(sys.n() * p);
    }
    else {
        evs.reserve(sys.n() * ipow(p, m - 2) * ((p / 2) + 1));
    }

    for (const auto &theta : get_grid_iterator<commensurate>(m, p)) {
        create_matrix(sys, A, r, 0., theta, p);
        eslvr.compute(A, false);
        std::ranges::copy(eslvr.eigenvalues(), std::back_inserter(evs));
    }
    const auto max_re = std::ranges::max(evs | std::views::transform([](const auto &e) -> double { return e.real(); }));
    const double Xi = 1.05 * std::sin(2.0 * std::numbers::pi / static_cast<double>(p)) * max_re;

    // TODO: Create spline at compile-time?
    auto a_spline = utils::CubicSpline(0., std::numbers::pi / 64., a);
    auto b_spline = utils::CubicSpline(0., std::numbers::pi / 64., b);
    auto compute_N = [&a_spline, &b_spline](const std::complex<double> &z) -> int {
        return std::ceil((std::abs(z) - b_spline(std::abs(std::arg(z)))) / a_spline(std::abs(std::arg(z))));
    };
    auto N_view =
        evs | std::views::filter([Xi](const std::complex<double> &el) { return el.real() >= 0 && el.real() <= Xi; }) |
        std::views::transform(compute_N);
    const auto max_it = std::ranges::max_element(N_view);
    const int N1 = (max_it != std::end(N_view)) ? *max_it : 0;

    int N2 = 0;

    for (const auto &theta : get_grid_iterator<commensurate>(m, p)) {
        create_matrix(sys, A, r, Xi, theta, p);
        eslvr.compute(A, false);
        for (const auto ev = eslvr.eigenvalues(); const auto &el : ev) {
            if (std::real(el) >= Xi) {
                N2 = std::max(compute_N(el), N2);
            }
        }
    }
    return std::max({N1, N2, 5});
}

std::size_t compute_N_rhp(const tds &sys, const double r, const std::size_t p)
{
    if (sys.mA() <= 4) {
        return compute_N_rhp_impl<false>(sys, r, p);
    }
    return compute_N_rhp_impl<true>(sys, r, p);
}
} // namespace tds